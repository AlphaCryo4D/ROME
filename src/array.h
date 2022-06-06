/***************************************************************************
 *
 * Authors: "Yongbei(Glow) Ma,Jiayi (Timmy) Wu, Youdong (Jack) Mao"
 * Dana-Farber Cancer Institute, Harvard Medical School and Peking University
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * This complete copyright notice must be included in any revised version of the
 * source code. Additional authorship citations may be added, but existing
 * author citations must be preserved.
 ***************************************************************************/

#ifndef ARRAY_H_
#define ARRAY_H_

#include <iostream>
#include "mkl.h"

#include "./memory.h"
//
#define ACCESS(Vol,k,i,j) (Vol.ptr[k*Vol.dimyx+i*Vol.dimx+j])
// 3D volume data structure
template<typename T>
class VolBase{
public:
    size_t dimx,dimy,dimz,dimyx,dimzyx;
    T* ptr;
    bool inHBM;
    VolBase():ptr(nullptr),dimx(0),dimy(0),dimz(0),dimyx(0),dimzyx(0),inHBM(false){}
    ~VolBase(){fini();}
    void fini(){
        if (ptr != nullptr) freeMemory(ptr, inHBM);ptr = nullptr;
        dimx = 0;dimy = 0;dimz = 0;dimyx = 0;dimzyx = 0;
    }
	T const* rptr() { return ptr; }
	T*       wptr() { return ptr; }
    virtual T& operator()(size_t k,size_t i,size_t j) = 0;
    virtual void fill(T t) = 0;
    virtual void zero()  = 0;
    void writeToDisk(std::string fn,bool do_fini = true){
        FILE* pFile;
        ERROR_CHECK((pFile = fopen((fn+".vol").c_str(), "wb")) == NULL,"open file fail.");
        ERROR_CHECK(fwrite(&dimx, sizeof(size_t), 1, pFile) != 1,"write file fail");
        ERROR_CHECK(fwrite(&dimy, sizeof(size_t), 1, pFile) != 1,"write file fail");
        ERROR_CHECK(fwrite(&dimz, sizeof(size_t), 1, pFile) != 1,"write file fail");
        ERROR_CHECK(fwrite(&dimyx, sizeof(size_t), 1, pFile) != 1,"write file fail");
        ERROR_CHECK(fwrite(&dimzyx, sizeof(size_t), 1, pFile) != 1,"write file fail");
        ERROR_CHECK(fwrite((char*)ptr, dimzyx*sizeof(T), 1, pFile) != 1,"write file fail");
        fclose(pFile);
        if(do_fini) fini();
    }
    void readFromDisk(std::string fn,bool do_init = true){
        FILE* pFile;
        ERROR_CHECK((pFile = fopen((fn+".vol").c_str(), "rb")) == NULL,"open file fail.");
        ERROR_CHECK(fread(&dimx, sizeof(size_t), 1, pFile) != 1,"write file fail");
        ERROR_CHECK(fread(&dimy, sizeof(size_t), 1, pFile) != 1,"write file fail");
        ERROR_CHECK(fread(&dimz, sizeof(size_t), 1, pFile) != 1,"write file fail");
        ERROR_CHECK(fread(&dimyx, sizeof(size_t), 1, pFile) != 1,"write file fail");
        ERROR_CHECK(fread(&dimzyx, sizeof(size_t), 1, pFile) != 1,"write file fail");
        /*assert(inHBM==false);inHBM==false;*/
        if(do_init) ptr = (T*)allocMemory(sizeof(T)*dimzyx, 64, inHBM, __FILE__, __LINE__);
        ERROR_CHECK(fread((char*)ptr, dimzyx*sizeof(T), 1, pFile) != 1,"write file fail");
        fclose(pFile);
    }
};

template<typename T>
class Vol : public VolBase<T> {
public:
    Vol(){}
    Vol(size_t _dimz,size_t _dimy,size_t _dimx,bool _inHBM = false){init(_dimz,_dimy,_dimx,_inHBM);}
    void init(size_t _dimz,size_t _dimy,size_t _dimx,bool _inHBM = false){
        this->dimx = _dimx;this->dimy = _dimy;this->dimz = _dimz;
        this->dimyx = _dimy*_dimx;this->dimzyx = _dimz*_dimy*_dimx;
        this->inHBM = _inHBM;
        this->ptr = (T*)allocMemory(sizeof(T)*_dimx*_dimy*_dimz, 64, this->inHBM, __FILE__, __LINE__);
        fill(0.);
    }
    void fill(T v){
        #pragma vector aligned
        for (size_t i = 0; i < this->dimzyx; i++)
            this->ptr[i] = v;
    }
    void zero(){fill(0.);}
    T& operator()(size_t k,size_t i,size_t j){
        return this->ptr[k*this->dimyx+i*this->dimx+j];
    }
};

template<>
class Vol<MKL_Complex16> : public VolBase<MKL_Complex16> {
public:
    Vol(){}
    Vol(size_t _dimz,size_t _dimy,size_t _dimx,bool _inHBM = false){init(_dimz,_dimy,_dimx,_inHBM);}
    void init(size_t _dimz,size_t _dimy,size_t _dimx,bool _inHBM = false){
        this->dimx = _dimx;this->dimy = _dimy;this->dimz = _dimz;
        this->dimyx = _dimy*_dimx;this->dimzyx = _dimz*_dimy*_dimx;
        this->inHBM = _inHBM;
        this->ptr = (MKL_Complex16*)allocMemory(sizeof(MKL_Complex16)*_dimx*_dimy*_dimz, 64, this->inHBM, __FILE__, __LINE__);
        MKL_Complex16 v;
        v.real = v.imag = 0;
        fill(v);
    }
    void fill(MKL_Complex16 v){
        for (size_t i = 0; i < this->dimzyx; i++){
            this->ptr[i].real = v.real;
            this->ptr[i].imag = v.imag;
        }
    }
    void zero(){MKL_Complex16 v;v.real = v.imag = 0;fill(v);}
    MKL_Complex16& operator()(size_t k,size_t i,size_t j){
        return this->ptr[k*this->dimyx+i*this->dimx+j];
    }
};

template<>
class Vol<MKL_Complex8> : public VolBase<MKL_Complex8> {
public:
    Vol(){}
    Vol(size_t _dimz,size_t _dimy,size_t _dimx,bool _inHBM = false){init(_dimz,_dimy,_dimx,_inHBM);}
    void init(size_t _dimz,size_t _dimy,size_t _dimx,bool _inHBM = false){
        this->dimx = _dimx;this->dimy = _dimy;this->dimz = _dimz;
        this->dimyx = _dimy*_dimx;this->dimzyx = _dimz*_dimy*_dimx;
        this->inHBM = _inHBM;
        this->ptr = (MKL_Complex8*)allocMemory(sizeof(MKL_Complex8)*_dimx*_dimy*_dimz, 64, this->inHBM, __FILE__, __LINE__);
        MKL_Complex8 v;
        v.real = v.imag = 0;
        fill(v);
    }
    void fill(MKL_Complex8 v){
        for (size_t i = 0; i < this->dimzyx; i++){
            this->ptr[i].real = v.real;
            this->ptr[i].imag = v.imag;
        }
    }
    void zero(){MKL_Complex8 v;v.real = v.imag = 0;fill(v);}
    MKL_Complex8& operator()(size_t k,size_t i,size_t j){
        return this->ptr[k*this->dimyx+i*this->dimx+j];
    }
};

// the sub volume same as Vol
// but it can save memory and but so same memory access method like Vol
template<typename T>
class VolSub : public VolBase<T>{
public:
    struct {// coordinate
        double startx,starty,startz,endx,endy,endz;
        size_t originx,originy,originz;
    } coords;
    struct {// cube maybe not the best choice
        size_t size,size2,size3,tilex,tiley,tilez,tileyx,tilezyx;
    } cube;
    VolSub(){}
    VolSub(size_t _dimz,size_t _dimy,size_t _dimx,size_t _cube_size,bool _inHBM = false){init(_dimz,_dimy,_dimx,_cube_size,_inHBM);}
    void init(size_t _dimz,size_t _dimy,size_t _dimx,size_t _cube_size,bool _inHBM = false){
        this->dimx = _dimx;this->dimy = _dimy;this->dimz = _dimz;
        this->dimyx = _dimy*_dimx;this->dimzyx = _dimz*_dimy*_dimx;
        cube.size = _cube_size;cube.size2 = _cube_size*_cube_size;cube.size3 = _cube_size*_cube_size*_cube_size;
        cube.tilex = ceil(double(_dimx)/_cube_size);cube.tiley = ceil(double(_dimy)/_cube_size);cube.tilez = ceil(double(_dimz)/_cube_size);
        cube.tileyx = cube.tilex*cube.tiley;cube.tilezyx = cube.tilex*cube.tiley*cube.tilez;
        inHBM = _inHBM;
        this->ptr = (T*)allocMemory(sizeof(T)*cube.size3, 64, inHBM);
        fill(0.);
    }
    void setCoords(size_t cube_index){
        assert(cube_index<cube.tilezyx);
        size_t cube_tilez_index = cube_index / cube.tileyx;
        size_t cube_tiley_index = (cube_index % cube.tileyx) / cube.tilex;
        size_t cube_tilex_index = (cube_index % cube.tileyx) % cube.tilex;
        coords.startz = cube_tilez_index*cube.size;
        coords.endz = (cube_tilez_index+1)*cube.size < this->dimz?(cube_tilez_index+1)*cube.size:this->dimz;
        coords.starty = cube_tiley_index*cube.size;
        coords.endy = (cube_tiley_index+1)*cube.size < this->dimy?(cube_tiley_index+1)*cube.size:this->dimy;
        coords.startx = cube_tilex_index*cube.size;
        coords.endx = (cube_tilex_index+1)*cube.size < this->dimx?(cube_tilex_index+1)*cube.size:this->dimx;
        coords.originx = coords.originy = coords.originz = 0;
    }
    void fill(T t){for (size_t i = 0; i < cube.size3; i++) this->ptr[i] = t;}
    void zero(){fill(0);}
    // TDOO
    inline bool inSub(double z,double y,double x){return !(z<coords.startz || z>=coords.endz ||
                                                           y<coords.starty || y>=coords.endy ||
                                                           x<coords.startx || x>=coords.endx);}
    //
    T& operator()(size_t k,size_t i,size_t j){
        k = k - coords.startz;i = i - coords.starty;j = j - coords.startx;
        return this->ptr[k*cube.size2+i*cube.size+j];
    }
    //
    void add(size_t k,size_t i,size_t j,T v){
        if (this->inSub(k,i,j)) {
            // std::cout<<"setnew "<<k<<" "<<i<<" "<<j<<std::endl;
            this->operator()(k,i,j) += v;
        }
    }
    
    void unitTestCorrection(){
        size_t dimz = 223,dimy = 117,dimx = 87;
        size_t originz = dimz/5,originy = dimy/2,originx = dimx/3;
        size_t cube_size = 40;
        size_t points_number = 200*200;
        size_t image_number = 100;
        //
        int maxthreads = omp_get_max_threads();
        std::cout<<"max number of threads : "<<maxthreads<<std::endl;
        VolSub<double> testVolSub[maxthreads];
        for (int thread = 0; thread < maxthreads; thread++) {
            testVolSub[thread].init(dimz, dimy, dimx, cube_size);
        }
        Vol<double> testVol,testVolBaseline;
        testVol.init(dimz, dimy, dimx);
        testVolBaseline.init(dimz, dimy, dimx);
        std::cout<<"number of cubes : "<<testVolSub[0].cube.tilezyx<<std::endl;
        
        std::vector<double> pointx(image_number*points_number),pointy(image_number*points_number),pointz(image_number*points_number);
        size_t x0,x1,y0,y1,z0,z1;
        for (size_t i = 0; i < image_number*points_number; i++) {
            pointx[i] = (double)rand() / RAND_MAX * dimz;
            pointy[i] = (double)rand() / RAND_MAX * dimy;
            pointz[i] = (double)rand() / RAND_MAX * dimx;
            // set some boandary
            if(pointx[i] < 5) pointx[i] = 0;if(pointx[i] > dimx-5) pointx[i] = dimx-1;
            if(pointy[i] < 5) pointy[i] = 0;if(pointy[i] > dimy-5) pointy[i] = dimy-1;
            if(pointz[i] < 5) pointz[i] = 0;if(pointz[i] > dimz-5) pointz[i] = dimz-1;
            pointx[i] -= originx;pointy[i] -= originy;pointz[i] -= originz;
            // std::cout<<pointx[i]<<" "<<pointy[i]<<" "<<pointz[i]<<std::endl;
        }
        auto notBoundary=[&](size_t k,size_t i,size_t j){
            if (k == dimz) return false;
            if (i == dimy) return false;
            if (j == dimx) return false;
            // std::cout<<"set "<<k<<" "<<i<<" "<<j<<std::endl;
            return true;
        };
        double start = dtime();
        // get the baseline
        for (size_t iimage = 0; iimage < image_number; iimage++){
            for (size_t i = 0; i < points_number; i++) {
                x0 = floor(pointx[iimage*points_number+i]) + originx;x1 = x0 + 1;
                y0 = floor(pointy[iimage*points_number+i]) + originy;y1 = y0 + 1;
                z0 = floor(pointz[iimage*points_number+i]) + originz;z1 = z0 + 1;
                if(notBoundary(z0,y0,x0)) testVolBaseline(z0,y0,x0) += 1.;
                if(notBoundary(z1,y0,x0)) testVolBaseline(z1,y0,x0) += 1.;
                if(notBoundary(z0,y1,x0)) testVolBaseline(z0,y1,x0) += 1.;
                if(notBoundary(z0,y0,x1)) testVolBaseline(z0,y0,x1) += 1.;
                if(notBoundary(z0,y1,x1)) testVolBaseline(z0,y1,x1) += 1.;
                if(notBoundary(z1,y0,x1)) testVolBaseline(z1,y0,x1) += 1.;
                if(notBoundary(z1,y1,x0)) testVolBaseline(z1,y1,x0) += 1.;
                if(notBoundary(z1,y1,x1)) testVolBaseline(z1,y1,x1) += 1.;
            }
        }
        double end = dtime();
        std::cout<<"NO multi-threads : "<<(end-start)<<" seconds"<<std::endl;
        
        start = dtime();
        // testing for sub vol
#pragma omp parallel for schedule(dynamic) private(x0,x1,y0,y1,z0,z1)
        for (size_t cube_index = 0; cube_index < testVolSub[0].cube.tilezyx; cube_index++)
        {
            int tid = omp_get_thread_num();
            testVolSub[tid].setCoords(cube_index);
            testVolSub[tid].fill(0);
            for (size_t iimage = 0; iimage < image_number; iimage++){
                for (size_t i = 0; i < points_number; i++) {
                    
                    x0 = floor(pointx[iimage*points_number+i]) + originx;x1 = x0 + 1;
                    y0 = floor(pointy[iimage*points_number+i]) + originy;y1 = y0 + 1;
                    z0 = floor(pointz[iimage*points_number+i]) + originz;z1 = z0 + 1;
                    testVolSub[tid].add(z0,y0,x0,1);testVolSub[tid].add(z1,y0,x0,1);
                    testVolSub[tid].add(z0,y1,x0,1);testVolSub[tid].add(z0,y0,x1,1);
                    testVolSub[tid].add(z0,y1,x1,1);testVolSub[tid].add(z1,y0,x1,1);
                    testVolSub[tid].add(z1,y1,x0,1);testVolSub[tid].add(z1,y1,x1,1);
                }
            }
            // add to full vol
            for (size_t k = 0; k < testVolSub[tid].dimz; k++)
                for (size_t i = 0; i < testVolSub[tid].dimy; i++)
                    for (size_t j = 0; j < testVolSub[tid].dimx; j++)
                        if(testVolSub[tid].inSub(k,i,j)) {
                            testVol(k,i,j) = testVolSub[tid](k,i,j);
                        }
            
            if (tid == 0) {
                std::cout<<"done cube_index ... "<<cube_index<<std::endl;
            }
        }
        end = dtime();
        std::cout<<"multi-threads and vol-sub : "<<(end-start)<<" seconds"<<std::endl;
        
        // compare with the baseline
        for (size_t k = 0; k < testVol.dimz; k++)
            for (size_t i = 0; i < testVol.dimy; i++)
                for (size_t j = 0; j < testVol.dimx; j++){
                    if(testVol(k,i,j) != testVolBaseline(k,i,j)){
                        std::cout<<"wrong index "<<k<<" "<<i<<" "<<j<<std::endl;
                        std::cout<<"wrong testVol["<<"]="<<testVol(k,i,j)<<__FILE__<<__LINE__<<std::endl;
                        std::cout<<"wrong testVolBaseline["<<"]="<<testVolBaseline(k,i,j)<<__FILE__<<__LINE__<<std::endl;
                        exit(1);
                    }
                }
        
    }
};


class VolSlice{
    
};

template<typename T>
class ArrayOperate{
public:
    T* ptr;size_t size;
    ArrayOperate(T* _ptr,size_t _size){ptr=_ptr;size = _size;}
    ~ArrayOperate(){}
    void operator=(ArrayOperate const & rhs) {
        assert(this->size == rhs.size);
#pragma vector aligned
        for (int i = 0; i < this->size; i++)
            this->ptr[i] = rhs.ptr[i];
    }
    void operator+=(ArrayOperate const & rhs) {
        assert(this->size == rhs.size);
#pragma vector aligned
        for (int i = 0; i < this->size; i++)
            this->ptr[i] += rhs.ptr[i];
    }
    void operator-=(ArrayOperate const & rhs) {
        assert(this->size == rhs.size);
#pragma vector aligned
        for (int i = 0; i < this->size; i++)
            this->ptr[i] += rhs.ptr[i];
    }
};

// 1D aligned data structure
template<typename T>
class Aligned1dArray{
public:
    Aligned1dArray():ptr(nullptr),dimy(0),length(0),inHBM(false){}
    Aligned1dArray(size_t _dimy,bool _inHBM = false):ptr(nullptr) {
        init(_dimy,_inHBM);
    }
    ~Aligned1dArray(){fini();}
    void init(size_t _dimy,bool _inHBM = false){
        fini();
        assert(_dimy > 0);
        dimy = _dimy;length = dimy;
        inHBM = _inHBM;
        ptr = (T*)allocMemory(sizeof(T)*length, 64, inHBM, __FILE__, __LINE__);
    }
    void fini(){
        if (ptr) freeMemory(ptr, inHBM);ptr = nullptr;inHBM = false;
        dimy = 0;length = 0;
    }
    void fill(T v){
        for (size_t i = 0; i < dimy; i++)
            ptr[i] = v;
    }
    void fill_with_first_touch(T v){
        assert(false);
    }
    // access and modify the data
    T* wptr(size_t i){
        assert(i < dimy);assert(size_t(ptr)%64==0);
        return ptr;
    }
    const T* rptr(size_t i) const{
        assert(i < dimy);assert(size_t(ptr)%64==0);
        return ptr;
    }
    T& operator[](int i) const {
        assert(0 <= i && i < dimy);
        return _ptr[i];
    }
    size_t size(){return length;}
private:
    size_t dimy;
    size_t length;
    T* ptr;
    bool inHBM;
};

// 2D padding-aligned data structure
template<typename T>
class Aligned2dArray{
public:
    Aligned2dArray():ptr(nullptr),dimy(0),dimx(0),dimx_pad(0),length(0),inHBM(false){}
    Aligned2dArray(size_t _dimy,size_t _dimx,bool _inHBM = false):ptr(nullptr) {
        init(_dimy,_dimx,_inHBM);
    }
    ~Aligned2dArray(){fini();}
    void init(size_t _dimy,size_t _dimx,bool _inHBM = false){
        fini();
        assert(_dimy > 0);assert(_dimx > 0);
        dimy = _dimy;dimx = _dimx;
        assert(64%sizeof(T)==0);
        dimx_pad = ( ( (_dimx*sizeof(T)+63) / 64 ) * 64 ) / sizeof(T);
        length = dimy*dimx_pad;
        inHBM = _inHBM;
        ptr = (T*)allocMemory(sizeof(T)*length, 64, inHBM, __FILE__, __LINE__);
    }
    void fini(){
        if (ptr) freeMemory(ptr, inHBM);ptr = nullptr;inHBM = false;
        dimy = 0;dimx = 0;dimx_pad = 0;length = 0;
    }
    void fill(T v){
        for (size_t i = 0; i < dimy; i++)
        {
            for (size_t j = 0; j < dimx_pad; j++)
                ptr[i*dimx_pad+j] = v;
        }
    }
    void fill_with_first_touch(T v){
#pragma omp parallel for
        for (size_t i = 0; i < dimy; i++)
        {
            for (size_t j = 0; j < dimx_pad; j++)
                ptr[i*dimx_pad+j] = v;
        }
    }
    size_t size(){return length;}
    T* wptr(size_t i,size_t _dimx){
        assert(i < dimy);assert(_dimx <= dimx);
        auto aligned_ptr = ptr+i*dimx_pad;
        assert(size_t(aligned_ptr)%64==0);
        return aligned_ptr;
    }
    const T* rptr(size_t i,size_t _dimx) const{
        assert(i < dimy);assert(_dimx <= dimx);
        auto aligned_ptr = ptr+i*dimx_pad;
        assert(size_t(aligned_ptr)%64==0);
        return aligned_ptr;
    }
private:
    size_t dimy,dimx,dimx_pad;
    size_t length;
    T* ptr;
    bool inHBM;
};

// 3D padding-aligned data structure
template<typename T>
class Aligned3dArray{
public:
    Aligned3dArray():ptr(nullptr),dimz(0),dimy(0),dimx(0),dimx_pad(0),length(0),inHBM(false){}
    Aligned3dArray(size_t _dimz,size_t _dimy,size_t _dimx,bool _inHBM = false):ptr(nullptr) {
        init(_dimz,_dimy,_dimx,_inHBM);
    }
    ~Aligned3dArray(){fini();}
    void init(size_t _dimz,size_t _dimy,size_t _dimx,bool _inHBM = false){
        fini();
        assert(_dimz > 0);assert(_dimy > 0);assert(_dimx > 0);
        dimz = _dimz;dimy = _dimy;dimx = _dimx;
        assert(64%sizeof(T)==0);
        dimx_pad = ( ( (_dimx*sizeof(T)+63) / 64) * 64 ) / sizeof(T);
        length = dimz*dimy*dimx_pad;
        inHBM = _inHBM;
        ptr = (T*)allocMemory(sizeof(T)*length, 64, inHBM, __FILE__, __LINE__);
    }
    void fini(){
        if (ptr) freeMemory(ptr, inHBM);ptr = nullptr;inHBM = false;
        dimz = 0;dimy = 0;dimx = 0;dimx_pad = 0;length = 0;
    }
    void fill(T v){
        for (size_t k = 0; k < dimz; k++)
        {
            for (size_t i = 0; i < dimy; i++)
                for (size_t j = 0; j < dimx_pad; j++)
                    ptr[k*dimy*dimx_pad+i*dimx_pad+j] = v;
        }
    }
    void fill_with_first_touch(T v){
#pragma omp parallel for
        for (size_t k = 0; k < dimz; k++)
        {
            for (size_t i = 0; i < dimy; i++)
                for (size_t j = 0; j < dimx_pad; j++)
                    ptr[k*dimy*dimx_pad+i*dimx_pad+j] = v;
        }
    }
    size_t size(){return length;}
    T* wptr(size_t k,size_t i,size_t _dimx){
        assert(k < dimz);assert(i < dimy);assert(_dimx<=dimx);
        auto aligned_ptr = ptr+k*dimy*dimx_pad+i*dimx_pad;
        assert(size_t(aligned_ptr)%64==0);
        return aligned_ptr;
    }
    const T* rptr(size_t k,size_t i,size_t _dimx) const{
        assert(k < dimz);assert(i < dimy);assert(_dimx<=dimx);
        auto aligned_ptr = ptr+k*dimy*dimx_pad+i*dimx_pad;
        assert(size_t(aligned_ptr)%64==0);
        return aligned_ptr;
    }
private:
    T* ptr;
    size_t dimz,dimy,dimx,dimx_pad;
    size_t length;
    bool inHBM;
};
#endif /* defined(ARRAY_H_) */
