/***************************************************************************
 *
 * Authors: "Yongbei(Glow) Ma,Jiayi (Timmy) Wu, Youdong (Jack) Mao"
 * Dana-Farber Cancer Institute, Harvard Medical School and Peking University
 * Authors: "Brett, Bevin"
 * Intel Corporation
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

#ifndef MATRIX_H_
#define MATRIX_H_

#include "./resmap_macros.h"
#include "./resmap_error.h"
//
// Forward declarations
template<typename T>
class Matrix1D;
template<typename T>
class Matrix2D;
//
template<typename T>
class Matrix1D{
private:
    std::vector<T> __data;
    int __dim;
    bool __isrow;// is a row vector
public:
    Matrix1D():__dim(0),__isrow(false){}
    Matrix1D(int dim):__isrow(false){resize(dim);}
    Matrix1D(const Matrix1D& rhs) = default;
    ~Matrix1D(){__data.resize(0);}
    void resize(int dim){
        //assert(__data.size()==0);
        __data.resize(dim);
        __dim=dim;
    }
    void initConstant(T v){
        for (auto& _data : __data) _data = v;
    }
    void initZeros(){
        for (auto& _data : __data) _data = 0;
    }
    // do nothing,change the 1D-Matrix between column and row
    Matrix1D<T> transpose(){
        __isrow=!__isrow;
        Matrix1D<T> result(*this);
        __isrow=!__isrow;
        return result;
    }
    // Normalize this vector, store the result here
    void selfNormalize(){
        FDOUBLE sum2 = 0;
        for (auto& _data : __data) sum2 += _data*_data;
        FDOUBLE m = sqrt(sum2);
        if (fabs(m) > XMIPP_EQUAL_ACCURACY)
        {
            T im=(T) (1.0/m);
            for (auto& _data : __data) _data *= im;
        }
        else
            initZeros();
    }
    bool isrow() 	const	{return __isrow;}
    int  size()		const	{return __dim;}
    T& operator()(int i) 		{assert(i<__dim);return __data[i];}
    T  operator()(int i) const 	{assert(i<__dim);return __data[i];}
    std::vector<T>& data(){return __data;}
    // matrix_1_d = matrix_1_d*matrix_2_d
    Matrix1D<T> operator*(Matrix2D<T> const & rhs) {
        assert(__isrow==true);
        assert(__dim == rhs.dimy());
        Matrix1D<T> result(rhs.dimx());result.__isrow = true;
        for (int j = 0; j < rhs.dimx(); j++){
            result(j) = 0;
            for (int i = 0; i < rhs.dimy(); i++)
                result(j) += (*this)(i) * rhs(i,j);
        }
        return result;
    }
    void operator=(Matrix1D<T> const & rhs) {
        __dim = rhs.__dim;__isrow = rhs.__isrow;__data = rhs.__data;
    }
    T& operator[](int index) {
        assert(index<__dim);
        return __data[index];
    }
};
#define XX(MAT) MAT(0)
#define YY(MAT) MAT(1)
#define ZZ(MAT) MAT(2)
template<typename T>
T dotProduct(const Matrix1D< T >& v1, const Matrix1D< T >& v2){
    if (v1.size()!=v2.size()) {
        ERROR_REPORT("Dot product: vectors of different size or shape");
    }
    T accumulate = 0;
    for (int j = 0; j < v1.size(); j++)
        accumulate += v1(j) * v2(j);
    return accumulate;
}
template<typename T>
Matrix1D<T> vectorR3(T x, T y, T z){
    Matrix1D<T> result(3);
    result(0) = x;
    result(1) = y;
    result(2) = z;
    return result;
}
//
template<typename T>
class Matrix2D{
private:
    T*  __data;			// std::vector<std::vector<T>> exposed a bug in the the Microsoft vector implementation 
    int __dimy,__dimx,__dim;
	T& data(int y, int x)       { return __data[y*__dimx + x]; }
	T  data(int y, int x) const { return __data[y*__dimx + x]; }
	void initUndefined(int Ydim, int Xdim) {
		if (__dimy == Ydim && __dimx == Xdim) return;
		vDelete(__data);
		__dimy = Ydim; __dimx = Xdim; __dim = __dimy*__dimx;
		if (__dim == 0) return;
		__data = vNew(T,__dim);
	}
public:
    Matrix2D(                  ) : __data(nullptr), __dimy(0), __dimx(0), __dim(0) { }
    Matrix2D(int Ydim, int Xdim) : __data(nullptr), __dimy(0), __dimx(0), __dim(0) { initZeros(Ydim,Xdim); }
    Matrix2D(T v[][3]          ) : __data(nullptr), __dimy(0), __dimx(0), __dim(0) {
		initUndefined(3,3);
        data(0,0) = v[0][0]; data(0,1) = v[0][1]; data(0,2) = v[0][2];
        data(1,0) = v[1][0]; data(1,1) = v[1][1]; data(1,2) = v[1][2];
        data(2,0) = v[2][0]; data(2,1) = v[2][1]; data(2,2) = v[2][2];
    }
    Matrix2D(const Matrix2D<T>& rhs){
        __dimy = rhs.__dimy;__dimx = rhs.__dimx;__dim = rhs.__dim;
        __data = vNew(T,__dim);
        for (int i = 0; i < __dim; i++) __data[i] = rhs.__data[i];
    }
    ~Matrix2D() {
		vDelete(__data);
		__dimy = __dimx = __dim = 0;
	}

    void resize(int Ydim, int Xdim){
		if (__dimy == Ydim && __dimx == Xdim) return;
        if (__dimy == 0 || __dimx == 0) {
            initZeros(Ydim, Xdim);
            return;
        }
        auto oldY = __dimy; auto oldX = __dimx; auto old = vNew(T,oldY*oldX);
        for (int i = 0; i < oldY*oldX; i++) old[i] = __data[i];
		initUndefined(Ydim, Xdim);
        for (int y = 0; y < __dimy; y++) {
            for (int x = 0; x < __dimx; x++) {
				data(y,x) = (y < oldY && x < oldX) ? old[y*oldX + x] : 0;
            }
        }
        vDelete(old);
    }

    void initIdentity(){
        assert(__dimx==__dimy);
		for (int y = 0; y < __dimy; y++) {
			for (int x = 0; x < __dimx; x++) {
				data(y, x) = T(y == x);
			}
		}
    }

    void initZeros(int Ydim, int Xdim){
		initUndefined(Ydim, Xdim);
		for (int i = 0; i < __dim; i++) __data[i] = T(0);
    }

    bool isIdentity(){
        assert(__dimx==__dimy);
		for (int y = 0; y < __dimy; y++) {
			for (int x = 0; x < __dimx; x++) {
				auto err = data(y, x) - T(y == x);
				if (fabs(err) > XMIPP_EQUAL_ACCURACY) return false;
			}
		}
        return true;
    }

    bool equal(Matrix2D<T> const & rhs,FDOUBLE accuracy = XMIPP_EQUAL_ACCURACY){
        if (__dimx!=rhs.dimx() || __dimy!=rhs.dimy()) return false;
		for (int y = 0; y < __dimy; y++) {
			for (int x = 0; x < __dimx; x++) {
				auto err = data(y, x) - rhs.__data[y*__dimx + x];
				if (fabs(err) > accuracy) return false;
			}
		}
        return true;
    }

	void setSmallValuesToZero(FDOUBLE accuracy = XMIPP_EQUAL_ACCURACY)
    {
		for (int i = 0; i < __dim; i++) if (fabs(__data[i]) < accuracy) __data[i] = T(0);
    }

	int dimx()	const 	{return __dimx;}
    int dimy()	const	{return __dimy;}
    int dim()	const	{return __dim;}
	
    const T* wptr()	{return __data;}
    T& operator()(int y,int x) 			{assert(0 <= y && y<__dimy);assert(0<=x && x<__dimx);return data(y,x);}
    T  operator()(int y,int x) const 	{assert(0 <= y && y<__dimy);assert(0<=x && x<__dimx);return data(y,x);}

    void operator=(Matrix2D<T> const & rhs) {
        assert(__dimx == rhs.dimx());assert(__dimy == rhs.dimy());
		for (int i = 0; i < __dim; i++) __data[i] = rhs.__data[i];
    }

    // matrix_2_d = matrix_2_d*matrix_2_d
    Matrix2D<T> operator*(Matrix2D<T> const & rhs) const {
        assert(__dimx==__dimy);
        assert(__dimx == rhs.dimx());assert(__dimy == rhs.dimy());
        Matrix2D<T> result(__dimx,__dimy);
        for (int i = 0; i < __dimy; i++)
            for (int j = 0; j < __dimx; j++){
                T tmp = 0;
                for (int k = 0; k < __dimx; k++)
                    tmp += (*this)(i, k) * rhs(k, j);
				result.data(i, j) = tmp;
            }
        return result;
    }

    // matrix_1_d = matrix_2_d*matrix_1_d
    Matrix1D<T> operator*(Matrix1D<T> const & rhs) const {
        assert(rhs.isrow()==false);// column 1d array
        assert(__dimx == rhs.size());
        Matrix1D<T> result(__dimy);
        for (int i = 0; i < __dimy; i++){
			T tmp = 0;
            for (int j = 0; j < __dimx; j++)
                tmp += data(i, j) * rhs(j);
			result(i) = tmp;
		}
        return result;
    }

    Matrix2D<T> transpose() const{
        assert(__dimx==__dimy);
        Matrix2D<T> result(__dimx, __dimy);
        for (int i = 0; i < __dimy; i++)
            for (int j = 0; j < __dimx; j++)
                result(i,j) = data(j,i);
        return result;
    }

    Matrix2D<T> inv() const{
        assert(__dimx==__dimy);
        Matrix2D<T> result(__dimx, __dimy);
        // TODO
        ERROR_REPORT("Matrix2D inv is not implemented yet.");
        return result;
    }

    void getRow(int i, Matrix1D<T>& v) const
    {
        assert(i >= 0);
        assert(i < __dimy);
        for (int j = 0; j < __dimx; j++)
            v[j] = data(i, j);
        //v.setRow();
    }
    
    friend std::ostream& operator<<(std::ostream& os, Matrix2D<T> const & rhs) {
        os << std::endl;
		for (int y = 0; y < rhs.dimy(); y++) {
			for (int x = 0; x < rhs.dimx(); x++) {
				os << std::setprecision(7) << std::setw(15) << rhs(y,x) << " ";
			}
            os << std::endl;
        }
        return os;
    }
};

/***************************************************************************
 *
 * Author: "Sjors H.W. Scheres"
 * MRC Laboratory of Molecular Biology
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
#if defined(_WIN32) || defined(_WIN64)
#define strtok_r strtok_s
#endif
// some utils
typedef std::string FileName;
#define textToDOUBLE(v) atof(v)
#define textToFloat(v) atof(v)
#define textToInteger(v) atoi(v)
//
inline char* firstToken(const std::string& str,char* &signal){
    return strtok_r((char*) str.c_str(), " \t\n", &signal);
}
inline char* nextToken(char* &signal){
    return strtok_r((char*) NULL, " \t\n", &signal);
}
/* from transformations.cpp ---------------------------------------------------- */
/* Rotation 3D around the system axes -------------------------------------- */
void rotation3DMatrix(FDOUBLE ang, char axis, Matrix2D< FDOUBLE > &result,
                      bool homogeneous = true);
/* Align a vector with Z axis */
void alignWithZ(const Matrix1D< FDOUBLE > &axis, Matrix2D< FDOUBLE >& result,
                bool homogeneous = true);
/* Rotation 3D around any axis -------------------------------------------- */
void rotation3DMatrix(FDOUBLE ang, const Matrix1D< FDOUBLE > &axis,
                      Matrix2D< FDOUBLE > &result, bool homogeneous = true);


#endif /* defined(MATRIX_H_) */
