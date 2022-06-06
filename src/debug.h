/***************************************************************************
 *
 * Authors: "Yongbei(Glow) Ma,Jiayi (Timmy) Wu, Youdong (Jack) Mao"
 * Dana-Farber Cancer Institute, Harvard Medical School and Peking University
 *			"Bevin R. Brett"
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

#ifndef _ROME_DEBUG
#define _ROME_DEBUG

#include "./util.h"


 //
// #define DATA_STREAM

class cacheInvestigator{
    //
    
};

#define NOT_IN_PARALLEL_REGION (omp_get_num_threads()==1)
#define IN_PARALLEL_REGION (omp_get_num_threads()!=1)
#define MASTER_THREAD (omp_get_thread_num()==0)
// write the intermediate data to disk and not print this out
// local_dataflow global_dataflow
class DataStream{
    template<typename T>
    int same(T* v1,T* v2,int len,double precision){
        int differentCounts = 0;
        for (int i = 0; i < len; i++){
            if(std::fabs(double(v1[i]-v2[i])) > precision) {
                if (differentCounts<10) {
                    std::cout<<"(i:"<<i<<",v1:"<<v1[i]<<",v2:"<<v2[i]<<")"<<" ";
                }
                differentCounts++;
            }
        }
        return differentCounts;
    }
private:
    // double,float,int pointer
    enum DataType{Double,Float,Int};
    class Ptr{
    protected:
        static const bool reduceDataLength = false;//true;
    public:
        void* data;int len;char* what;char* file;int line;
        DataType dataType;
        size_t dataTypeSize;
        Ptr(){}
        Ptr(void* _data,int _len,char* _what,char* _file,int _line){
            data = _data;len = _len;what = _what;file = _file;line = _line;
        }
        ~Ptr(){}
        DataType type(){return dataType;}
        size_t sizeoftype(){return dataTypeSize;}
        template<typename T>
        void stat(T* v,int len,T& min,T& max,T& average){
            min = std::numeric_limits<T>::max();
            max = std::numeric_limits<T>::min();
            average = 0;
            for (int i = 0; i < len; i++) {
                if (v[i] < min && v[i] != 0) min = v[i];
                if (v[i] > max) max = v[i];
                average += (T)((double)v[i]/(double)len);
            }
        }
    };
    class DoublePtr : public Ptr {
    public:
        double data_reduce[3];//minimum,average,maximum
		static DoublePtr* make(double* _data, int _len, char* _what, char* _file, int _line) {
#include "./util_heap_undefs.h"
			return sNewA(DoublePtr,(_data, _len, _what, _file, _line));
#include "./util_heap_defs.h"
		}
		DoublePtr(double* _data,int _len,char* _what,char* _file,int _line) 
		  : Ptr(_data,_len,_what,_file,_line) {
            if (reduceDataLength) {
                stat(_data,len,data_reduce[0],data_reduce[1],data_reduce[2]);
                data = data_reduce;len = 3;
            }
            dataType = Double;
            dataTypeSize = sizeof(double);
        }
        ~DoublePtr(){}
    };
    // for AOS,complex like [real,imag,real,imag,real,imag....]
    // convert to SOA like [real,real,real,....imag,imag,imag....]
    // NOTE : this will be trouble for check AOS-data
    class DoubleComplexPtr : public Ptr{
    public:
        double* complex_SOA_data;
        double data_reduce[6];
		static DoubleComplexPtr* make(const std::complex<double>* _data,int _len,char* _what,char* _file,int _line) {
#include "./util_heap_undefs.h"
			return sNewA(DoubleComplexPtr,(_data,_len,_what,_file,_line));
#include "./util_heap_defs.h"
		}
		DoubleComplexPtr(const std::complex<double>* _data,int _len,char* _what,char* _file,int _line)
        {
            len = 2*_len;what = _what;file = _file;line = _line;
            complex_SOA_data = (double*)aMalloc(sizeof(double)*len,64);
            // copy the real part
            for (int i = 0; i < _len; i++) complex_SOA_data[i] = _data[i].real();
            // copy the imag part
            for (int i = 0; i < _len; i++) complex_SOA_data[_len+i] = _data[i].imag();
            data = complex_SOA_data;
            if (reduceDataLength) {
                stat(complex_SOA_data,_len,data_reduce[0],data_reduce[1],data_reduce[2]);
                stat(complex_SOA_data+_len,_len,data_reduce[3],data_reduce[4],data_reduce[5]);
                data = data_reduce;len = 6;
            }
            dataType = Double;
            dataTypeSize = sizeof(double);
        }
        ~DoubleComplexPtr(){aFree(complex_SOA_data);}
    };
    
    class FloatPtr : public Ptr {
    public:
        float data_reduce[3];//minimum,average,maximum
		static FloatPtr* make(float* _data,int _len,char* _what,char* _file,int _line) {
#include "./util_heap_undefs.h"
			return sNewA(FloatPtr,(_data,_len,_what,_file,_line));
#include "./util_heap_defs.h"
		}
		FloatPtr(float* _data,int _len,char* _what,char* _file,int _line) :
        Ptr(_data,_len,_what,_file,_line) {
            if (reduceDataLength) {
                stat(_data,len,data_reduce[0],data_reduce[1],data_reduce[2]);
                data = data_reduce;len = 3;
            }
        	dataType = Float;
            dataTypeSize = sizeof(float);
        }
        ~FloatPtr(){}
    };
    class IntPtr : public Ptr {
    public:
        int data_reduce[3];//minimum,average,maximum
		static IntPtr* make(int* _data,int _len,char* _what,char* _file,int _line) {
#include "./util_heap_undefs.h"
			return sNewA(IntPtr,(_data,_len,_what,_file,_line));
#include "./util_heap_defs.h"
		}
		IntPtr(int* _data,int _len,char* _what,char* _file,int _line) :
        Ptr(_data,_len,_what,_file,_line) {
            if (reduceDataLength) {
                stat(_data,len,data_reduce[0],data_reduce[1],data_reduce[2]);
                data = data_reduce;len = 3;
            }
        	dataType = Int;
            dataTypeSize = sizeof(int);
        }
        ~IntPtr(){}
    };
    class IntVar : public Ptr{
    public:
        int var;
		static IntVar* make(int _data,int _len,char* _what,char* _file,int _line) {
#include "./util_heap_undefs.h"
			return sNewA(IntVar,(_data,_len,_what,_file,_line));
#include "./util_heap_defs.h"
		}
		IntVar(int _data,int _len,char* _what,char* _file,int _line)
        {
            var = _data;
            data = &var;
            len = _len;what = _what;file = _file;line = _line;
            dataType = Int;
            dataTypeSize = sizeof(int);
        }
        ~IntVar(){}
    };
    class DoubleVar : public Ptr{
    public:
        double var;
		static DoubleVar* make(double _data,int _len,char* _what,char* _file,int _line) {
#include "./util_heap_undefs.h"
			return sNewA(DoubleVar,(_data,_len,_what,_file,_line));
#include "./util_heap_defs.h"
		}
		DoubleVar(double _data,int _len,char* _what,char* _file,int _line)
        {
            var = _data;
            data = &var;
            len = _len;what = _what;file = _file;line = _line;
            dataType = Double;
            dataTypeSize = sizeof(double);
        }
        ~DoubleVar(){}
    };
    std::string instream_fn,outstream_fn;
    FILE *instream_file,*outstream_file;
    std::vector<Ptr*> cachedPtr;
    size_t stream_size;
public:
    bool doDataStream = false;
public:
    DataStream():instream_fn("NULL"),outstream_fn("NULL"){}
    ~DataStream(){}
    void init(std::string out_fn,std::string in_fn){
        stream_size = 0;
        instream_fn = in_fn;
        outstream_fn = out_fn;
        // open file
        if(instream_fn!="NULL") {instream_file = fopen(instream_fn.c_str(),"rb");doDataStream = true;}
        if(outstream_fn!="NULL") {outstream_file = fopen(outstream_fn.c_str(),"wb");doDataStream = true;}
        std::cout.precision(20);
    }
    void fini(){
        if(instream_fn!="NULL") fclose(instream_file);
        if(outstream_fn!="NULL") fclose(outstream_file);
    }
    void turnOn(){if(instream_fn!="NULL" || outstream_fn!="NULL") doDataStream = true;}
    void turnOff(){if(instream_fn!="NULL" || outstream_fn!="NULL") doDataStream = false;}
    inline void foutDouble(double* _data,int _len,char* _what,char* _file,int _line){
        if(doDataStream&&NOT_IN_PARALLEL_REGION) cachedPtr.push_back(DoublePtr::make(_data,_len,_what,_file,_line));
    }
    inline void foutDoubleComplex(const std::complex<double>* _data,int _len,char* _what,char* _file,int _line){
        if(doDataStream&&NOT_IN_PARALLEL_REGION) cachedPtr.push_back(DoubleComplexPtr::make(_data,_len,_what,_file,_line));
    }
    inline void foutDouble(double _data,char* _what,char* _file,int _line){
        if(doDataStream&&NOT_IN_PARALLEL_REGION) cachedPtr.push_back(DoubleVar::make(_data,1,_what,_file,_line));
    }
    inline void foutFloat(float* _data,int _len,char* _what,char* _file,int _line){
        if(doDataStream&&NOT_IN_PARALLEL_REGION) cachedPtr.push_back(FloatPtr::make(_data,_len,_what,_file,_line));
    }
    inline void foutInt(int* _data,int _len,char* _what,char* _file,int _line){
        if(doDataStream&&NOT_IN_PARALLEL_REGION) cachedPtr.push_back(IntPtr::make(_data,_len,_what,_file,_line));
    }
    inline void foutInt(int _data,char* _what,char* _file,int _line){
        if(doDataStream&&NOT_IN_PARALLEL_REGION) cachedPtr.push_back(IntVar::make(_data,1,_what,_file,_line));
    }
    inline void foutInt(long _data,char* _what,char* _file,int _line){
        if(doDataStream&&NOT_IN_PARALLEL_REGION) cachedPtr.push_back(IntVar::make(int(_data),1,_what,_file,_line));
    }
    inline void foutInt(size_t _data,char* _what,char* _file,int _line){
        if(doDataStream&&NOT_IN_PARALLEL_REGION) cachedPtr.push_back(IntVar::make(int(_data),1,_what,_file,_line));
    }
    void flush(){
        // write the file
        if(doDataStream&&outstream_fn!="NULL")
        {
            if (IN_PARALLEL_REGION && MASTER_THREAD)
                std::cerr<<"Warning : write data in parallel....it will only write out thread 0 data...."<<std::endl;
            if (NOT_IN_PARALLEL_REGION)
            {
                for (auto& ptr : cachedPtr){
                    fwrite((char*)ptr->data,ptr->sizeoftype()*ptr->len,1,outstream_file);
                    stream_size += ptr->sizeoftype()*ptr->len;
                    sDelete(ptr);
                }
                if (stream_size%(1024)==0) {
                    double size = stream_size/(1024.*1024.*1024.);
                    // std::cout<<"write to "<<outstream_fn<<",file size : "<<size<<" GB."<<std::endl;
                }
            }
        }
        cachedPtr.resize(0);
    }
    void check(){
        if(doDataStream&&instream_fn!="NULL")
        {
            if (IN_PARALLEL_REGION && MASTER_THREAD)
                std::cerr<<"Warning : check data in parallel....it will only check thread 0 data...."<<std::endl;
            if (NOT_IN_PARALLEL_REGION)
            {
                double precision = 1e-9;
                assert(cachedPtr.size()>0);
                for (auto& ptr : cachedPtr){
                    void* data = (void*)aMalloc(ptr->sizeoftype()*ptr->len,64);
                    int diffcounts;
                    fread((char*)data,ptr->sizeoftype()*ptr->len,1,instream_file);
                    switch (ptr->type()) {
                        case Double:
                            diffcounts = same((double*)data,(double*)ptr->data,ptr->len,precision);
                            break;
                        case Float:
                            diffcounts = same((float*)data,(float*)ptr->data,ptr->len,precision);
                            break;
                        case Int:
                            diffcounts = same((int*)data,(int*)ptr->data,ptr->len,precision);
                            break;
                        default:
                            break;
                    }
                    aFree(data);
                    if (diffcounts!=0) {
                        std::cerr<<"!!!diff array..."<<ptr->what<<" in : "<<ptr->file<<" line : "<<ptr->line
                        		<<" diff percentage : "<< (double(diffcounts)/double(ptr->len))*100<<"%"<<std::endl;
                    }
                }
            }
        }
    }
};

#if 0
//-----------------------------------------------
// to define debug which function......
//-----------------------------------------------
// debug getFourierTransformsAndCtfs()
//#define YBDEBUG_GFTAC

// debug precalculateShiftedImagesCtfsAndInvSigma2s()
//#define YBDEBUG_PRESHIFT

//debug getAllSquaredDifferences()
//#define YBDEBUG_GASD

// debug convertAllSquaredDifferencesToWeights()
// #define YBDEBUG_CASDTW

// debug updateOtherParams()
// #define YBDEBUG_PARAM // a little diff

// debug backProjection()
// #define YBDEBUG_BACK

// debug storeWeightedSums()
// #define YBDEBUG_SWS

//debug maximization()
//#define YBDEBUG_MAXI

//-----------------------------------------------
#define debug_iter 0
#define YBDEBUG_STOP(_iter,_count) {static int count = 0;if(iter==_iter) count++;if(iter==_iter&&count>_count) exit(1);}
//#define YBDEBUG_PASS(what,flag) std::cerr<<what<<" "<<flag<<" "<<__FILE__<<" "<<__LINE__<<std::endl;
#define YBDEBUG_PASS(what,flag)
//-----------------------------------------------
// debug for input parameter

#if (defined(DEBUG_ROME) || defined(DEBUG_RELION))

// define something
#define YBDEBUG_SCALE_VAR(NAME,VAR) if(iter >= debug_iter) debug_output_file<<"~~~ "<<#NAME<<" , "<<VAR<<std::endl;

#define YBDEBUG_SCALE_VARS_START(WHAT) if(iter >= debug_iter) { debug_output_file<<"~~~ "<<WHAT<<" , ";
#define PV(VAR) debug_output_file<<" "<<VAR<<" ";
#define YBDEBUG_SCALE_VARS_END debug_output_file<<std::endl; }

#define YBDEBUG_VEC_VAR(NAME,VEC,LEN) if(iter >= debug_iter){ \
debug_output_file.precision(12); \
if((LEN) < 3 ) \
debug_output_file<<"~~~ "<<#NAME<<" , "<<(VEC)[0]<<" "<<sumVec(VEC,LEN)<<std::endl; \
else \
debug_output_file<<"~~~ "<<#NAME<<" , "<<(VEC)[0]<<" "<<(VEC)[1]<<" "<<(VEC)[2]<<" "<<sumVec(VEC,LEN)<<" "<<(VEC)[LEN-3]<<" "<<(VEC)[LEN-2]<<" "<<(VEC)[LEN-1]<<std::endl; \
}

#define YBDEBUG_VEC_VAR_REAL(NAME,VEC,LEN) if(iter >= debug_iter){ \
debug_output_file.precision(12);\
debug_output_file<<"~~~ "<<#NAME<<" , "<<(VEC)[0].real<<" "<<(VEC)[1].real<<" "<<(VEC)[2].real<<" "<<(sumVec(VEC,LEN)).real<<" "<<(VEC)[LEN-3].real<<" "<<(VEC)[LEN-2].real<<" "<<(VEC)[LEN-1].real<<std::endl; \
}

#define YBDEBUG_VEC_VAR_IMAG(NAME,VEC,LEN) if(iter >= debug_iter) {\
debug_output_file.precision(12);\
debug_output_file<<"~~~ "<<#NAME<<" , "<<(VEC)[0].imag<<" "<<(VEC)[1].imag<<" "<<(VEC)[2].imag<<" "<<(sumVec(VEC,LEN)).imag<<" "<<(VEC)[LEN-3].imag<<" "<<(VEC)[LEN-2].imag<<" "<<(VEC)[LEN-1].imag<<std::endl; \
}

#else

// define nothing
#define YBDEBUG_SCALE_VAR(NAME,VAR)
#define YBDEBUG_SCALE_VARS_START(WHAT)
#define PV(VAR)
#define YBDEBUG_SCALE_VARS_END
#define YBDEBUG_VEC_VAR(NAME,VEC,LEN)
#define YBDEBUG_VEC_VAR_REAL(NAME,VEC,LEN)
#define YBDEBUG_VEC_VAR_IMAG(NAME,VEC,LEN)

#endif

#endif // 0

template<typename T>
T sumVec(T* data,int length){
    T sum = 0;
    for (int i = 0; i < length; i++)
        sum += data[i];
    return sum;
}

template<typename T>
T sumVec(const T* data,int length){
    T sum = 0;
    for (int i = 0; i < length; i++)
        sum += data[i];
    return sum;
}
//template<class T1,typename T2>
//void printList(T1& outstream,T2 t){outstream<<" "<<t<<" ";}
//
//template<class T1,typename T2,typename... Args>
//void printList(T1& outstream,T2 t,Args... args)
//{
//    outstream<<" "<<t<<" ";
//    printList(outstream,args...);
//}

//template<class T>
//void printList(std::initializer_list<T> list){
//    for (auto elem : list) {
//        std::cout<<" "<<elem<<" ";
//    }
//}

#endif /* defined(_DEBUG_H_) */
