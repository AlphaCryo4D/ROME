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

#ifndef FFT_BASE_H_
#define FFT_BASE_H_

#include "fftw/fftw3.h"

#include "./error.h"
#include "./macros.h"

// this fftw3 transformer do basic job of forward and backward data
// without additional memory alloc ,it can save memory but cannot get SOA of complex data
// notice both the input and output data will be changed when execute plan
// NOTE : fftwf_(float),fftw_(double),fftwl_(long double)
#define CLASS_NAME FFTWBase
#define DATA_TYPE double
#define FFT_TYPE(v) fftw_##v
class CLASS_NAME
{
protected:
    // The dimension of 2D/3D data
    size_t zsize,ysize,xsize;
    size_t size3,Fsize3;
    
    // 2D data in Real space
    DATA_TYPE* real;
    
    // 2D data in Fourier space
    FFT_TYPE(complex) *complex;
    
    // FFTW forward and backward plan
    FFT_TYPE(plan) forward_plan;
    FFT_TYPE(plan) backward_plan;
    
public:
    //
    int nr_threads;
    // Constructor
    CLASS_NAME() : real(nullptr),complex(nullptr),zsize(0),ysize(0),xsize(0),size3(0),Fsize3(0) {}
    void init(DATA_TYPE* _real,FFT_TYPE(complex)* _complex,size_t _xsize,size_t _ysize,size_t _zsize = 1,int _nr_threads = 1);
    // Destructor
    ~CLASS_NAME(){fini();}
    void fini();
    
    // FFT forward 2D/3D data
    void FourierTransform(bool normalize = true);
    
    // FFT backward 2D/3D data
    void inverseFourierTransform();
    
    // plan manager
    void reset_plan(DATA_TYPE* _real,FFT_TYPE(complex)* _complex,size_t _xsize,size_t _ysize,size_t _zsize = 1);
    
protected:
    void create_plan();
    
    void delete_plan();
    
    void set(DATA_TYPE* _real,FFT_TYPE(complex)* _complex,size_t _xsize,size_t _ysize,size_t _zsize);
    
    void normalizeComplex();
    
    // execute the plan
    void forward_execute();
    void backward_execute();
};
#undef CLASS_NAME
#undef DATA_TYPE
#undef FFT_TYPE

//    ----------------------------------------------------------------------------   //
//    ----------   duplicate above code and change MACRO #define   ---------------   //
//    ----------------------------------------------------------------------------   //

#define CLASS_NAME FFTWFBase
#define DATA_TYPE float
#define FFT_TYPE(v) fftwf_##v
class CLASS_NAME
{
protected:
    // The dimension of 2D/3D data
    size_t zsize,ysize,xsize;
    size_t size3,Fsize3;
    
    // 2D data in Real space
    DATA_TYPE* real;
    
    // 2D data in Fourier space
    FFT_TYPE(complex) *complex;
    
    // FFTW forward and backward plan
    FFT_TYPE(plan) forward_plan;
    FFT_TYPE(plan) backward_plan;
    
public:
    //
    int nr_threads;
    // Constructor
    CLASS_NAME():real(nullptr),complex(nullptr),zsize(0),ysize(0),xsize(0),size3(0),Fsize3(0){}
    void init(DATA_TYPE* _real,FFT_TYPE(complex)* _complex,size_t _xsize,size_t _ysize,size_t _zsize = 1,int _nr_threads = 1);
    // Destructor
    ~CLASS_NAME(){fini();}
    void fini();
    
    // FFT forward 2D/3D data
    void FourierTransform(bool normalize = true);
    
    // FFT backward 2D/3D data
    void inverseFourierTransform();
    
    // plan manager
    void reset_plan(DATA_TYPE* _real,FFT_TYPE(complex)* _complex,size_t _xsize,size_t _ysize,size_t _zsize = 1);
    
protected:
    void create_plan();
    
    void delete_plan();
    
    void set(DATA_TYPE* _real,FFT_TYPE(complex)* _complex,size_t _xsize,size_t _ysize,size_t _zsize);
    
    void normalizeComplex();
    
    // execute the plan
    void forward_execute();
    void backward_execute();
};
#undef CLASS_NAME
#undef DATA_TYPE
#undef FFT_TYPE

#if defined(FLOAT_PRECISION)
typedef FFTWFBase FourierTransformerBase;
#else
typedef FFTWBase FourierTransformerBase;
#endif

#endif /* defined(FFT_BASE_H_) */
