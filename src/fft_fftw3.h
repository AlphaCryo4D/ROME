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

#ifndef FFT_FFTW3_H
#define FFT_FFTW3_H

#include "./fft_base.h"


template<typename T>
struct SOAComplexBase{
    T* real;
    T* imag;
};

template<>
class SOAComplexBase<double> {
public:
    double* real;
    double* imag;
	SOAComplexBase() : real(nullptr), imag(nullptr) {}
    SOAComplexBase(double* real, double* imag) { this->real = real; this->imag = imag; }
};
typedef SOAComplexBase<double> SOAComplexDouble;
typedef SOAComplexBase<float> SOAComplexFloat;

template<>
class SOAComplexBase<double const> {
public:
    const double * real;
    const double * imag;
	SOAComplexBase() : real(nullptr), imag(nullptr) {}
    SOAComplexBase(double const * real, double const * imag) { this->real = real; this->imag = imag; }
	SOAComplexBase(SOAComplexDouble& src) : real(src.real), imag(src.imag) {}
};
typedef SOAComplexBase<double const> SOAComplexReadonly;

template<typename T>
struct SOAComplexArray : public SOAComplexBase<T> {
    size_t length;
    SOAComplexArray(size_t _length) : length(_length) {
        this->real = (T*)aMalloc(sizeof(T)*length,64);
        this->imag = (T*)aMalloc(sizeof(T)*length,64);
    }
    SOAComplexArray(const SOAComplexArray& from) : length(from.length) {
        this->real = (T*)aMalloc(sizeof(T)*length,64);
        this->imag = (T*)aMalloc(sizeof(T)*length,64);
        for(size_t i = 0;i < length;i++){
            this->real[i] = from.real[i];
            this->imag[i] = from.imag[i];
        }
    }
    void operator=(const SOAComplexArray& from) {
        assert(length == from.length);
        for(size_t i = 0;i < length;i++){
            this->real[i] = from.real[i];
            this->imag[i] = from.imag[i];
        }
    }
    ~SOAComplexArray() {
        if(this->real) aFree(this->real); this->real = NULL;
        if(this->imag) aFree(this->imag); this->imag = NULL;
    }
};

#if defined(FLOAT_PRECISION)
typedef MKL_Complex8 MKL_Complex;
#else
typedef MKL_Complex16 MKL_Complex;
#endif

// NOTE : fftwf_(float),fftw_(double),fftwl_(long double)
#define CLASS_NAME(v) FFTW##v
#define DATA_TYPE double
#define FFT_TYPE(v) fftw_##v
class CLASS_NAME(Transformer) : public CLASS_NAME(Base)
{
    
private:
    
    // 2D data in Real space
    DATA_TYPE* realData;
    
    // 2D data in Fourier space
    FFT_TYPE(complex)* complexData;
    
public:
    // Constructor
    CLASS_NAME(Transformer)             (size_t _xsize, size_t _ysize, size_t _zsize = 1);
	static CLASS_NAME(Transformer)* make(size_t _xsize, size_t _ysize, size_t _zsize = 1) {
#include "./util_heap_undefs.h"
		typedef CLASS_NAME(Transformer) T;
		return sNewA(T,(_xsize, _ysize, _zsize));
#include "./util_heap_defs.h"
	}

    // Destructor
    ~CLASS_NAME(Transformer)();
    
    // FFT forward 2D/3D data
    void FourierTransform(const DATA_TYPE* data, SOAComplexBase<DATA_TYPE>& Fdata,bool normalize = true);
    
    void FourierTransform(const DATA_TYPE* data, DATA_TYPE* Fdata_real,DATA_TYPE* Fdata_imag,bool normalize = true);
    
    // FFT backward 2D/3D data
    void inverseFourierTransform(const SOAComplexBase<DATA_TYPE>& Fdata,DATA_TYPE* data);
    
    void inverseFourierTransform(const DATA_TYPE* Fdata_real,const DATA_TYPE* Fdata_imag,DATA_TYPE* data);
    
    void testCorrection();
};
#undef CLASS_NAME
#undef DATA_TYPE
#undef FFT_TYPE

//    ----------------------------------------------------------------------------   //
//    ----------   duplicate above code and change MACRO #define   ---------------   //
//    ----------------------------------------------------------------------------   //

#define CLASS_NAME(v) FFTWF##v
#define DATA_TYPE float
#define FFT_TYPE(v) fftwf_##v
class CLASS_NAME(Transformer) : public CLASS_NAME(Base)
{
    
private:
    
    // 2D data in Real space
    DATA_TYPE* realData;
    
    // 2D data in Fourier space
    FFT_TYPE(complex)* complexData;
    
public:
    // Constructor
    CLASS_NAME(Transformer)             (size_t _xsize, size_t _ysize, size_t _zsize = 1);
	static CLASS_NAME(Transformer)* make(size_t _xsize, size_t _ysize, size_t _zsize = 1) {
#include "./util_heap_undefs.h"
		typedef CLASS_NAME(Transformer) T;
		return sNewA(T,(_xsize, _ysize, _zsize));
#include "./util_heap_defs.h"
	}

    // Destructor
    ~CLASS_NAME(Transformer)();
    
    // FFT forward 2D/3D data
    void FourierTransform(const DATA_TYPE* data, SOAComplexBase<DATA_TYPE>& Fdata,bool normalize = true);
    
    void FourierTransform(const DATA_TYPE* data, DATA_TYPE* Fdata_real,DATA_TYPE* Fdata_imag,bool normalize = true);
    
    // FFT backward 2D/3D data
    void inverseFourierTransform(const SOAComplexBase<DATA_TYPE>& Fdata,DATA_TYPE* data);
    
    void inverseFourierTransform(const DATA_TYPE* Fdata_real,const DATA_TYPE* Fdata_imag,DATA_TYPE* data);
    
};
#undef CLASS_NAME
#undef DATA_TYPE
#undef FFT_TYPE
//
#if defined(FLOAT_PRECISION)
typedef FFTWFTransformer FourierTransformer;
typedef fftwf_complex FourierComplex;
#else
typedef FFTWTransformer FourierTransformer;
typedef fftw_complex FourierComplex;
#endif

#endif
