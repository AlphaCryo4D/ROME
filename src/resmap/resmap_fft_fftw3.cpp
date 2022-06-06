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
 
#include "resmap_fft_fftw3.h"

#define CLASS_NAME(v) FFTW##v
#define DATA_TYPE double
#define FFT_TYPE(v) fftw_##v
// Constructor
CLASS_NAME(Transformer)::CLASS_NAME(Transformer)(size_t _xsize,size_t _ysize,size_t _zsize /*= 1*/)
{
    xsize = _xsize;ysize = _ysize;zsize = _zsize;
    size_t size2 = ysize*xsize;size3 = zsize*ysize*xsize;
    size_t Fsize2 = ysize*(xsize/2+1);Fsize3 = zsize*Fsize2;
    realData = (DATA_TYPE*)FFT_TYPE(malloc)(sizeof(DATA_TYPE)*size3);
    complexData = (FFT_TYPE(complex)*) FFT_TYPE(malloc)(sizeof(FFT_TYPE(complex)) * Fsize3);
    init(realData, complexData, _xsize, _ysize, _zsize);
}

// Destructor
CLASS_NAME(Transformer)::~CLASS_NAME(Transformer)(){
    FFT_TYPE(free)(realData);
    FFT_TYPE(free)(complexData);
}

// FFT forward 2D/3D data
void CLASS_NAME(Transformer)::FourierTransform(const DATA_TYPE* data, SOAComplexBase<DATA_TYPE>& Fdata,bool normalize /*= true*/)
{
    DATA_TYPE* Fdata_real = Fdata.real;
    DATA_TYPE* Fdata_imag = Fdata.imag;
    FourierTransform(data, Fdata_real,Fdata_imag,normalize);
}

void CLASS_NAME(Transformer)::FourierTransform(const DATA_TYPE* data, DATA_TYPE* Fdata_real,DATA_TYPE* Fdata_imag,bool normalize /*= true*/)
{
    // copy data
    for (size_t i = 0; i < size3; i++)
        realData[i] = data[i];
    
    CLASS_NAME(Base)::FourierTransform(normalize);
    
    // copy data
    for (size_t i = 0; i < Fsize3; i++){
        Fdata_real[i] = complexData[i][0];
        Fdata_imag[i] = complexData[i][1];
    }
}

// FFT backward 2D/3D data
void CLASS_NAME(Transformer)::inverseFourierTransform(const SOAComplexBase<DATA_TYPE>& Fdata,DATA_TYPE* data)
{
    DATA_TYPE* Fdata_real = Fdata.real;
    DATA_TYPE* Fdata_imag = Fdata.imag;
    inverseFourierTransform(Fdata_real,Fdata_imag,data);
}

void CLASS_NAME(Transformer)::inverseFourierTransform(const DATA_TYPE* Fdata_real,const DATA_TYPE* Fdata_imag,DATA_TYPE* data)
{
    // copy data
    for (size_t i = 0; i < Fsize3; i++){
        complexData[i][0] = Fdata_real[i];
        complexData[i][1] = Fdata_imag[i];
    }
    
    CLASS_NAME(Base)::inverseFourierTransform();
    
    // copy data
    for (size_t i = 0; i < size3; i++)
        data[i] = realData[i];
}
#undef CLASS_NAME
#undef DATA_TYPE
#undef FFT_TYPE

//    ----------------------------------------------------------------------------   //
//    ----------   duplicate above code and change MACRO #define   ---------------   //
//    ----------------------------------------------------------------------------   //

#define CLASS_NAME(v) FFTWF##v
#define DATA_TYPE float
#define FFT_TYPE(v) fftwf_##v
// Constructor
CLASS_NAME(Transformer)::CLASS_NAME(Transformer)(size_t _xsize,size_t _ysize,size_t _zsize /*= 1*/)
{
    xsize = _xsize;ysize = _ysize;zsize = _zsize;
    size_t size2 = ysize*xsize;size3 = zsize*ysize*xsize;
    size_t Fsize2 = ysize*(xsize/2+1);Fsize3 = zsize*Fsize2;
    realData = (DATA_TYPE*)FFT_TYPE(malloc)(sizeof(DATA_TYPE)*size3);
    complexData = (FFT_TYPE(complex)*) FFT_TYPE(malloc)(sizeof(FFT_TYPE(complex)) * Fsize3);
    init(realData, complexData, _xsize, _ysize, _zsize);
}

// Destructor
CLASS_NAME(Transformer)::~CLASS_NAME(Transformer)(){
    FFT_TYPE(free)(realData);
    FFT_TYPE(free)(complexData);
}

// FFT forward 2D/3D data
void CLASS_NAME(Transformer)::FourierTransform(const DATA_TYPE* data, SOAComplexBase<DATA_TYPE>& Fdata,bool normalize /*= true*/)
{
    DATA_TYPE* Fdata_real = Fdata.real;
    DATA_TYPE* Fdata_imag = Fdata.imag;
    FourierTransform(data, Fdata_real,Fdata_imag,normalize);
}

void CLASS_NAME(Transformer)::FourierTransform(const DATA_TYPE* data, DATA_TYPE* Fdata_real,DATA_TYPE* Fdata_imag,bool normalize /*= true*/)
{
    // copy data
    for (size_t i = 0; i < size3; i++)
        realData[i] = data[i];
    
    CLASS_NAME(Base)::FourierTransform(normalize);
    
    // copy data
    for (size_t i = 0; i < Fsize3; i++){
        Fdata_real[i] = complexData[i][0];
        Fdata_imag[i] = complexData[i][1];
    }
}

// FFT backward 2D/3D data
void CLASS_NAME(Transformer)::inverseFourierTransform(const SOAComplexBase<DATA_TYPE>& Fdata,DATA_TYPE* data)
{
    DATA_TYPE* Fdata_real = Fdata.real;
    DATA_TYPE* Fdata_imag = Fdata.imag;
    inverseFourierTransform(Fdata_real,Fdata_imag,data);
}

void CLASS_NAME(Transformer)::inverseFourierTransform(const DATA_TYPE* Fdata_real,const DATA_TYPE* Fdata_imag,DATA_TYPE* data)
{
    // copy data
    for (size_t i = 0; i < Fsize3; i++){
        complexData[i][0] = Fdata_real[i];
        complexData[i][1] = Fdata_imag[i];
    }
    
    CLASS_NAME(Base)::inverseFourierTransform();
    
    // copy data
    for (size_t i = 0; i < size3; i++)
        data[i] = realData[i];
}
#undef CLASS_NAME
#undef DATA_TYPE
#undef FFT_TYPE

// ---------------------------------   //

#include "resmap_time.h"
#include "mkl.h"
#include "resmap_array.h"
void FFTWTransformer::testCorrection()
{
    //    std::cout.precision(30);
    // FFTWBase newtransformer;

    double t1 = dtime();
    for (int size = 100; size < 1100; size+=100)
    {
        double init_time,set_time,forward_time,backward_time;
        double start = dtime();

        // old
        FFTWTransformer transformer(size,size,size);
        Vol<double> Mconv,Fconv_real,Fconv_imag;
        Mconv.init_zero(size,size,size);Fconv_real.init_zero(size,size,size/2+1);Fconv_imag.init_zero(size,size,size/2+1);

        // new
        double* realData = (double*)aMalloc(sizeof(double)*size*size*size,64);
        MKL_Complex16* complexData = (MKL_Complex16*) aMalloc(sizeof(MKL_Complex16) * size*size*(size/2+1),64);
        
        // if(size == 100) newtransformer.init(realData, (fftw_complex*)complexData, size, size, size, 72);
        // else newtransformer.reset_plan(realData, (fftw_complex*)complexData, size, size, size);
        
        FFTWBase newtransformer;
        newtransformer.init(realData, (fftw_complex*)complexData, size, size, size, 1);
        
        double end = dtime();
        init_time = end - start;
        
        
        start = dtime();
        // set
        for (int i = 0; i < size*size*size; i++) {
            (Mconv.wptr())[i] = (double)rand() / RAND_MAX;
            realData[i] = (Mconv.wptr())[i];
        }
        end = dtime();
        set_time = end - start;

        start = dtime();
        // new
        newtransformer.FourierTransform();
        // old
        transformer.FourierTransform(Mconv.wptr(), Fconv_real.wptr(), Fconv_imag.wptr());
        for (int i = 0; i < size*size*(size/2+1); i++) {
            if(fabs(complexData[i].real - (Fconv_real.wptr())[i]) > 1e-6 ||
               fabs(complexData[i].imag - (Fconv_imag.wptr())[i]) > 1e-6){
                std::cout<<complexData[i].real<<" "<<(Fconv_real.wptr())[i]<<std::endl;
                std::cout<<complexData[i].imag<<" "<<(Fconv_imag.wptr())[i]<<std::endl;
            }
        }
        end = dtime();
        forward_time = end - start;

        start = dtime();
        // new
        newtransformer.inverseFourierTransform();

        // old
        transformer.inverseFourierTransform(Fconv_real.wptr(), Fconv_imag.wptr(), Mconv.wptr());
        for (int i = 0; i < size*size*size; i++) {
            if (fabs(realData[i] - (Mconv.wptr())[i]) > 1e-6) {
                std::cout<<realData[i]<<" "<<(Mconv.wptr())[i]<<std::endl;
            }
        }

        end = dtime();
        backward_time = end - start;

        std::cout<<"size : "<<size<<" init : "<<init_time<<" seconds, setdata : "<<set_time<<" seconds, FFT-forward : "<<forward_time<<" seconds, FFT_backward : "<<backward_time<<" seconds."<<" ";
        std::cout<<"memory : "<<size*size*size*8./1024./1024./1024.*2<<" GB."<<std::endl;

        Mconv.fini();Fconv_real.fini();Fconv_imag.fini();
        aFree(realData);fftw_free(complexData);
    }

    double t2 = dtime();
    std::cout<<(t2-t1)<<std::endl;

}