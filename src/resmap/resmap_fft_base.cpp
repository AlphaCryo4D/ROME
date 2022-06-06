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

#include "resmap_util.h"		// used for building precompiled headers on Windows

#include "resmap_fft_base.h"

#define CLASS_NAME FFTWBase
#define DATA_TYPE double
#define FFT_TYPE(v) fftw_##v
void CLASS_NAME::init(DATA_TYPE* _real,FFT_TYPE(complex)* _complex,size_t _xsize,size_t _ysize,size_t _zsize /*= 1*/,int _nr_threads /*= 1*/)
{
    //
    nr_threads = _nr_threads;
    int val = FFT_TYPE(init_threads)(); //val returned non zero value
    if(!val) ERROR_REPORT("Connot init multi-threads for FFTW3.");
    FFT_TYPE(plan_with_nthreads)(nr_threads);
    
    set(_real, _complex, _xsize, _ysize, _zsize);
    
    create_plan();
}

void CLASS_NAME::fini()
{
    if (zsize) {
        delete_plan();
        FFT_TYPE(cleanup_threads)();
        zsize = 0;ysize = 0;xsize = 0;
        size3 = 0;Fsize3 = 0;
    }
}
    
// FFT forward 2D/3D data
void CLASS_NAME::FourierTransform(bool normalize /*= true*/)
{
    forward_execute();
    // normalize the data
    if (normalize == true) normalizeComplex();
}
    
// FFT backward 2D/3D data
void CLASS_NAME::inverseFourierTransform(){
    backward_execute();
}
    
// fftw plan manager
void CLASS_NAME::reset_plan(DATA_TYPE* _real,FFT_TYPE(complex)* _complex,size_t _xsize,size_t _ysize,size_t _zsize /*= 1*/)
{
    delete_plan();
    
    set(_real, _complex, _xsize, _ysize, _zsize);
    
    create_plan();
}
    

void CLASS_NAME::create_plan()
{
    if (zsize == 1) { // 2D dimension
        forward_plan = FFT_TYPE(plan_dft_r2c_2d)(ysize,xsize,real,complex,FFTW_ESTIMATE);
        backward_plan = FFT_TYPE(plan_dft_c2r_2d)(ysize,xsize,complex,real,FFTW_ESTIMATE);
    }
    else { // 3D dimension....NOTE : FFTW_ESTIMATE may change to FFTWF_ESTIMATE for float datatype
        forward_plan = FFT_TYPE(plan_dft_r2c_3d)(zsize,ysize,xsize,real,complex,FFTW_ESTIMATE);
        backward_plan = FFT_TYPE(plan_dft_c2r_3d)(zsize,ysize,xsize,complex,real,FFTW_ESTIMATE);
    }
}
    
void CLASS_NAME::delete_plan()
{
    FFT_TYPE(destroy_plan)(forward_plan);
    FFT_TYPE(destroy_plan)(backward_plan);
    FFT_TYPE(cleanup)();
}
    
void CLASS_NAME::set(DATA_TYPE* _real,FFT_TYPE(complex)* _complex,size_t _xsize,size_t _ysize,size_t _zsize)
{
    xsize = _xsize;ysize = _ysize;zsize = _zsize;
    size_t size2 = ysize*xsize;size3 = zsize*ysize*xsize;
    size_t Fsize2 = ysize*(xsize/2+1);Fsize3 = zsize*Fsize2;
    real = _real;
    complex = _complex;
}
    
void CLASS_NAME::normalizeComplex()
{
    auto complex_data = (DATA_TYPE*)complex;
    // for (int i = 0; i < Fsize3; i++){
    //     complex[i][0] = complex[i][0]/size3;
    //     complex[i][1] = complex[i][1]/size3;
    // }
    for (int i = 0; i < 2*Fsize3; i++) {
        complex_data[i] = complex_data[i]/size3;
    }
}
    
// execute the fftw plan
void CLASS_NAME::forward_execute(){FFT_TYPE(execute)(forward_plan);}
void CLASS_NAME::backward_execute(){FFT_TYPE(execute)(backward_plan);};

#undef CLASS_NAME
#undef DATA_TYPE
#undef FFT_TYPE

//    ----------------------------------------------------------------------------   //
//    ----------   duplicate above code and change MACRO #define   ---------------   //
//    ----------------------------------------------------------------------------   //

#define CLASS_NAME FFTWFBase
#define DATA_TYPE float
#define FFT_TYPE(v) fftwf_##v
void CLASS_NAME::init(DATA_TYPE* _real,FFT_TYPE(complex)* _complex,size_t _xsize,size_t _ysize,size_t _zsize /*= 1*/,int _nr_threads /*= 1*/)
{
    //
    nr_threads = _nr_threads;
    int val = FFT_TYPE(init_threads)(); //val returned non zero value
    if(!val) ERROR_REPORT("Connot init multi-threads for FFTW3.");
    FFT_TYPE(plan_with_nthreads)(nr_threads);
    
    set(_real, _complex, _xsize, _ysize, _zsize);
    
    create_plan();
}

void CLASS_NAME::fini()
{
    if (zsize) {
        delete_plan();
        FFT_TYPE(cleanup_threads)();
        zsize = 0;ysize = 0;xsize = 0;
        size3 = 0;Fsize3 = 0;
    }
}

// FFT forward 2D/3D data
void CLASS_NAME::FourierTransform(bool normalize /*= true*/)
{
    forward_execute();
    // normalize the data
    if (normalize == true) normalizeComplex();
}

// FFT backward 2D/3D data
void CLASS_NAME::inverseFourierTransform(){
    backward_execute();
}

// fftw plan manager
void CLASS_NAME::reset_plan(DATA_TYPE* _real,FFT_TYPE(complex)* _complex,size_t _xsize,size_t _ysize,size_t _zsize /*= 1*/)
{
    delete_plan();
    
    set(_real, _complex, _xsize, _ysize, _zsize);
    
    create_plan();
}


void CLASS_NAME::create_plan()
{
    if (zsize == 1) { // 2D dimension
        forward_plan = FFT_TYPE(plan_dft_r2c_2d)(ysize,xsize,real,complex,FFTW_ESTIMATE);
        backward_plan = FFT_TYPE(plan_dft_c2r_2d)(ysize,xsize,complex,real,FFTW_ESTIMATE);
    }
    else { // 3D dimension...
        forward_plan = FFT_TYPE(plan_dft_r2c_3d)(zsize,ysize,xsize,real,complex,FFTW_ESTIMATE);
        backward_plan = FFT_TYPE(plan_dft_c2r_3d)(zsize,ysize,xsize,complex,real,FFTW_ESTIMATE);
    }
}

void CLASS_NAME::delete_plan()
{
    FFT_TYPE(destroy_plan)(forward_plan);
    FFT_TYPE(destroy_plan)(backward_plan);
    FFT_TYPE(cleanup)();
}

void CLASS_NAME::set(DATA_TYPE* _real,FFT_TYPE(complex)* _complex,size_t _xsize,size_t _ysize,size_t _zsize)
{
    xsize = _xsize;ysize = _ysize;zsize = _zsize;
    size_t size2 = ysize*xsize;size3 = zsize*ysize*xsize;
    size_t Fsize2 = ysize*(xsize/2+1);Fsize3 = zsize*Fsize2;
    real = _real;
    complex = _complex;
}

void CLASS_NAME::normalizeComplex()
{
    for (size_t i = 0; i < Fsize3; i++){
        complex[i][0] = complex[i][0]/size3;
        complex[i][1] = complex[i][1]/size3;
    }
}

// execute the fftw plan
void CLASS_NAME::forward_execute(){FFT_TYPE(execute)(forward_plan);}
void CLASS_NAME::backward_execute(){FFT_TYPE(execute)(backward_plan);};

#undef CLASS_NAME
#undef DATA_TYPE
#undef FFT_TYPE