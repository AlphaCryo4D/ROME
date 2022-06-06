/***************************************************************************
 *
 * Authors: "Jiayi (Timmy) Wu, Yongbei(Glow) Ma, Youdong (Jack) Mao"
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

#ifndef FFT_DFTI_H_
#define FFT_DFTI_H_

#include "./resmap_util.h"


class FourierTransformerDFTI
{
public:
    // The dimension of 2D data
    int ysize,xsize;
    
    // 2D data in Real space
	float* realData;
    
    // 2D data in Fourier space
	MKL_Complex8 *complexData;
    
    // FFTW forward and backward plan
	MKL_LONG status;
	DFTI_DESCRIPTOR_HANDLE forward_plan;
	DFTI_DESCRIPTOR_HANDLE backward_plan;
    
public:
    // Constructor
	FourierTransformerDFTI(int _ysize,int _xsize){
        
		ysize = _ysize;xsize=_xsize;
        
		realData = (float*)aMalloc(sizeof(float)*ysize*xsize,64);
		complexData = (MKL_Complex8*)aMalloc(ysize*(xsize/2+1)*sizeof(MKL_Complex8),64);

		// Setup for forward FFT plan
		MKL_LONG dim[2];dim[0] = ysize;dim[1] = xsize;
		MKL_LONG is[3]; is[0] = 0; is[1] = xsize; is[2] = 1;
		MKL_LONG os[3]; os[0] = 0; os[1] = xsize/2+1; os[2] = 1;
		status = DftiCreateDescriptor(&forward_plan, DFTI_SINGLE,DFTI_REAL, 2, dim);
		status = DftiSetValue(forward_plan, DFTI_CONJUGATE_EVEN_STORAGE, DFTI_COMPLEX_COMPLEX);
		status = DftiSetValue(forward_plan, DFTI_PACKED_FORMAT, DFTI_CCE_FORMAT);
		status = DftiSetValue(forward_plan, DFTI_PLACEMENT, DFTI_NOT_INPLACE);
		status = DftiSetValue(forward_plan, DFTI_INPUT_STRIDES, is);
		status = DftiSetValue(forward_plan, DFTI_OUTPUT_STRIDES, os);
		status = DftiCommitDescriptor(forward_plan);

		// Setup for backword FFT plan
		is[1] = xsize/2+1;os[1] = xsize;
		status = DftiCreateDescriptor(&backward_plan, DFTI_SINGLE,DFTI_REAL, 2, dim);
		status = DftiSetValue(backward_plan, DFTI_CONJUGATE_EVEN_STORAGE, DFTI_COMPLEX_COMPLEX);
		status = DftiSetValue(backward_plan, DFTI_PACKED_FORMAT, DFTI_CCE_FORMAT);
		status = DftiSetValue(backward_plan, DFTI_PLACEMENT, DFTI_NOT_INPLACE);
		status = DftiSetValue(backward_plan, DFTI_INPUT_STRIDES, is);
		status = DftiSetValue(backward_plan, DFTI_OUTPUT_STRIDES, os);
		status = DftiCommitDescriptor(backward_plan);
	}
    
    // Destructor
	~FourierTransformerDFTI(){
		aFree(realData);
		aFree(complexData);
		status = DftiFreeDescriptor(&forward_plan);
		status = DftiFreeDescriptor(&backward_plan);
	}


	// FFT forward 2D data
	void FourierTransform(const float* data, MKL_Complex8* Fdata){
        
        int size2 = ysize*xsize;
        realData[0:size2] = data[0:size2];

		status = DftiComputeForward(forward_plan,realData,Fdata);

	}
	
	//FFT backward 2D data
	void inverseFourierTransform(const float* FdataReal,const float* FdataImag,float* data){

		int Fsize2 = ysize*(xsize/2+1);
		for(int i = 0;i < Fsize2;i++){
			complexData[i].real = FdataReal[i];
			complexData[i].imag = FdataImag[i];
		}

		status = DftiComputeBackward(backward_plan,complexData,data);

	}

};



#endif