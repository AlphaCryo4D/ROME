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

#ifndef INITIALIZE_H_
#define INITIALIZE_H_

#include "mkl.h"

#include <assert.h>
#include <cmath>
#include <iostream>
#include <cstring> /*memcpy*/

#include "./macros.h"
#include "./error.h"
#include "./mrcs.h"
#include "./fft_fftw3.h"
#include "./array_vector.h"

class ListOfImages {
public:
    // TBD : why this virtual deconstructor cannot work on my Mac
    //virtual ~ListOfImages() {}
    virtual int nr_images() = 0;
    virtual int imageSide() = 0;				// images are assumed to be square
    virtual const double* rImage(int i) = 0;	// Only reading the elements
    virtual		  double* wImage(int i) = 0;	// Only writing the elements (may read written data)
    virtual		  double* mImage(int i) = 0;	// Reading elements before writing them
};

class ListOfImagesFromDoubleVector : public ::ListOfImages {
	double * const _dv;
	int _nr_images;
	int _imageSide;
	int _imageStride;
public:
	ListOfImagesFromDoubleVector(double * dv, int nr_images, int imageSide, int eltsPerImage) 
		: _dv(dv), _nr_images(nr_images), _imageSide(imageSide), _imageStride(eltsPerImage) {
		assert(imageSide*imageSide == eltsPerImage);
	}
	virtual ~ListOfImagesFromDoubleVector() {}
	virtual int nr_images() { return _nr_images; } 
	virtual int imageSide() { return _imageSide; }
    virtual const double* rImage(int i) { return image(i); }	// Only reading the elements
    virtual		  double* wImage(int i) { return image(i); }	// Only writing the elements (may read written data)
    virtual		  double* mImage(int i) { return image(i); }	// Reading elements before writing them
private:
	double* image(int i) { assert(0 <= i && i < _nr_images); return _dv + i*_imageStride; }
};

// generate one gauss circle image reference
void genGaussCircleRef(float *Iref,int size);

// check data normalize or not
void checkNormalize(const Images& images,int size,double particle_diameter,double pixel_size);
void checkNormalize(double const *images,int size,int nr_images,double particle_diameter,double pixel_size);
bool checkNormalizedImage(const Image& img,int size,double particle_diameter,double pixel_size, double& sum, double& sum2);

// normalize the data
void normalizeData(float *Vec,int length);

//  solve W for GTM   --------------------------------------
void solveW(double *PHI,int K,int M,float *W,float *Y,int D);
void solveW(double *PHI,int K,int M,double *W,double *Y,int D);
void solveW(double *PHI,int K,int M,double *W,double *Y,int D,int cfm_num);

//  solve A[N][M]*X = B[N][D],result in B
//    Uses QR or LQ factorization to solve a overdetermined or underdetermined linear system with full rank matrix.
//    notice this will overwrite A and B!

void mkl_solveNotdetermined(double *A,int N,int M,double *B,int D);

void mkl_solveNotdetermined(float *A,int N,int M,float *B,int D);


#endif