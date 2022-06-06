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

 /***@date 2014.08.28
 /***@version 1.00
/@description implement of principal component analysis.
this code can be compiled fine for g++ or intel cimpuler.
*/
/***@date 2014.09.30
1)change PCA constructor parameter order
2)modify initializePHI()
***/
/***@date 2016.08.07
 re-add this to rome
 ***/
#ifndef _PCA_H_
#define _PCA_H_

#include "./resmap_util.h"

#define DEBUG_PCA

class Pca{
public:
	double *X;  // K*L
	int K;
	int L;

	double *PHI;// K*M

	double *T;// N*D
	int N;
	int D;

	double *W;// M*D
	int M;

	double *MU;
	int Mnl;
	double *infos;

	double *Belta;

	double sigma;

public:
	// constructor
	// input parameter,data space : input:T[N][D],latent space : X[K][L],MU[Mnl][L],MU_infos
	// output parameter,biases : PHI[K][M],W[M=Mnl+L+1][D],Belta
	Pca(double *_T,int _N,int _D,double *_X,int _K,int _L,double *_PHI,double *_W,int _M,double* _Belta);

	~Pca(){}

	// using principal components,initialize W,Belta
	void initializeW();
	void initializeBelta();

	// run the Pca algorithm
	void run_pca();

private:

    // -------- define some basic function using mkl library

	// Vector Statistical Library (VSL) Statistical Functions///////////////////
	// covariation for each column,X is NxM,cols base!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	int mkl_cov(double *X,int N,int M,double *cov,double *mean);

	// compute the sample mean for a one-dimensional dataset, initialize a variable for the mean value
	// from <<intel nath kernel library reference>>-2834 example
	double mkl_meanVec(double *X,int N);

	// compute the sample mean for a two-dimensional dataset, initialize a variable for the mean value
	double *mkl_mean(double *X,int N,int M);


	// compute the sample minimum for a two-dimensional dataset, initialize a variable for the mean value
	double *mkl_min(double *X,int N,int M);


	// computes all the eigenvalues of a real symmetric tridiagonal matrix A
	// DESCENDING order!!!!!!!!!!!!!!!!!!!!!!
	// for 16cpus core,N = 10000,it costs 214.688 s..................
	double *mkl_eig(double *A,int N);

	// Uses QR or LQ factorization to solve a overdetermined or underdetermined linear system with full rank matrix.
	// this silly differents from my own version for e-13 small number!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	void mkl_solveOverdetermined(double *A,int N,int M,double *B,int D);

	// BLAS
	// calculate the A[N][M] * B[M][D]
	double *mkl_multiplyMat(double *X,int N,int M,double *Y,int M2,int D);

	// not found suitable mkl function for this
	double *mkl_DIST(double *X,int N,int M);

	// not found standard deviation in mkl
	void mkl_stddev();
    
    
    // -----------  define some basic mathematical function
    // get Matrix's column or row means
    double *mean(double *X,int rows,int cols,bool doColumn = true);
    // standard deviation
    double *stddev(double *X,int rows,int cols,bool doColumn = true);
    
};

#undef DEBUG_PCA

#endif
