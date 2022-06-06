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

#include "util.h"		// used for building precompiled headers on Windows

#include "pca_optimizer.h"

// change the parater order@date 2014.09.30
Pca::Pca(double *_T,int _N,int _D,double *_X,int _K,int _L,double *_PHI,double *_W,int _M,double* _Belta){

    T = _T;N = _N;D = _D;
    X = _X;K = _K;L = _L;

    PHI = _PHI;

    W = _W;M = _M;
    Belta = _Belta;

    // std::cout<<"construct PCA!\n"<<std::endl;
}


// using principal components,initialize W,Belta
void Pca::initializeW(){

	double *meanT = (double*)aMalloc(sizeof(double)*D,64);
	double *covs = (double*)aMalloc(sizeof(double)*D*D,64);
	//it will get covs and meanT
	mkl_cov(T,N,D,covs,meanT);  //symmetry
	
	// std::cout<<"covs = "<<std::endl;
	// prtMatrix(covs,D,D);
	// std::cout<<"covs complete."<<std::endl;
	std::cout<<"starting soving a "<<D<<"x"<<D<<" symmetrical real matrixs eigenvalue and eigenvectors."<<std::endl;

	double *eigVal = mkl_eig(covs,D);  //DESCENDING order

	std::cout<<"soving eig completed."<<std::endl;

	// covs is D*D demension
	// std::cout<<"eigvector = (opposite with matlab)"<<std::endl;
	// prtMatrix(covs,D,D);
	// std::cout<<"eigvalue = "<<std::endl;
	// prtMatrix(eigVal,1,D);

	// do first L element
	double *A = (double*)aMalloc(sizeof(double)*D*L,64);

#pragma simd   //////we can put sqrt() out !!!!!!!!!!!!!!!!!!!!!!!!
	for(int i = 0;i < D;i++)
		for(int j = 0;j < L;j++){
		  A[i*L+j] = covs[i*D+j]*sqrt(eigVal[j]);  //need check again
	}

	// initializeBelta
	if(L < D)   //if X demension is smaller than T demension
		*Belta = 1./eigVal[L];
	else
		*Belta = std::numeric_limits<double>::max();

	std::cout<<"belta(in initializeW) = "<<*Belta<<std::endl;

	// delete [] eigVal;
	aFree(eigVal);
	aFree(covs);

	// std::cout<<"A = "<<std::endl;
	// prtMatrix(A,D,L);

	double *normX = (double*)aMalloc(sizeof(double)*K*L,64);
	double *meanX = mkl_mean(X,K,L);
	double *stdX = stddev(X,K,L);  //not found suitable mkl function for this
#pragma simd
	for(int i = 0;i < K;i++)
		for(int j = 0;j < L;j++)
		  normX[i*L+j] = (X[i*L+j] - meanX[j])/stdX[j];

	// std::cout<<"normX = "<<std::endl;
	// prtMatrix(normX,K,L);
	// std::cout<<"stdx = "<<std::endl;
	// prtMatrix(stdX,1,L);
	// std::cout<<"mean = "<<mean(X,K)<<std::endl;  //this is not equal to matlab,but it is normal


	// FI*W = normX*A',solve FI*W(d) = (normX*A')(d)
	// here,K must larger than M,for Overdetermined system
	double *B = (double*)aMalloc(sizeof(double)*K*D,64);//normX*A'
	double *AA = (double*)aMalloc(sizeof(double)*K*M,64);

#pragma omp parallel for
	for(int k = 0;k < K;k++){
		for(int d = 0;d < D;d++){
		  B[k*D+d] = 0.;
		  for(int l = 0;l < L;l++)
		    B[k*D+d] += normX[k*L+l]*A[d*L+l];
		}
	}

	// rotation image
	// for(int k = 1;k < K;k++){
	// 	double angle = k/360.*2*3.1415926535897;
	// 	memcpy(B+k*D,B,sizeof(double)*D);
	// 	image_rotation(B+k*D,sqrt(D),angle);
	// }

	// delete [] A;
	aFree(A);
	aFree(normX);

	memcpy(AA,PHI,sizeof(double)*K*M);   //K*M matrix PHI

	// std::cout<<"A = "<<std::endl;
	// prtMatrix(AA,K,M);

	// std::cout<<"b = "<<std::endl;
	// prtMatrix(B,K,D);

	// this result silly different from own PBML version!!!!!!!
	//////!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	mkl_solveOverdetermined(AA,K,M,B,D);  //K > M acually not sure overdetermined or undetermined system 
	// std::cout<<"x = "<<std::endl;
	// prtMatrix(B,M,D);

	memcpy(W,B,sizeof(double)*M*D);

	memcpy(W+(M-1)*D,meanT,sizeof(double)*D);

	aFree(meanT);
	aFree(B);
	aFree(AA);

	// std::cout<<"W = "<<std::endl;
	// prtMatrix(W,M,1);

}

void Pca::initializeBelta(){

	// prtMatrix(PHI,1,M);
	double *Ydot = mkl_multiplyMat(PHI,K,M,W,M,D);//K*D

	double *DIST = mkl_DIST(Ydot,K,D);


	// std::cout<<"Ydot = "<<std::endl;
	// prtMatrix(Ydot,K,D);

	aFree(Ydot);

	for(int i = 0;i < K;i++)
	  DIST[i*K+i] += std::numeric_limits<double>::max();


	// std::cout<<"DIST = "<<std::endl;
	// prtMatrix(DIST,K,K);

	double BeltaTemp = 2./mkl_meanVec(mkl_min(DIST,K,K),K);
	std::cout<<"BeltaTemp(in initializeBelta) = "<<BeltaTemp<<std::endl;

	*Belta = (*Belta < BeltaTemp)?*Belta:BeltaTemp;

	std::cout<<"Belta(final) = "<<*Belta<<std::endl;

}

void Pca::run_pca(){
 	initializeW();
 	initializeBelta();
}


// Vector Statistical Library (VSL) Statistical Functions///////////////////
// covariation for each column,X is NxM,cols base!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
int Pca::mkl_cov(double *X,int N,int M,double *cov,double *mean){

	VSLSSTaskPtr task;

	// double *cov = (double*)aMalloc(sizeof(double)*M*M,64);
	// double *mean = (double*)aMalloc(sizeof(double)*M,64);

	MKL_INT p, n, xstorage;
	int status;
	/* initialize variables used in the computations of sample mean */
	p = M;
	n = N;
	xstorage = VSL_SS_MATRIX_STORAGE_COLS;   //can change this
	// mean = 0.0;
	/* create task */
	status = vsldSSNewTask( &task, &p, &n, &xstorage, X, 0, 0 );
	/* initialize task parameters */
	// status = vsldSSEditTask( task, VSL_SS_ED_MEAN, mean );
	status = vsldSSEditCovCor( task , mean , cov , &xstorage , 0 , 0 );

	/* compute mean using SS fast method */
	// status = vsldSSCompute(task, VSL_SS_MEAN, VSL_SS_METHOD_FAST );
	status = vsldSSCompute(task, VSL_SS_COV, VSL_SS_METHOD_FAST );
	/* deallocate task resources */
	status = vslSSDeleteTask( &task );

	// aFree(mean);
	// return cov;
	return status;
}

// compute the sample mean for a one-dimensional dataset, initialize a variable for the mean value
// from <<intel nath kernel library reference>>-2834 example
double Pca::mkl_meanVec(double *X,int N){
	VSLSSTaskPtr task;
	double mean;
	MKL_INT p, n, xstorage;
	int status;
	/* initialize variables used in the computations of sample mean */
	p = 1;//Dimension of the task, number of variables;
	n = N;//Number of observations
	xstorage = VSL_SS_MATRIX_STORAGE_ROWS;
	mean = 0.0;
	/* create task */
	status = vsldSSNewTask( &task, &p, &n, &xstorage, X, 0, 0 );
	/* initialize task parameters */
	status = vsldSSEditTask( task, VSL_SS_ED_MEAN, &mean );
	/* compute mean using SS fast method */
	status = vsldSSCompute(task, VSL_SS_MEAN, VSL_SS_METHOD_FAST );
	/* deallocate task resources */
	status = vslSSDeleteTask( &task );

	// std::cout<<"it is ok for mkl_meanVec.mean = "<<mean<<std::endl;

	// for(int n = 0;n < N;n++)
	// 	if(std::isnan(X[n]))
	// 	    std::cout<<"n = "<<n<<std::endl;
	// std::cout<<"in mkl_meanvec:isnan??"<<std::isnan(mean)<<std::endl;

	return mean;
}

// compute the sample mean for a two-dimensional dataset, initialize a variable for the mean value
double *Pca::mkl_mean(double *X,int N,int M){
	VSLSSTaskPtr task;

	double *mean = (double*)aMalloc(sizeof(double)*M,64);

	MKL_INT p, n, xstorage;
	int status;
	/* initialize variables used in the computations of sample mean */
	p = M;//Dimension of the task, number of variables;
	n = N;//Number of observations
	xstorage = VSL_SS_MATRIX_STORAGE_COLS;
	/* create task */
	status = vsldSSNewTask( &task, &p, &n, &xstorage, X, 0, 0 );
	/* initialize task parameters */
	status = vsldSSEditTask( task, VSL_SS_ED_MEAN, mean );
	/* compute mean using SS fast method */
	status = vsldSSCompute(task, VSL_SS_MEAN, VSL_SS_METHOD_FAST );
	/* deallocate task resources */
	status = vslSSDeleteTask( &task );

	// std::cout<<"it is ok for mkl_mean()."<<std::endl;

	// prtMatrix(mean,1,10);
	// for(int m = 0;m < M;m++)
	// 	if(std::isnan(X[m]))
	// 	    std::cout<<"m = "<<m<<std::endl;

	return mean;
}


// compute the sample minimum for a two-dimensional dataset, initialize a variable for the mean value
double *Pca::mkl_min(double *X,int N,int M){
	VSLSSTaskPtr task;

	double *min = (double*)aMalloc(sizeof(double)*M,64);

	for(int m = 0;m < M;m++) /////////!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
		min[m] = X[m];      //must do this initializition!!!!!!!!!!!!!!!!!!!!!!

	MKL_INT p, n, xstorage;
	int status;
	/* initialize variables used in the computations of sample min */
	p = M;//Dimension of the task, number of variables;
	n = N;//Number of observations
	xstorage = VSL_SS_MATRIX_STORAGE_COLS;
	/* create task */
	status = vsldSSNewTask( &task, &p, &n, &xstorage, X, 0, 0 );
	/* initialize task parameters */
	status = vsldSSEditTask( task, VSL_SS_ED_MIN, min );
	/* compute min using SS fast method */
	status = vsldSSCompute(task, VSL_SS_MIN, VSL_SS_METHOD_FAST );
	/* deallocate task resources */
	status = vslSSDeleteTask( &task );

	for(int m = 0;m < M;m++)
		if(std::isnan(min[m]))
		    std::cout<<"m = "<<m<<std::endl;

	return min;
}

// computes all the eigenvalues of a real symmetric tridiagonal matrix A
// DESCENDING order!!!!!!!!!!!!!!!!!!!!!!
// for 16cpus core,N = 10000,it costs 214.688 s..................
double *Pca::mkl_eig(double *A,int N){

	double* d = (double*)aMalloc(sizeof(double)*N,64);
	double* e = (double*)aMalloc(sizeof(double)*N,64);
	double* tau = (double*)aMalloc(sizeof(double)*N,64);
	
	// mkl_set_num_threads(32);

	//Reduces a real symmetric matrix to tridiagonal form.
	//lapack_int LAPACKE_<?>sytrd ( int matrix_layout , char uplo , lapack_int n , <datatype>* a ,
	//lapack_int lda , <datatype>* d , <datatype>* e , <datatype>* tau );
	LAPACKE_dsytrd (LAPACK_ROW_MAJOR,'U',N,A,N,d,e,tau);

	//lapack_int LAPACKE_<?>sterf ( lapack_int n , <datatype>* d , <datatype>* e );
	// int info = LAPACKE_dsterf(N,d,e);

	//Generates the real orthogonal matrix Q determined by ?sytrd
	//lapack_int LAPACKE_<?>orgtr ( int matrix_layout , char uplo , lapack_int n , <datatype>* a ,
	//lapack_int lda , const <datatype>* tau );
    int info = LAPACKE_dorgtr (LAPACK_ROW_MAJOR,'U',N,A,N,tau);

    if(info < 0)
    	std::cerr<<"For dorgtr,"<<(-info)<<"-th parameter had an illegal value."<<std::endl;

	//Computes all eigenvalues and eigenvectors of asymmetric or Hermitian matrix reduced to tridiagonal
	//form (QR algorithm).
	//lapack_int LAPACKE_dsteqr ( int matrix_layout , char compz , lapack_int n , double* d ,
	//double* e , double* z , lapack_int ldz );
	info = LAPACKE_dsteqr (LAPACK_ROW_MAJOR,'V',N,d,e,A,N);

	if(info > 0)
		std::cerr<<"For dorgtr,"<<info<<" off-diagonal elements have not converged to zero."<<std::endl;
	else if(info < 0)
		std::cerr<<"For dorgtr,"<<(-info)<<"-th parameter had an illegal value."<<std::endl;


// exchenge the eigvalue and eigvector
#pragma simd
	for(int i = 0;i < N;i++){
		e[i] = d[N-1-i];
	}

	for(int i = 0;i < N;i++){
		//d as temp
#pragma simd
		for(int j = 0;j < N;j++)
			d[j] = A[i*N+N-1-j];
#pragma simd
		for(int j = 0;j < N;j++)
			A[i*N+j] = d[j];
	}


	aFree(d);
	aFree(tau);
	return e;
}


// Uses QR or LQ factorization to solve a overdetermined or underdetermined linear system with full rank matrix.
// this silly differents from my own version for e-13 small number!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
void Pca::mkl_solveOverdetermined(double *A,int N,int M,double *B,int D){

	// Uses QR or LQ factorization to solve a overdetermined or underdetermined linear system with full rank matrix.
	// lapack_int LAPACKE_<?>gels ( int matrix_layout , char trans , lapack_int m , lapack_int n ,
	// lapack_int nrhs , <datatype>* a , lapack_int lda , <datatype>* b , lapack_int ldb );

	int info = LAPACKE_dgels ( LAPACK_ROW_MAJOR , 'N' , N , M ,D , A , M , B , D );

	if(info > 0)
		std::cerr<<"For dgels,the "<<info<<"-th diagonal element of the triangular factor of A is zero, so\
					that A does not have full rank; the least squares solution could not be\
					computed"<<std::endl;
	else if(info < 0)
		std::cerr<<"For dgels,the "<<(-info)<<"-th parameter had an illegal value."<<std::endl;

}

// BLAS
// calculate the A[N][M] * B[M][D]
double *Pca::mkl_multiplyMat(double *X,int N,int M,double *Y,int M2,int D){

	double *C = (double*)aMalloc(sizeof(double)*N*D,64);
	// this fucntion calculate [C := alpha*op(A)*op(B) + beta*C]
	// void cblas_dgemm ( const CBLAS_LAYOUT Layout , const CBLAS_TRANSPOSE transa , const
	// CBLAS_TRANSPOSE transb , const MKL_INT m , const MKL_INT n , const MKL_INT k , const double
	// alpha , const double a , const MKL_INT lda , const double *b , const MKL_INT ldb , const
	// double beta , double *c , const MKL_INT ldc );
	cblas_dgemm (CblasRowMajor,CblasNoTrans,CblasNoTrans,N,D,M,1.,X,M,Y,D,0.,C,D);

	return C;
}

// not found suitable mkl function for this
double *Pca::mkl_DIST(double *X,int N,int M){

	std::cout<<"calling DIST,it may take lots of time."<<std::endl;

	double *DIST = (double*)aMalloc(sizeof(double)*N*N,64);
	
	double distance;
	int index1,index2;
	bool flag = true;
#pragma omp parallel for collapse(2) private(distance,index1,index2)
	for(int i = 0;i < N;i++){
	    for(int j = 0;j < N;j++){//DIST is symmetry
	    	if(j < i) continue;//cannot use j = i index
	    	distance = 0.;
	    	index1 = i*M;  //change here!!!!!!!!!!!!!!!!!!!!!!!!!!
	    	index2 = j*M;   //from N ==> M !!!!!!!!!!!!!!!!!!!!!!!
	#pragma simd reduction(+:distance) //may been slowwd down the speed
	    	for (int m = 0; m < M; m++)
	    		distance += (X[index1+m]-X[index2+m])*(X[index1+m]-X[index2+m]);
	    	// std::cout<<"it is ok for mkl_DIST."<<std::endl;
	    	// DIST[i*M+j] = DIST[j*M+i] = distance;
	    	// prtMatrix(X,1,1);
	    	// exit(EXIT_FAILURE);
	    	DIST[i*N+j] = DIST[j*N+i] = distance;

	    	// if(std::isnan(distance)) {
	    	// 	std::cout<<"i = "<<i<<",j = "<<j<<std::endl;
	    	// }
	    	// if((i*N + j)% 10000000 == 0){
	    	// 	std::cout<<"i*N+j = "<<(i*N + j)<<std::endl;
	    	// 	std::cout<<"distance = "<<distance<<std::endl;
	    	// }

	    }
	}
	// std::cout<<"it is ok for mkl_DIST."<<std::endl;
		return DIST;
}

// not found standard deviation in mkl
void Pca::mkl_stddev(){
	std::cout<<"not found."<<std::endl;
}

// get Matrix's column or row means
double *Pca::mean(double *X,int rows,int cols,bool doColumn/** = true**/){
    double *means;
    if(doColumn == true){ //calculate each column's mean
        means = vNew(double,cols);
        double sum;
        for(int j = 0;j < cols;j++){
            sum = 0;
            for(int i = 0;i < rows;i++)
                sum += X[i*cols+j];
            means[j] = sum/rows;
        }
    }
    else{
        means = vNew(double,rows);
        double sum;
        for(int i = 0;i < rows;i++){
            sum = 0;
            for(int j = 0;j < cols;j++)
                sum += X[i*cols+j];
            means[i] = sum/cols;
        }
    }
    return means;
}

// standard deviation
double *Pca::stddev(double *X,int rows,int cols,bool doColumn/** = true**/){
    double *stdMeanX = mean(X,rows,cols,doColumn);
    double sum;
    if(doColumn == true){
        if(rows < 2) throw("not existed stdandard deviation.");
        sum = 0;
        for(int j = 0;j < cols;j++){
            sum = 0;
            for(int i = 0;i < rows;i++)  ///index by column,not effective
                sum += (X[i*cols+j] - stdMeanX[j])*(X[i*cols+j] - stdMeanX[j]);
            stdMeanX[j] = sqrt(sum/(rows-1));
        }
    }
    else{
        sum = 0;
        for(int i = 0;i < rows;i++){
            if(cols < 2) throw("not existed stdandard deviation.");
            sum = 0;
            for(int j = 0;j < cols;j++)
                sum += (X[i*cols+j] - stdMeanX[i])*(X[i*cols+j] - stdMeanX[i]);
            stdMeanX[i] = sqrt(sum/(cols-1));
        }
    }
    return stdMeanX;
}
