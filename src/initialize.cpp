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

#include "initialize.h"

void genGaussCircleRef(float *Iref,int size)
{
    int NX = size;
    int NY = size;
    int DX = (NX/2) + 1;
    int DY = (NY/2) + 1;
    int SX = 0.25*NX;
    int SY = 0.25*NY;
    int norder = 1;
    
    double QUADPI = 3.1415926535897932;
    double GNM = 1.0 / SX / SY / 2.0 / QUADPI;
    //TNM   = ALOG(1.0 / TINY(GNM));
    
    double SXSQ  = SX * SX;
    double SYSQ  = SY * SY;
    
    for (int i = 1; i <= NY; i++) {
        for (int k = 1; k <= NX; k++) {
            double EEE = 0.5 * ((k-DX)*(k-DX) / SXSQ + (i-DY)*(i-DY) / SYSQ);
            if (EEE >= 87.)
                Iref[(i-1)*NX+k-1] = 0.0;
            else{
                EEE   = 0.5 * pow((2*EEE),norder);
                Iref[(i-1)*NX+k-1] = GNM * exp(-EEE);
            }
            //std::cout<<EEE<<std::endl;
        }
    }
}

void checkNormalize(const Images& images,int size,double particle_diameter,double pixel_size)
{
    int nr_images = images.size();
    assert(nr_images>0);
    for (int iimage = 0; iimage < nr_images; iimage++) {
        // Read image from disc
        double sum, sum2;
        if (checkNormalizedImage(images[iimage], size, particle_diameter, pixel_size, sum, sum2)) continue;
        std::cerr << " fn_img= " << iimage << " bg_avg= " << sum << " bg_stddev= " << sum2 << std::endl;
        ERROR_REPORT("ERROR: It appears that these images have not been normalised?");
    }
}

void checkNormalize(double const *images,int size,int nr_images,double particle_diameter,double pixel_size)
{
    Image tmpImg(size*size);auto tmpImg_data = tmpImg.wptrAll();
	for (int iimage = 0; iimage < nr_images; iimage++) {
        // Read image from disc
		double const* img = images+iimage*size*size;
        for (int i = 0; i < size*size; i++) tmpImg_data[i] = img[i];
		double sum, sum2;
		if (checkNormalizedImage(tmpImg, size, particle_diameter, pixel_size, sum, sum2)) continue;
        std::cerr << " fn_img= " << iimage << " bg_avg= " << sum << " bg_stddev= " << sum2 << std::endl;
        ERROR_REPORT("ERROR: It appears that these images have not been normalised?");
	}
}
        
bool checkNormalizedImage(const Image& img,int size,double particle_diameter,double pixel_size, double& sum, double& sum2)
{
    auto img_data = img.rptrAll();
    // Check that the average in the noise area is approximately zero and the stddev is one
    int bg_radius2 = round(particle_diameter / (2. * pixel_size));
    bg_radius2 *= bg_radius2;
	sum = 0.;
	sum2 = 0.;
    double nn = 0.;
    int ori_size_shift =  -(long int)((float) (size) / 2.0);
    for (int i = 0; i < size; i++){
        for (int j = 0; j < size; j++)
        {
            int x = i + ori_size_shift;
            int y = j + ori_size_shift;
            if (x*x+y*y > bg_radius2)
            {
                double bit = img_data[i*size+j];
                sum += bit;
                sum2 += bit*bit;
                nn += 1.;
            }
        }
    }
    // stddev
    sum2 -= sum*sum/nn;
    sum2 = sqrt(sum2/nn);
    // average
    sum /= nn;
    
    // Average should be close to zero, i.e. max +/-50% of stddev...
    // Stddev should be close to one, i.e. larger than 0.5 and smaller than 2)
	if (fabs(sum/sum2) > 0.5 || sum2 < 0.5 || sum2 > 2.0) return false;
	
	return true;
}

void normalizeData(float *image,int size2){
    float mean = 0.;
    for(int i = 0;i < size2;i++)
        mean += image[i];
    mean /= size2;
    
    float norm = 0.;
    for(int i = 0;i < size2;i++)
        norm += (image[i]-mean)*(image[i]-mean);
    norm = sqrt(norm/(size2));
    if (norm == 0) {
        return;
    }
    // normalizeData
    for(int i = 0;i < size2;i++)
        image[i] = (image[i]-mean)/norm;
}


void solveW(double *PHI,int K,int M,float *W,float *Y,int D){

	int maxKM = K > M?K:M;
	float *B = (float*)aMalloc(sizeof(float)*maxKM*D,64);
	float *A = (float*)aMalloc(sizeof(float)*K*M,64);


	// memcpy(B,Y,sizeof(float)*K*D);
    for (int i = 0; i < K*D; i++) {
        B[i] = Y[i];
    }
	// memcpy(A,PHI,sizeof(float)*K*M);
	for(int i = 0;i < K*M;i++)
		A[i] = PHI[i];

	mkl_solveNotdetermined(A,K,M,B,D);

	// memcpy(W,B,sizeof(float)*M*D);
    for (int i = 0; i < M*D; i++) {
        W[i] = B[i];
    }
    
	aFree(B);
	aFree(A);
}

void solveW(double *PHI,int K,int M,double *W,double *Y,int D){

	int maxKM = K > M?K:M;
	double *B = (double*)aMalloc(sizeof(double)*maxKM*D,64);
	double *A = (double*)aMalloc(sizeof(double)*K*M,64);

	// memcpy(B,Y,sizeof(double)*K*D);
    for (int i = 0; i < K*D; i++) {
        B[i] = Y[i];
    }
	// memcpy(A,PHI,sizeof(double)*K*M);
    for (int i = 0; i < K*M; i++) {
        A[i] = PHI[i];
    }
    
	mkl_solveNotdetermined(A,K,M,B,D);

	// memcpy(W,B,sizeof(double)*M*D);
    for (int i = 0; i < M*D; i++) {
        W[i] = B[i];
    }
    
	aFree(B);
	aFree(A);
}
	
void solveW(double *PHI,int K,int M,double *W,double *Y,int D,int cfm_num){

	int cfm_K = K/cfm_num;
	int maxKM = cfm_K > M?cfm_K:M;
	double *A = (double*)aMalloc(sizeof(double)*cfm_K*M,64);
	double *B = (double*)aMalloc(sizeof(double)*maxKM*D,64);


	// memcpy(B,Y,sizeof(double)*K*D);
	for(int cfm_index = 0;cfm_index < cfm_num;cfm_index++){

        // std::cout<<cfm_index<<std::endl;

		memcpy(A,PHI+cfm_index*cfm_K*M,sizeof(double)*cfm_K*M);

		// memcpy(B,Y+cfm_index*cfm_K*D,sizeof(double)*cfm_K*D);
		memcpy(B,Y+cfm_index*cfm_K*D,sizeof(double)*cfm_K*D);

		mkl_solveNotdetermined(A,cfm_K,M,B,D);

		memcpy(W+cfm_index*M*D,B,sizeof(double)*M*D);

	}

	aFree(B);
	aFree(A);
}


void mkl_solveNotdetermined(double *A,int N,int M,double *B,int D){

	// Uses QR or LQ factorization to solve a overdetermined or underdetermined linear system with full rank matrix.
	// lapack_int LAPACKE_<?>gels ( int matrix_layout , char trans , lapack_int m , lapack_int n ,
	// lapack_int nrhs , <datatype>* a , lapack_int lda , <datatype>* b , lapack_int ldb );
	// std::cout<<"ok???"<<std::endl;
	int info = LAPACKE_dgels ( LAPACK_ROW_MAJOR , 'N' , N , M ,D , A , M , B , D );

    if(info > 0)
    {
		std::cerr<<"For dgels,the "<<info<<"-th diagonal element of the triangular factor of A is zero, so\
					that A does not have full rank; the least squares solution could not be\
					computed"<<std::endl;
        ERROR_REPORT("");
    }
    else if(info < 0)
    {
		std::cerr<<"For dgels,the "<<(-info)<<"-th parameter had an illegal value."<<std::endl;
        ERROR_REPORT("");
    }

}

void mkl_solveNotdetermined(float *A,int N,int M,float *B,int D){

	// Uses QR or LQ factorization to solve a overdetermined or underdetermined linear system with full rank matrix.
	// lapack_int LAPACKE_<?>gels ( int matrix_layout , char trans , lapack_int m , lapack_int n ,
	// lapack_int nrhs , <datatype>* a , lapack_int lda , <datatype>* b , lapack_int ldb );
	// std::cout<<"ok???"<<std::endl;
	int info = LAPACKE_sgels ( LAPACK_ROW_MAJOR , 'N' , N , M ,D , A , M , B , D );

    if(info > 0)
    {
		std::cerr<<"For dgels,the "<<info<<"-th diagonal element of the triangular factor of A is zero, so\
					that A does not have full rank; the least squares solution could not be\
					computed"<<std::endl;
        ERROR_REPORT("");
    }
    else if(info < 0)
    {
		std::cerr<<"For dgels,the "<<(-info)<<"-th parameter had an illegal value."<<std::endl;
        ERROR_REPORT("");
    }

}

