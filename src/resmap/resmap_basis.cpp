
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

#include "resmap_util.h"		// used for building precompiled headers on Windows

#include "./resmap_basis.h"

namespace GTMBasis
{
// debug flag
// #define DEBUG
    static double* alignedDoubleMalloc(size_t size,int align){return (double*)aMalloc(sizeof(double)*size,align);}
    static void alignedDoubleFree(double* &ptr){aFree(ptr);}
    
    static void prtMatrix(const double* A,int M,int N)
    {
        for (int i = 0; i < M; i++) {
            for (int j = 0; j < N; j++) {
                std::cout<<A[i*N+j]<<" ";
            }
            std::cout<<std::endl;
        }
    }
        
    // initialize PHI[K][M],M = Mnl+L+1
    void set_shift_basis(const double *X,int K,int L,
                         const double *MU,const double *MU_infos,int Mnl,
                         double *PHI,int M)
    {
        static double s = 1;
        double sigma;

        double *expPart = alignedDoubleMalloc(L,64);
        
        // calculate sigma for basic_sampling case
        for(int l = 0;l < L;l++){
            // sigma = (infos[l*3+1] - infos[l*3+0])/((1 << (int)infos[l*3+2]) - 1);
            sigma = (MU_infos[l*3+1] - MU_infos[l*3+0])/( (int)MU_infos[l*3+2] - 1);
            sigma *= s;
            expPart[l] = -1./(2.*sigma*sigma);
        }

    #ifdef DEBUG
        std::cout<<"sigma = "<<sigma<<std::endl;
    #endif

    #pragma omp parallel for collapse(2)
        for (int k = 0; k < K; k++)
            for (int m = 0; m < M; m++){   // for each PHI_km
                if(m+1 <= Mnl){
                    PHI[k*M+m] = 0.;
                    for(int l = 0;l < L;l++)
                        PHI[k*M+m] += expPart[l]*(X[k*L+l]-MU[m*L+l])*(X[k*L+l]-MU[m*L+l]);
                    PHI[k*M+m] = exp(PHI[k*M+m]);
                }
                // else if(m+1 <= Mnl + L)
                //     PHI[k*M+m] = X[k*L + m-Mnl];
                else
                    PHI[k*M+m] = 1;
        }

        alignedDoubleFree(expPart);

     #ifdef DEBUG
        std::cout<<"shift,initialize PHI!"<<std::endl;
        prtMatrix(PHI,K,M);
    #endif
    }

    // initialize PHI[K][M],M = Mnl+L+1
    void set_shift_basis2(const double *X,int K,int L,
                          const double *MU,const double *MU_infos,int Mnl,
                          double *PHI,int M)
    {
        static double s = 1;
        double sigma;

        double *expPart = alignedDoubleMalloc(L,64);

        // calculate sigma for basic_sampling case
        for(int l = 0;l < L;l++){
            // sigma = (infos[l*3+1] - infos[l*3+0])/((1 << (int)infos[l*3+2]) - 1);
            sigma = (MU_infos[l*3+1] - MU_infos[l*3+0])/( (int)MU_infos[l*3+2] - 1);
            sigma *= s;
            expPart[l] = -1./(2.*sigma*sigma);
        }

    #ifdef DEBUG
        std::cout<<"sigma = "<<sigma<<std::endl;
    #endif

    #pragma omp parallel for collapse(2)
        for (int k = 0; k < K; k++)
            for (int m = 0; m < M; m++){   // for each PHI_km
                if(m+1 <= Mnl){
                    PHI[k*M+m] = 0.;
                    int l = 0;
                    PHI[k*M+m] += expPart[l]*(X[k*L+l]-MU[m*L+l])*(X[k*L+l]-MU[m*L+l]);

                    PHI[k*M+m] = exp(PHI[k*M+m]);
                }
                else if(m+1 <= 2*Mnl){
                    PHI[k*M+m] = 0.;
                    int l = 1;
                    PHI[k*M+m] += expPart[l]*(X[k*L+l]-MU[m*L+l])*(X[k*L+l]-MU[m*L+l]);

                    PHI[k*M+m] = exp(PHI[k*M+m]);
                }
                // else if(m+1 <= Mnl + L)
                // PHI[k*M+m] = X[k*L + m-Mnl];
                else
                    PHI[k*M+m] = 1;
        }

        alignedDoubleFree(expPart);

     #ifdef DEBUG
        std::cout<<"initialize PHI!"<<std::endl;
        prtMatrix(PHI,K,M);
     #endif
    }

    // initialize PHI[K][M],M = Mnl,L = 1
    void set_rotation_basis(const double *X,int K,int L,
                            const double *MU,const double *MU_infos,int Mnl,
                            double *PHI,int M)
    {
        static double s = 2;
        double temp;
        
        double sigma = (MU_infos[1] - MU_infos[0])/( (int)MU_infos[2] - 1);
        s = 1;
        sigma *= s;
        // double sigma = 2*(MU[1]-MU[0]);
        // double sigma = pi/8;
        // double sigma = 2*pi*(MU[1]-MU[0]);
        // double sigma = 32*(MU[1]-MU[0]);
        // std::cout<<"sigma = "<<sigma<<std::endl;
        double expPart = -1./(2.*sigma*sigma);
        // double expPart = -1./(sigma*sigma);

    #pragma omp parallel for collapse(2) private(temp)
        for (int k = 0; k < K; k++)
            for (int m = 0; m < M; m++){   // for each PHI_km

                if(m+1 <= Mnl){

                    // for period case
                    // temp = 2*pi - fabs(X[k]-MU[m]);
                    // temp = fabs(X[k]-MU[m]) < temp?fabs(X[k]-MU[m]):temp;

                    // for un-period case
                    temp = fabs(X[k]-MU[m]);

                    PHI[k*M+m] = expPart*temp*temp;
                    PHI[k*M+m] = exp(PHI[k*M+m]);

                }
                // else if(m+1 <= Mnl + L)
                    // PHI[k*M+m] = X[k*L + m-Mnl];
                else
                    PHI[k*M+m] = 1;
        }


     #ifdef DEBUG
         std::cout<<"initialize PHI!"<<std::endl;
         prtMatrix(PHI,K,M);
     #endif
    }


    // initialize PHI[K][M],M = Mnl,L = 1,set mutil conformation rotation basis
    void set_mutilcfm_rotation_basis(const double *X,int K,int L,
                                     const double *MU,const double *MU_infos,int Mnl,
                                     double *PHI,int M)
    {
        static double s = 2;
        double temp;
        
        double sigma = (MU_infos[1] - MU_infos[0])/( (int)MU_infos[2] - 1);
        s = 1;
        sigma *= s;
        // double sigma = 2*(MU[1]-MU[0]);
        // double sigma = 3.1415*2*51/49;
        // double sigma = 2*pi*(MU[1]-MU[0]);
        std::cout<<"sigma = "<<sigma<<std::endl;
        double expPart = -1./(2.*sigma*sigma);
        
    #pragma omp parallel for collapse(2) private(temp)
        for (int k = 0; k < K; k++)
            for (int m = 0; m < M; m++){   // for each PHI_km

                if(m+1 <= Mnl){

                    // for peroid case
                    // temp = 2*pi - fabs(X[k*L]-MU[m]);
                    // temp = fabs(X[k*L]-MU[m]) < temp?fabs(X[k*L]-MU[m]):temp;

                    // for un-peroid case
                    temp = fabs(X[k*L]-MU[m]);
                    
                    PHI[k*M+m] = expPart*temp*temp;
                    PHI[k*M+m] = exp(PHI[k*M+m]);

                }
                // else if(m+1 <= Mnl + L-1) //L-1 ignore conformation item
                // PHI[k*M+m] = X[k*L + m-Mnl+1]/4;
                else
                    PHI[k*M+m] = 1;
        }


    #ifdef DEBUG
        std::cout<<"initialize PHI!"<<std::endl;
        prtMatrix(PHI,K,M);
    #endif
    }

    // initialize PHI[K][M],M = Mnl+L+1,notice the MU_infos last three element is rotation info.
    void set_shift_rotation_basis(const double *X,int K,int L,
                                  const double *MU,double *MU_infos,int Mnl,
                                  double *PHI,int M)
    {
        static double s = 1;
        double sigma;
        double temp;

        double *expPart = alignedDoubleMalloc(L,64);

        // calculate sigma for shift case
        for(int l = 0;l < L-1;l++){
            // sigma = (infos[l*3+1] - infos[l*3+0])/((1 << (int)infos[l*3+2]) - 1);
            sigma = (MU_infos[l*3+1] - MU_infos[l*3+0])/( (int)MU_infos[l*3+2] - 1);
            sigma *= s;
            expPart[l] = -1./(2.*sigma*sigma);
        }
        // calculate sigma for period rotation case
        sigma = (MU_infos[2*3+1] - MU_infos[2*3+0])/( (int)MU_infos[2*3+2] - 1);
        sigma *= s;
        expPart[2] = -1./(2.*sigma*sigma);

    #ifdef DEBUG
        std::cout<<"sigma = "<<sigma<<std::endl;
    #endif
        
    #pragma omp parallel for collapse(2)
        for (int k = 0; k < K; k++)
            for (int m = 0; m < M; m++){   // for each PHI_km
                if(m+1 <= Mnl){
                    PHI[k*M+m] = 0.;
                    // for shift
                    for(int l = 0;l < L-1;l++)
                        PHI[k*M+m] += expPart[l]*(X[k*L+l]-MU[m*L+l])*(X[k*L+l]-MU[m*L+l]);
                    // for period rotation
                    // temp = 2*pi - fabs(X[k*L+2]-MU[m*L+2]);
                    // temp = fabs(X[k*L+2]-MU[m*L+2]) < temp?fabs(X[k*L+2]-MU[m*L+2]):temp;

                    //for unperiod rotation
                    temp = fabs(X[k*L+2]-MU[m*L+2]);


                    PHI[k*M+m] += expPart[2]*temp*temp;

                    PHI[k*M+m] = exp(PHI[k*M+m]);
                }
                // else if(m+1 <= Mnl + L)
                // PHI[k*M+m] = X[k*L + m-Mnl];
                else
                    PHI[k*M+m] = 1;
        }

        alignedDoubleFree(expPart);

    #ifdef DEBUG
        std::cout<<"initialize PHI!"<<std::endl;
        prtMatrix(PHI,K,M);
        prtMatrix(PHI,100,M);
    #endif
    }


    // initialize PHI[K][M],M = Mnl+L+1,notice the MU_infos last three element is rotation info!
    void set_oneshift_rotation_basis(const double *X,int K,int L,
                                     const double *MU,const double *MU_infos,int Mnl,
                                     double *PHI,int M)
    {
        static double s = 1;
        double sigma;
        double temp;

        double *expPart = alignedDoubleMalloc(L,64);

        // calculate sigma for shift case

        // sigma = (infos[l*3+1] - infos[l*3+0])/((1 << (int)infos[l*3+2]) - 1);
        sigma = (MU_infos[0*3+1] - MU_infos[0*3+0])/( (int)MU_infos[0*3+2] - 1);
        sigma *= s;
        expPart[0] = -1./(2.*sigma*sigma);

        // calculate sigma for period rotation case
        sigma = (MU_infos[1*3+1] - MU_infos[1*3+0])/( (int)MU_infos[1*3+2] - 1);
        sigma *= s;
        expPart[1] = -1./(2.*sigma*sigma);

    #ifdef DEBUG
        std::cout<<"sigma = "<<sigma<<std::endl;
    #endif

    #pragma omp parallel for collapse(2)
        for (int k = 0; k < K; k++)
            for (int m = 0; m < M; m++){   // for each PHI_km
                if(m+1 <= Mnl){
                    PHI[k*M+m] = 0.;
                    // for shift
                    PHI[k*M+m] += expPart[0]*(X[k*L+0]-MU[m*L+0])*(X[k*L+0]-MU[m*L+0]);

                    // for period rotation
                    temp = 2*PI - fabs(X[k*L+1]-MU[m*L+1]);
                    temp = fabs(X[k*L+1]-MU[m*L+1]) < temp?fabs(X[k*L+1]-MU[m*L+1]):temp;

                    // for unperiod rotation
                    // temp = fabs(X[k*L+1]-MU[m*L+1]);


                    PHI[k*M+m] += expPart[2]*temp*temp;

                    PHI[k*M+m] = exp(PHI[k*M+m]);
                }
                // else if(m+1 <= Mnl + L)
                // PHI[k*M+m] = X[k*L + m-Mnl];
                else
                    PHI[k*M+m] = 1;
        }

        alignedDoubleFree(expPart);

    #ifdef DEBUG
        std::cout<<"initialize PHI!"<<std::endl;
        prtMatrix(PHI,K,M);
        prtMatrix(PHI,100,M);
    #endif
    }


    void set_2agnle_healpix_basis(const double *X,int K,int L,
                                  const double *MU,double sigma1,double sigma2,int Mnl,
                                  double *PHI,int M)
    {

        static double s = 1;
        if(L != 2){
            std::cout<<"L is not 2."<<std::endl;
            EXIT_ABNORMALLY;
        }

        double *expPart = alignedDoubleMalloc(L,64);

        // calculate sigma for basic_sampling case
        sigma1 *= s;
        expPart[0] = -1./(2.*sigma1*sigma1);
        sigma2 *= s;
        expPart[1] = -1./(2.*sigma2*sigma2);

    #pragma omp parallel for collapse(2)
        for (int k = 0; k < K; k++)
            for (int m = 0; m < M; m++){   // for each PHI_km
                if(m+1 <= Mnl){
                    PHI[k*M+m] = 0.;
                    for(int l = 0;l < L;l++)
                        PHI[k*M+m] += expPart[l]*(X[k*L+l]-MU[m*L+l])*(X[k*L+l]-MU[m*L+l]);
                    PHI[k*M+m] = exp(PHI[k*M+m]);
                }
                // else if(m+1 <= Mnl + L)
                // PHI[k*M+m] = X[k*L + m-Mnl];
                else
                    PHI[k*M+m] = 1;
        }

        alignedDoubleFree(expPart);

    #ifdef DEBUG
        std::cout<<"initialize PHI!"<<std::endl;
        prtMatrix(PHI,K,M);
    #endif
    }

    void set_2angle_basis(const double *X,int K,int L,
                          const double *MU,const double *MU_infos,int Mnl,
                          double *PHI,int M,bool period)
    {
        static double s = 1;
        double sigma;
        double temp;
        
        double *expPart = alignedDoubleMalloc(L,64);

        // calculate sigma for basic_sampling case
        for(int l = 0;l < L;l++){
            // sigma = (infos[l*3+1] - infos[l*3+0])/((1 << (int)infos[l*3+2]) - 1);
            sigma = (MU_infos[l*3+1] - MU_infos[l*3+0])/( (int)MU_infos[l*3+2] - 1);
            sigma *= s;
            expPart[l] = -1./(2.*sigma*sigma);
        }

    #ifdef DEBUG
        std::cout<<"sigma = "<<sigma<<std::endl;
    #endif

    #pragma omp parallel for collapse(2)
        for (int k = 0; k < K; k++)
            for (int m = 0; m < M; m++){   // for each PHI_km
                if(m+1 <= Mnl){
                    PHI[k*M+m] = 0.;
                    for(int l = 0;l < L;l++){

                        // for period rotation
                        if(period){
                            temp = 2*PI - fabs(X[k*L+l]-MU[m*L+l]);
                            temp = fabs(X[k*L+l]-MU[m*L+l]) < temp?fabs(X[k*L+l]-MU[m*L+l]):temp;
                        }
                        else
                            temp = fabs(X[k*L+l]-MU[m*L+l]); 	// for unperiod rotation

                        PHI[k*M+m] += expPart[l]*temp*temp;
                    }
                        
                    PHI[k*M+m] = exp(PHI[k*M+m]);
                }
                // else if(m+1 <= Mnl + L)
                // PHI[k*M+m] = X[k*L + m-Mnl];
                else
                    PHI[k*M+m] = 1;
        }

        alignedDoubleFree(expPart);

    #ifdef DEBUG
        std::cout<<"initialize PHI!"<<std::endl;
        prtMatrix(PHI,K,M);
    #endif
    }


    void set_3angle_basis(const double *X,int K,int L,
                          const double *MU,const double *MU_infos,int Mnl,
                          double *PHI,int M)
    {
        static double s = 1;
        double sigma;

        double *expPart = alignedDoubleMalloc(L,64);

        // calculate sigma for basic_sampling case
        for(int l = 0;l < L;l++){
            // sigma = (infos[l*3+1] - infos[l*3+0])/((1 << (int)infos[l*3+2]) - 1);
            sigma = (MU_infos[l*3+1] - MU_infos[l*3+0])/( (int)MU_infos[l*3+2] - 1);
            sigma *= s;
            expPart[l] = -1./(2.*sigma*sigma);
        }

    #ifdef DEBUG
        std::cout<<"sigma = "<<sigma<<std::endl;
    #endif

    #pragma omp parallel for collapse(2)
        for (int k = 0; k < K; k++)
            for (int m = 0; m < M; m++){   // for each PHI_km
                if(m+1 <= Mnl){
                    PHI[k*M+m] = 0.;
                    for(int l = 0;l < L;l++)
                        PHI[k*M+m] += expPart[l]*(X[k*L+l]-MU[m*L+l])*(X[k*L+l]-MU[m*L+l]);
                    PHI[k*M+m] = exp(PHI[k*M+m]);
                }
                // else if(m+1 <= Mnl + L)
                // PHI[k*M+m] = X[k*L + m-Mnl];
                else
                    PHI[k*M+m] = 1;
        }

        alignedDoubleFree(expPart);

    #ifdef DEBUG
        std::cout<<"initialize PHI!"<<std::endl;
        prtMatrix(PHI,K,M);
    #endif
    }

}