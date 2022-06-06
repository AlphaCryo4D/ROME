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

// #include "util.h"		// used for building precompiled headers on Windows

#include "spider_dep.h"

void spider_fq_q(float *data,int NX,int NY,double filter)
{
//    int NX = 80;
//    int NY = 80;
    int N2X = NX*2;
    int N2Y = NY*2;
    
    int LSD = N2X + 2;
    
    bool spider_sign = true;
    
    fftwf_plan forward_plan,backward_plan;
    
    float *B = (float*)aMalloc(sizeof(float)*N2Y*N2X,64);
    memset(B, 0, sizeof(float)*N2Y*N2X);
    //read image in B(NX*NY)
    for (int j = 0; j < NY; j++) {
        for (int i = 0; i < NX; i++) {
            B[j*N2X+i] = data[j*NX+i];
        }
    }
    double AVE = 0;
    for (int i = 0; i < NX; i++) {
        AVE = AVE + B[i*N2Y+0] + B[i*N2Y+NY-1];
    }
    for (int j = 1; j < NY-1; j++) {
        AVE = AVE + B[0*N2Y+j] + B[(NX-1)*N2Y+j];
    }
    AVE /= (2*(NX+NY)-4);
    
//    std::cout<<AVE<<std::endl;
    for (int j = 0; j < NY; j++)
        for (int i = NX; i < N2X; i++)
            B[j*N2X+i] = float(AVE);
    
    for (int j = NX; j < N2X; j++)
        for (int i = 0; i < N2Y; i++)
            B[j*N2X+i] = float(AVE);
    
    fftwf_complex* Fdata = (fftwf_complex*) fftwf_malloc(sizeof(fftwf_complex)*N2Y*(N2X/2+1));
    forward_plan = fftwf_plan_dft_r2c_2d(N2Y,N2X,B,Fdata,FFTW_ESTIMATE);
    fftwf_execute(forward_plan);
    fftwf_destroy_plan(forward_plan);
    
    if (spider_sign) {
        for (int i = 0; i < N2Y*(N2X/2+1); i++) {
            Fdata[i][1] = -1*Fdata[i][1];
        }
    }
    
//    for (int j = 0; j < N2Y; j++) {
//        std::cout<<(j+1)<<std::endl;
//        for (int i = 0; i < (N2X/2+1); i++) {
//            std::cout<<Fdata[j*(N2X/2+1)+i][0]<<" "<<Fdata[j*(N2X/2+1)+i][1]<<" ";
//        }
//        std::cout<<std::endl;
//    }
    

    double PARM1 = filter;
    double PARM2 = PARM1;
    if (filter < 0.0 || filter > 1.0) PARM1 = 0.5*PARM1/(NX/2);

    
    int NR2    = N2Y / 2;
    double X1    = (N2X/2)*(N2X/2);
    double Y1     = NR2*NR2;
    double PARM   = PARM1*PARM1;
    double PARM22 = PARM2*PARM2;
    
    for (int j = 1; j <= N2Y; j++) {
        
        int IY = j-1;
        if(IY > NR2) IY = IY - N2Y;
        
        for (int i = 1; i <= (N2X/2+1); i++) {
            
            int IX = i-1;
            if(0.25*((IX*IX)/X1/PARM + (IY*IY)/Y1/PARM22) > 1.0){
                Fdata[(j-1)*(N2X/2+1)+i-1][0] = 0.0;
                Fdata[(j-1)*(N2X/2+1)+i-1][1] = 0.0;
            }
            
        }
    }
    
    if (spider_sign) {
        for (int i = 0; i < N2Y*(N2X/2+1); i++) {
            Fdata[i][1] = -1*Fdata[i][1];
        }
    }
    
    backward_plan = fftwf_plan_dft_c2r_2d(N2Y,N2X,Fdata,B,FFTW_ESTIMATE);
    fftwf_execute(backward_plan);
    fftwf_destroy_plan(backward_plan);
    
    for (int j = 0; j < NY; j++) {
        for (int i = 0; i < NX; i++) {
            data[j*NX+i] = B[j*N2X+i];
        }
    }
    
    aFree(B);
    fftwf_free(Fdata);
    
}

float spider_fbs2(float x,float y,int nx,int ny,const float* data,const float* X1,const float* Y1,const float* XY2)
{
    int i,j,i2,j2,i3,j3;
    float A0,A1,A2,A3,ADX,BDX,DADX,DBDX;
    
    i = int(floor(x));
    j = int(floor(y));
    
    float dx    = x - i;
    float dy    = y - j;
    float dxsqr = dx*dx;
    float dxcub = dx*dx*dx;
    
    bool chkbound = true;
    if (chkbound) {
        i2 = modulo(i-1,nx);
        j2 = modulo(j-1,ny);
        i3 = modulo(i,nx);
        j3 = modulo(j,ny);
    }
    else{
        i2 = i - 1;
        j2 = j - 1;
        i3 = i;
        j3 = j;
    }
    
    A0   = data[j2*nx+i2];
    A1 = X1[j2*nx+i2];
    A2   = 3*( data[j2*nx+i3]-A0) -2*A1 - X1[j2*nx+i3];
    A3   = 2*(-data[j2*nx+i3]+A0)  + A1 + X1[j2*nx+i3];
    ADX  = A0 + A1*dx + A2*dxsqr + A3*dxcub;
    
    A0   = data[j3*nx+i2];
    A1   = X1[j3*nx+i2];
    A2   = 3*( data[j3*nx+i3]-A0) - 2*A1 - X1[j3*nx+i3];
    A3   = 2*(-data[j3*nx+i3]+A0) +   A1 + X1[j3*nx+i3];
    BDX  = A0 + A1*dx + A2*dxsqr + A3*dxcub;
    
    A0   = Y1[j2*nx+i2];
    A1   = XY2[j2*nx+i2];
    A2   = 3*( Y1[j2*nx+i3]-A0) -2*A1 - XY2[j2*nx+i3];
    A3   = 2*(-Y1[j2*nx+i3]+A0) +  A1 + XY2[j2*nx+i3];
    DADX = A0 + A1*dx + A2*dxsqr + A3*dxcub;
    
    A0   = Y1[j3*nx+i2];
    A1   = XY2[j3*nx+i2];
    A2   = 3*( Y1[j3*nx+i3]-A0) - 2*A1 - XY2[j3*nx+i3];
    A3   = 2*(-Y1[j3*nx+i3]+A0) +   A1 + XY2[j3*nx+i3];
    DBDX = A0 + A1*dx + A2*dxsqr + A3*dxcub;
    
    A2   = 3*(BDX - ADX) - 2*DADX - DBDX;
    A3   = 2*(ADX - BDX) +   DADX + DBDX;
    
    float FBS2 = ADX + DADX * dy + A2 * dy*dy + A3 * dy*dy*dy;
    
    //    std::cout<<ADX<<" "<<BDX<<" "<<DADX<<" "<<DBDX<<" "<<FBS2<<" ";
    //    exit(1);
    return FBS2;
}

void spider_rtsf(const float* data,float *out_data,int nx,int ny,float angle,float shiftx,float shifty)
{
    bool spider_sign = true;
    bool spider_scale = true;
    float scale = 1.;
    float spider_PI2 = 6.28318530717958647692f;
    float spider_PI  = 3.14159265358979323846f;
    
    /******  FBS2_PREP  ******/
    
    fftwf_plan forward_plan,backward_plan;
    
    float *data2 = (float*)aMalloc(sizeof(float)*ny*nx,64);
    memcpy(data2, data, sizeof(float)*ny*nx);
    fftwf_complex* Fdata = (fftwf_complex*) fftwf_malloc(sizeof(fftwf_complex)*ny*(nx/2+1));
    fftwf_complex* Ftemp = (fftwf_complex*) fftwf_malloc(sizeof(fftwf_complex)*ny*(nx/2+1));
    float *WX = (float*)aMalloc(sizeof(float)*(nx/2+1),64);
    float *WY = (float*)aMalloc(sizeof(float)*ny,64);
    float* X1 = (float*)aMalloc(sizeof(float)*ny*nx,64);
    float *Y1 = (float*)aMalloc(sizeof(float)*ny*nx,64);
    float *XY2 = (float*)aMalloc(sizeof(float)*ny*nx,64);
    
    forward_plan = fftwf_plan_dft_r2c_2d(ny,nx,data2,Fdata,FFTW_ESTIMATE);
    fftwf_execute(forward_plan);
    fftwf_destroy_plan(forward_plan);
    
    if (spider_sign) {
        for (int i = 0; i < ny*(nx/2+1); i++) {
            Fdata[i][1] = -1*Fdata[i][1];
        }
    }
    
    float A4 = spider_PI2/nx;
    for (int i = 0; i < (nx/2+1); i++) {
        WX[i] = i*A4;
    }

    A4 = spider_PI2/ny;
    for (int i = 0; i < ny/2+1; i++) {
        WY[i] = i*A4;
    }
    for (int i = ny/2+1; i < ny; i++) {
        WY[i] = (i-ny)*A4;
    }
    
    
    for (int j = 0; j < ny; j++) {
        for (int i = 0; i < nx/2+1; i++) {
            Ftemp[j*(nx/2+1)+i][0] = Fdata[j*(nx/2+1)+i][1]*WX[i];
            Ftemp[j*(nx/2+1)+i][1] = -Fdata[j*(nx/2+1)+i][0]*WX[i];
        }
    }
    if (spider_sign) {
        for (int i = 0; i < ny*(nx/2+1); i++) {
            Ftemp[i][1] = -1*Ftemp[i][1];
        }
    }
    
    backward_plan = fftwf_plan_dft_c2r_2d(ny,nx,Ftemp,X1,FFTW_ESTIMATE);
    fftwf_execute(backward_plan);
    fftwf_destroy_plan(backward_plan);
    
    if (spider_scale) {
        float pix = 1.0f/(ny*nx);
        for (int i = 0; i < ny*nx; i++) {
            X1[i] *= pix;
        }
    }
    
    
    for (int j = 0; j < ny; j++) {
        for (int i = 0; i < nx/2+1; i++) {
            Ftemp[j*(nx/2+1)+i][0] = Fdata[j*(nx/2+1)+i][1]*WY[j];
            Ftemp[j*(nx/2+1)+i][1] = -Fdata[j*(nx/2+1)+i][0]*WY[j];
        }
    }
    
    if (spider_sign) {
        for (int i = 0; i < ny*(nx/2+1); i++) {
            Ftemp[i][1] = -1*Ftemp[i][1];
        }
    }
    
    backward_plan = fftwf_plan_dft_c2r_2d(ny,nx,Ftemp,Y1,FFTW_ESTIMATE);
    fftwf_execute(backward_plan);
    fftwf_destroy_plan(backward_plan);
    
    if (spider_scale) {
        float pix = 1.0f/(ny*nx);
        for (int i = 0; i < ny*nx; i++) {
            Y1[i] *= pix;
        }
    }
    
    
    for (int j = 0; j < ny; j++) {
        for (int i = 0; i < nx/2+1; i++) {
            Ftemp[j*(nx/2+1)+i][0] = Fdata[j*(nx/2+1)+i][1]*WY[j];
            Ftemp[j*(nx/2+1)+i][1] = -Fdata[j*(nx/2+1)+i][0]*WY[j];
        }
    }
    
    if (spider_sign) {
        for (int i = 0; i < ny*(nx/2+1); i++) {
            Ftemp[i][1] = -1*Ftemp[i][1];
        }
    }
    
    backward_plan = fftwf_plan_dft_c2r_2d(ny,nx,Ftemp,XY2,FFTW_ESTIMATE);
    fftwf_execute(backward_plan);
    fftwf_destroy_plan(backward_plan);
    
    if (spider_scale) {
        float pix = 1.0f/(ny*nx);
        for (int i = 0; i < ny*nx; i++) {
            XY2[i] *= pix;
        }
    }
    
    fftwf_free(Fdata);
    fftwf_free(Ftemp);
    aFree(data2);
    aFree(WX);
    aFree(WY);
    /******   RTSF ********/
    
    float theta = angle * spider_PI / 180;
    float costh = cos(theta);
    float sinth = sin(theta);
    
    int cx    = nx / 2 + 1;
    int cy    = ny / 2 + 1;
    
    float shx = modulo(shiftx,(float)nx);
    float shy = modulo(shifty,(float)ny);
//    std::cout<<shx<<" "<<shiftx<<" "<<shy<<" "<<shifty<<std::endl;
    
    if (scale == 1.) {
        
        float fy0    = - shy - cy;
        float fy1    = - shy + ny - cy;
        float fy2    = - shy - ny - cy;
        
        float fx0    = - shx - cx;
        float fx1    = - shx + nx - cx;
        float fx2    = - shx - nx - cx;
        
        float shypny = shy + ny;
        float shxpnx = shx + nx;
        
        for (int iy = 1; iy <= ny; iy++) {
            float fy = iy + fy0;
            if((iy-1) < shy) fy = iy + fy1;
            if((iy-1) >= shypny) fy = iy + fy2;
            float ycod =  costh*fy+cy;
            float ysid = -sinth*fy+cx;
            
            for (int ix = 1; ix <= nx; ix++) {
                float fx = ix+fx0;
                if ((ix-1) < shx) fx = ix + fx1;
                if ((ix-1) >= shxpnx) fx = ix + fx2;
                
                float xold = costh*fx + ysid;
                float yold = sinth*fx + ycod;
                
                out_data[(iy-1)*nx+(ix-1)] = spider_fbs2(xold,yold,nx,ny,data,X1,Y1,XY2);
            }
        }
    }
    
    aFree(X1);
    aFree(Y1);
    aFree(XY2);
}

