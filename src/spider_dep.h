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

#ifndef SPIDER_OLD_H_
#define SPIDER_OLD_H_

#include "fftw/fftw3.h"

#include "./util.h"


template<typename T>
inline T modulo(T A,T P){
    return T((double)A - floor((double)A/P) * P);
}

void spider_fq_q(float *data,int NX,int NY,double filter);


float spider_fbs2(float x,float y,int nx,int ny,const float* data,const float* X1,const float* Y1,const float* XY2);

//anticlockwise rotation(positive),shift from left(negative) to righ(positive)t,top(negative) to bottom(positive)
void spider_rtsf(const float* data,float *out_data,int nx,int ny,float angle,float shiftx,float shifty);

#endif 