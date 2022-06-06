/***************************************************************************
 *
 * Intel® Parallel Computing Center for Structural Biology
 * Principal Investigator : Youdong (Jack) Mao (Youdong_Mao@dfci.harvard.edu)
 * Dana-Farber Cancer Institute, Harvard Medical School and Peking University
 *
 * Authors: "Yong Bei Ma(galowma@gmail.com)"
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
/***************************************************************************
 *
 * Author: "Sjors H.W. Scheres"
 * MRC Laboratory of Molecular Biology
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

#ifndef RECONSTRUCT_H_
#define RECONSTRUCT_H_

#include "./resmap_image.h"
#include "./resmap_mrcs.h"
#include "./resmap_math.h"
#include "./resmap_time.h"
#include "./resmap_array.h"
#include "./resmap_symmetries.h"
#include "./resmap_mpi.h"

#ifndef NEAREST_NEIGHBOUR
#define NEAREST_NEIGHBOUR 0
#endif
#ifndef TRILINEAR
#define TRILINEAR 1
#endif

//#define VECTORIZED_RECONS
#define RADIUS_LE(radius2,kr,ir) (radius2-kr*kr-ir*ir) < 0 ? 0 : floor(sqrt(double(radius2-kr*kr-ir*ir))) + 1
#define RADIUS_L(radius2,kr,ir)  (radius2-kr*kr-ir*ir) < 0 ? 0 : ceil (sqrt(double(radius2-kr*kr-ir*ir)));

#define XMIPP_ORIGIN(x) -(int)((float) (x) / 2.0)


// - - -- - - - -- - - - -- - -debug - - -- - - -
//#define PRINT_RESULT
inline void printSum(FDOUBLE* A,int L,char* what){
#ifdef PRINT_RESULT
    double sum = 0;
    std::cout.precision(50);
    for(int i = 0;i < L;i++) sum += A[i];
    std::cout<<what<<"1~3 elements : "<<A[0]<<" "<<A[1]<<" "<<A[2]<<",sum = "<<sum<<std::endl;
#endif
};
inline void printSumDouble(double* A,int L,char* what){
#ifdef PRINT_RESULT
    double sum = 0;
    std::cout.precision(50);
    for(int i = 0;i < L;i++) sum += A[i];
    std::cout<<what<<"1~3 elements : "<<A[0]<<" "<<A[1]<<" "<<A[2]<<",sum = "<<sum<<std::endl;
#endif
};
inline void printSum_real(MKL_Complex* A,int L,char* what){
#ifdef PRINT_RESULT
    double sum = 0;
    std::cout.precision(50);
    for(int i = 0;i < L;i++) sum += A[i].real;
    std::cout<<what<<"1~3 elements : "<<A[0].real<<" "<<A[1].real<<" "<<A[2].real<<",sum = "<<sum<<std::endl;
#endif
};
inline void printSum_imag(MKL_Complex* A,int L,char* what){
#ifdef PRINT_RESULT
    double sum = 0;
    std::cout.precision(50);
    for(int i = 0;i < L;i++) sum += A[i].imag;
    std::cout<<what<<"1~3 elements : "<<A[0].imag<<" "<<A[1].imag<<" "<<A[2].imag<<",sum = "<<sum<<std::endl;
#endif
};
// - - -- - - - -


template<typename T>
void centerFFT3DBackword(const T* v_in, T* v_out, int dim)
{
    int dim2 = dim*dim;
    // 3D
    int l, shift;
    
    // Shift in the X direction
    l = dim;
    shift = (int)(l / 2);
    
    bool forward = false;
    if (!forward)
        shift = -shift;
    
    for (int k = 0; k < dim; k++)
    {
        int kp = k + shift;
        if (kp < 0) kp += l;
        else if (kp >= l) kp -= l;
        
        for (int i = 0; i < dim; i++)
        {
            int ip = i + shift;
            if (ip < 0) ip += l;
            else if (ip >= l) ip -= l;
            
            // Shift the input in an auxiliar vector
            for (int j = 0; j < l; j++)
            {
                int jp = j + shift;
                if (jp < 0) jp += l;
                else if (jp >= l) jp -= l;
                
                v_out[kp*dim2+ip*dim+jp] = v_in[k*dim2+i*dim+j];
                
            }// end loop j
        }// end loop i
    } // end loop k
}

template<typename T>
void centerFFT3DForward(const T* v_in, T* v_out, int dim)
{
    int dim2 = dim*dim;
    // 3D
    int l, shift;
    
    // Shift in the X direction
    l = dim;
    assert(l%2==0);
    shift = (int)(l / 2);
    
    // bool forward = true;
    // if (!forward)
    //     shift = -shift;
    assert(shift>0);
    
    for (int k = 0; k < dim; k++)
    {
        int kp = k + shift;
        if (kp < 0) kp += l;
        else if (kp >= l) kp -= l;
        
        for (int i = 0; i < dim; i++)
        {
            int ip = i + shift;
            if (ip < 0) ip += l;
            else if (ip >= l) ip -= l;
            
            auto v_in_part = v_in+k*dim2+i*dim;
            // Shift the input in an auxiliar vector
            auto v_out_part = v_out+kp*dim2+ip*dim+shift;
            #pragma ivdep
            for (int j = 0; j < l/2; j++)
            {
                v_out_part[j] = v_in_part[j];
            }// end loop j
            v_out_part = v_out+kp*dim2+ip*dim+shift-l;
            #pragma ivdep
            for (int j = l/2; j < l; j++)
            {
                v_out_part[j] = v_in_part[j];
            }// end loop j
        }// end loop i
    } // end loop k
}

template <typename T>
void centerFFT2D(const T* v_in, T* v_out, int dim,bool forward)
{
    for (int i = 0; i < dim*dim; i++)
        v_out[i] = v_in[i];
    
    // 2D
    // Shift in the X,Y direction
    int l, shift;
    l = dim;
    
    shift = (int)(l / 2);
    
    //
    if (!forward)
        shift = -shift;
    
    for (int i = 0; i < dim; i++)
    {
        int ip = i + shift;
        
        if (ip < 0)
            ip += l;
        else if (ip >= l)
            ip -= l;
        
        // Shift the input in an auxiliar vector
        for (int j = 0; j < l; j++)
        {
            int jp = j + shift;
            
            if (jp < 0)
                jp += l;
            else if (jp >= l)
                jp -= l;
            
            v_out[ip*dim+jp] = v_in[i*dim+j];
        }
    }
}

template<typename T>
void window(Vol<T>& vol_in,Vol<T>& vol_out,T init_value = 0)
{
    // TODO : Check size of the input array
    // vol_in.dimy > vol_out.dimy
    int offset_z = XMIPP_ORIGIN(vol_out.dimz) - XMIPP_ORIGIN(vol_in.dimz);
    int offset_y = XMIPP_ORIGIN(vol_out.dimy) - XMIPP_ORIGIN(vol_in.dimy);
    int offset_x = XMIPP_ORIGIN(vol_out.dimx) - XMIPP_ORIGIN(vol_in.dimx);
    for (int k = 0; k < vol_out.dimz; k++) {
        int k2 = k + offset_z;
        for (int i = 0; i < vol_out.dimy; i++) {
            int i2 = i + offset_y;
            for (int j = 0; j < vol_out.dimx; j++) {
                int j2 = j + offset_x;
                if (k2 > -1 && k2 < vol_in.dimz &&
                    i2 > -1 && i2 < vol_in.dimy &&
                    j2 > -1 && j2 < vol_in.dimx)
                    ACCESS(vol_out, k, i, j) = ACCESS(vol_in, k2, i2, j2);
                else
                    ACCESS(vol_out, k, i, j) = init_value;
            }
        }
    }
}

// Mask out corners outside sphere (replace by average value)
// Apply a soft mask (raised cosine with cosine_width pixels width)
template<typename T>
void softMaskOutsideMap(Vol<T> &vol, double radius = -1., double cosine_width = 3, Vol<T> *Mnoise = nullptr)
{
    // TODO check size
    // TODO Mnoise equal to vol
    
    double r_2, radius_p, radius_2, raisedcos, radius_p_2, sum_bg = 0., sum = 0.;
    if (radius < 0)
        radius = (double)vol.dimy/2.;
    radius_p = radius + cosine_width;
    
    radius_2 = radius*radius;
    radius_p_2 = radius_p*radius_p;
    
    int vol_origin_offset_z = XMIPP_ORIGIN(vol.dimz);
    int vol_origin_offset_y = XMIPP_ORIGIN(vol.dimy);
    int vol_origin_offset_x = XMIPP_ORIGIN(vol.dimx);
    
    if (Mnoise == nullptr)
    {
        // Calculate average background value
        for (int k = 0; k < vol.dimz; k++){
            int kp = k + vol_origin_offset_z;
            for (int i = 0; i < vol.dimy; i++){
                int ip = i + vol_origin_offset_y;
                int kp_ip_2 = kp*kp + ip*ip;
                for (int j = 0; j < vol.dimx; j++)
                {
                    int jp = j + vol_origin_offset_x;
                    // NOTE : only calculate pow(r,2) then compare with radius_2
                    // r = sqrt((double)(kp_ip_2 + jp*jp));
                    int r_2 = kp_ip_2 + jp*jp;
                    if (r_2 < radius_2)
                        continue;
                    else if (r_2 > radius_p_2)
                    {
                        sum    += 1.;
                        sum_bg += ACCESS(vol, k, i, j);
                    }
                    else
                    {
                        raisedcos = 0.5 + 0.5 * cos(rome_pi * (radius_p - sqrt(r_2)) / cosine_width );
                        sum += raisedcos;
                        sum_bg += raisedcos * ACCESS(vol, k, i, j);
                    }
                }
            }
        }
        sum_bg /= sum;
    }
    
    // Apply noisy or average background value
    for (int k = 0; k < vol.dimz; k++){
        int kp = k + vol_origin_offset_z;
        for (int i = 0; i < vol.dimy; i++){
            int ip = i + vol_origin_offset_y;
            int kp_ip_2 = kp*kp + ip*ip;
            for (int j = 0; j < vol.dimx; j++)
            {
                int jp = j + vol_origin_offset_x;
                // r = sqrt((double)(kp_ip_2 + jp*jp));
                r_2 = kp_ip_2 + jp*jp;
                if (r_2 < radius_2)
                {
                    continue;
                }
                else if (r_2 > radius_p_2)
                {
                    vol(k, i, j) = (Mnoise == nullptr) ? sum_bg : (*Mnoise)(k, i, j);
                }
                else
                {
                    raisedcos = 0.5 + 0.5 * cos(rome_pi* (radius_p - sqrt(r_2)) / cosine_width );
                    double add = (Mnoise == nullptr) ?  sum_bg : (*Mnoise)(k, i, j);
                    vol(k, i, j) = (1 - raisedcos) * vol(k, i, j) + raisedcos * add;
                }
            }
        }
    }
}

// Go from the Projector-centered fourier transform back to FFTW-uncentered one
template <typename T1,typename T2>
void decenter(Vol<T1> &Fin, Vol<T2> &Fout, int my_rmax2)
{
    // TODO check size
    // TOOD 3D case slower than original version
    // NOTE : Fin size should larger than Fout....
    assert(Fout.dimy <= Fin.dimy);
    int Fin_origin_offset_z = XMIPP_ORIGIN(Fin.dimz);
    int Fin_origin_offset_y = XMIPP_ORIGIN(Fin.dimy);
    // Mout should already have the right size
    for (int k = 0; k < Fout.dimz; k++){
        for (int i = 0;i < Fout.dimy; i++){
            int kp = k < Fout.dimx?k:(k-Fout.dimz);
            int ip = i < Fout.dimx?i:(i-Fout.dimy);
            // set xmipp origin
            int k2 = kp - Fin_origin_offset_z;
            int i2 = ip - Fin_origin_offset_y;
#ifdef VECTORIZED_RECONS // vectorized
            int j_max = RADIUS_LE(my_rmax2,kp,ip);
            j_max = (Fout.dimx)<j_max?Fout.dimx:j_max;
            // jp = j2 = j
            #pragma simd
            for (int j = 0; j < j_max; j++) {
                ACCESS(Fout, k, i, j) = ACCESS(Fin, k2, i2 ,j);
            }
            #pragma simd
            for (int j = j_max; j < Fout.dimx; j++) {
                ACCESS(Fout, k, i, j) = 0;
            }
#else
            for (int j = 0; j < Fout.dimx; j++)
            {
                int jp = j;
                if (kp*kp + ip*ip + jp*jp <= my_rmax2){
                    int j2 = jp;
                    assert(i2*Fin.dimx+j2 < Fin.dimzyx);
                    Fout(k, i, j) = Fin(k2, i2, j2);
                }
                else// Initialize to zero
                    Fout(k, i, j) = 0;
            }
#endif
        }
    }
}

template <>
void decenter(Vol<MKL_Complex16> &Fin, Vol<MKL_Complex16> &Fout, int my_rmax2)
{
    // TODO check size
    // TOOD 3D case slower than original version
    // NOTE : Fin size should larger than Fout....
    assert(Fout.dimy <= Fin.dimy);
    int Fin_origin_offset_z = XMIPP_ORIGIN(Fin.dimz);
    int Fin_origin_offset_y = XMIPP_ORIGIN(Fin.dimy);
    // Mout should already have the right size
    for (int k = 0; k < Fout.dimz; k++){
        for (int i = 0;i < Fout.dimy; i++){
            int kp = k < Fout.dimx?k:(k-Fout.dimz);
            int ip = i < Fout.dimx?i:(i-Fout.dimy);
            // set xmipp origin
            int k2 = kp - Fin_origin_offset_z;
            int i2 = ip - Fin_origin_offset_y;
#ifdef VECTORIZED_RECONS // vectorized
            int j_max = RADIUS_LE(my_rmax2, kp, ip);
            j_max = j_max < Fout.dimx ? j_max : Fout.dimx;
            #pragma simd
            for (int j = 0; j < j_max; j++)
            {
                // jp = j
                ACCESS(Fout, k, i, j) = ACCESS(Fin,k2, i2, j);
            }
            #pragma simd
            for (int j = j_max; j < Fout.dimx; j++)
            {
                // jp = j
                ACCESS(Fout, k, i, j).real = Fout(k, i, j).imag = 0;
            }
#else
            for (int j = 0,jp = 0; j < Fout.dimx; j++,jp = j)
            {
                if (kp*kp + ip*ip + jp*jp <= my_rmax2){
                    int j2 = jp;
                    assert(i2*Fin.dimx+j2 < Fin.dimzyx);
                    Fout(k, i, j) = Fin(k2, i2, j2);
                }
                else// Initialize to zero
                    Fout(k, i, j).real = Fout(k, i, j).imag = 0;
            }
#endif
        }
    }
}

template <>
void decenter(Vol<MKL_Complex8> &Fin, Vol<MKL_Complex8> &Fout, int my_rmax2)
{
    // TODO check size
    // TOOD 3D case slower than original version
    // NOTE : Fin size should larger than Fout....
    assert(Fout.dimy <= Fin.dimy);
    int Fin_origin_offset_z = XMIPP_ORIGIN(Fin.dimz);
    int Fin_origin_offset_y = XMIPP_ORIGIN(Fin.dimy);
    // Mout should already have the right size
    for (int k = 0; k < Fout.dimz; k++){
        for (int i = 0;i < Fout.dimy; i++){
            int kp = k < Fout.dimx?k:(k-Fout.dimz);
            int ip = i < Fout.dimx?i:(i-Fout.dimy);
            // set xmipp origin
            int k2 = kp - Fin_origin_offset_z;
            int i2 = ip - Fin_origin_offset_y;
#ifdef VECTORIZED_RECONS // vectorized
            int j_max = RADIUS_LE(my_rmax2, kp, ip);
            j_max = j_max < Fout.dimx ? j_max : Fout.dimx;
#pragma simd
            for (int j = 0; j < j_max; j++)
            {
                // jp = j
                ACCESS(Fout, k, i, j) = ACCESS(Fin,k2, i2, j);
            }
#pragma simd
            for (int j = j_max; j < Fout.dimx; j++)
            {
                // jp = j
                ACCESS(Fout, k, i, j).real = Fout(k, i, j).imag = 0;
            }
#else
            for (int j = 0,jp = 0; j < Fout.dimx; j++,jp = j)
            {
                if (kp*kp + ip*ip + jp*jp <= my_rmax2){
                    int j2 = jp;
                    assert(i2*Fin.dimx+j2 < Fin.dimzyx);
                    Fout(k, i, j) = Fin(k2, i2, j2);
                }
                else// Initialize to zero
                    Fout(k, i, j).real = Fout(k, i, j).imag = 0;
            }
#endif
        }
    }
}

template<typename T>
void getSpectrum(Vol<T> &Min,FDOUBLE* spectrum,int spectrum_type)
{
    assert(Min.dimy == Min.dimx);
    Vol<FDOUBLE> Faux_real,Faux_imag;
    int xsize = Min.dimx;
    int Fxsize = Min.dimx/2+1;
    FourierTransformer *transformer;
    if (Min.dimz == 1) {
        transformer = FourierTransformer::make(Min.dimy,Min.dimy);
        Faux_real.init_nonzero(1, Min.dimy, Min.dimx/2+1);
        Faux_imag.init_nonzero(1, Min.dimy, Min.dimx/2+1);
    }
    else{
        assert(Min.dimz == Min.dimy);
        transformer = FourierTransformer::make(Min.dimy,Min.dimy,Min.dimy);
        Faux_real.init_nonzero(Min.dimz, Min.dimy, Min.dimx/2+1);
        Faux_imag.init_nonzero(Min.dimz, Min.dimy, Min.dimx/2+1);
    }
    
    std::vector<FDOUBLE> count(xsize,0);
    
    for (int i = 0; i < xsize; i++) spectrum[i] = 0;
    transformer->FourierTransform(Min.wptr(), Faux_real.wptr(), Faux_imag.wptr());
    
    if (spectrum_type == AMPLITUDE_SPECTRUM)
    {
        // FOR_ALL_ELEMENTS_IN_FFTW_TRANSFORM(Faux)
        for (int k = 0; k < Min.dimz; k++){
            int kp = (k < Fxsize) ? k : k - Min.dimz;
            for (int i = 0; i < Min.dimy; i++){
                int ip = (i < Fxsize) ? i : i - Min.dimy;
                int kp_ip_2 = kp*kp + ip*ip;
                for (int j = 0; j < Fxsize; j++)
                {
                    int jp = j;
                    int idx = round(sqrt(kp_ip_2 + jp*jp));
                    spectrum[idx] += sqrt(ACCESS(Faux_real, k, i, j)*ACCESS(Faux_real, k, i, j) + \
                                          ACCESS(Faux_imag, k, i, j)*ACCESS(Faux_imag, k, i, j));
                    count[idx] += 1.;
                }
            }
        }
    }
	else
    {
        for (int k = 0; k < Min.dimz; k++){
            int kp = (k < Fxsize) ? k : k - Min.dimz;
            for (int i = 0; i < Min.dimy; i++){
                int ip = (i < Fxsize) ? i : i - Min.dimy;
                int kp_ip_2 = kp*kp + ip*ip;
                for (int j = 0; j < Fxsize; j++)
                {
                    int jp = j;
                    int idx = round(sqrt(kp_ip_2 + jp*jp));
                    spectrum[idx] += ACCESS(Faux_real, k, i, j)*ACCESS(Faux_real, k, i, j) + \
                    					ACCESS(Faux_imag, k, i, j)*ACCESS(Faux_imag, k, i, j);
                    count[idx] += 1.;
                }
            }
        }
    }
    for (int i = 0; i < xsize; i++)
        if (count[i] > 0.)
            spectrum[i] /= count[i];
    sDelete(transformer);
    Faux_real.fini();Faux_imag.fini();
}

// Window an FFTW-centered Fourier-transform to a given size
template<typename T>
void windowFourierTransform(Vol<T>& Fin, Vol<T>& Fout)
{
    // TODO : Check size of the input array
    // if (YSIZE(in) > 1 && YSIZE(in)/2 + 1 != XSIZE(in))
    //     REPORT_ERROR("windowFourierTransform ERROR: the Fourier transform should be of an image with equal sizes in all dimensions!");
    // long int newhdim = newdim/2 + 1;
    
    // If same size, just return input
    //
    if (Fin.dimx == Fout.dimx)
    {
        if (Fin.dimy != Fout.dimy){
            std::cerr<<"warning : windowFourierTransform different dimension y."<<std::endl;
            // TODO
            assert(false);
			EXIT_ABNORMALLY;
        }
        for (int n = 0; n < Fin.dimzyx; n++) {
            Fout(0,0,n) = Fin(0,0,n);
        }
        return;
    }
    
    if (Fout.dimx > Fin.dimx)
    {
        // Otherwise apply a windowing operation
        // Initialise output array
        Fout.fill_zero();
        
        int max_r2 = (Fin.dimx -1) * (Fin.dimx - 1);
        // FOR_ALL_ELEMENTS_IN_FFTW_TRANSFORM(in)
        for (int k = 0; k<Fin.dimz; k++){
            for (int i = 0; i<Fin.dimy; i++){
                int kp = (k < Fin.dimx) ? k : k - Fin.dimz;
                int ip = (i < Fin.dimx) ? i : i - Fin.dimy;
                int k2 = (kp<0)? kp+Fout.dimz : kp;
                int i2 = (ip<0)? ip+Fout.dimy : ip;
                for (int j = 0, jp = 0; j<Fin.dimx; j++, jp = j)
                {
                    // Make sure windowed FT has nothing in the corners, otherwise we end up with an asymmetric FT!
                    if (kp*kp + ip*ip + jp*jp <= max_r2){
                        int j2 = jp;
                        Fout(k2, i2, j2) = Fin(k, i, j);
                    }
                }
        }}
    }
    else
    {
        // FOR_ALL_ELEMENTS_IN_FFTW_TRANSFORM(out)
        for (int k = 0; k<Fout.dimz; k++){
            for (int i = 0; i<Fout.dimy; i++){
                int kp = (k < Fout.dimx) ? k : k - Fout.dimz;
                int ip = (i < Fout.dimx) ? i : i - Fout.dimy;
                int k2 = (kp<0)? kp+Fin.dimz : kp;
                int i2 = (ip<0)? ip+Fin.dimy : ip;
                for (int j = 0, jp = 0; j<Fout.dimx; j++, jp = j)
                {
                    int j2 = jp;
                    Fout(k, i, j) = Fin(k2, i2, j2);
                }
        }}
    }
}

/** Linear interpolation
 *
 * From low (when a=0) to high (when a=1). The following value is returned
 * (equal to (a*h)+((1-a)*l)
 */
//template<typename T>
//    T lin_interp(T a, T l, T h){return l + (h - l) * (a);}

#define lin_interp(a,l,h) ((l) + ((h) - (l)) * (a))

class ProjectorBase{
public:
    // The Fourier-space image data array
    Vol<MKL_Complex> data;
    
    // Only points within this many pixels from the origin (in the original size) will be interpolated
    int r_max;
    
    // Radius of sphere within TRILINEAR interpolation will be used in NEAREST_NEIGHBOUR interpolator
    int r_min_nn;
    
    // Original size of the real-space map
    int ori_size,ori_Fsize;
    
    // Padded size of the map in Fourier-space
    int pad_size,pad_Fsize;
    
    // Interpolation scheme (TRILINEAR or NEAREST_NEIGHBOUR, for BackProjector also CONVOLUTE_BLOB)
    int interpolator;
    
    // Oversample FT by padding in real space
    int padding_factor;
    
    // Dimension of the reference (currently allowed 2 or 3)
    int ref_dim;
    ProjectorBase():padding_factor(2),interpolator(TRILINEAR),r_min_nn(10){}
    ~ProjectorBase(){finalize();}
    void initialize(int _ori_size,int _current_size,int _ref_dim,int _padding_factor = 2,int _interpolator = TRILINEAR,int _r_min_nn = 10);
    void finalize(){data.fini();}
    //
    void griddingCorrect(Vol<FDOUBLE>& vol);
};

class MyProjector : public ProjectorBase
{
public:
    
    MyProjector(){}
    ~MyProjector(){}
    
    // Fill data array with oversampled Fourier transform, and calculate its power spectrum
    void computeFourierTransformMap(Vol<FDOUBLE>& ref_in,int _ori_size,int _current_size,FDOUBLE* power_spectrum,bool do_gridding = true);
    
    //
    void project(FDOUBLE* f2d_real,FDOUBLE* f2d_imag,int f2d_size, const FDOUBLE A[][3],bool inv);
    //
    void projectOneTile(FDOUBLE* f2d_real,FDOUBLE* f2d_imag,int n_start,int n_end,int f2d_size,const FDOUBLE A[][3],bool inv);
    //
    void projectOneTileByShell(FDOUBLE* f2d_real,FDOUBLE* f2d_imag,int shell_n_start,int shell_n_end,
                               int f2d_size,const FDOUBLE A[][3],bool inv,const int* nIndex);
    //
    void rotate2D(FDOUBLE* f2d_real,FDOUBLE* f2d_imag,int f2d_size,const FDOUBLE A[][3],bool inv);
    //
};


class MyBackProjector : public ProjectorBase
{
    class TabFtBlob{
    private:
        double radius;
        double alpha;
        int order;
        double sampling;
        std::vector<double> tabulatedValues;
    public:
        // Empty constructor
        TabFtBlob() {}
        ~TabFtBlob(){tabulatedValues.resize(0);}
        // Constructor (with parameters)
        void init(double _radius, double _alpha, int _order, const int _nr_elem = 10000){
            radius = _radius;
            alpha = _alpha;
            order = _order;
            sampling = 0.5 / (double)_nr_elem;
            tabulatedValues.resize(_nr_elem);
            for (int i = 0; i < _nr_elem; i++)
            {
                double xx = (double) i * sampling;
                tabulatedValues[i] = kaiser_Fourier_value(xx, radius, alpha, order);
            }
        }
        // Value access
        double operator()(double val) const {
            int idx = (int)( fabs(val) / sampling);
            if (idx >= tabulatedValues.size())
                return 0.;
            else
                return tabulatedValues[idx];
        }
    };
    
public:
    // For backward projection: sum of weights
    Vol<FDOUBLE> weight;
    
    // Tabulated blob values
    TabFtBlob tab_ftblob;
    
    // Symmetry object
    std::string fn_sym;
    
    MyBackProjector(){}
    ~MyBackProjector(){finalize();}
    
    void initialize(int _ori_size,int _current_size,int _ref_dim,std::string _fn_sym = "C1",
                    int _padding_factor = 2,int _interpolator = TRILINEAR,int _r_min_nn = 10,
                    int _blob_order = 0, double _blob_radius = 1.9, double _blob_alpha = 15);
    
    void finalize(){ProjectorBase::finalize();weight.fini();}
    
    void backproject(const FDOUBLE* f2d_real,const FDOUBLE* f2d_imag,int f2d_size,const FDOUBLE A[][3], bool inv,
                     Vol<MKL_Complex>& data_td,Vol<FDOUBLE>& weight_td,
                     const FDOUBLE* Mweight = nullptr,int z_start = 0,int z_end = 100000);
    
    void backprojectOneTileByShell(const FDOUBLE* f2d_real,const FDOUBLE* f2d_imag,int shell_n_start,int shell_n_end,
                                   int f2d_size,const FDOUBLE A[][3],bool inv, Vol<MKL_Complex>& data_td,
                                   Vol<FDOUBLE>& weight_td,const FDOUBLE* Mweight,const int* nIndex);
    
    void backrotate2D(const FDOUBLE* f2d_real,const FDOUBLE* f2d_imag,int f2d_size,const FDOUBLE A[][3], bool inv,
                      Vol<MKL_Complex>& data_td,Vol<FDOUBLE>& weight_td,const FDOUBLE* Mweight = nullptr);
    
    // Enforce hermitian symmetry on data and on weight (all points in the x==0 plane)
    // Because the interpolations are numerical, hermitian symmetry may be broken.
    // Repairing it here gives like a 2-fold averaging correction for interpolation errors...
    void enforceHermitianSymmetry(Vol<MKL_Complex>& my_data,Vol<FDOUBLE>& my_weight);
    
    // Applies the symmetry from the SymList object to the weight and the data array
    void symmetrise(Vol<MKL_Complex>& my_data,Vol<FDOUBLE>& my_weight, int my_rmax2);
    
    // Convolute in Fourier-space with the blob by multiplication in real-space
    // Note the convlution is done on the complex array inside the transformer object!!
    void convoluteBlobRealSpace(Vol<MKL_Complex>& Fconv, FourierTransformerBase& transformer, bool do_mask = false);
    
    // Calculate the inverse FFT of Fin and windows the result to ori_size
    // Also pass the transformer, to prevent making and clearing a new one before clearing the one in reconstruct()
    void windowToOridimRealSpace(Vol<MKL_Complex> &Fconv, Vol<FDOUBLE> &Vol_out, FourierTransformerBase& transformer);
    
    //
    void reconstruct(Vol<FDOUBLE> &vol_out,
                     int max_iter_preweight,
                     bool do_map,
                     double tau2_fudge,
                     FDOUBLE* tau2,
                     FDOUBLE* sigma2,
                     FDOUBLE* data_vs_prior,
                     const FDOUBLE* fsc, // only input
                     double normalise = 1.,
                     bool update_tau2_with_fsc = false,
                     bool is_whole_instead_of_half = false,
                     int minres_map = -1,
                     int nr_threads = 1,
                     std::string tmp_folder = "NULL");

    // Get only the lowest resolution components from the data and weight array
    // (to be joined together for two independent halves in order to force convergence in the same orientation)
    void getLowResDataAndWeight(Vol<MKL_Complex > &lowres_data, Vol<FDOUBLE> &lowres_weight,int lowres_r_max);
    
    // Set only the lowest resolution components from the data and weight array
    // (to be joined together for two independent halves in order to force convergence in the same orientation)
    void setLowResDataAndWeight(Vol<MKL_Complex > &lowres_data, Vol<FDOUBLE> &lowres_weight,int lowres_r_max);
    
    //  Get complex array at the original size as the straightforward average
    //  padding_factor*padding_factor*padding_factor voxels
    //  This will then be used for FSC calculation between two random halves
    void getDownsampledAverage(Vol<MKL_Complex > &avg);
    
    // From two of the straightforward downsampled averages, calculate an FSC curve
    void calculateDownSampledFourierShellCorrelation(Vol<MKL_Complex > &avg1,
                                                     Vol<MKL_Complex > &avg2,
                                                     FDOUBLE* fsc);
};




#endif /* defined(RECONSTRUCT_H_) */
