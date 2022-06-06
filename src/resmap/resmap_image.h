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
#ifndef IMAGE_H_
#define IMAGE_H_

#include "./resmap_util.h"

#include "./resmap_fft_fftw3.h"
#include "./resmap_macros.h"
#include "./resmap_bessel.h"
#include "./resmap_error.h"
#include "./resmap_memory.h"
#include "./resmap_time.h"
#include "./resmap_macros.h"

#include "../../src_tools/l2cachemodel.h"

#define IS_INV true
#define IS_NOT_INV false
#define DONT_WRAP false
#define WRAP true
#ifndef DBL_EPSILON
#define DBL_EPSILON 1e-50
#endif

#ifndef NEAREST_NEIGHBOUR
#define NEAREST_NEIGHBOUR 0
#endif
#ifndef TRILINEAR
#define TRILINEAR 1
#endif

#undef LIN_INTERP
#define LIN_INTERP(a, l, h) ((l) + ((h) - (l)) * (a))


template<typename T>
void applyGeometry(const T *V1,T* V2,int dim,const double A[][3],bool inv,bool wrap,T outside = 0)
{
//#define DEBUG_APPLYGEO 1
    
    if (inv == IS_INV) {
        std::cerr<<"not implement inv.";
        throw("");
    }
    double Aref[3][3];
    for (int i = 0; i < 3; i++)
        for (int j = 0; j < 3; j++)
            Aref[i][j] = A[i][j];

    if (!inv)
    {
        Aref[0][0] =   A[2][2]*A[1][1]-A[2][1]*A[1][2];
        Aref[0][1] = -(A[2][2]*A[0][1]-A[2][1]*A[0][2]);
        Aref[0][2] =   A[1][2]*A[0][1]-A[1][1]*A[0][2];
        Aref[1][0] = -(A[2][2]*A[1][0]-A[2][0]*A[1][2]);
        Aref[1][1] =   A[2][2]*A[0][0]-A[2][0]*A[0][2];
        Aref[1][2] = -(A[1][2]*A[0][0]-A[1][0]*A[0][2]);
        Aref[2][0] =   A[2][1]*A[1][0]-A[2][0]*A[1][1];
        Aref[2][1] = -(A[2][1]*A[0][0]-A[2][0]*A[0][1]);
        Aref[2][2] =   A[1][1]*A[0][0]-A[1][0]*A[0][1];
        double tmp = A[0][0] * Aref[0][0]+A[1][0] * Aref[0][1] + A[2][0] * Aref[0][2];
        for (int i = 0; i < 3; i++)
            for (int j = 0; j < 3; j++)
                Aref[i][j] = Aref[i][j]/tmp;

    }
    
    // For scalings the output matrix is resized outside to the final
    // size instead of being resized inside the routine with the
    // same size as the input matrix

    // 2D transformation
    
    int m1, n1, m2, n2;
    double x, y, xp, yp;
    double minxp, minyp, maxxp, maxyp;
    int cen_x, cen_y, cen_xp, cen_yp;
    double wx, wy;
    int Xdim, Ydim;
    
    // Find center and limits of image
    cen_y  = (int)(dim / 2);
    cen_x  = (int)(dim / 2);
    cen_yp = (int)(dim / 2);
    cen_xp = (int)(dim / 2);
    minxp  = -cen_xp;
    minyp  = -cen_yp;
    maxxp  = dim - cen_xp - 1;
    maxyp  = dim - cen_xp - 1;
    Xdim   = dim;
    Ydim   = dim;
    
    // Now we go from the output image to the input image, ie, for any pixel
    // in the output image we calculate which are the corresponding ones in
    // the original image, make an interpolation with them and put this value
    // at the output pixel
    
#ifdef DEBUG_APPLYGEO
    std::cout<<Aref[0][0]<<" "<<Aref[0][1]<<" "<<Aref[0][2]<<std::endl;
    std::cout<<Aref[1][0]<<" "<<Aref[1][1]<<" "<<Aref[1][2]<<std::endl;
    std::cout<<Aref[2][0]<<" "<<Aref[2][1]<<" "<<Aref[2][2]<<std::endl;
    std::cout<< "(cen_x ,cen_y )=(" << cen_x  << "," << cen_y  << ")\n"
    << "(cen_xp,cen_yp)=(" << cen_xp << "," << cen_yp << ")\n"
    << "(min_xp,min_yp)=(" << minxp  << "," << minyp  << ")\n"
    << "(max_xp,max_yp)=(" << maxxp  << "," << maxyp  << ")\n";
#endif
    
    for (int i = 0; i < dim; i++)
    {
        // Calculate position of the beginning of the row in the output image
        x = -cen_x;
        y = i - cen_y;
        
        // Calculate this position in the input image according to the
        // geometrical transformation
        // they are related by
        // coords_output(=x,y) = A * coords_input (=xp,yp)

        xp = x * Aref[0][0] + y * Aref[0][1] + Aref[0][2];
        yp = x * Aref[1][0] + y * Aref[1][1] + Aref[1][2];

        for (int j = 0; j < dim; j++)
        {
            bool interp;
            T tmp;
            
#ifdef DEBUG_APPLYGEO
            
            std::cout << "Computing (" << i << "," << j << ")\n";
            std::cout << "   (y, x) =(" << y << "," << x << ")\n"
            << "   before wrapping (y',x')=(" << yp << "," << xp << ") "
            << std::endl;
#endif
            // If the point is outside the image, apply a periodic extension
            // of the image, what exits by one side enters by the other
            interp = true;
            if (wrap)
            {
                if (xp < minxp - ROME_EQUAL_ACCURACY ||
                    xp > maxxp + ROME_EQUAL_ACCURACY)
                    xp = realWRAP(xp, minxp - 0.5, maxxp + 0.5);
                
                if (yp < minyp - ROME_EQUAL_ACCURACY ||
                    yp > maxyp + ROME_EQUAL_ACCURACY)
                    
                    yp = realWRAP(yp, minyp - 0.5, maxyp + 0.5);
            }
            else
            {
                if (xp < minxp - ROME_EQUAL_ACCURACY ||
                    xp > maxxp + ROME_EQUAL_ACCURACY)
                    interp = false;
                
                if (yp < minyp - ROME_EQUAL_ACCURACY ||
                    yp > maxyp + ROME_EQUAL_ACCURACY)
                    interp = false;
            }
            
#ifdef DEBUG_APPLYGEO
            std::cout << "   after wrapping (y',x')=(" << yp << "," << xp << ") "
            << std::endl;
            std::cout << "   Interp = " << interp << std::endl;
            // The following line sounds dangerous...
            // x++;
#endif
            
            if (interp)
            {
                // Linear interpolation
                
                // Calculate the integer position in input image, be careful
                // that it is not the nearest but the one at the top left corner
                // of the interpolation square. Ie, (0.7,0.7) would give (0,0)
                // Calculate also weights for point m1+1,n1+1
                wx = xp + cen_xp;
                m1 = (int) wx;
                wx = wx - m1;
                m2 = m1 + 1;
                wy = yp + cen_yp;
                n1 = (int) wy;
                wy = wy - n1;
                n2 = n1 + 1;
                
                // m2 and n2 can be out by 1 so wrap must be checked here
                if (wrap)
                {
                    if (m2 >= Xdim)
                        m2 = 0;
                    if (n2 >= Ydim)
                        n2 = 0;
                }
                
#ifdef DEBUG_APPLYGEO
                std::cout << "   From (" << n1 << "," << m1 << ") and ("
                << n2 << "," << m2 << ")\n";
                std::cout << "   wx= " << wx << " wy= " << wy << std::endl;
#endif
                
                // Perform interpolation
                // if wx == 0 means that the rightest point is useless for this
                // interpolation, and even it might not be defined if m1=xdim-1
                // The same can be said for wy.
                tmp  = (T)((1 - wy) * (1 - wx) * V1[n1*dim+m1]);
                
                if (wx != 0 && m2 < dim)
                    tmp += (T)((1 - wy) * wx * V1[n1*dim+m2]);
                
                if (wy != 0 && n2 < dim)
                {
                    tmp += (T)(wy * (1 - wx) * V1[n2*dim+m1]);
                    
                    if (wx != 0 && m2 < dim)
                        tmp += (T)(wy * wx * V1[n2*dim+m2]);
                }
                
                V2[i*dim+j] = tmp;

            } // if interp
            else
                V2[i*dim+j] = outside;
            
            // Compute new point inside input image
            xp += Aref[0][0];
            yp += Aref[1][0];
        }
    }
}

template<typename T>
void translate(const T* V1,T *V2,int dim,double offset_x,double offset_y,bool wrap = WRAP, T outside = 0)
{
    double transMat[3][3];
    
    transMat[0][0] = 1;transMat[0][1] = 0;transMat[0][2] = offset_x;
    
    transMat[1][0] = 0;transMat[1][1] = 1;transMat[1][2] = offset_y;
    
    transMat[2][0] = 0;transMat[2][1] = 0;transMat[2][2] = 1;
    
    applyGeometry(V1, V2, dim,transMat, IS_NOT_INV, wrap, outside);
}

template<typename T>
void selfTranslate(T* V1,int dim,double offset_x,double offset_y,bool wrap = WRAP, T outside = 0)
{
    int dim2 = dim*dim;
    
    T* aux = (T*)aMalloc(sizeof(T)*dim2,64);
    
    for (int i = 0; i < dim2; i++) {
        aux[i] = V1[i];
    }
    
    translate(aux, V1, dim, offset_x,offset_y, wrap, outside);
    
    aFree(aux);
}

// Center an array, to have its origin at the origin of the FFTW
template <typename T>
void CenterFFT(T* v, int dim,bool forward)
{
    // 2D
    int l, shift;
    // Shift in the X direction
    l = dim;

    T* aux = (T*)aMalloc(sizeof(T)*l,64);
    
    shift = (int)(l / 2);
    
    if (!forward)
        shift = -shift;
    
    for (int i = 0; i < dim; i++)
    {
        // Shift the input in an auxiliar vector
        for (int j = 0; j < l; j++)
        {
            int jp = j + shift;
            
            if (jp < 0)
                jp += l;
            else if (jp >= l)
                jp -= l;
            
            aux[jp] = v[i*dim+j];
        }
        
        // Copy the vector
        for (int j = 0; j < l; j++)
            v[i*dim+j] = aux[j];
    }
    
    // Shift in the Y direction
    l = dim;
    shift = (int)(l / 2);
    
    if (!forward)
        shift = -shift;
    
    for (int j = 0; j < dim; j++)
    {
        // Shift the input in an auxiliar vector
        for (int i = 0; i < l; i++)
        {
            int ip = i + shift;
            
            if (ip < 0)
                ip += l;
            else if (ip >= l)
                ip -= l;
            
            aux[ip] = v[i*dim+j];
        }
        
        // Copy the vector
        for (int i = 0; i < l; i++)
            v[i*dim+j] = aux[i];
    }
    
    aFree(aux);
}


// Center an array,to have its origin at the origin of the FFTW
template <typename T>
void CenterFFT(const T& v_in, T& v_out, int dim,bool forward)
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


// Go from the Projector-centered fourier transform back to FFTW-uncentered one
template <typename T>
void decenter(const T *Min,int Min_size, T *Mout, int Mout_size,int my_rmax2)
{
    
    // Mout should already have the right size
    // NOTE : Initialize to zero,because we not sure the data type,so initialize zero outside this function.....
    int Min_size_shift = -(int)((float) (Min_size) / 2.0);
    int Min_Fsize = Min_size/2+1;
    int Mout_Fsize = Mout_size/2+1;
    
    for (int i = 0; i<Mout_size; i++){
        
        int ip = (i < Mout_Fsize) ? i : i - Mout_size;
        // my_rmax2 should be larger than ip*ip
        int jp2 = my_rmax2-ip*ip;
        if (jp2 < 0) continue;
        int j_max = sqrt(jp2)+1;
        j_max = j_max<Mout_Fsize?j_max:Mout_Fsize;
        ip = ip - Min_size_shift;
        
#pragma ivdep
        for (int j = 0; j<j_max; j++)
        {
            int jp = j;
            // if (ip*ip + jp*jp <= my_rmax2)
            Mout[i*Mout_Fsize+j] = Min[ip*Min_Fsize+jp];//only y dim shift
        }
    }
}

inline void decenter(const SOAComplexDouble& Min,int Min_size, SOAComplexDouble& Mout, int Mout_size,int my_rmax2);

// windows image,in_size or out_size should be even number
template<class T>
void windowTransform(const T *in,int in_size,T *out,int out_size)
{
    int in_Fsize = in_size/2 + 1;
    int out_Fsize = out_size/2 + 1;
    int out_Fsize2 = out_Fsize*out_size;
    // Check size of the input array
    if (in_size == out_size)
    {
        for (int i = 0; i < out_Fsize2; i++) {
            out[i] = in[i];
        }
        return;
    }
    
    
    for (int i = 0; i < out_Fsize2; i++) {
        out[i] = 0;
    }
    
    if (out_Fsize > in_Fsize)
    {
        int max_r2 = (in_Fsize -1) * (in_Fsize - 1);
        
        for (int i = 0; i<in_size; i++){
            
            int ip = (i < in_Fsize) ? i : i - in_size;
            int in_i = ((ip < 0) ? (ip + in_size) : (ip));
            int out_i = ((ip < 0) ? (ip + out_size) : (ip));
            int jp2 = max_r2-ip*ip;
            if (jp2 < 0) continue;
            int j_max = sqrt(jp2)+1;
            j_max = j_max < in_Fsize?j_max:in_Fsize;
            
#pragma ivdep
            for (int j = 0; j<j_max; j++)
            {
                // Make sure windowed FT has nothing in the corners, otherwise we end up with an asymmetric FT!
                int jp = j;
                int in_j = jp;
                int out_j = jp;
                
                //if (ip*ip + jp*jp <= max_r2)
                out[out_i*out_Fsize+out_j] = in[in_i*in_Fsize+in_j];
            }
        }
    }
    else
    {
        for ( int i = 0; i<out_size; i++){
            
            int ip = (i < out_Fsize) ? i : i - out_size;
            int in_i = ((ip < 0) ? (ip + in_size) : (ip));
            int out_i = ((ip < 0) ? (ip + out_size) : (ip));

#pragma ivdep
            for ( int j = 0; j<out_Fsize; j++)
            {
                
                int jp = j;
                int in_j = jp;
                int out_j = jp;
                
                out[out_i*out_Fsize+out_j] = in[in_i*in_Fsize+in_j];
            }
        }
    }
}

void windowFourierTransform(SOAComplexReadonly& in,int in_size,SOAComplexDouble& out,int out_size);

inline void windowFourierTransform(SOAComplexDouble& in,int in_size,SOAComplexDouble& out,int out_size) {
    SOAComplexReadonly inReadonly(in);
	windowFourierTransform(inReadonly, in_size, out, out_size);
}

template<typename T1,typename T2>
void windowFourierTransform(const T1& in_real,const T1& in_imag,int in_size,T2& out_real,T2& out_imag,int out_size)
{
    int in_Fsize = in_size/2 + 1;
    int out_Fsize = out_size/2 + 1;
    int out_Fsize2 = out_Fsize*out_size;
    // Check size of the input array
    if (in_size == out_size)
    {
        for (int i = 0; i < out_Fsize2; i++) {
            out_real[i] = in_real[i];
            out_imag[i] = in_imag[i];
        }
        return;
    }
    
    
    for (int i = 0; i < out_Fsize2; i++) {
        out_real[i] = 0;
        out_imag[i] = 0;
    }
    
    if (out_Fsize > in_Fsize)
    {
        int max_r2 = (in_Fsize -1) * (in_Fsize - 1);
        
        for (int i = 0; i<in_size; i++){
            
            int ip = (i < in_Fsize) ? i : i - in_size;
            int in_i = ((ip < 0) ? (ip + in_size) : (ip));
            int out_i = ((ip < 0) ? (ip + out_size) : (ip));
            // max_r2 should be larger than ip*ip
            int jp2 = max_r2-ip*ip;
            if (jp2 < 0) continue;
            int j_max = sqrt(jp2)+1;
            j_max = j_max < in_Fsize?j_max:in_Fsize;
            
#pragma ivdep
            for (int j = 0; j<j_max; j++)
            {
                // Make sure windowed FT has nothing in the corners, otherwise we end up with an asymmetric FT!
                int jp = j;
                int in_j = jp;
                int out_j = jp;
                
                // if (ip*ip + jp*jp <= max_r2)
                out_real[out_i*out_Fsize+out_j] = in_real[in_i*in_Fsize+in_j];
                out_imag[out_i*out_Fsize+out_j] = in_imag[in_i*in_Fsize+in_j];
            }
        }
    }
    else
    {
        for ( int i = 0; i<out_size; i++){
            
            int ip = (i < out_Fsize) ? i : i - out_size;
            int in_i = ((ip < 0) ? (ip + in_size) : (ip));
            int out_i = ((ip < 0) ? (ip + out_size) : (ip));
            
#pragma ivdep
            for ( int j = 0; j<out_Fsize; j++)
            {
                
                int jp = j;
                int in_j = jp;
                int out_j = jp;
                
                out_real[out_i*out_Fsize+out_j] = in_real[in_i*in_Fsize+in_j];
                out_imag[out_i*out_Fsize+out_j] = in_imag[in_i*in_Fsize+in_j];
            }
        }
    }
}

// Mask out corners outside sphere (replace by average value) Apply a soft mask (raised cosine with cosine_width pixels width)
template<typename T>
void softMaskOutsideMap(T* vol, int vol_size,double radius = -1., double cosine_width = 3, const T *Mnoise = nullptr)
{
    double r, radius_p, raisedcos, sum_bg = 0., sum = 0.;
    if (radius < 0)
        radius = (double)vol_size/2.;
    radius_p = radius + cosine_width;
    
    int vol_size_origin = -(int)((float) (vol_size) / 2.0);
    
    if (Mnoise == nullptr)
    {
        // Calculate average background value
        for (int i = 0; i<vol_size; i++)
            for (int j = 0; j < vol_size; j++)
            {
                r = sqrt((double)((i+vol_size_origin)*(i+vol_size_origin) + (j+vol_size_origin)*(j+vol_size_origin)));
                if (r < radius)
                    continue;
                else if (r > radius_p)
                {
                    sum    += 1.;
                    sum_bg += vol[i*vol_size+j];
                }
                else
                {
                    raisedcos = 0.5 + 0.5 * cos(PI * (radius_p - r) / cosine_width );
                    sum += raisedcos;
                    
                    sum_bg += raisedcos * vol[i*vol_size+j];
                }
            }
        sum_bg /= sum;
    }
    // Apply noisy or average background value
    for (int i = 0; i < vol_size; i++){
        
        for (int j = 0; j < vol_size; j++)
        {
            
            r = sqrt((double)((i+vol_size_origin)*(i+vol_size_origin) + (j+vol_size_origin)*(j+vol_size_origin)));
            
            if (r < radius)
            {
                continue;
            }
            
            if (r <= radius_p)
            {
                raisedcos = 0.5 + 0.5 * cos(PI * (radius_p - r) / cosine_width );
                
                double add = (Mnoise == nullptr) ?  sum_bg : Mnoise[i*vol_size+j];
                
                vol[i*vol_size+j] = (1 - raisedcos) * vol[i*vol_size+j] + raisedcos * add;
                
            }
            
            if (r > radius_p)
            {
                vol[i*vol_size+j] = (Mnoise == nullptr) ? sum_bg : Mnoise[i*vol_size+j];
            }
        }
        
    }
}

template<typename T>
void softMaskOutsideMap(const T& vol_in,T& vol_out,int vol_size,double radius, double cosine_width, const double *Mnoise)
{
    for (int i = 0; i < vol_size*vol_size; i++)
        vol_out[i] = vol_in[i];
    
    softMaskOutsideMap(&vol_out[0],vol_size,radius,cosine_width,Mnoise);
}

#define FIX_TAB_ROUND_ERROR // Yongbei debug,keep the complete same as old one.....
// shift image in fourier space
namespace TabulatedSinCos {
    static const int    nr_element      = 5000;//100000;//5000;
	static const double radiansPerEntry = 2 * PI / double(nr_element);
	class Table {
	public:
		Table();
		inline void sinCos(double & sin, double & cos, double dotpDividedBy2Pi) const {
#if !defined(FIX_TAB_ROUND_ERROR)
			auto i = int( fabs(dotpDividedBy2Pi)*nr_element) % nr_element;
#else
            int    nr_element      = 5000;//100000;//5000;
            double sampling_cos_sin = 2*PI/(double)nr_element;
            double dotpDivided = 2*PI*dotpDividedBy2Pi;
            int i = ((int)(fabs(dotpDivided)/sampling_cos_sin)) % nr_element;
#endif
#ifdef L2_CACHE_MODELING
			mostRecentAccess[i].store(epoch);
#endif
			cos = tabulated[i].cos;
			sin = tabulated[i].sin;
			sin = (dotpDividedBy2Pi < 0)?-sin:sin;
		}
		int initNotingSeqAcc() const;
		void finiNotingSeqAcc(int initResult) const;

	private:
		struct {double sin; double cos;} tabulated[nr_element];
#ifdef L2_CACHE_MODELING
		mutable std::atomic<int> epoch;
		mutable std::atomic<int> mostRecentAccess[nr_element];
#endif
	};
	extern const Table table;	// yuck but it keeps the tables const...
};

extern IntPerformanceCounter ShiftImageInFourierTransformNew_performanceCounter;
extern IntPerformanceCounter ShiftImageInFourierTransformNew_init_performanceCounter;

// Translation/ Time-Shifting : f(x-a)=e^(-iaw)F(w)=(cos(aw)+isin(aw))F(w)
extern bool ShiftImageInFourierTransformNew_useSinCosTable;

template<
typename T1,	// some indexable type
typename T2>	// some indexable type
class ShiftImageInFourierTransformNew {
private:
    int    inout_size;
    double shift_x;
    double shift_y;
    int    ori_size;
    
    int	 inout_Fsize;
    int	 inout_Fsize2;
    
    double xshift;
    double yshift;
    bool   justMove;
    
    int	 size1;
    int	 size2;
    
    double* aTable; // cos
    double* bTable; // sin
    int tableCapacity;
public:
    ShiftImageInFourierTransformNew()
    : 	inout_size(0), shift_x(0), shift_y(0), ori_size(0),
    	inout_Fsize (0),inout_Fsize2(0),
    	xshift(0),yshift(0),
        justMove(false),
        size1(0),size2(0),tableCapacity(0),
    	aTable(nullptr),bTable(nullptr)
    {
    }
    ~ShiftImageInFourierTransformNew() { 
		fini(); 
		finiTable(); 
	}
    void setNeededTableCapacity(int needed){
		if (needed <= tableCapacity) return;
        finiTable();
        tableCapacity = needed;
        aTable = (double*)aMalloc(sizeof(double)*tableCapacity,64);
        bTable = (double*)aMalloc(sizeof(double)*tableCapacity,64);
    }
    void finiTable(){
        tableCapacity = 0;
        aFree(aTable);
        aFree(bTable);
    }
    void init(int inout_size, double shift_x, double shift_y, int ori_size)
    {
        fini();
		ShiftImageInFourierTransformNew_init_performanceCounter.count.v++;
        inout_size = inout_size; shift_x = shift_x; shift_y = shift_y; ori_size = ori_size;
        inout_Fsize = inout_size / 2 + 1;
        inout_Fsize2 = inout_size*inout_Fsize;
        xshift = double(shift_x) / double(-ori_size);
        yshift = double(shift_y) / double(-ori_size);
        justMove = (fabs(xshift) < ROME_EQUAL_ACCURACY) && (fabs(yshift) < ROME_EQUAL_ACCURACY);
        size1  = (             inout_Fsize)*inout_Fsize;
        size2  = (inout_size - inout_Fsize)*inout_Fsize;
        assert(inout_Fsize2 == size1 + size2);
        assert(aTable);assert(bTable);assert(inout_Fsize2<=tableCapacity);
		L2CacheModel::seqAcc("ShiftImageInFourierTransformNew a_table", 0, &aTable[0], inout_Fsize2);
		L2CacheModel::seqAcc("ShiftImageInFourierTransformNew b_table", 0, &bTable[0], inout_Fsize2);
		const auto sinCosInit = TabulatedSinCos::table.initNotingSeqAcc();

		// top-half ring
        for (int i = 0; i < inout_Fsize; i++) {
            const double y = i;
            for (int j = 0; j < inout_Fsize; j++) {
                const double x = j;
                const double dotpDividedBy2Pi = (x * xshift + y * yshift);
                const int index = i*inout_Fsize + j;
                assert(0 <= index && index < size1);
                if (ShiftImageInFourierTransformNew_useSinCosTable) {
                    TabulatedSinCos::table.sinCos(bTable[index],
                                                  aTable[index],
                                                  dotpDividedBy2Pi);
                }
                else{
                    aTable[index] = cos(2*PI*dotpDividedBy2Pi);
                    bTable[index] = sin(2*PI*dotpDividedBy2Pi);
                }
            }
        }
        // bottom-half ring
        for (int i = inout_size - 1; i >= inout_Fsize; i--) {
            const double y = i - inout_size;
            for (int j = 0; j < inout_Fsize; j++) {
                const int index = (i - inout_Fsize)*inout_Fsize+j;
                assert(0 <= index && index < size2);
                const double x = j;
                const double dotpDividedBy2Pi = (x * xshift + y * yshift);
                if (ShiftImageInFourierTransformNew_useSinCosTable) {
                    TabulatedSinCos::table.sinCos(bTable[size1+index],
                                                  aTable[size1+index],
                                                  dotpDividedBy2Pi);
                }
                else{
                    aTable[size1+index] = cos(2*PI*dotpDividedBy2Pi);
                    bTable[size1+index] = sin(2*PI*dotpDividedBy2Pi);
                }
            }
        }
		TabulatedSinCos::table.finiNotingSeqAcc(sinCosInit);
	}
    void fini()
    {
        inout_size = 0; shift_x = 0; shift_y = 0; ori_size = 0;
        inout_Fsize = 0; inout_Fsize2 = 0;
        xshift = 0; yshift = 0;
        justMove = false;
        size1  = 0; size2  = 0;
    }
    inline       double* aTable_wptr()       { return aTable; }
    inline       double* bTable_wptr()       { return bTable; }
	inline const double* aTable_rptr() const { return aTable; }
	inline const double* bTable_rptr() const { return bTable; }
	inline       int     tableSize  () const { return inout_Fsize2; }
    
    void transform(const T1* fin_real, const T1* fin_imag, T2* fout_real, T2* fout_imag) const {
        transformOneTile(0, inout_Fsize2,
                         fin_real, fin_imag, fout_real, fout_imag);
    }
    void transformOneTile(const int nLo, const int nHi,
                          const T1* fin_real, const T1* fin_imag, T2* fout_real, T2* fout_imag) const {
        
		ShiftImageInFourierTransformNew_performanceCounter.count.v++;

		L2CacheModel::seqAcc("fout_real", -1, &fout_real[nLo], nHi - nLo);
		L2CacheModel::seqAcc("fout_imag", -1, &fout_imag[nLo], nHi - nLo);
		L2CacheModel::seqAcc("fin_real", -1,  &fin_real [nLo], nHi - nLo);
		L2CacheModel::seqAcc("fin_imag", -1,  &fin_imag [nLo], nHi - nLo);

		assert(0 <= nLo); assert(nHi <= inout_Fsize2);
        if (justMove) {
            for (int n = nLo; n < nHi; n++) {
                fout_real[n] = CHECK_NOT_IND(fin_real[n]);
                fout_imag[n] = CHECK_NOT_IND(fin_imag[n]);
            }
            return;
        }

		L2CacheModel::seqAcc("aTable", -1, &aTable[nLo], nHi - nLo);
		L2CacheModel::seqAcc("bTable", -1, &bTable[nLo], nHi - nLo);

#define LOOP \
        for (int n = nLo; n < nHi; n++) {									\
            const double a = aTable	 [n];									\
            const double b = bTable	 [n];									\
            const double c = fin_real[n];									\
            const double d = fin_imag[n];									\
            const double ac = a * c;										\
            const double bd = b * d;										\
            const double ab_cd = (a + b) * (c + d);							\
            fout_real[n] = CHECK_NOT_IND(ac - bd        	        );		\
            fout_imag[n] = CHECK_NOT_IND(ab_cd - ac - bd			);		\
        }																	\
// end of macro
        if (isVectorAligned(&aTable[nLo])
            && isVectorAligned(&bTable[nLo])
            && isVectorAligned(&fin_real[nLo])
            && isVectorAligned(&fin_imag[nLo])
            && isVectorAligned(&fout_real[nLo])
            && isVectorAligned(&fout_imag[nLo])
            ) {
#pragma vector aligned
#pragma ivdep
            LOOP
        } else {
#pragma ivdep
            LOOP
        }
#undef LOOP
    }
    //
    // four cases : shiftx+shifty,shiftx,shifty,none shift
    // Yongbei: this is a template class and use NOINLINE
    // warning #2196: routine is both "inline" and "noinline"
    void /*NOINLINE*/ transformOneTileByShiftXY(const int nLo, const int nHi,
                                   const T1* fin_real, const T1* fin_imag, T2* fout_real, T2* fout_imag,
                                   const T1* aTableX, const T1* bTableX, int XSize, double Xstatus,
                                   const T1* aTableY, const T1* bTableY, int YSize, double Ystatus)
    {
		ShiftImageInFourierTransformNew_performanceCounter.count.v++;
		
		L2CacheModel::seqAcc("fout_real", -1, &fout_real[nLo], nHi - nLo);
        L2CacheModel::seqAcc("fout_imag", -1, &fout_imag[nLo], nHi - nLo);
        L2CacheModel::seqAcc("fin_real", -1, &fin_real[nLo], nHi - nLo);
        L2CacheModel::seqAcc("fin_imag", -1, &fin_imag[nLo], nHi - nLo);
        
        assert(XSize == YSize);
        assert(0 <= nLo);assert(nHi <= XSize);assert(nLo <= nHi);
        
        L2CacheModel::seqAcc("aTable", -1, &aTableX[nLo], nHi - nLo);
        L2CacheModel::seqAcc("bTable", -1, &bTableX[nLo], nHi - nLo);
        L2CacheModel::seqAcc("aTable", -1, &aTableY[nLo], nHi - nLo);
        L2CacheModel::seqAcc("bTable", -1, &bTableY[nLo], nHi - nLo);
        
#define LOOP_SHIFT_XY \
        for (int n = nLo; n < nHi; n++) {									\
            const double a1 = aTableX [n];									\
            const double b1 = Xstatus*bTableX [n];							\
            const double c1 = aTableY [n];									\
            const double d1 = Ystatus*bTableY [n];							\
            const double ac1 = a1 * c1;										\
            const double bd1 = b1 * d1;										\
            const double ab_cd1 = (a1 + b1) * (c1 + d1);					\
            const double a = CHECK_NOT_IND(ac1 - bd1         );				\
            const double b = CHECK_NOT_IND(ab_cd1 - ac1 - bd1);				\
            const double c = fin_real[n];									\
            const double d = fin_imag[n];									\
            const double ac = a * c;										\
            const double bd = b * d;										\
            const double ab_cd = (a + b) * (c + d);							\
            fout_real[n] = CHECK_NOT_IND(ac - bd        	        );		\
            fout_imag[n] = CHECK_NOT_IND(ab_cd - ac - bd			);		\
        }														// end of macro
        if (isVectorAligned(&aTableX[nLo])
            && isVectorAligned(&bTableX[nLo])
            && isVectorAligned(&aTableY[nLo])
            && isVectorAligned(&bTableY[nLo])
            && isVectorAligned(&fin_real[nLo])
            && isVectorAligned(&fin_imag[nLo])
            && isVectorAligned(&fout_real[nLo])
            && isVectorAligned(&fout_imag[nLo])
            ) {
#pragma vector aligned
#pragma ivdep
            LOOP_SHIFT_XY
        } else {
#pragma ivdep
            LOOP_SHIFT_XY
        }
#undef LOOP_SHIFT_XY
    }

    void /*NOINLINE*/ transformOneTileByShiftX(const int nLo, const int nHi,
                                   const T1* fin_real, const T1* fin_imag, T2* fout_real, T2* fout_imag,
                                   const T1* aTableX, const T1* bTableX, int XSize, double Xstatus)
    {
		ShiftImageInFourierTransformNew_performanceCounter.count.v++;

		L2CacheModel::seqAcc("fout_real", -1, &fout_real[nLo], nHi - nLo);
        L2CacheModel::seqAcc("fout_imag", -1, &fout_imag[nLo], nHi - nLo);
        L2CacheModel::seqAcc("fin_real", -1, &fin_real[nLo], nHi - nLo);
        L2CacheModel::seqAcc("fin_imag", -1, &fin_imag[nLo], nHi - nLo);
        
        assert(0 <= nLo);assert(nHi <= XSize);assert(nLo <= nHi);
        
        L2CacheModel::seqAcc("aTable", -1, &aTableX[nLo], nHi - nLo);
        L2CacheModel::seqAcc("bTable", -1, &bTableX[nLo], nHi - nLo);
        
#define LOOP_SHIFT_X \
        for (int n = nLo; n < nHi; n++) {									\
            const double a = aTableX [n];									\
            const double b = Xstatus*bTableX [n];							\
            const double c = fin_real[n];									\
            const double d = fin_imag[n];									\
            const double ac = a * c;										\
            const double bd = b * d;										\
            const double ab_cd = (a + b) * (c + d);							\
            fout_real[n] = CHECK_NOT_IND(ac - bd        	        );		\
            fout_imag[n] = CHECK_NOT_IND(ab_cd - ac - bd			);		\
        }														// end of macro
        
        if (isVectorAligned(&aTableX[nLo])
            && isVectorAligned(&bTableX[nLo])
            && isVectorAligned(&fin_real[nLo])
            && isVectorAligned(&fin_imag[nLo])
            && isVectorAligned(&fout_real[nLo])
            && isVectorAligned(&fout_imag[nLo])
            ) {
#pragma vector aligned
#pragma ivdep
            LOOP_SHIFT_X
        } else {
#pragma ivdep
            LOOP_SHIFT_X
        }
#undef LOOP_SHIFT_X
    }

    void /*NOINLINE*/ transformOneTileByShiftY(const int nLo, const int nHi,
                                   const T1* fin_real, const T1* fin_imag, T2* fout_real, T2* fout_imag,
                                   const T1* aTableY, const T1* bTableY, int YSize, double Ystatus)
    {
		ShiftImageInFourierTransformNew_performanceCounter.count.v++;

		L2CacheModel::seqAcc("fout_real", -1, &fout_real[nLo], nHi - nLo);
        L2CacheModel::seqAcc("fout_imag", -1, &fout_imag[nLo], nHi - nLo);
        L2CacheModel::seqAcc("fin_real", -1, &fin_real[nLo], nHi - nLo);
        L2CacheModel::seqAcc("fin_imag", -1, &fin_imag[nLo], nHi - nLo);
        
        assert(0 <= nLo);assert(nHi <= YSize);assert(nLo <= nHi);
        
        L2CacheModel::seqAcc("aTable", -1, &aTableY[nLo], nHi - nLo);
        L2CacheModel::seqAcc("bTable", -1, &bTableY[nLo], nHi - nLo);
        
#define LOOP_SHIFT_Y \
        for (int n = nLo; n < nHi; n++) {									\
            const double a = aTableY [n];									\
            const double b = Ystatus*bTableY [n];							\
            const double c = fin_real[n];									\
            const double d = fin_imag[n];									\
            const double ac = a * c;										\
            const double bd = b * d;										\
            const double ab_cd = (a + b) * (c + d);							\
            fout_real[n] = CHECK_NOT_IND(ac - bd        	        );		\
            fout_imag[n] = CHECK_NOT_IND(ab_cd - ac - bd			);		\
        }														// end of macro
        
        if (   isVectorAligned(&aTableY  [nLo])
            && isVectorAligned(&bTableY  [nLo])
            && isVectorAligned(&fin_real [nLo])
            && isVectorAligned(&fin_imag [nLo])
            && isVectorAligned(&fout_real[nLo])
            && isVectorAligned(&fout_imag[nLo])
            ) {
#pragma vector aligned
#pragma ivdep
            LOOP_SHIFT_Y
        } else {
#pragma ivdep
            LOOP_SHIFT_Y
        }
#undef LOOP_SHIFT_Y
    }

    void /*NOINLINE*/ transformOneTileByCopy(const int nLo, const int nHi,
                                const T1* fin_real, const T1* fin_imag, T2* fout_real, T2* fout_imag)
    {
		ShiftImageInFourierTransformNew_performanceCounter.count.v++;

		L2CacheModel::seqAcc("fout_real", -1, &fout_real[nLo], nHi - nLo);
        L2CacheModel::seqAcc("fout_imag", -1, &fout_imag[nLo], nHi - nLo);
        L2CacheModel::seqAcc("fin_real",  -1, &fin_real [nLo], nHi - nLo);
        L2CacheModel::seqAcc("fin_imag",  -1, &fin_imag [nLo], nHi - nLo);
        
#define LOOP_COPY \
        for (int n = nLo; n < nHi; n++) {									\
            fout_real[n] = fin_real[n];										\
            fout_imag[n] = fin_imag[n];										\
        }														// end of macro
        
        if (isVectorAligned(&fin_real[nLo])
            && isVectorAligned(&fin_imag[nLo])
            && isVectorAligned(&fout_real[nLo])
            && isVectorAligned(&fout_imag[nLo])
            ) {
#pragma vector aligned
#pragma ivdep
            LOOP_COPY
        } else {
#pragma ivdep
            LOOP_COPY
        }
#undef LOOP_COPY
    }
};

extern IntPerformanceCounter ShiftImageInFourierTransform1_performanceCounter;
extern IntPerformanceCounter ShiftImageInFourierTransform2_performanceCounter;
extern IntPerformanceCounter ShiftImageInFourierTransform4_performanceCounter;

template<
	typename T1,	// some indexable type
	typename T2>	// some indexable type
class ShiftImageInFourierTransform {
protected:
	const int    inout_size;
	const double shift_x;
	const double shift_y;
	const int    ori_size;

	const int	 inout_Fsize;
	const int	 inout_Fsize2;

	const double xshift;
	const double yshift;
	const bool   justMove;

	const int	 size1;
	const int	 size2;
	const int    _tileSize;
	const int    _tileCount;

	double* aTable1;
	double* bTable1;
	double* aTable2;
	double* bTable2;

	bool isVectorAligned(void const* p) { return (size_t(p) % 64 == 0); }

public:
	ShiftImageInFourierTransform(
		int inout_size, double shift_x, double shift_y, int ori_size, bool useTabulation = true)
	  : inout_size(inout_size), shift_x(shift_x), shift_y(shift_y), ori_size(ori_size),
	  	inout_Fsize (inout_size / 2 + 1),
		inout_Fsize2(inout_size*inout_Fsize),
		xshift(double(shift_x) / double(-ori_size)),
		yshift(double(shift_y) / double(-ori_size)),
		justMove(fabs(xshift) < ROME_EQUAL_ACCURACY && fabs(yshift) < ROME_EQUAL_ACCURACY),

		size1 ((             inout_Fsize)*inout_Fsize),
		size2 ((inout_size - inout_Fsize)*inout_Fsize),

		_tileSize(256),		// fine grained tiling to get transform4 to fit - maybe :-)
		_tileCount(std::max(1,(inout_Fsize2 + _tileSize-1) / _tileSize))
	{
		assert(inout_Fsize2 == size1 + size2);

		aTable1 = (double*)aMalloc(sizeof(double)*size1,64); // new double[size1];
		bTable1 = (double*)aMalloc(sizeof(double)*size1,64); // new double[size1];
		for (int i = 0; i < inout_Fsize; i++) {
			const double y = i;
			for (int j = 0; j < inout_Fsize; j++) {
				const double x = j;
				const double dotpDividedBy2Pi = (x * xshift + y * yshift);
				const int index = i*inout_Fsize + j;
				assert(0 <= index && index < size1);
                if (useTabulation) {
                    TabulatedSinCos::table.sinCos(bTable1[index],
                                                  aTable1[index], 
                                                  dotpDividedBy2Pi);
                }
                else{
                    aTable1[index] = cos(2*PI*dotpDividedBy2Pi);
                    bTable1[index] = sin(2*PI*dotpDividedBy2Pi);
                }
			}
		}

		aTable2 = (double*)aMalloc(sizeof(double)*size2,64); // new double[size2];
		bTable2 = (double*)aMalloc(sizeof(double)*size2,64); // new double[size2];
		for (int i = inout_size - 1; i >= inout_Fsize; i--) {
			const double y = i - inout_size;
			for (int j = 0; j < inout_Fsize; j++) {
				const int index = (i - inout_Fsize)*inout_Fsize+j;
				assert(0 <= index && index < size2);
				const double x = j;
				const double dotpDividedBy2Pi = (x * xshift + y * yshift);
                if (useTabulation) {
                    TabulatedSinCos::table.sinCos(bTable2[index],
                                                  aTable2[index], 
                                                  dotpDividedBy2Pi);
                }
                else{
                    aTable2[index] = cos(2*PI*dotpDividedBy2Pi);
                    bTable2[index] = sin(2*PI*dotpDividedBy2Pi);
                }
			}
		}
	}

	~ShiftImageInFourierTransform() {
        if(aTable1) aFree(aTable1);
		if(bTable1) aFree(bTable1);
		if(aTable2) aFree(aTable2);
		if(bTable2) aFree(bTable2);
	}
    
	int tileCount() {
		return _tileCount;
	}

	void transform(
		const T1* fin_real, const T1* fin_imag, T2* fout_real, T2* fout_imag) {
		transformOneTile(
			0, inout_Fsize2, 
			fin_real, fin_imag, fout_real, fout_imag);
	}

	void transformTiled(int tile,
		const T1* fin_real, const T1* fin_imag, T2* fout_real, T2* fout_imag) {
		transformOneTile(
			tile*_tileSize, std::min(tile*_tileSize+_tileSize,inout_Fsize2),
			fin_real, fin_imag, fout_real, fout_imag);
	}

	void transformOneTile(
		const int xLo, const int xHi, 
		const T1* fin_real, const T1* fin_imag, T2* fout_real, T2* fout_imag) {

		ShiftImageInFourierTransform1_performanceCounter.count.v++;

		assert(0 <= xLo); assert(xHi <= inout_Fsize2);
		if (justMove) {
			for (int x = xLo; x < xHi; x++) {
				fout_real[x] = CHECK_NOT_IND(fin_real[x]);
				fout_imag[x] = CHECK_NOT_IND(fin_imag[x]);
			}
			return;
		}

		const int xLo1 = xLo; 
		const int xHi1 = std::min(xHi,size1); 
#define LOOP \
		for (int x = xLo1; x < xHi1; x++) {									\
			const double a = aTable1 [x];									\
			const double b = bTable1 [x];									\
			const double c = fin_real[x];									\
			const double d = fin_imag[x];									\
			const double ac = a * c;										\
			const double bd = b * d;										\
			const double ab_cd = (a + b) * (c + d);							\
			fout_real[x] = CHECK_NOT_IND(ac - bd        	        );		\
			fout_imag[x] = CHECK_NOT_IND(ab_cd - ac - bd			);		\
		}																	\
		// end of macro
		if (isVectorAligned(&aTable1[xLo1])
            && isVectorAligned(&bTable1[xLo1])
            && isVectorAligned(&fin_real[xLo1])
            && isVectorAligned(&fin_imag[xLo1])
            && isVectorAligned(&fout_real[xLo1])
            && isVectorAligned(&fout_imag[xLo1])
			) {
			#pragma vector aligned
			#pragma ivdep
			LOOP
		} else {
			#pragma ivdep
			LOOP
		}
#undef LOOP

		const int indexLo = std::max(size1, xLo); 
		const int indexHi = std::min(inout_Fsize2, xHi); 
#define LOOP															\
		for (int index = indexLo; index < indexHi; index++) {			\
			const int x = index - size1;								\
			const double a = aTable2 [x];								\
			const double b = bTable2 [x];								\
			assert(index < (inout_size-1)*inout_Fsize + inout_Fsize);	\
			const double c = fin_real[index];							\
			const double d = fin_imag[index];							\
			const double ac = a * c;									\
			const double bd = b * d;									\
			const double ab_cd = (a + b) * (c + d);						\
			fout_real[index] = CHECK_NOT_IND(ac - bd        	);		\
			fout_imag[index] = CHECK_NOT_IND(ab_cd - ac - bd	);		\
		}																\
		// end of macro
		if (isVectorAligned(&aTable2[indexLo-size1])
            && isVectorAligned(&bTable2[indexLo-size1])
            && isVectorAligned(&fin_real[indexLo])
            && isVectorAligned(&fin_imag[indexLo])
            && isVectorAligned(&fout_real[indexLo])
            && isVectorAligned(&fout_imag[indexLo])
			) {
			#pragma vector aligned
			#pragma ivdep
			LOOP
		} else {
			#pragma ivdep
			LOOP
		}
#undef LOOP
	}

	void transform2(
		const T1* fin0_real, const T1* fin0_imag, T2* fout0_real, T2* fout0_imag,
		const T1* fin1_real, const T1* fin1_imag, T2* fout1_real, T2* fout1_imag
		) {
		transform2OneTile(
			0, inout_Fsize2, 
			fin0_real, fin0_imag, fout0_real, fout0_imag, 
			fin1_real, fin1_imag, fout1_real, fout1_imag);
	}
	void transform2Tiled(int tile,
		const T1* fin0_real, const T1* fin0_imag, T2* fout0_real, T2* fout0_imag,
		const T1* fin1_real, const T1* fin1_imag, T2* fout1_real, T2* fout1_imag
		) {
		transform2OneTile(
			tile*_tileSize, std::min(tile*_tileSize+_tileSize,inout_Fsize2),
			fin0_real, fin0_imag, fout0_real, fout0_imag, 
			fin1_real, fin1_imag, fout1_real, fout1_imag);
	}
	void transform2OneTile(
		const int xLo, const int xHi,
		const T1* fin0_real, const T1* fin0_imag, T2* fout0_real, T2* fout0_imag,
		const T1* fin1_real, const T1* fin1_imag, T2* fout1_real, T2* fout1_imag
		) {

		ShiftImageInFourierTransform2_performanceCounter.count.v++;

		assert(0 <= xLo); assert(xHi <= inout_Fsize2);

		if (justMove) {
			for (int x = xLo; x < xHi; x++)  {
#define M(N) \
				fout##N##_real[x] = CHECK_NOT_IND(fin##N##_real[x]); \
				fout##N##_imag[x] = CHECK_NOT_IND(fin##N##_imag[x]); \
				// end of macro
				M(0) M(1)
#undef M
			}
			return;
		}

#define M(N) {\
			const double c = fin##N##_real[x];								\
			const double d = fin##N##_imag[x];								\
			const double ac = a * c;										\
			const double bd = b * d;										\
			const double ab_cd = (a + b) * (c + d);							\
			fout##N##_real[x] = CHECK_NOT_IND(ac - bd        			);	\
			fout##N##_imag[x] = CHECK_NOT_IND(ab_cd - ac - bd			);	\
			}	// end of macro
		const int xLo1 = xLo; 
		const int xHi1 = std::min(xHi,size1); 
#define LOOP \
		for (int x = xLo1; x < xHi1; x++) {					\
			const double a = aTable1 [x];					\
			const double b = bTable1 [x];					\
			M(0) M(1)										\
		}													\
		// end of macro
		if (isVectorAligned(fin0_real) && isVectorAligned(fin0_imag) && isVectorAligned(fout0_real) && isVectorAligned(fout0_imag) &&
		    isVectorAligned(fin1_real) && isVectorAligned(fin1_imag) && isVectorAligned(fout1_real) && isVectorAligned(fout1_imag)
			) {
			#pragma vector aligned
			#pragma ivdep
			LOOP
		} else {
			#pragma ivdep
			LOOP
		}
#undef LOOP
#undef M

#define M(N) {\
			const double c = fin##N##_real[index];			\
			const double d = fin##N##_imag[index];			\
			const double ac = a * c;						\
			const double bd = b * d;						\
			const double ab_cd = (a + b) * (c + d);			\
			fout##N##_real[index] = CHECK_NOT_IND(ac - bd		 );	\
			fout##N##_imag[index] = CHECK_NOT_IND(ab_cd - ac - bd);	\
			}	// end of macro

		const int indexLo = std::max(size1, xLo); 
		const int indexHi = std::min(inout_Fsize2, xHi); 
#define LOOP															\
		for (int index = indexLo; index < indexHi; index++) {			\
			const int x = index - size1;								\
			const double a = aTable2 [x];								\
			const double b = bTable2 [x];								\
			M(0) M(1)													\
		}																\
		// end of macro
		if (isVectorAligned(fin0_real) && isVectorAligned(fin0_imag) && isVectorAligned(fout0_real) && isVectorAligned(fout0_imag) &&
		    isVectorAligned(fin1_real) && isVectorAligned(fin1_imag) && isVectorAligned(fout1_real) && isVectorAligned(fout1_imag)
			) {
			#pragma vector aligned
			#pragma ivdep
			LOOP
		} else {
			#pragma ivdep
			LOOP
		}
#undef LOOP
#undef M
	}

	void transform4(
		const T1* fin0_real, const T1* fin0_imag, T2* fout0_real, T2* fout0_imag,
		const T1* fin1_real, const T1* fin1_imag, T2* fout1_real, T2* fout1_imag,
		const T1* fin2_real, const T1* fin2_imag, T2* fout2_real, T2* fout2_imag,
		const T1* fin3_real, const T1* fin3_imag, T2* fout3_real, T2* fout3_imag
		) {
		transform4OneTile(
			0, inout_Fsize2, 
			fin0_real, fin0_imag, fout0_real, fout0_imag, 
			fin1_real, fin1_imag, fout1_real, fout1_imag, 
			fin2_real, fin2_imag, fout2_real, fout2_imag, 
			fin3_real, fin3_imag, fout3_real, fout3_imag);
	}
	void transform4Tiled(int tile,
		const T1* fin0_real, const T1* fin0_imag, T2* fout0_real, T2* fout0_imag,
		const T1* fin1_real, const T1* fin1_imag, T2* fout1_real, T2* fout1_imag,
		const T1* fin2_real, const T1* fin2_imag, T2* fout2_real, T2* fout2_imag,
		const T1* fin3_real, const T1* fin3_imag, T2* fout3_real, T2* fout3_imag
		) {
		transform4OneTile(
			tile*_tileSize, std::min(tile*_tileSize+_tileSize,inout_Fsize2),
			fin0_real, fin0_imag, fout0_real, fout0_imag, 
			fin1_real, fin1_imag, fout1_real, fout1_imag, 
			fin2_real, fin2_imag, fout2_real, fout2_imag, 
			fin3_real, fin3_imag, fout3_real, fout3_imag);
	}
	void transform4OneTile(
		const int xLo, const int xHi,
		const T1* fin0_real, const T1* fin0_imag, T2* fout0_real, T2* fout0_imag,
		const T1* fin1_real, const T1* fin1_imag, T2* fout1_real, T2* fout1_imag,
		const T1* fin2_real, const T1* fin2_imag, T2* fout2_real, T2* fout2_imag,
		const T1* fin3_real, const T1* fin3_imag, T2* fout3_real, T2* fout3_imag
		) {
		assert(0 <= xLo); assert(xHi <= inout_Fsize2);

		ShiftImageInFourierTransform4_performanceCounter.count.v++;

		if (justMove) {
			for (int x = xLo; x < xHi; x++) {
#define M(N) \
				fout##N##_real[x] = CHECK_NOT_IND(fin##N##_real[x]);	\
				fout##N##_imag[x] = CHECK_NOT_IND(fin##N##_imag[x]);	\
				// end of macro
				M(0) M(1) M(2) M(3)
#undef M
			}
			return;
		}

#define M(N) {\
			const double c = fin##N##_real[x];						\
			const double d = fin##N##_imag[x];						\
			const double ac = a * c;								\
			const double bd = b * d;								\
			const double ab_cd = (a + b) * (c + d);					\
			fout##N##_real[x] = CHECK_NOT_IND(ac - bd        );		\
			fout##N##_imag[x] = CHECK_NOT_IND(ab_cd - ac - bd);		\
			}	// end of macro

		const int xLo1 = xLo; 
		const int xHi1 = std::min(xHi,size1); 
#define LOOP \
		for (int x = xLo1; x < xHi1; x++) {					\
				const double a = aTable1 [x];				\
				const double b = bTable1 [x];				\
				M(0) M(1) M(2) M(3)							\
			}												\
			// end of macro
		if (   isVectorAligned( fin0_real) && isVectorAligned( fin1_real) && isVectorAligned( fin2_real) && isVectorAligned( fin3_real)
			&& isVectorAligned( fin0_imag) && isVectorAligned( fin1_imag) && isVectorAligned( fin2_imag) && isVectorAligned( fin3_imag)
			&& isVectorAligned(fout0_real) && isVectorAligned(fout1_real) && isVectorAligned(fout2_real) && isVectorAligned(fout3_real)
			&& isVectorAligned(fout0_imag) && isVectorAligned(fout1_imag) && isVectorAligned(fout2_imag) && isVectorAligned(fout3_imag)
			) {
			#pragma vector aligned
			#pragma ivdep
			LOOP
		} else {
			#pragma ivdep
			LOOP
		}

#undef LOOP
#undef M

#define M(N) {\
			const double c = fin##N##_real[index];			\
			const double d = fin##N##_imag[index];			\
			const double ac = a * c;						\
			const double bd = b * d;						\
			const double ab_cd = (a + b) * (c + d);			\
			fout##N##_real[index] = CHECK_NOT_IND(ac - bd		 );	\
			fout##N##_imag[index] = CHECK_NOT_IND(ab_cd - ac - bd);	\
			}	// end of macro

		const int indexLo = std::max(size1, xLo); 
		const int indexHi = std::min(inout_Fsize2, xHi); 
#define LOOP															\
		for (int index = indexLo; index < indexHi; index++) {			\
			const int x = index - size1;								\
			const double a = aTable2 [x];								\
			const double b = bTable2 [x];								\
			M(0) M(1) M(2) M(3)											\
		}																\
		// end of macro

		if (   isVectorAligned( fin0_real) && isVectorAligned( fin1_real) && isVectorAligned( fin2_real) && isVectorAligned( fin3_real)
			&& isVectorAligned( fin0_imag) && isVectorAligned( fin1_imag) && isVectorAligned( fin2_imag) && isVectorAligned( fin3_imag)
			&& isVectorAligned(fout0_real) && isVectorAligned(fout1_real) && isVectorAligned(fout2_real) && isVectorAligned(fout3_real)
			&& isVectorAligned(fout0_imag) && isVectorAligned(fout1_imag) && isVectorAligned(fout2_imag) && isVectorAligned(fout3_imag)
			) {
			#pragma vector aligned
			#pragma ivdep
			LOOP
		} else {
			#pragma ivdep
			LOOP
		}
#undef LOOP
#undef M
	}
};


extern IntPerformanceCounter ShiftImageInFourierTransformSimple_performanceCounter;

template<typename T1,typename T2>
void shiftImageInFourierTransform(const T1& fin_real,const T1& fin_imag,T2& fout_real,T2& fout_imag,int inout_size,
                                  double shift_x,double shift_y,int ori_size)
{
	ShiftImageInFourierTransformSimple_performanceCounter.count.v++;

    static int nr_element = 5000;
    static double tabulated_sin[5000];
    static double tabulated_cos[5000];
    static bool reCalculate = true;
    
    double dotp, a, b, c, d, ac, bd, ab_cd, x, y, xshift, yshift;// z,zshift;
    int inout_Fsize = inout_size/2+1;
    int inout_Fsize2 = inout_size*(inout_size/2+1);
    double sampling_sin = 2 * PI / (double) nr_element;
    double sampling_cos = 2 * PI / (double) nr_element;
    
    xshift = (double)shift_x/(-ori_size);
    yshift = (double)shift_y/(-ori_size);
    
    // initialize for sin element value and cos element value element
    if (reCalculate) {
        
        for (int i = 0; i < nr_element; i++)
        {
            double xx = (double) i * sampling_sin;
            // positive number
            tabulated_sin[i] = sin(xx);
        }
        
        for (int i = 0; i < nr_element; i++)
        {
            double xx = (double) i * sampling_cos;
            tabulated_cos[i] = cos(xx);
        }
        // std::cout<<"calculate tabulated_sin."<<std::endl;
        reCalculate = false;
    }
    
    if (fabs(xshift) < ROME_EQUAL_ACCURACY && fabs(yshift) < ROME_EQUAL_ACCURACY)
    {
        for (int i = 0; i < inout_Fsize2; i++) {
            fout_real[i] = CHECK_NOT_IND(fin_real[i]);
            fout_imag[i] = CHECK_NOT_IND(fin_imag[i]);
        }
        return;
    }
    
    for (int i=0; i<inout_Fsize; i++)
    {
        y = i;
        
        for (int j=0; j<inout_Fsize; j++)
        {
            x = j;
            
            dotp = 2 * PI * (x * xshift + y * yshift);
            a = tabulated_cos[((int)( fabs(dotp) / sampling_cos)) % nr_element];
            b = tabulated_sin[((int)( fabs(dotp) / sampling_sin)) % nr_element];
            b = (dotp < 0)?-b:b;
            
            c = fin_real[i*inout_Fsize+j];
            d = fin_imag[i*inout_Fsize+j];
            
            ac = a * c;
            bd = b * d;
            ab_cd = (a + b) * (c + d);
            fout_real[i*inout_Fsize+j] = CHECK_NOT_IND(ac - bd);
            fout_imag[i*inout_Fsize+j] = CHECK_NOT_IND(ab_cd - ac - bd);
        }
    }
    
    for (int i=inout_size-1; i>=inout_Fsize; i--)
    {
        y = i - inout_size;
        
        for (int j=0; j<inout_Fsize; j++)
        {
            x = j;
            dotp = 2 * PI * (x * xshift + y * yshift);
            
            a = tabulated_cos[((int)( fabs(dotp) / sampling_cos)) % nr_element];
            b = tabulated_sin[((int)( fabs(dotp) / sampling_sin)) % nr_element];
            b = (dotp < 0)?-b:b;
            
            c = fin_real[i*inout_Fsize+j];
            d = fin_imag[i*inout_Fsize+j];
            
            ac = a * c;
            bd = b * d;
            ab_cd = (a + b) * (c + d);
            fout_real[i*inout_Fsize+j] = CHECK_NOT_IND(ac - bd);
            fout_imag[i*inout_Fsize+j] = CHECK_NOT_IND(ab_cd - ac - bd);
        }
    }
    
}


inline void shiftImageInFourierTransform(
	const SOAComplexReadonly& fin, const SOAComplexDouble& fout, int inout_size, double shift_x, double shift_y, int ori_size) {

	shiftImageInFourierTransform(
		fin.real, fin.imag, fout.real, fout.imag, inout_size, shift_x, shift_y, ori_size);
}

inline void shiftImageInFourierTransform(
	const SOAComplexDouble& fin, const SOAComplexDouble& fout, int inout_size, double shift_x, double shift_y, int ori_size) {

	shiftImageInFourierTransform(
		fin.real, fin.imag, fout.real, fout.imag, inout_size, shift_x, shift_y, ori_size);
}

template<typename T1,typename T2>
void shiftImageInFourierTransformNoTab(const T1* fin_real,const T1* fin_imag,T2* fout_real,T2* fout_imag,int inout_size,
                                       double shift_x,double shift_y,int ori_size)
{
    int inout_Fsize = inout_size/2+1;
    int inout_Fsize2 = inout_size*(inout_size/2+1);
    double dotp, a, b, c, d, ac, bd, ab_cd, x, y, z;
    shift_x /= -ori_size;
    shift_y /= -ori_size;
    if (fabs(shift_x) < ROME_EQUAL_ACCURACY && fabs(shift_y) < ROME_EQUAL_ACCURACY)
    {
        for (int i = 0; i < inout_Fsize2; i++) {
            fout_real[i] = fin_real[i];
            fout_imag[i] = fin_imag[i];
        }
        return;
    }
    for (int i=0; i<inout_Fsize; i++)
        for (int j=0; j<inout_Fsize; j++)
        {
            x = j;
            y = i;
            dotp = 2 * PI * (x * shift_x + y * shift_y);
            a = cos(dotp);
            b = sin(dotp);
            c = fin_real[i*inout_Fsize+j];
            d = fin_imag[i*inout_Fsize+j];
            ac = a * c;
            bd = b * d;
            ab_cd = (a + b) * (c + d);
            fout_real[i*inout_Fsize+j] = ac - bd;
            fout_imag[i*inout_Fsize+j] = ab_cd - ac - bd;
        }
    for (int i=inout_size-1; i>=inout_Fsize; i--)
    {
        y = i - inout_size;
        for (int j=0; j<inout_Fsize; j++)
        {
            x = j;
            dotp = 2 * PI * (x * shift_x + y * shift_y);
            a = cos(dotp);
            b = sin(dotp);
            c = fin_real[i*inout_Fsize+j];
            d = fin_imag[i*inout_Fsize+j];
            ac = a * c;
            bd = b * d;
            ab_cd = (a + b) * (c + d);
            fout_real[i*inout_Fsize+j] = ac - bd;
            fout_imag[i*inout_Fsize+j] = ab_cd - ac - bd;
        }
    }
}

// get spectrum for image
#define POWER_SPECTRUM 0
#define AMPLITUDE_SPECTRUM 1
void getSpectrum(const double* Min,int Min_size,double* spectrum,int spectrum_type);
void getSpectrum(const float* Min,int Min_size,float* spectrum,int spectrum_type);
// inverse 3x3 matrix
void inverse(double A[][3]);

void shiftImageInFourierTransformUnitTestCorrectness();

void shiftImageInFourierTransformUnitTestPerformance();

void testImageModule();

#endif /* defined(IMAGE_H_) */
