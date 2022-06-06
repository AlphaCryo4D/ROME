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
#include "resmap_util.h"		// used for building precompiled headers on Windows

#include "resmap_reconstruct.h"

void ProjectorBase::initialize(int _ori_size,int _current_size,int _ref_dim,
                             int _padding_factor/* = 2*/,int _interpolator/* = TRILINEAR*/,int _r_min_nn/* = 10*/)
{
    padding_factor = _padding_factor;
    interpolator = _interpolator;
    r_min_nn = _r_min_nn;
    ori_size = _ori_size;ori_Fsize = _ori_size/2+1;
    r_max = std::min(_current_size / 2, ori_size / 2);
    pad_size = 2 * (padding_factor * r_max + 1) + 1;
    pad_Fsize = pad_size/2 + 1;
    ref_dim = _ref_dim;
    //
    int dimz = ref_dim==3?pad_size:1;
    data.fini();
    data.init_zero(dimz, pad_size, pad_Fsize, true);
}

void ProjectorBase::griddingCorrect(Vol<FDOUBLE>& vol)
{
    assert(vol.dimx == vol.dimy);
    // assert(vol.dimx == ori_size);
    // Correct real-space map by dividing it by the Fourier transform of the interpolator(s)
    int vol_shift_z = -(int)((float) (vol.dimz) / 2.0);
    int vol_shift_y = -(int)((float) (vol.dimy) / 2.0);
    int vol_shift_x = -(int)((float) (vol.dimx) / 2.0);
    int x,y,z;
    
    
    // double ftblob = blob_Fourier_val(rval, blob) / blob_Fourier_val(0., blob);
    // Interpolation (goes with "interpolator") to go from arbitrary to fine grid
    if (interpolator==NEAREST_NEIGHBOUR && r_min_nn == 0)
    {
        for (int k = 0; k < vol.dimz; k++){
            z = k + vol_shift_z;
            for (int i = 0; i < vol.dimy; i++){
                y = i + vol_shift_y;
                auto vol_data = vol.wptr() + k*vol.dimy*vol.dimx + i*vol.dimx;
                for (int j = 0; j < vol.dimx; j++)
                {
                    x = j + vol_shift_x;
                    FDOUBLE r = sqrt((FDOUBLE)(x*x+y*y+z*z));
                    // if r==0: do nothing (i.e. divide by 1)
                    if (r > 0.)
                    {
                        // vol.dimx == vol.dimy
                        FDOUBLE rval = r / (vol.dimx * padding_factor);
                        FDOUBLE sinc = sin(rome_pi * rval) / ( rome_pi * rval);
                        // NN interpolation is convolution with a rectangular pulse, which FT is a sinc function
                        vol_data[j] /= sinc;
                    }// if r > 0.
                }// end loop j
            }// end loop i
        }// end loop k
    }
    else if (interpolator==TRILINEAR || (interpolator==NEAREST_NEIGHBOUR && r_min_nn > 0) )
    {
        for (int k = 0; k < vol.dimz; k++){
            z = k + vol_shift_z;
            for (int i = 0; i < vol.dimy; i++){
                y = i + vol_shift_y;
                auto vol_data = vol.wptr() + k*vol.dimy*vol.dimx + i*vol.dimx;
                if (z!=0 || y!=0) // only z=0 and x=0 the r maybe zero!!!
                {
                    #pragma ivdep
                    for (int j = 0; j < vol.dimx; j++)
                    {
                        x = j + vol_shift_x;
                        FDOUBLE r = sqrt((FDOUBLE)(x*x+y*y+z*z));
                        // if r==0: do nothing (i.e. divide by 1)
                        // vol.dimx == vol.dimy
                        FDOUBLE rval = r / (vol.dimx * padding_factor);
                        FDOUBLE sinc = sin(rome_pi * rval) / ( rome_pi * rval);
                        // trilinear interpolation is convolution with a triangular pulse, which FT is a sinc^2 function
                        vol_data[j] /= sinc * sinc;
                    }// end loop j
                }
				else
                {
                    for (int j = 0; j < vol.dimx; j++)
                    {
                        x = j + vol_shift_x;
                        FDOUBLE r = sqrt((FDOUBLE)(x*x+y*y+z*z));
                        // if r==0: do nothing (i.e. divide by 1)
                        if (r > 0.)
                        {
                            // vol.dimx == vol.dimy
                            FDOUBLE rval = r / (vol.dimx * padding_factor);
                            FDOUBLE sinc = sin(rome_pi * rval) / ( rome_pi * rval);
                            // trilinear interpolation is convolution with a triangular pulse, which FT is a sinc^2 function
                            vol_data[j] /= sinc * sinc;
                        }// if r > 0.
                    }// end loop j
                }
            }// end loop i
        }// end loop k
    }
}

// Fill data array with oversampled Fourier transform, and calculate its power spectrum
void MyProjector::computeFourierTransformMap(Vol<FDOUBLE>& ref_in,int _ori_size,int _current_size,FDOUBLE* power_spectrum,bool do_gridding/* = true*/)
{
    assert(ref_in.dimx==ref_in.dimy);
    assert(ref_in.dimx==_ori_size);
    // initialize
    int _ref_dim = ref_in.dimz!=1?3:2;
    initialize(_ori_size, _current_size, _ref_dim);
    
    std::vector<FDOUBLE> counter(ori_Fsize);
    
    int padoridim = padding_factor * ori_size;
    int padForidim = padoridim/2 + 1;
    FDOUBLE normfft;
    
    for (int i = 0; i < ori_Fsize; i++)
        power_spectrum[i] = counter[i] = 0;
    
    if (ref_dim == 3)
        normfft = (FDOUBLE)(padding_factor * padding_factor*padding_factor*ori_size);
    else
        normfft = (FDOUBLE)(padding_factor * padding_factor);
    
#ifdef DEBUG_CLASSPROJECTOR
    // write out ref_in
#endif
    
    if (do_gridding)
        griddingCorrect(ref_in);
    
#ifdef DEBUG_CLASSPROJECTOR
    // write out ref_in
#endif

    Vol<FDOUBLE> Mpad,MpadCenter;
    Mpad      .init_zero(_ref_dim==3?padoridim:1, padoridim, padoridim);
    MpadCenter.init_nonzero(_ref_dim==3?padoridim:1, padoridim, padoridim);
    //
    int ori_size_shift = -(int)((float) (ori_size) / 2.0);
    int padoridim_shift = -(int)((float) (padoridim) / 2.0);
    for (int k = 0; k < ref_in.dimz; k++) {
        int z = ref_dim==3?(k + ori_size_shift - padoridim_shift):0;
        for (int i = 0; i < ref_in.dimy; i++) {
            int y = i + ori_size_shift - padoridim_shift;
            #pragma ivdep
            for (int j = 0; j < ref_in.dimx; j++) {
                int x = j + ori_size_shift - padoridim_shift;
                ACCESS(Mpad, z, y, x) = ACCESS(ref_in, k, i, j);
            }
        }
    }
    
    // Translate padded map to put origin of FT in the center
    if(ref_dim==3)
        centerFFT3DForward(Mpad.rptr(),MpadCenter.wptr(),padoridim);
    else
    	centerFFT2D(Mpad.rptr(),MpadCenter.wptr(),padoridim,true);
    
    Mpad.fini();
    Vol<MKL_Complex> Faux;
    Faux.init_nonzero(_ref_dim==3?padoridim:1, padoridim, padoridim/2+1);
    //
    FourierTransformerBase transformer;
    // assert(omp_get_num_threads()==1);// not in parallel region
    transformer.init(MpadCenter.wptr(), (FourierComplex*)Faux.wptr(), padoridim, padoridim, (ref_dim==3?padoridim:1), 1/*omp_get_max_threads()*/);
    // also normalize
    transformer.FourierTransform(true);
    
    MpadCenter.fini();
    // By default r_max is half ori_size
    // int r_max = (current_size < 0)?(ori_size/2):(current_size/2);
    // Never allow r_max beyond Nyquist...
    // r_max = std::min(r_max, ori_size / 2);
    int max_r2 = r_max * r_max * padding_factor * padding_factor;
    int refPad_size_shift = -(int)((float) (pad_size) / 2.0);
    // This will also work for 2D
    for (int k = 0; k < Faux.dimz; k++) {
        int kp = ref_dim==3?((k < padForidim) ? k : k - padoridim) : 0;
        for (int i = 0; i < Faux.dimy; i++) {
            int ip = (i < padForidim) ? i : i - padoridim;
            int kp_ip_2 = kp*kp + ip*ip;
            for (int j = 0; j < Faux.dimx; j++)
            {
                int jp = j;
                int r2 = kp_ip_2 + jp*jp;
                // The Fourier Transforms are all "normalised" for 2D transforms of size = ori_size x ori_size
                if (r2 <= max_r2)
                {
                    // Set data array
                    int z = ref_dim==3?(kp - refPad_size_shift) : 0;
                    int y = ip - refPad_size_shift;
                    int x = jp;
                    
                    ACCESS(data,z,y,x).real = ACCESS(Faux,k,i,j).real * normfft;
                    ACCESS(data,z,y,x).imag = ACCESS(Faux,k,i,j).imag * normfft;
                    
                    // Calculate power spectrum
                    int ires = round( sqrt((FDOUBLE)r2) / padding_factor );
                    // Factor two because of two-dimensionality of the complex plane
                    power_spectrum[ires] += (ACCESS(data,z,y,x).real*ACCESS(data,z,y,x).real + ACCESS(data,z,y,x).imag*ACCESS(data,z,y,x).imag) / 2.;
                    counter[ires] += 1.;
                }
            }
        }
    }
    Faux.fini();
#ifdef DEBUG_CLASSPROJECTOR
    // FFT back to real space
    // write out refPad
#endif
    // Calculate radial average of power spectrum
    for (int i=0; i < ori_Fsize; i++)
    {
        if (counter[i] < 1.)
            power_spectrum[i] = 0.;
        else
            power_spectrum[i] /= counter[i];
    }
    
}


void MyProjector::project(FDOUBLE* f2d_real,FDOUBLE* f2d_imag,int f2d_size, const FDOUBLE A[][3],bool inv)
{
    // TODO: it need to test other padding_factor,interpolator,r_min_nn can be work right!
    
    FDOUBLE fx, fy, fz, xp, yp,zp;
    int x0, x1, y0, y1, z0, z1, y, y2, r2;
    bool is_neg_x;
    
#define COMPLEX(N) FDOUBLE N##_real,N##_imag;
    
    COMPLEX(d000) COMPLEX(d001) COMPLEX(d010)
    COMPLEX(d011) COMPLEX(d100) COMPLEX(d101)
    COMPLEX(d110) COMPLEX(d111)
    
    COMPLEX(dx00) COMPLEX(dx01) COMPLEX(dx10)
    COMPLEX(dx11) COMPLEX(dxy0) COMPLEX(dxy1)
    
#undef COMPLEX
    
    FDOUBLE Ainv[3][3];
    
    assert(f2d_size/2 <= r_max);
    int my_r_max = std::min(r_max,f2d_size/2);
    // Go from the 2D slice coordinates to the map coordinates
    int max_r2 = my_r_max * my_r_max;
    int min_r2_nn = r_min_nn * r_min_nn;
    
#ifdef DEBUG_PROJECT
    std::cout<<" r_max = "<<r_max<<std::endl;
    std::cout<<" r_min_nn = "<<r_min_nn<<std::endl;
    std::cout<<" max_r2 = "<<max_r2<<std::endl;
    std::cout<<" min_r2_nn = "<<min_r2_nn<<std::endl;
    std::cout<<" pad_Fsize = "<<pad_Fsize<<std::endl;
    std::cout<<" Fref_size = "<<Fref_size<<std::endl;
#endif
    
    int f2d_Fsize = f2d_size/2 + 1;
    int pad_size_shift = -(int)((float)(pad_size)/2.0);
    for (int i = 0; i < f2d_size*f2d_Fsize; i++)
        f2d_real[i] = f2d_imag[i] = 0;
    
    // Use the inverse matrix
    if (inv){
        for (int i = 0; i < 3; i++)
            for (int j = 0; j < 3; j++)
                Ainv[i][j] = A[i][j]*(FDOUBLE)padding_factor;
    }
    else{
        for (int i = 0; i < 3; i++)
            for (int j = 0; j < 3; j++)
                Ainv[i][j] = A[j][i]*(FDOUBLE)padding_factor;
    }
    
    for (int i = 0; i < f2d_size; i++)
    {
        // Don't search beyond square with side max_r
        if (i <= my_r_max)
        {
            y = i;
        }
        else if (i >= f2d_size - my_r_max)
        {
            y = i - f2d_size;
        }
        else
            continue;
        
        y2 = y * y;
        
        // TODO,remove if branch in for loop.
        // int my_r_min = sqrt(max_r2-y2)
        
        for (int x = 0; x <= my_r_max; x++)
        {
            // Only include points with radius < max_r (exclude points outside circle in square)
            r2 = x * x + y2;
            if (r2 > max_r2)
                continue;
            
            // Get logical coordinates in the 3D map
            xp = Ainv[0][0] * x + Ainv[0][1] * y;
            yp = Ainv[1][0] * x + Ainv[1][1] * y;
            zp = Ainv[2][0] * x + Ainv[2][1] * y;
            
#ifdef DEBUG_PROJECT
            std::cout<<xp<<" "<<yp<<" "<<zp<<std::endl;
#endif
            if (interpolator == TRILINEAR || r2 < min_r2_nn)
            {
                // Only asymmetric half is stored
                if (xp < 0)
                {
                    // Get complex conjugated hermitian symmetry pair
                    xp = -xp;
                    yp = -yp;
                    zp = -zp;
                    is_neg_x = true;
                }
                else
                {
                    is_neg_x = false;
                }
                
                // Trilinear interpolation (with physical coords)
                // Subtract STARTINGY to accelerate access to data (STARTINGX=0)
                x0 = floor(xp);
                fx = xp - x0;
                x1 = x0 + 1;
                
                y0 = floor(yp);
                fy = yp - y0;
                y0 -= pad_size_shift;
                y1 = y0 + 1;
                
                z0 = floor(zp);
                fz = zp - z0;
                z0 -= pad_size_shift;
                z1 = z0 + 1;
                
#ifdef DEBUG_PROJECT
                std::cout<<x0<<" "<<x1<<" "<<y0<<" "<<y1<<" "<<z0<<" "<<z1<<std::endl;
#endif
                // Matrix access can be accelerated through pre-calculation of z0*xydim etc.
                const auto data000 = &data(z0, y0, x0).real;
                d000_real = data000[0];d000_imag = data000[1];
                d001_real = data000[2];d001_imag = data000[3];
                
				const auto data010 = &data(z0, y1, x0).real;
                d010_real = data010[0];d010_imag = data010[1];
                d011_real = data010[2];d011_imag = data010[3];
                
				const auto data100 = &data(z1, y0, x0).real;
                d100_real = data100[0];d100_imag = data100[1];
                d101_real = data100[2];d101_imag = data100[3];
                
                const auto data110 = &data(z1,y1,x0).real;
                d110_real = data110[0];d110_imag = data110[1];
                d111_real = data110[2];d111_imag = data110[3];
                
                // Set the interpolated value in the 2D output array
                
                dx00_real = lin_interp(fx, d000_real, d001_real);
                dx00_imag = lin_interp(fx, d000_imag, d001_imag);
                dx01_real = lin_interp(fx, d100_real, d101_real);
                dx01_imag = lin_interp(fx, d100_imag, d101_imag);
                dx10_real = lin_interp(fx, d010_real, d011_real);
                dx10_imag = lin_interp(fx, d010_imag, d011_imag);
                dx11_real = lin_interp(fx, d110_real, d111_real);
                dx11_imag = lin_interp(fx, d110_imag, d111_imag);
                
                dxy0_real = lin_interp(fy, dx00_real, dx10_real);
                dxy0_imag = lin_interp(fy, dx00_imag, dx10_imag);
                dxy1_real = lin_interp(fy, dx01_real, dx11_real);
                dxy1_imag = lin_interp(fy, dx01_imag, dx11_imag);
                
                f2d_real[i*f2d_Fsize+x] = lin_interp(fz, dxy0_real, dxy1_real);
                f2d_imag[i*f2d_Fsize+x] = lin_interp(fz, dxy0_imag, dxy1_imag);
                //std::cout<<Fref_real[i*Fref_Fsize+x]<<" "<<Fref_imag[i*Fref_Fsize+x]<<std::endl;
                // Take complex conjugated for half with negative x
                if (is_neg_x){
                    f2d_real[i*f2d_Fsize+x] = f2d_real[i*f2d_Fsize+x];
                    f2d_imag[i*f2d_Fsize+x] = -1.0*f2d_imag[i*f2d_Fsize+x];
                }
            } // endif TRILINEAR
            else if (interpolator == NEAREST_NEIGHBOUR )
            {
                x0 = round(xp);
                y0 = round(yp);
                z0 = round(zp);
                
                if (x0 < 0){
                    f2d_real[i*f2d_Fsize+x] = data(-y0-pad_size_shift,-x0,-z0-pad_size_shift).real;
                    f2d_imag[i*f2d_Fsize+x] = -1*data(-y0-pad_size_shift,-x0,-z0-pad_size_shift).imag;
                }
                else{
                    f2d_real[i*f2d_Fsize+x] = data(y0-pad_size_shift,x0,z0-pad_size_shift).real;
                    f2d_imag[i*f2d_Fsize+x] = data(y0-pad_size_shift,x0,z0-pad_size_shift).imag;
                }
            } // endif NEAREST_NEIGHBOUR
        } // endif x-loop
    } // endif y-loop
#ifdef DEBUG_CLASSPROJECTOR
    // write out 2D ref
#endif
}

void MyProjector::projectOneTile(FDOUBLE* f2d_real,FDOUBLE* f2d_imag,int n_start,int n_end,int f2d_size,const FDOUBLE A[][3],bool inv)
{
    // TODO: it need to test other padding_factor,interpolator,r_min_nn can be work right!
    assert(0<=n_start);assert(n_start<=n_end);
    n_end = std::min(n_end, f2d_size*(f2d_size/2+1));
    assert(n_end<=f2d_size*(f2d_size/2+1));
    FDOUBLE fx, fy, fz, xp, yp,zp;
    int x0, x1, y0, y1, z0, z1, y, y2, r2;
    bool is_neg_x;
    
#define COMPLEX(N) FDOUBLE N##_real,N##_imag;
    
    COMPLEX(d000) COMPLEX(d001) COMPLEX(d010)
    COMPLEX(d011) COMPLEX(d100) COMPLEX(d101)
    COMPLEX(d110) COMPLEX(d111)
    
    COMPLEX(dx00) COMPLEX(dx01) COMPLEX(dx10)
    COMPLEX(dx11) COMPLEX(dxy0) COMPLEX(dxy1)
    
#undef COMPLEX
    
    FDOUBLE Ainv[3][3];
    
    assert(f2d_size/2 <= r_max);
    int my_r_max = std::min(r_max,f2d_size/2);
    // Go from the 2D slice coordinates to the map coordinates
    int max_r2 = my_r_max * my_r_max;
    int min_r2_nn = r_min_nn * r_min_nn;
    
#ifdef DEBUG_PROJECT
    std::cout<<" r_max = "<<r_max<<std::endl;
    std::cout<<" r_min_nn = "<<r_min_nn<<std::endl;
    std::cout<<" max_r2 = "<<max_r2<<std::endl;
    std::cout<<" min_r2_nn = "<<min_r2_nn<<std::endl;
    std::cout<<" pad_Fsize = "<<pad_Fsize<<std::endl;
    std::cout<<" Fref_size = "<<Fref_size<<std::endl;
#endif
    
    int f2d_Fsize = f2d_size/2 + 1;
    //int Fref_pad_size_shift = -(int)((float) (fref_pad_size) / 2.0);
    int pad_size_shift = -(int)((float)(pad_size)/2.0);
    // for (int i = 0; i < f2d_size*f2d_Fsize; i++)
    //     f2d_real[i] = f2d_imag[i] = 0;
    
    // Use the inverse matrix
    if (inv){
        for (int i = 0; i < 3; i++)
            for (int j = 0; j < 3; j++)
                Ainv[i][j] = A[i][j]*(FDOUBLE)padding_factor;
    }
    else{
        for (int i = 0; i < 3; i++)
            for (int j = 0; j < 3; j++)
                Ainv[i][j] = A[j][i]*(FDOUBLE)padding_factor;
    }
    
	// Crudely approximate the accessed cells.  Crude for f2d. Even more so for data since doesn't include the transform
	// but since our intent is to do ones that use the same cells, it will do for the current analysis...
	L2CacheModel::seqAcc("projectOneTile f2d_real", -1, &f2d_real[0], n_end - n_start + 1, sizeof(f2d_real[0]));
	L2CacheModel::seqAcc("projectOneTile f2d_imag", -1, &f2d_imag[0], n_end - n_start + 1, sizeof(f2d_imag[0]));
	L2CacheModel::seqAcc("projectOneTile data",     -1, &data(0, 0, 0), size_t(my_r_max*my_r_max*3.1415), sizeof(data(0, 0, 0)));

	// Do it
	int i_start = n_start/(f2d_size/2+1);
    int i_end = (n_end-1)/(f2d_size/2+1);
    int x_start,x_end;
    // initialize to zero
    for (int i = i_start; i <= i_end; i++)
    {
        if (i == i_start) x_start = n_start%(f2d_size/2+1);
        else x_start = 0;
        if (i == i_end) x_end = (n_end-1)%(f2d_size/2+1);
        else x_end = my_r_max;
        for (int x = x_start; x <= x_end; x++)
            f2d_real[i*f2d_Fsize+x] = f2d_imag[i*f2d_Fsize+x] = 0.;
    }
    
    for (int i = i_start; i <= i_end; i++)
    {
        // Don't search beyond square with side max_r
        if (i <= my_r_max)
        {
            y = i;
        }
        else if (i >= f2d_size - my_r_max)
        {
            y = i - f2d_size;
        }
        else
            continue;
        
        y2 = y * y;
        
        // TODO,remove if branch in for loop.
        // int my_r_min = sqrt(max_r2-y2)
        if (i == i_start) x_start = n_start%(f2d_size/2+1);
        else x_start = 0;
        if (i == i_end) x_end = (n_end-1)%(f2d_size/2+1);
        else x_end = my_r_max;
        for (int x = x_start; x <= x_end; x++)
        {
            // Only include points with radius < max_r (exclude points outside circle in square)
            r2 = x * x + y2;
            if (r2 > max_r2)
                continue;
            
            // Get logical coordinates in the 3D map
            xp = Ainv[0][0] * x + Ainv[0][1] * y;
            yp = Ainv[1][0] * x + Ainv[1][1] * y;
            zp = Ainv[2][0] * x + Ainv[2][1] * y;
            
#ifdef DEBUG_PROJECT
            std::cout<<xp<<" "<<yp<<" "<<zp<<std::endl;
#endif
            if (interpolator == TRILINEAR || r2 < min_r2_nn)
            {
                // Only asymmetric half is stored
                if (xp < 0)
                {
                    // Get complex conjugated hermitian symmetry pair
                    xp = -xp;
                    yp = -yp;
                    zp = -zp;
                    is_neg_x = true;
                }
                else
                {
                    is_neg_x = false;
                }
                
                // Trilinear interpolation (with physical coords)
                // Subtract STARTINGY to accelerate access to data (STARTINGX=0)
                x0 = floor(xp);
                fx = xp - x0;
                x1 = x0 + 1;
                
                y0 = floor(yp);
                fy = yp - y0;
                y0 -= pad_size_shift;
                y1 = y0 + 1;
                
                z0 = floor(zp);
                fz = zp - z0;
                z0 -= pad_size_shift;
                z1 = z0 + 1;
                
#ifdef DEBUG_PROJECT
                std::cout<<x0<<" "<<x1<<" "<<y0<<" "<<y1<<" "<<z0<<" "<<z1<<std::endl;
#endif
                // Matrix access can be accelerated through pre-calculation of z0*xydim etc.
                const auto data000 = &data(z0, y0, x0).real;
                d000_real = data000[0];d000_imag = data000[1];
                d001_real = data000[2];d001_imag = data000[3];
                
                const auto data010 = &data(z0, y1, x0).real;
                d010_real = data010[0];d010_imag = data010[1];
                d011_real = data010[2];d011_imag = data010[3];
                
                const auto data100 = &data(z1, y0, x0).real;
                d100_real = data100[0];d100_imag = data100[1];
                d101_real = data100[2];d101_imag = data100[3];
                
                const auto data110 = &data(z1,y1,x0).real;
                d110_real = data110[0];d110_imag = data110[1];
                d111_real = data110[2];d111_imag = data110[3];
                
                // Set the interpolated value in the 2D output array
                
                dx00_real = lin_interp(fx, d000_real, d001_real);
                dx00_imag = lin_interp(fx, d000_imag, d001_imag);
                dx01_real = lin_interp(fx, d100_real, d101_real);
                dx01_imag = lin_interp(fx, d100_imag, d101_imag);
                dx10_real = lin_interp(fx, d010_real, d011_real);
                dx10_imag = lin_interp(fx, d010_imag, d011_imag);
                dx11_real = lin_interp(fx, d110_real, d111_real);
                dx11_imag = lin_interp(fx, d110_imag, d111_imag);
                
                dxy0_real = lin_interp(fy, dx00_real, dx10_real);
                dxy0_imag = lin_interp(fy, dx00_imag, dx10_imag);
                dxy1_real = lin_interp(fy, dx01_real, dx11_real);
                dxy1_imag = lin_interp(fy, dx01_imag, dx11_imag);
                
                f2d_real[i*f2d_Fsize+x] = lin_interp(fz, dxy0_real, dxy1_real);
                f2d_imag[i*f2d_Fsize+x] = lin_interp(fz, dxy0_imag, dxy1_imag);
                //std::cout<<Fref_real[i*Fref_Fsize+x]<<" "<<Fref_imag[i*Fref_Fsize+x]<<std::endl;
                // Take complex conjugated for half with negative x
                if (is_neg_x){
                    f2d_real[i*f2d_Fsize+x] = f2d_real[i*f2d_Fsize+x];
                    f2d_imag[i*f2d_Fsize+x] = -1.0*f2d_imag[i*f2d_Fsize+x];
                }
            } // endif TRILINEAR
            else if (interpolator == NEAREST_NEIGHBOUR )
            {
                x0 = round(xp);
                y0 = round(yp);
                z0 = round(zp);
                
                if (x0 < 0){
                    f2d_real[i*f2d_Fsize+x] = data(-y0-pad_size_shift,-x0,-z0-pad_size_shift).real;
                    f2d_imag[i*f2d_Fsize+x] = -1*data(-y0-pad_size_shift,-x0,-z0-pad_size_shift).imag;
                }
                else{
                    f2d_real[i*f2d_Fsize+x] = data(y0-pad_size_shift,x0,z0-pad_size_shift).real;
                    f2d_imag[i*f2d_Fsize+x] = data(y0-pad_size_shift,x0,z0-pad_size_shift).imag;
                }
            } // endif NEAREST_NEIGHBOUR
        } // endif x-loop
    } // endif y-loop
#ifdef DEBUG_CLASSPROJECTOR
    // write out 2D ref
#endif
}

void MyProjector::projectOneTileByShell(FDOUBLE* f2d_real,FDOUBLE* f2d_imag,int shell_n_start,int shell_n_end,
                                        int f2d_size,const FDOUBLE A[][3],bool inv,const int* nIndex)
{
    // TODO: it need to test other padding_factor,interpolator,r_min_nn can be work right!
    assert(0<=shell_n_start);assert(shell_n_start<=shell_n_end);
    assert(shell_n_end<=f2d_size*(f2d_size/2+1));
    FDOUBLE fx, fy, fz, xp, yp,zp;
    int x0, x1, y0, y1, z0, z1, y, y2, r2;
    bool is_neg_x;
    
#define COMPLEX(N) FDOUBLE N##_real,N##_imag;
    
    COMPLEX(d000) COMPLEX(d001) COMPLEX(d010)
    COMPLEX(d011) COMPLEX(d100) COMPLEX(d101)
    COMPLEX(d110) COMPLEX(d111)
    
    COMPLEX(dx00) COMPLEX(dx01) COMPLEX(dx10)
    COMPLEX(dx11) COMPLEX(dxy0) COMPLEX(dxy1)
    
#undef COMPLEX
    
    FDOUBLE Ainv[3][3];
    
    assert(f2d_size/2 <= r_max);
    int my_r_max = std::min(r_max,f2d_size/2);
    // Go from the 2D slice coordinates to the map coordinates
    int max_r2 = my_r_max * my_r_max;
    int min_r2_nn = r_min_nn * r_min_nn;
    
#ifdef DEBUG_PROJECT
    std::cout<<" r_max = "<<r_max<<std::endl;
    std::cout<<" r_min_nn = "<<r_min_nn<<std::endl;
    std::cout<<" max_r2 = "<<max_r2<<std::endl;
    std::cout<<" min_r2_nn = "<<min_r2_nn<<std::endl;
    std::cout<<" pad_Fsize = "<<pad_Fsize<<std::endl;
    std::cout<<" Fref_size = "<<Fref_size<<std::endl;
#endif
    
    int f2d_Fsize = f2d_size/2 + 1;
    //int Fref_pad_size_shift = -(int)((float) (fref_pad_size) / 2.0);
    int pad_size_shift = -(int)((float)(pad_size)/2.0);
    // for (int i = 0; i < f2d_size*f2d_Fsize; i++)
    //     f2d_real[i] = f2d_imag[i] = 0;
    
    // Use the inverse matrix
    if (inv){
        for (int i = 0; i < 3; i++)
            for (int j = 0; j < 3; j++)
                Ainv[i][j] = A[i][j]*(FDOUBLE)padding_factor;
    }
    else{
        for (int i = 0; i < 3; i++)
            for (int j = 0; j < 3; j++)
                Ainv[i][j] = A[j][i]*(FDOUBLE)padding_factor;
    }

	// Crudely approximate the accessed cells.  Ok for f2d but for data doesn't include the transform
	// but since our intent is to do ones that use the same cells, it will do for the current analysis...
	L2CacheModel::seqAcc("projectOneTileByShell f2d_real", -1, &f2d_real[shell_n_start], shell_n_end - shell_n_start,	 sizeof(f2d_real[shell_n_start]));
	L2CacheModel::seqAcc("projectOneTileByShell f2d_imag", -1, &f2d_imag[shell_n_start], shell_n_end - shell_n_start,	 sizeof(f2d_imag[shell_n_start]));
	L2CacheModel::seqAcc("projectOneTileByShell data",     -1, &data(0, 0, 0),		   size_t(my_r_max*my_r_max*3.1415), sizeof(data(0, 0, 0)));

    // initialize to zero
    for (int shell_n = shell_n_start; shell_n < shell_n_end; shell_n++) {
        f2d_real[shell_n] = f2d_imag[shell_n] = 0.;
    }
    
    //
    for (int shell_n = shell_n_start; shell_n < shell_n_end; shell_n++)
    {
        int n = nIndex[shell_n];
        int i = n / (f2d_size/2+1);
        // Don't search beyond square with side max_r
        if (i <= my_r_max)
        {
            y = i;
        }
        else if (i >= f2d_size - my_r_max)
        {
            y = i - f2d_size;
        }
        else
            continue;
        //
        int x = n % (f2d_size/2+1);
        y2 = y * y;
        // Only include points with radius < max_r (exclude points outside circle in square)
        r2 = x * x + y2;
        if (r2 > max_r2) // always not meet??
            continue;
            
        // Get logical coordinates in the 3D map
        xp = Ainv[0][0] * x + Ainv[0][1] * y;
        yp = Ainv[1][0] * x + Ainv[1][1] * y;
        zp = Ainv[2][0] * x + Ainv[2][1] * y;
        
#ifdef DEBUG_PROJECT
        std::cout<<xp<<" "<<yp<<" "<<zp<<std::endl;
#endif
        if (interpolator == TRILINEAR || r2 < min_r2_nn)
        {
            // Only asymmetric half is stored
            if (xp < 0)
            {
                // Get complex conjugated hermitian symmetry pair
                xp = -xp;
                yp = -yp;
                zp = -zp;
                is_neg_x = true;
            }
            else
            {
                is_neg_x = false;
            }
            
            // Trilinear interpolation (with physical coords)
            // Subtract STARTINGY to accelerate access to data (STARTINGX=0)
            x0 = floor(xp);
            fx = xp - x0;
            x1 = x0 + 1;
            
            y0 = floor(yp);
            fy = yp - y0;
            y0 -= pad_size_shift;
            y1 = y0 + 1;
            
            z0 = floor(zp);
            fz = zp - z0;
            z0 -= pad_size_shift;
            z1 = z0 + 1;
            
#ifdef DEBUG_PROJECT
            std::cout<<x0<<" "<<x1<<" "<<y0<<" "<<y1<<" "<<z0<<" "<<z1<<std::endl;
#endif
            // Matrix access can be accelerated through pre-calculation of z0*xydim etc.
            const auto data000 = &data(z0, y0, x0).real;
            d000_real = data000[0];d000_imag = data000[1];
            d001_real = data000[2];d001_imag = data000[3];
            
            const auto data010 = &data(z0, y1, x0).real;
            d010_real = data010[0];d010_imag = data010[1];
            d011_real = data010[2];d011_imag = data010[3];
            
            const auto data100 = &data(z1, y0, x0).real;
            d100_real = data100[0];d100_imag = data100[1];
            d101_real = data100[2];d101_imag = data100[3];
            
            const auto data110 = &data(z1,y1,x0).real;
            d110_real = data110[0];d110_imag = data110[1];
            d111_real = data110[2];d111_imag = data110[3];
            
            // Set the interpolated value in the 2D output array
            
            dx00_real = lin_interp(fx, d000_real, d001_real);
            dx00_imag = lin_interp(fx, d000_imag, d001_imag);
            dx01_real = lin_interp(fx, d100_real, d101_real);
            dx01_imag = lin_interp(fx, d100_imag, d101_imag);
            dx10_real = lin_interp(fx, d010_real, d011_real);
            dx10_imag = lin_interp(fx, d010_imag, d011_imag);
            dx11_real = lin_interp(fx, d110_real, d111_real);
            dx11_imag = lin_interp(fx, d110_imag, d111_imag);
            
            dxy0_real = lin_interp(fy, dx00_real, dx10_real);
            dxy0_imag = lin_interp(fy, dx00_imag, dx10_imag);
            dxy1_real = lin_interp(fy, dx01_real, dx11_real);
            dxy1_imag = lin_interp(fy, dx01_imag, dx11_imag);
            
            f2d_real[shell_n] = lin_interp(fz, dxy0_real, dxy1_real);
            f2d_imag[shell_n] = lin_interp(fz, dxy0_imag, dxy1_imag);
            //std::cout<<Fref_real[i*Fref_Fsize+x]<<" "<<Fref_imag[i*Fref_Fsize+x]<<std::endl;
            // Take complex conjugated for half with negative x
            if (is_neg_x){
                f2d_real[shell_n] = f2d_real[shell_n];
                f2d_imag[shell_n] = -1.0*f2d_imag[shell_n];
            }
        } // endif TRILINEAR
        else if (interpolator == NEAREST_NEIGHBOUR )
        {
            x0 = round(xp);
            y0 = round(yp);
            z0 = round(zp);
            
            if (x0 < 0){
                f2d_real[shell_n] = data(-y0-pad_size_shift,-x0,-z0-pad_size_shift).real;
                f2d_imag[shell_n] = -1*data(-y0-pad_size_shift,-x0,-z0-pad_size_shift).imag;
            }
            else{
                f2d_real[shell_n] = data(y0-pad_size_shift,x0,z0-pad_size_shift).real;
                f2d_imag[shell_n] = data(y0-pad_size_shift,x0,z0-pad_size_shift).imag;
            }
        } // endif NEAREST_NEIGHBOUR
    } // endif shell_n
#ifdef DEBUG_CLASSPROJECTOR
    // write out 2D ref
#endif
}


void MyProjector::rotate2D(FDOUBLE* f2d_real,FDOUBLE* f2d_imag,int f2d_size, const FDOUBLE A[][3],bool inv)
{
    // TODO: it need to test other padding_factor,interpolator,r_min_nn can be work right!
    
    FDOUBLE fx, fy, xp, yp;
    int x0, x1, y0, y1, y, y2, r2;
    bool is_neg_x;
    
#define COMPLEX(N) FDOUBLE N##_real,N##_imag;
    
    COMPLEX(d00) COMPLEX(d01)
    COMPLEX(d10) COMPLEX(d11)
    
    COMPLEX(dx0) COMPLEX(dx1)
    
#undef COMPLEX
    
    FDOUBLE Ainv[3][3];
    // The f2d image may be smaller than r_max, in that case also make sure not to fill the corners!
    assert(f2d_size/2 <= r_max);
    int my_r_max = std::min(r_max,f2d_size/2);
    // Go from the 2D slice coordinates to the map coordinates
    int max_r2 = my_r_max * my_r_max;
    int min_r2_nn = r_min_nn * r_min_nn;
    
    int f2d_Fsize = f2d_size/2 + 1;
    int pad_size_shift = -(int)((float)(pad_size)/2.0);
    for (int i = 0; i < f2d_size*f2d_Fsize; i++)
        f2d_real[i] = f2d_imag[i] = 0;
    
    // Use the inverse matrix
    if (inv){
        for (int i = 0; i < 3; i++)
            for (int j = 0; j < 3; j++)
                Ainv[i][j] = A[i][j]*(FDOUBLE)padding_factor;
    }
    else{
        for (int i = 0; i < 3; i++)
            for (int j = 0; j < 3; j++)
                Ainv[i][j] = A[j][i]*(FDOUBLE)padding_factor;
    }
    
    for (int i=0; i < f2d_size; i++)
    {
        // Don't search beyond square with side max_r
        if (i <= my_r_max)
        {
            y = i;
        }
        else if (i >= f2d_size - my_r_max)
        {
            y = i - f2d_size;
        }
        else
            continue;
        
        y2 = y * y;
        
        for (int x=0; x <= my_r_max; x++)
        {
            // Only include points with radius < max_r (exclude points outside circle in square)
            r2 = x * x + y2;
            if (r2 > max_r2)
                continue;
            
            // Get logical coordinates in the 3D map
            xp = Ainv[0][0] * x + Ainv[0][1] * y;
            yp = Ainv[1][0] * x + Ainv[1][1] * y;
            if (interpolator == TRILINEAR || r2 < min_r2_nn)
            {
                // Only asymmetric half is stored
                if (xp < 0)
                {
                    // Get complex conjugated hermitian symmetry pair
                    xp = -xp;
                    yp = -yp;
                    is_neg_x = true;
                }
                else
                {
                    is_neg_x = false;
                }
                
                // Trilinear interpolation (with physical coords)
                // Subtract STARTINGY to accelerate access to data (STARTINGX=0)
                x0 = floor(xp);
                fx = xp - x0;
                x1 = x0 + 1;
                
                y0 = floor(yp);
                fy = yp - y0;
                y0 -= pad_size_shift;
                y1 = y0 + 1;
                
                // Matrix access can be accelerated through pre-calculation of z0*xydim etc.
                const auto data00 = &data(0, y0, x0).real;
                d00_real = data00[0]; d00_imag = data00[1];
                d01_real = data00[2]; d01_imag = data00[3];
                
                const auto data10 = &data(0, y1, x0).real;
                d10_real = data10[0]; d10_imag = data10[1];
                d11_real = data10[2]; d11_imag = data10[3];
                
                // Set the interpolated value in the 2D output array
                dx0_real = LIN_INTERP(fx, d00_real, d01_real);
                dx0_imag = LIN_INTERP(fx, d00_imag, d01_imag);
                dx1_real = LIN_INTERP(fx, d10_real, d11_real);
                dx1_imag = LIN_INTERP(fx, d10_imag, d11_imag);
                
                f2d_real[i*f2d_Fsize+x] = LIN_INTERP(fy, dx0_real, dx1_real);
                f2d_imag[i*f2d_Fsize+x] = LIN_INTERP(fy, dx0_imag, dx1_imag);
                
                // Take complex conjugated for half with negative x
                if (is_neg_x){
                    f2d_real[i*f2d_Fsize+x] = f2d_real[i*f2d_Fsize+x];
                    f2d_imag[i*f2d_Fsize+x] = -1.0*f2d_imag[i*f2d_Fsize+x];
                }
            } // endif TRILINEAR
            
            else if (interpolator == NEAREST_NEIGHBOUR )
            {
                x0 = round(xp);
                y0 = round(yp);
                if (x0 < 0){
                    f2d_real[i*f2d_Fsize+x] = data(0,-y0-pad_size_shift,-x0).real;
                    f2d_imag[i*f2d_Fsize+x] = -1*data(0,-y0-pad_size_shift,-x0).imag;
                }
                else{
                    f2d_real[i*f2d_Fsize+x] = data(0,y0-pad_size_shift,x0).real;
                    f2d_imag[i*f2d_Fsize+x] = data(0,y0-pad_size_shift,x0).imag;
                }
            } // endif NEAREST_NEIGHBOUR
            
        } // endif x-loop
    } // endif y-loop
}


void MyBackProjector::initialize(int _ori_size,int _current_size,int _ref_dim,std::string _fn_sym/* = "C1"*/,
                                 int _padding_factor/* = 2*/,int _interpolator/* = TRILINEAR*/,int _r_min_nn/* = 10*/,
                                 int _blob_order/* = 0*/, double _blob_radius/* = 1.9*/, double _blob_alpha/* = 15*/)
{
    ProjectorBase::initialize(_ori_size,_current_size,_ref_dim,_padding_factor,_interpolator,_r_min_nn);
    //
    int dimz = ref_dim==3?pad_size:1;
    weight.fini();
    weight.init_zero(dimz, pad_size, pad_Fsize, true);
    // Precalculate tabulated ftblob values
    tab_ftblob.init(_blob_radius * padding_factor, _blob_alpha, _blob_order, 10000);
    //
    fn_sym = _fn_sym;
}


void MyBackProjector::backproject(const FDOUBLE* f2d_real,const FDOUBLE* f2d_imag,int f2d_size,
                                  const FDOUBLE A[][3], bool inv,
                                  Vol<MKL_Complex>& data_td,Vol<FDOUBLE>& weight_td,
                                  const FDOUBLE* Mweight/* = nullptr*/,int z_start /*= 0*/,int z_end /*= 100000*/)
{
    FDOUBLE fx, fy, fz, mfx, mfy, mfz, xp, yp, zp;
    int first_x, x0, x1, y0, y1, z0, z1, y, y2, r2;
    bool is_neg_x;
    FDOUBLE dd000, dd001, dd010, dd011, dd100, dd101, dd110, dd111;
    FDOUBLE my_val_real,my_val_imag;
    FDOUBLE Ainv[3][3];
    FDOUBLE my_weight = 1.;
    
    int f2d_Fsize = f2d_size/2 + 1;
    int pad_size_shift = -(int)((float)(pad_size)/2.0);
    // f2d should already be in the right size (ori_size,orihalfdim)
    // AND the points outside max_r should already be zero...
    
    // Use the inverse matrix
    if (inv){
        for (int i = 0; i < 3; i++)
            for (int j = 0; j < 3; j++)
                Ainv[i][j] = A[i][j]*(FDOUBLE)padding_factor;
    }
    else{
        for (int i = 0; i < 3; i++)
            for (int j = 0; j < 3; j++)
                Ainv[i][j] = A[j][i]*(FDOUBLE)padding_factor;
    }
    
    // Go from the 2D slice coordinates to the 3D coordinates
    int max_r2 = r_max * r_max;
    int min_r2_nn = r_min_nn * r_min_nn;
    
    //#define DEBUG_BACKP
#ifdef DEBUG_BACKP
    //
#endif
    
    for (int i=0; i < f2d_size; i++)
    {
        // Dont search beyond square with side max_r
        if (i <= r_max)
        {
            y = i;
            first_x = 0;
        }
        else if (i >= f2d_size - r_max)
        {
            y = i - f2d_size;
            // x==0 plane is stored twice in the FFTW format. Dont set it twice in BACKPROJECTION!
            first_x = 1;
        }
        else
            continue;
        
        y2 = y * y;
        for (int x=first_x; x <= r_max; x++)
        {
            // Only include points with radius < max_r (exclude points outside circle in square)
            r2 = x * x + y2;
            if (r2 > max_r2)
                continue;
            
            // Get the relevant value in the input image
            my_val_real = f2d_real[i*f2d_Fsize+x];
            my_val_imag = f2d_imag[i*f2d_Fsize+x];
            // Get the weight
            if (Mweight != NULL)
                my_weight = Mweight[i*f2d_Fsize+x];
            // else: my_weight was already initialised to 1.
            
            if (my_weight > 0.)
            {
                
                // Get logical coordinates in the 3D map
                xp = Ainv[0][0] * x + Ainv[0][1] * y;
                yp = Ainv[1][0] * x + Ainv[1][1] * y;
                zp = Ainv[2][0] * x + Ainv[2][1] * y;
                
                if (interpolator == TRILINEAR || r2 < min_r2_nn)
                {
                    
                    // Only asymmetric half is stored
                    if (xp < 0)
                    {
                        // Get complex conjugated hermitian symmetry pair
                        xp = -xp;
                        yp = -yp;
                        zp = -zp;
                        is_neg_x = true;
                    }
                    else
                    {
                        is_neg_x = false;
                    }
                    
                    // Trilinear interpolation (with physical coords)
                    // Subtract STARTINGY and STARTINGZ to accelerate access to data (STARTINGX=0)
                    // In that way use DIRECT_A3D_ELEM, rather than A3D_ELEM
                    x0 = floor(xp);
                    fx = xp - x0;
                    x1 = x0 + 1;
                    
                    y0 = floor(yp);
                    fy = yp - y0;
                    y0 -=  pad_size_shift;
                    y1 = y0 + 1;
                    
                    z0 = floor(zp);
                    fz = zp - z0;
                    z0 -= pad_size_shift;
                    z1 = z0 + 1;
                    
                    // z_start <= zp < z_end
                    if (z0 < z_start || z0 >= z_end) continue;
                    
                    mfx = 1. - fx;
                    mfy = 1. - fy;
                    mfz = 1. - fz;
                    
                    dd000 = mfz * mfy * mfx;
                    dd001 = mfz * mfy *  fx;
                    dd010 = mfz *  fy * mfx;
                    dd011 = mfz *  fy *  fx;
                    dd100 =  fz * mfy * mfx;
                    dd101 =  fz * mfy *  fx;
                    dd110 =  fz *  fy * mfx;
                    dd111 =  fz *  fy *  fx;
                    
                    if (is_neg_x){
                        my_val_real = my_val_real;
                        my_val_imag = -my_val_imag;
                    }
                    // Store slice in 3D weighted sum
                    assert(x0>=0);assert(y0>=0);assert(z0>=0);
                    // TODO. merge (real,imag,weight) to better memory fetch
                    auto data_td000 = &data_td(z0, y0, x0).real;
                    data_td000[0] += dd000 * my_val_real;data_td000[1] += dd000 * my_val_imag;
                    data_td000[2] += dd001 * my_val_real;data_td000[3] += dd001 * my_val_imag;
                    auto data_td010 = &data_td(z0, y1, x0).real;
                    data_td010[0] += dd010 * my_val_real;data_td010[1] += dd010 * my_val_imag;
                    data_td010[2] += dd011 * my_val_real;data_td010[3] += dd011 * my_val_imag;
                    auto data_td100 = &data_td(z1, y0, x0).real;
                    data_td100[0] += dd100 * my_val_real;data_td100[1] += dd100 * my_val_imag;
                    data_td100[2] += dd101 * my_val_real;data_td100[3] += dd101 * my_val_imag;
                    auto data_td110 = &data_td(z1, y1, x0).real;
                    data_td110[0] += dd110 * my_val_real;data_td110[1] += dd110 * my_val_imag;
                    data_td110[2] += dd111 * my_val_real;data_td110[3] += dd111 * my_val_imag;
                    
                    auto weight_td000 = &weight_td(z0, y0, x0);
                    weight_td000[0] += dd000 * my_weight;weight_td000[1] += dd001 * my_weight;
                    auto weight_td010 = &weight_td(z0, y1, x0);
                    weight_td010[0] += dd010 * my_weight;weight_td010[1] += dd011 * my_weight;
                    auto weight_td100 = &weight_td(z1, y0, x0);
                    weight_td100[0] += dd100 * my_weight;weight_td100[1] += dd101 * my_weight;
                    auto weight_td110 = &weight_td(z1, y1, x0);
                    weight_td110[0] += dd110 * my_weight;weight_td110[1] += dd111 * my_weight;
                } // endif TRILINEAR
                else if (interpolator == NEAREST_NEIGHBOUR )
                {
                    
                    x0 = round(xp);
                    y0 = round(yp);
                    z0 = round(zp);
                    
                    if (x0 < 0)
                    {
                        ACCESS(data_td, -z0, -y0, -x0).real += my_val_real;
                        ACCESS(data_td, -z0, -y0, -x0).imag += (-my_val_imag);
                        ACCESS(weight_td, -z0, -y0, -x0) += my_weight;
                    }
                    else
                    {
                        ACCESS(data_td, z0, y0, x0).real += my_val_real;
                        ACCESS(data_td, z0, y0, x0).imag += my_val_imag;
                        ACCESS(weight_td, z0, y0, x0) += my_weight;
                    }
                    
                } // endif NEAREST_NEIGHBOUR
                else
                {
                    ERROR_REPORT("FourierInterpolator::backproject%%ERROR: unrecognized interpolator ");
                }
            } // endif weight>0.
        } // endif x-loop
    } // endif y-loop
}

void MyBackProjector::backprojectOneTileByShell(const FDOUBLE* f2d_real,const FDOUBLE* f2d_imag,int shell_n_start,int shell_n_end,
                                                int f2d_size,const FDOUBLE A[][3], bool inv,Vol<MKL_Complex>& data_td,
                                                Vol<FDOUBLE>& weight_td,const FDOUBLE* Mweight,const int* nIndex)
{
    FDOUBLE fx, fy, fz, mfx, mfy, mfz, xp, yp, zp;
    int first_x, x0, x1, y0, y1, z0, z1, y, y2, r2;
    bool is_neg_x;
    FDOUBLE dd000, dd001, dd010, dd011, dd100, dd101, dd110, dd111;
    FDOUBLE my_val_real,my_val_imag;
    FDOUBLE Ainv[3][3];
    FDOUBLE my_weight = 1.;
    
    int f2d_Fsize = f2d_size/2 + 1;
    int pad_size_shift = -(int)((float)(pad_size)/2.0);
    // f2d should already be in the right size (ori_size,orihalfdim)
    // AND the points outside max_r should already be zero...
    
    // Use the inverse matrix
    if (inv){
        for (int i = 0; i < 3; i++)
            for (int j = 0; j < 3; j++)
                Ainv[i][j] = A[i][j]*(FDOUBLE)padding_factor;
    }
    else{
        for (int i = 0; i < 3; i++)
            for (int j = 0; j < 3; j++)
                Ainv[i][j] = A[j][i]*(FDOUBLE)padding_factor;
    }
    
    // Go from the 2D slice coordinates to the 3D coordinates
    int max_r2 = r_max * r_max;
    int min_r2_nn = r_min_nn * r_min_nn;
    
    //#define DEBUG_BACKP
#ifdef DEBUG_BACKP
    //
#endif
    
    for (int shell_n = shell_n_start; shell_n < shell_n_end; shell_n++)
    {
        int n = nIndex[shell_n];
        int i = n / (f2d_size/2+1);
        // Dont search beyond square with side max_r
        if (i <= r_max)
        {
            y = i;
            first_x = 0;
        }
        else if (i >= f2d_size - r_max)
        {
            y = i - f2d_size;
            // x==0 plane is stored twice in the FFTW format. Dont set it twice in BACKPROJECTION!
            first_x = 1;
        }
        else
            continue;
        //
        int x = n % (f2d_size/2+1);
        if (x < first_x)
            continue;
        y2 = y * y;
        // Only include points with radius < max_r (exclude points outside circle in square)
        r2 = x * x + y2;
        if (r2 > max_r2) // always not meet??
            continue;
        
        // Get the relevant value in the input image
        my_val_real = f2d_real[shell_n];
        my_val_imag = f2d_imag[shell_n];
        // Get the weight
        if (Mweight != NULL)
            my_weight = Mweight[shell_n];
        // else: my_weight was already initialised to 1.
        // std::cout<<shell_n<<" "<<y<<" "<<x<<" "<<my_val_real<<" "<<my_val_imag<<" "<<my_weight<<std::endl;
        if (my_weight > 0.)
        {
            
            // Get logical coordinates in the 3D map
            xp = Ainv[0][0] * x + Ainv[0][1] * y;
            yp = Ainv[1][0] * x + Ainv[1][1] * y;
            zp = Ainv[2][0] * x + Ainv[2][1] * y;
            
            if (interpolator == TRILINEAR || r2 < min_r2_nn)
            {
                
                // Only asymmetric half is stored
                if (xp < 0)
                {
                    // Get complex conjugated hermitian symmetry pair
                    xp = -xp;
                    yp = -yp;
                    zp = -zp;
                    is_neg_x = true;
                }
                else
                {
                    is_neg_x = false;
                }
                
                // Trilinear interpolation (with physical coords)
                // Subtract STARTINGY and STARTINGZ to accelerate access to data (STARTINGX=0)
                // In that way use DIRECT_A3D_ELEM, rather than A3D_ELEM
                x0 = floor(xp);
                fx = xp - x0;
                x1 = x0 + 1;
                
                y0 = floor(yp);
                fy = yp - y0;
                y0 -=  pad_size_shift;
                y1 = y0 + 1;
                
                z0 = floor(zp);
                fz = zp - z0;
                z0 -= pad_size_shift;
                z1 = z0 + 1;
                
                mfx = 1. - fx;
                mfy = 1. - fy;
                mfz = 1. - fz;
                
                dd000 = mfz * mfy * mfx;
                dd001 = mfz * mfy *  fx;
                dd010 = mfz *  fy * mfx;
                dd011 = mfz *  fy *  fx;
                dd100 =  fz * mfy * mfx;
                dd101 =  fz * mfy *  fx;
                dd110 =  fz *  fy * mfx;
                dd111 =  fz *  fy *  fx;
                
                if (is_neg_x){
                    my_val_real = my_val_real;
                    my_val_imag = -my_val_imag;
                }
                // Store slice in 3D weighted sum
                assert(x0>=0);assert(y0>=0);assert(z0>=0);
                // TODO. merge (real,imag,weight) to better memory fetch
                auto data_td000 = &data_td(z0, y0, x0).real;
                data_td000[0] += dd000 * my_val_real;data_td000[1] += dd000 * my_val_imag;
                data_td000[2] += dd001 * my_val_real;data_td000[3] += dd001 * my_val_imag;
                auto data_td010 = &data_td(z0, y1, x0).real;
                data_td010[0] += dd010 * my_val_real;data_td010[1] += dd010 * my_val_imag;
                data_td010[2] += dd011 * my_val_real;data_td010[3] += dd011 * my_val_imag;
                auto data_td100 = &data_td(z1, y0, x0).real;
                data_td100[0] += dd100 * my_val_real;data_td100[1] += dd100 * my_val_imag;
                data_td100[2] += dd101 * my_val_real;data_td100[3] += dd101 * my_val_imag;
                auto data_td110 = &data_td(z1, y1, x0).real;
                data_td110[0] += dd110 * my_val_real;data_td110[1] += dd110 * my_val_imag;
                data_td110[2] += dd111 * my_val_real;data_td110[3] += dd111 * my_val_imag;
                
                auto weight_td000 = &weight_td(z0, y0, x0);
                weight_td000[0] += dd000 * my_weight;weight_td000[1] += dd001 * my_weight;
                auto weight_td010 = &weight_td(z0, y1, x0);
                weight_td010[0] += dd010 * my_weight;weight_td010[1] += dd011 * my_weight;
                auto weight_td100 = &weight_td(z1, y0, x0);
                weight_td100[0] += dd100 * my_weight;weight_td100[1] += dd101 * my_weight;
                auto weight_td110 = &weight_td(z1, y1, x0);
                weight_td110[0] += dd110 * my_weight;weight_td110[1] += dd111 * my_weight;
            } // endif TRILINEAR
            else if (interpolator == NEAREST_NEIGHBOUR )
            {
                
                x0 = round(xp);
                y0 = round(yp);
                z0 = round(zp);
                
                if (x0 < 0)
                {
                    ACCESS(data_td, -z0, -y0, -x0).real += my_val_real;
                    ACCESS(data_td, -z0, -y0, -x0).imag += (-my_val_imag);
                    ACCESS(weight_td, -z0, -y0, -x0) += my_weight;
                }
                else
                {
                    ACCESS(data_td, z0, y0, x0).real += my_val_real;
                    ACCESS(data_td, z0, y0, x0).imag += my_val_imag;
                    ACCESS(weight_td, z0, y0, x0) += my_weight;
                }
                
            } // endif NEAREST_NEIGHBOUR
            else
            {
                ERROR_REPORT("FourierInterpolator::backproject%%ERROR: unrecognized interpolator ");
            }
        } // endif weight>0.
    } // endif shell-n-loop
}

void MyBackProjector::backrotate2D(const FDOUBLE* f2d_real,const FDOUBLE* f2d_imag,int f2d_size,const FDOUBLE A[][3], bool inv,
                  					Vol<MKL_Complex>& data_td,Vol<FDOUBLE>& weight_td,const FDOUBLE* Mweight /*= nullptr*/)
{
    FDOUBLE fx, fy, mfx, mfy, xp, yp;
    int first_x, x0, x1, y0, y1, y, y2, r2;
    bool is_neg_x;
    FDOUBLE dd00, dd01, dd10, dd11;
    FDOUBLE my_val_real,my_val_imag;
    FDOUBLE Ainv[3][3];
    FDOUBLE my_weight = 1.;
    
    // Use the inverse matrix
    if (inv){
        for (int i = 0; i < 3; i++)
            for (int j = 0; j < 3; j++)
                Ainv[i][j] = A[i][j]*(FDOUBLE)padding_factor;
    }
    else{
        for (int i = 0; i < 3; i++)
            for (int j = 0; j < 3; j++)
                Ainv[i][j] = A[j][i]*(FDOUBLE)padding_factor;
    }
    
    // Go from the 2D slice coordinates to the data-array coordinates
    int f2d_Fsize = f2d_size/2 + 1;
    int pad_size_shift = -(int)((float)(pad_size)/2.0);
    int max_r2 = r_max * r_max;
    int min_r2_nn = r_min_nn * r_min_nn;
    
    for (int i=0; i < f2d_size; i++)
    {
        // Don't search beyond square with side max_r
        if (i <= r_max)
        {
            y = i;
            first_x = 0;
        }
        else if (i >= f2d_size - r_max)
        {
            y = i - f2d_size;
            // x==0 plane is stored twice in the FFTW format. Dont set it twice in BACKPROJECTION!
            first_x = 1;
        }
        else
            continue;
        
        y2 = y * y;
        for (int x=first_x; x <= r_max; x++)
        {
            // Only include points with radius < max_r (exclude points outside circle in square)
            r2 = x * x + y2;
            if (r2 > max_r2)
                continue;
            
            // Get the relevant value in the input image
            my_val_real = f2d_real[i*f2d_Fsize+x];
            my_val_imag = f2d_imag[i*f2d_Fsize+x];
            
            // Get the weight
            if (Mweight != nullptr)
                my_weight = Mweight[i*f2d_Fsize+x];
            // else: my_weight was already initialised to 1.
            
            if (my_weight > 0.)
            {
                // Get logical coordinates in the 3D map
                xp = Ainv[0][0] * x + Ainv[0][1] * y;
                yp = Ainv[1][0] * x + Ainv[1][1] * y;
                
                if (interpolator == TRILINEAR || r2 < min_r2_nn)
                {
                    // Only asymmetric half is stored
                    if (xp < 0)
                    {
                        // Get complex conjugated hermitian symmetry pair
                        xp = -xp;
                        yp = -yp;
                        is_neg_x = true;
                    }
                    else
                    {
                        is_neg_x = false;
                    }
                    
                    // Trilinear interpolation (with physical coords)
                    // Subtract STARTINGY to accelerate access to data (STARTINGX=0)
                    // In that way use DIRECT_A2D_ELEM, rather than A2D_ELEM
                    x0 = floor(xp);
                    fx = xp - x0;
                    x1 = x0 + 1;
                    
                    y0 = floor(yp);
                    fy = yp - y0;
                    y0 -= pad_size_shift;
                    y1 = y0 + 1;
                    
                    mfx = 1. - fx;
                    mfy = 1. - fy;
                    
                    dd00 = mfy * mfx;
                    dd01 = mfy *  fx;
                    dd10 =  fy * mfx;
                    dd11 =  fy *  fx;
                    
                    if (is_neg_x){
                        my_val_imag = -1*my_val_imag;
                    }
                    // Store slice in 3D weighted sum
                    auto data_td00 = &data_td(0, y0, x0).real;
                    data_td00[0] += dd00 * my_val_real; data_td00[1] += dd00 * my_val_imag;
                    data_td00[2] += dd01 * my_val_real; data_td00[3] += dd01 * my_val_imag;
                    auto data_td10 = &data_td(0, y1, x0).real;
                    data_td10[0] += dd10 * my_val_real; data_td10[1] += dd10 * my_val_imag;
                    data_td10[2] += dd11 * my_val_real; data_td10[3] += dd11 * my_val_imag;
                    
                    // Store corresponding weights
                    auto weight_td00 = &weight_td(0, y0, x0);
                    weight_td00[0] += dd00 * my_weight; weight_td00[1] += dd01 * my_weight;
                    auto weight_td10 = &weight_td(0, y1, x0);
                    weight_td10[0] += dd10 * my_weight; weight_td10[1] += dd11 * my_weight;
                } // endif TRILINEAR
                
                else if (interpolator == NEAREST_NEIGHBOUR )
                {
                    std::cerr<<"interpolator == NEAREST_NEIGHBOUR has not test yet."<<std::endl;
                    x0 = round(xp);
                    y0 = round(yp);
                    if (x0 < 0)
                    {
                        data_td(0,-y0,-x0).real += my_val_real;
                        data_td(0,-y0,-x0).imag += -1*my_val_imag;
                        
                        weight_td(0,-y0,-x0) += my_weight;
                    }
                    else
                    {
                        data_td(0,y0,x0).real += my_val_real;
                        data_td(0,y0,x0).imag += my_val_imag;
                        
                        weight_td(0,y0,x0) += my_weight;
                    }
                } // endif NEAREST_NEIGHBOUR
                
            } // endif weight > 0.
        } // endif x-loop
    } // endif y-loop
}

void MyBackProjector::enforceHermitianSymmetry(Vol<MKL_Complex>& my_data,Vol<FDOUBLE>& my_weight)
{
    assert(my_data.dimz%2==1);
    // TOOD check size
    // my_data_real equal my_weight
    int my_data_origin_z = XMIPP_ORIGIN(my_data.dimz);
    int my_data_origin_y = XMIPP_ORIGIN(my_data.dimy);
    
    // set xmipp origin,iz from (-dimz/2) ~ (dimz/2)
    int startz = my_data_origin_z;
    int endz = my_data.dimz + my_data_origin_z;
    for (int iz = startz; iz < endz; iz++)
    {
        int z_left = -iz - my_data_origin_z;
        int z_right = iz - my_data_origin_z;
        // Make sure all points are only included once.
        // iy from 0(or 1) ~ (dimy/2)
        int starty = (iz < 0) ? 0 : 1;
        int endy = my_data.dimy + my_data_origin_y;
        for (int iy = starty; iy < endy; iy++)
        {
            // do symmetry in x=0 plane
            int y_left = -iy - my_data_origin_y;
            int y_right = iy - my_data_origin_y;
            // I just need to sum the two points, not divide by 2!
            FDOUBLE fsum_real = my_data(z_right,y_right,0).real + my_data(z_left,y_left,0).real;
            FDOUBLE fsum_imag = my_data(z_right,y_right,0).imag - my_data(z_left,y_left,0).imag;
            my_data(z_right,y_right,0).real = fsum_real;
            my_data(z_right,y_right,0).imag = fsum_imag;
            my_data(z_left,y_left,0).real = fsum_real;
            my_data(z_left,y_left,0).imag = -fsum_imag;
            FDOUBLE sum = my_weight(z_right,y_right,0) + my_weight(z_left,y_left,0);
            my_weight(z_right,y_right,0) = sum;
            my_weight(z_left,y_left,0) = sum;
        }
    }
}

void MyBackProjector::symmetrise(Vol<MKL_Complex>& my_data,
                               	 Vol<FDOUBLE>& my_weight, int my_rmax2)
{
    //#define DEBUG_SYMM
#ifdef DEBUG_SYMM
    std::cerr << " SL.SymsNo()= " << SL.SymsNo() << std::endl;
    std::cerr << " SL.true_symNo= " << SL.true_symNo << std::endl;
#endif
    SymList SL(fn_sym);
    if (SL.SymsNo() > 0 && ref_dim == 3)
    {
        Matrix2D<FDOUBLE> L(4, 4), R(4, 4); // A matrix from the list
        FDOUBLE x, y, z, fx, fy, fz, xp, yp, zp, r2;
        bool is_neg_x;
        int x0, x1, y0, y1, z0, z1;
        
#define COMPLEX(N) FDOUBLE N##_real,N##_imag;
        COMPLEX(d000) COMPLEX(d001) COMPLEX(d010) COMPLEX(d011)
        COMPLEX(d100) COMPLEX(d101) COMPLEX(d110) COMPLEX(d111)
        COMPLEX(dx00) COMPLEX(dx01) COMPLEX(dx10) COMPLEX(dx11) COMPLEX(dxy0) COMPLEX(dxy1)
#undef COMPLEX
        
        FDOUBLE dd000, dd001, dd010, dd011, dd100, dd101, dd110, dd111;
        FDOUBLE ddx00, ddx01, ddx10, ddx11, ddxy0, ddxy1;
        
        int xmipp_origin_z = XMIPP_ORIGIN(my_data.dimz);
        int xmipp_origin_y = XMIPP_ORIGIN(my_data.dimy);
        
        Vol<FDOUBLE> sum_weight;sum_weight.init_nonzero(my_weight.dimz,my_weight.dimy,my_weight.dimx);
        Vol<MKL_Complex> sum_data;sum_data.init_nonzero(my_data.dimz,my_data.dimy,my_data.dimx);
        // First symmetry operator (not stored in SL) is the identity matrix
        for (int i = 0; i < my_data.dimzyx; i++) {
            sum_data.wptr()[i] = my_data.wptr()[i];
            sum_weight.wptr()[i] = my_weight.wptr()[i];
        }
        // Loop over all other symmetry operators
        for (int isym = 0; isym < SL.SymsNo(); isym++)
        {
            SL.get_matrices(isym, L, R);
#ifdef DEBUG_SYMM
            std::cerr << " isym= " << isym << " R= " << R << std::endl;
#endif
            
            // Loop over all points in the output (i.e. rotated, or summed) array
            // FOR_ALL_ELEMENTS_IN_ARRAY3D(sum_weight)
            for (int k=0; k<my_data.dimz; k++)
            for (int i=0; i<my_data.dimy; i++)
            for (int j=0; j<my_data.dimx; j++)
            {
                x = (FDOUBLE)j; // STARTINGX(sum_weight) is zero!
                y = (FDOUBLE)i+xmipp_origin_y;
                z = (FDOUBLE)k+xmipp_origin_z;
                r2 = x*x + y*y + z*z;
                if (r2 <= my_rmax2)
                {
                    // coords_output(x,y) = A * coords_input (xp,yp)
                    xp = x * R(0, 0) + y * R(0, 1) + z * R(0, 2);
                    yp = x * R(1, 0) + y * R(1, 1) + z * R(1, 2);
                    zp = x * R(2, 0) + y * R(2, 1) + z * R(2, 2);
                    
                    // Only asymmetric half is stored
                    if (xp < 0)
                    {
                        // Get complex conjugated hermitian symmetry pair
                        xp = -xp;
                        yp = -yp;
                        zp = -zp;
                        is_neg_x = true;
                    }
                    else
                    {
                        is_neg_x = false;
                    }
                    
                    // Trilinear interpolation (with physical coords)
                    // Subtract STARTINGY and STARTINGZ to accelerate access to data (STARTINGX=0)
                    // In that way use DIRECT_A3D_ELEM, rather than A3D_ELEM
                    x0 = floor(xp);
                    fx = xp - x0;
                    x1 = x0 + 1;
                    
                    y0 = floor(yp);
                    fy = yp - y0;
                    y0 -=  xmipp_origin_y;//STARTINGY(my_data)
                    y1 = y0 + 1;
                    
                    z0 = floor(zp);
                    fz = zp - z0;
                    z0 -= xmipp_origin_z;//STARTINGZ(my_data);
                    z1 = z0 + 1;
                    
#ifdef CHECK_SIZE
                    if (x0 < 0 || y0 < 0 || z0 < 0 ||
                        x1 < 0 || y1 < 0 || z1 < 0 ||
                        x0 >= XSIZE(my_data) || y0  >= YSIZE(my_data) || z0 >= ZSIZE(my_data) ||
                        x1 >= XSIZE(my_data) || y1  >= YSIZE(my_data)  || z1 >= ZSIZE(my_data) 	)
                    {
                        std::cerr << " x0= " << x0 << " y0= " << y0 << " z0= " << z0 << std::endl;
                        std::cerr << " x1= " << x1 << " y1= " << y1 << " z1= " << z1 << std::endl;
                        my_data.printShape();
                        REPORT_ERROR("BackProjector::symmetrise: checksize!!!");
                    }
#endif
                    // First interpolate (complex) data
                    auto my_data000 = &my_data(z0, y0, x0).real;
                    d000_real = my_data000[0];d000_imag = my_data000[1];
                    d001_real = my_data000[2];d001_imag = my_data000[3];
                    auto my_data010 = &my_data(z0, y1, x0).real;
                    d010_real = my_data010[0];d010_imag = my_data010[1];
                    d011_real = my_data010[2];d011_imag = my_data010[3];
                    auto my_data100 = &my_data(z1, y0, x0).real;
                    d100_real = my_data100[0];d100_imag = my_data100[1];
                    d101_real = my_data100[2];d101_imag = my_data100[3];
					auto my_data110 = &my_data(z1, y1, x0).real;
                    d110_real = my_data110[0];d110_imag = my_data110[1];
                    d111_real = my_data110[2];d111_imag = my_data110[3];
                    
                    dx00_real = LIN_INTERP(fx, d000_real, d001_real);
                    dx00_imag = LIN_INTERP(fx, d000_imag, d001_imag);
                    dx01_real = LIN_INTERP(fx, d100_real, d101_real);
                    dx01_imag = LIN_INTERP(fx, d100_imag, d101_imag);
                    dx10_real = LIN_INTERP(fx, d010_real, d011_real);
                    dx10_imag = LIN_INTERP(fx, d010_imag, d011_imag);
                    dx11_real = LIN_INTERP(fx, d110_real, d111_real);
                    dx11_imag = LIN_INTERP(fx, d110_imag, d111_imag);
                    dxy0_real = LIN_INTERP(fy, dx00_real, dx10_real);
                    dxy0_imag = LIN_INTERP(fy, dx00_imag, dx10_imag);
                    dxy1_real = LIN_INTERP(fy, dx01_real, dx11_real);
                    dxy1_imag = LIN_INTERP(fy, dx01_imag, dx11_imag);
                    
                    // Take complex conjugated for half with negative x
                    if (is_neg_x){
                        sum_data(k,i,j).real += LIN_INTERP(fz, dxy0_real, dxy1_real);
                        sum_data(k,i,j).imag += -1*( LIN_INTERP(fz, dxy0_imag, dxy1_imag) );
                    }
                    else{
                        sum_data(k,i,j).real += LIN_INTERP(fz, dxy0_real, dxy1_real);
                        sum_data(k,i,j).imag += LIN_INTERP(fz, dxy0_imag, dxy1_imag);
                    }
                    // Then interpolate (real) weight
                    auto my_weight000 = &my_weight(z0, y0, x0);
                    dd000 = my_weight000[0];dd001 = my_weight000[1];
                    auto my_weight010 = &my_weight(z0, y1, x0);
                    dd010 = my_weight010[0];dd011 = my_weight010[1];
                    auto my_weight100 = &my_weight(z1, y0, x0);
                    dd100 = my_weight100[0];dd101 = my_weight100[1];
                    auto my_weight110 = &my_weight(z1, y1, x0);
                    dd110 = my_weight110[0];dd111 = my_weight110[1];
                    
                    ddx00 = LIN_INTERP(fx, dd000, dd001);
                    ddx01 = LIN_INTERP(fx, dd100, dd101);
                    ddx10 = LIN_INTERP(fx, dd010, dd011);
                    ddx11 = LIN_INTERP(fx, dd110, dd111);
                    ddxy0 = LIN_INTERP(fy, ddx00, ddx10);
                    ddxy1 = LIN_INTERP(fy, ddx01, ddx11);
                    
                    sum_weight(k, i, j) +=  LIN_INTERP(fz, ddxy0, ddxy1);
                    
                } // end if r2 <= my_rmax2
                
            } // end loop over all elements of sum_weight
            
        } // end loop over symmetry operators
        
        for (int i = 0; i < my_data.dimzyx; i++) {
            my_data.wptr()[i] = sum_data.wptr()[i];
            my_weight.wptr()[i] = sum_weight.wptr()[i];
        }
        sum_data.fini();sum_weight.fini();
        // Average
        // The division should only be done if we would search all (C1) directions, not if we restrict the angular search!
        /*
         FOR_ALL_DIRECT_ELEMENTS_IN_MULTIDIMARRAY(data)
         {
         DIRECT_MULTIDIM_ELEM(data, n) = DIRECT_MULTIDIM_ELEM(sum_data, n) / (FDOUBLE)(SL.SymsNo() + 1);
         DIRECT_MULTIDIM_ELEM(weight, n) = DIRECT_MULTIDIM_ELEM(sum_weight, n) / (FDOUBLE)(SL.SymsNo() + 1);
         }
         */
    }
    
}

void MyBackProjector::convoluteBlobRealSpace(Vol<MKL_Complex>& Fconv, FourierTransformerBase& transformer, bool do_mask /*= false*/)
{
    TIMEPOINT_INIT
    TIMEPOINT
    // TODO,check size
    assert(pad_size%2 == 1);
    assert(Fconv.dimy == pad_size);

    Vol<FDOUBLE> Mconv;
    int padhdim = pad_size / 2;
    
    // Set up right dimension of real-space array
    // TODO: resize this according to r_max!!!
    Mconv.init_nonzero(Fconv.dimz,Fconv.dimy, Fconv.dimy);
    transformer.reset_plan(Mconv.wptr(),(FourierComplex*)Fconv.wptr(),Mconv.dimx,Mconv.dimy,Mconv.dimz);
    // inverse FFT
    transformer.inverseFourierTransform();
    TIMEPOINT
    // Blob normalisation in Fourier space
    FDOUBLE normftblob = tab_ftblob(0.);
    
    // TMP DEBUGGING
    //struct blobtype blob;
    //blob.order = 0;
    //blob.radius = 1.9 * padding_factor;
    //blob.alpha = 15;
    
    // Multiply with FT of the blob kernel
	if (1) {
		static int previouslyReported = 0;
		if (previouslyReported != transformer.nr_threads) {
			#pragma omp critical 
			if (previouslyReported != transformer.nr_threads) {
				previouslyReported = transformer.nr_threads;
				MASTERNODE std::cerr << __FILE__ << ":" << __LINE__ << " convoluteBlobRealSpace() transformer.nr_threads:" << transformer.nr_threads << std::endl;
			}
		}
	}

#pragma omp parallel for collapse(2) if(transformer.nr_threads > 1) num_threads(transformer.nr_threads)
    for (int k = 0;k < Mconv.dimz;k++){
        for (int i = 0;i < Mconv.dimy;i++){
            int kp = (k < padhdim) ? k : k - pad_size;
            int ip = (i < padhdim) ? i : i - pad_size;
            for (int j = 0;j < Mconv.dimx;j++)
            {
                int jp = (j < padhdim) ? j : j - pad_size;
                FDOUBLE rval = sqrt ( (FDOUBLE)(kp * kp + ip * ip + jp * jp) ) / (ori_size * padding_factor);
                //if (kp==0 && ip==0 && jp > 0)
                //	std::cerr << " jp= " << jp << " rval= " << rval << " tab_ftblob(rval) / normftblob= " << tab_ftblob(rval) / normftblob << " ori_size/2= " << ori_size/2 << std::endl;
                // In the final reconstruction: mask the real-space map beyond its original size to prevent aliasing ghosts
                // Note that rval goes until 1/2 in the oversampled map
                if (do_mask && rval > 1./(2. * padding_factor))
                    Mconv(k, i, j) = 0.;
                else
                    Mconv(k, i, j) *= (tab_ftblob(rval) / normftblob);
            }
        }
    }
    TIMEPOINT
    // forward FFT to go back to Fourier-space
    transformer.FourierTransform();
    
    Mconv.fini();
    TIMEPOINT
    TIMEPOINT_FINA
    
}

void MyBackProjector::windowToOridimRealSpace(Vol<MKL_Complex> &Fconv, Vol<FDOUBLE> &Vol_out, FourierTransformerBase& transformer)
{
    TIMEPOINT_INIT
    TIMEPOINT
    int padoridim = padding_factor * ori_size;
    Vol<MKL_Complex> Ftmp;
    Vol<FDOUBLE> Mout,MoutCenter;
    
    if (ref_dim == 2) {
        Ftmp.init_nonzero(1, padoridim, padoridim/2+1);
    }
    else{
        Ftmp.init_nonzero(padoridim, padoridim, padoridim/2+1);
    }
    
    TIMEPOINT
    FDOUBLE normfft;
    
    //#define DEBUG_WINDOWORIDIMREALSPACE
#ifdef DEBUG_WINDOWORIDIMREALSPACE
    // Image<FDOUBLE> tt;
    // tt().resize(ZSIZE(Fin), YSIZE(Fin), XSIZE(Fin));
    // FOR_ALL_DIRECT_ELEMENTS_IN_MULTIDIMARRAY(Fin)
    // {
    //     DIRECT_MULTIDIM_ELEM(tt(), n) = abs(DIRECT_MULTIDIM_ELEM(Fin, n));
    // }
    // tt.write("windoworidim_Fin.spi");
#endif
    
    // the Vol_out should be in right size
    if (ref_dim == 2)
    {
        normfft = (FDOUBLE)(padding_factor * padding_factor);
    }
    else
    {
        normfft = (FDOUBLE)(padding_factor * padding_factor * padding_factor * ori_size);
    }
    
    // Resize incoming complex array to the correct size
    windowFourierTransform(Fconv, Ftmp);
    
    Fconv.fini();
    if (ref_dim == 2) {
        Mout.init_nonzero(1, padoridim, padoridim);
    }
    else{
        Mout.init_nonzero(padoridim, padoridim, padoridim);
    }
    
    transformer.reset_plan(Mout.wptr(), (FourierComplex*)Ftmp.wptr(),Mout.dimx,Mout.dimy,Mout.dimz);
    
    TIMEPOINT
#ifdef DEBUG_WINDOWORIDIMREALSPACE
    // tt().resize(ZSIZE(Ftmp), YSIZE(Ftmp), XSIZE(Ftmp));
    // FOR_ALL_DIRECT_ELEMENTS_IN_MULTIDIMARRAY(Ftmp)
    // {
    //     DIRECT_MULTIDIM_ELEM(tt(), n) = abs(DIRECT_MULTIDIM_ELEM(Ftmp, n));
    // }
    // tt.write("windoworidim_Fresized.spi");
#endif
    
    // Do the inverse FFT
    transformer.inverseFourierTransform();
    
    // clear the memory
    Ftmp.fini();
    if (ref_dim == 2)
    {
        MoutCenter.init_nonzero(1, padoridim, padoridim);
    }
    else
    {
        MoutCenter.init_nonzero(padoridim, padoridim, padoridim);
    }
    
    TIMEPOINT
    // Shift the map back to its origin
    if (ref_dim == 2)
    {
        centerFFT2D(Mout.wptr(),MoutCenter.wptr(),padoridim,true);
    }
    else
    {
        centerFFT3DForward(Mout.wptr(),MoutCenter.wptr(),padoridim);
    }
    TIMEPOINT
#ifdef DEBUG_WINDOWORIDIMREALSPACE
    // tt()=Mout;
    // tt.write("windoworidim_Munwindowed.spi");
#endif
    
    // clear the memory
    Mout.fini();
    
	TIMEPOINT
    // Window in real-space
    if (ref_dim==2)
    {
        window(MoutCenter,Vol_out);
    }
    else
    {
        window(MoutCenter,Vol_out);
    }
    
    // clear the memory
    MoutCenter.fini();
    
    TIMEPOINT
    // Normalisation factor of FFTW
    // The Fourier Transforms are all "normalised" for 2D transforms of size = ori_size x ori_size
    for (int n = 0; n < Vol_out.dimzyx; n++)
        Vol_out.wptr()[n] /= normfft;
    
#ifdef DEBUG_WINDOWORIDIMREALSPACE
    // tt()=Mout;
    // tt.write("windoworidim_Mwindowed.spi");
#endif
    TIMEPOINT
    // Mask out corners to prevent aliasing artefacts
    softMaskOutsideMap(Vol_out);
	TIMEPOINT
    
#ifdef DEBUG_WINDOWORIDIMREALSPACE
    // tt()=Mout;
    // tt.write("windoworidim_Mwindowed_masked.spi");
    // FourierTransformer ttf;
    // ttf.FourierTransform(Mout, Ftmp);
    // tt().resize(ZSIZE(Ftmp), YSIZE(Ftmp), XSIZE(Ftmp));
    // FOR_ALL_DIRECT_ELEMENTS_IN_MULTIDIMARRAY(Ftmp)
    // {
    //     DIRECT_MULTIDIM_ELEM(tt(), n) = abs(DIRECT_MULTIDIM_ELEM(Ftmp, n));
    // }
    // tt.write("windoworidim_Fnew.spi");
#endif
    
}

void MyBackProjector::reconstruct(	Vol<FDOUBLE> &vol_out,
                                  	int max_iter_preweight,
                                  	bool do_map,
                                  	double tau2_fudge,
                                  	FDOUBLE* tau2,
                                  	FDOUBLE* sigma2,
                                  	FDOUBLE* data_vs_prior,
                                  	const FDOUBLE* fsc, // only input
                                  	double normalise/* = 1.*/,
                                  	bool update_tau2_with_fsc/* = false*/,
                                  	bool is_whole_instead_of_half/* = false*/,
                                  	int minres_map/* = -1*/,
									int nr_threads/* = 1*/,
                                    std::string tmp_folder/* = "NULL"*/)
{
    // the Vol_out,taus,sigma2,data_vs_prior and fsc should be already in right size
    FourierTransformerBase transformer;
    Vol<MKL_Complex> Fconv;
    Vol<FDOUBLE> Fweight;
    // NOTE : Fnewweight can become too large for a float: always keep this one in double-precision
    Vol<double> Fnewweight;
    
    TIMEPOINT_INIT
    TIMEPOINT
    if(ref_dim == 3){
        Fconv.init_nonzero(pad_size,pad_size,pad_Fsize);
        Fweight.init_nonzero(pad_size,pad_size,pad_Fsize);
        Fnewweight.init_nonzero(pad_size, pad_size, pad_Fsize);
    }
    else{
        Fconv.init_nonzero(1,pad_size,pad_Fsize);
        Fweight.init_nonzero(1,pad_size,pad_Fsize);
        Fnewweight.init_nonzero(1, pad_size, pad_Fsize);
    }
    TIMEPOINT
    transformer.init(vol_out.wptr(), (FourierComplex*)Fconv.wptr(), vol_out.dimy, vol_out.dimy, vol_out.dimz, nr_threads);
    
    int max_r2 = r_max * r_max * padding_factor * padding_factor;
    TIMEPOINT
    //#define DEBUG_RECONSTRUCT
#ifdef DEBUG_RECONSTRUCT
    // Image<FDOUBLE> ttt;
    // FileName fnttt;
    // ttt()=weight;
    // ttt.write("reconstruct_initial_weight.spi");
#endif
    
    printSum_real(data.wptr(),data.dimzyx,"data_real before enforceHermitianSymmetry");
    printSum_imag(data.wptr(),data.dimzyx,"data_imag before enforceHermitianSymmetry");
    
    // At the x=0 line, we have collected either the positive y-z coordinate, or its negative Friedel pair.
    // Sum these two together for both the data and the weight arrays
    enforceHermitianSymmetry(data, weight);
    printSum_real(data.wptr(),data.dimzyx,"data_real after enforceHermitianSymmetry");
    printSum_imag(data.wptr(),data.dimzyx,"data_imag after enforceHermitianSymmetry");
    printSum(weight.wptr(),weight.dimzyx,"weight after enforceHermitianSymmetry");

#ifdef DEBUG_RECONSTRUCT
    // ttt()=weight;
    // ttt.write("reconstruct_hermitian_weight.spi");
#endif
    TIMEPOINT
    // First enforce Hermitian symmetry, then symmetry!
    // This way the redundancy at the x=0 plane is handled correctly
    symmetrise(data, weight, max_r2);
#ifdef DEBUG_RECONSTRUCT
    // ttt()=weight;
    // ttt.write("reconstruct_symmetrised_weight.spi");
    // FOR_ALL_DIRECT_ELEMENTS_IN_MULTIDIMARRAY(data)
    // {
    //     DIRECT_MULTIDIM_ELEM(ttt(), n) = DIRECT_MULTIDIM_ELEM(data, n).real;
    // }
    // ttt.write("reconstruct_symmetrised_data_real.spi");
    //
    // std::cerr << " pad_size= " << pad_size << " padding_factor= " << padding_factor << " max_r2= " << max_r2 << std::endl;
#endif
    TIMEPOINT
    // Go from projector-centered to FFTW-uncentered
    decenter(weight, Fweight, max_r2);
    printSum(weight.wptr(),weight.dimzyx,"weight after decenter");

    // Take oversampling into account
    FDOUBLE oversampling_correction = (ref_dim == 3) ? (padding_factor * padding_factor * padding_factor) : (padding_factor * padding_factor);

    std::vector<FDOUBLE> counter(ori_Fsize);
    
    // First calculate the radial average of the (inverse of the) power of the noise in the reconstruction
    // This is the left-hand side term in the nominator of the Wiener-filter-like update formula
    // and it is stored inside the weight vector
    // Then, if (do_map) add the inverse of tau2-spectrum values to the weight
    for (int i = 0; i < ori_Fsize; i++)
        sigma2[i] = counter[i] = 0.;
	TIMEPOINT
    
    for (int k = 0; k < Fweight.dimz; k++){
        for (int i = 0; i < Fweight.dimy; i++){
            int kp = (k < Fweight.dimx) ? k : k - Fweight.dimz;
            int ip = (i < Fweight.dimx) ? i : i - Fweight.dimy;
            for (int j = 0, jp = 0; j < Fweight.dimx; j++, jp = j)
            {
                int r2 = kp * kp + ip * ip + jp * jp;
                if (r2 < max_r2)
                {
                    int ires = round( sqrt((FDOUBLE)r2) / padding_factor );
                    FDOUBLE invw = oversampling_correction * Fweight(k, i, j);
                    sigma2[ires] += invw;
                    counter[ires] += 1.;
                }
            }
    }}
    
    printSum(&sigma2[0],ori_Fsize,"sigma2 after init");
    printSum(&counter[0],ori_Fsize,"counter after init");
	TIMEPOINT
    // Average (inverse of) sigma2 in reconstruction
    for (int i = 0; i < ori_Fsize; i++)
    {
        if (sigma2[i] > 1e-10)
            sigma2[i] = counter[i] / sigma2[i];
        else if (sigma2[i] == 0)
            sigma2[i] = 0.;
        else
        {
            std::cerr << " DIRECT_A1D_ELEM(sigma2, i)= " << sigma2[i] << std::endl;
            ERROR_REPORT("BackProjector::reconstruct: ERROR: unexpectedly small, yet non-zero sigma2 value, this should not happen...a");
        }
    }
    printSum(&sigma2[0],ori_Fsize,"sigma2 after init 2");
    if (update_tau2_with_fsc)
    {
        for (int i = 0; i < ori_Fsize; i++)
            data_vs_prior[i] = 0.;
        // Then calculate new tau2 values, based on the FSC
        for (int i = 0; i < ori_Fsize; i++)
        {
            // FSC cannot be negative or zero for conversion into tau2
            FDOUBLE myfsc = std::max((FDOUBLE)0.001, fsc[i]);
            if (is_whole_instead_of_half)
            {
                // Factor two because of twice as many particles
                // Sqrt-term to get 60-degree phase errors....
                myfsc = sqrt(2. * myfsc / (myfsc + 1.));
            }
            myfsc = std::min((FDOUBLE)0.999, myfsc);
            FDOUBLE myssnr = myfsc / (1. - myfsc);
            FDOUBLE fsc_based_tau = myssnr * sigma2[i];
            tau2[i] = fsc_based_tau;
            // data_vs_prior is merely for reporting: it is not used for anything in the reconstruction
            data_vs_prior[i] = myssnr;
        }
    }
    TIMEPOINT
    // Apply MAP-additional term to the Fnewweight array
    // This will regularise the actual reconstruction
    if (do_map)
    {
        // Then, add the inverse of tau2-spectrum values to the weight
        // and also calculate spherical average of data_vs_prior ratios
        if (!update_tau2_with_fsc)
            for (int i = 0; i < ori_Fsize; i++)
                data_vs_prior[i] = 0.;
        for (int i = 0; i < ori_Fsize; i++)
            counter[i] = 0.;
        //FOR_ALL_ELEMENTS_IN_FFTW_TRANSFORM(Fconv)
        for (int k = 0; k < Fweight.dimz; k++){
            for (int i = 0; i < Fweight.dimy; i++){
                int kp = (k < Fweight.dimx) ? k : k - Fweight.dimz;
                int ip = (i < Fweight.dimx) ? i : i - Fweight.dimy;
                for (int j = 0; j < Fweight.dimx; j++)
                {
                    int jp = j;
                    int r2 = kp * kp + ip * ip + jp * jp;
                    if (r2 < max_r2)
                    {
                        int ires = round( sqrt((FDOUBLE)r2) / padding_factor );
                        FDOUBLE invw = Fweight(k, i, j);
                        
                        FDOUBLE invtau2;
                        if (tau2[ires] > 0.)
                        {
                            // Calculate inverse of tau2
                            invtau2 = 1. / (oversampling_correction * tau2_fudge * tau2[ires]);
                        }
                        else if (tau2[ires] == 0.)
                        {
                            // If tau2 is zero, use small value instead
                            invtau2 = 1./ ( 0.001 * invw);
                        }
                        else
                        {
                            std::cerr << " ires= "<<ires<<std::endl;
                            std::cerr << " sigma2= " << sigma2[ires] << std::endl;
                            std::cerr << " fsc= " << fsc[ires] << std::endl;
                            std::cerr << " tau2= " << tau2[ires] << std::endl;
                            ERROR_REPORT("ERROR BackProjector::reconstruct: Negative or zero values encountered for tau2 spectrum!");
                        }
                        
                        // Keep track of spectral evidence-to-prior ratio and remaining noise in the reconstruction
                        if (!update_tau2_with_fsc)
                            data_vs_prior[ires] += invw / invtau2;
                        counter[ires] += 1.;
                        
                        // Only for (ires >= minres_map) add Wiener-filter like term
                        if (ires >= minres_map)
                        {
                            // Now add the inverse-of-tau2_class term
                            invw += invtau2;
                            // Store the new weight again in Fweight
                            Fweight(k, i, j) = invw;
                        }
                    }
                }
        }}
        
        // Average data_vs_prior
        if (!update_tau2_with_fsc)
        {
            //FOR_ALL_DIRECT_ELEMENTS_IN_ARRAY1D(data_vs_prior)
            for (int i = 0; i < ori_Fsize; i++)
            {
                if (i > r_max)
                    data_vs_prior[i] = 0.;
                else if (counter[i] < 0.001)
                    data_vs_prior[i] = 999.;
                else
                    data_vs_prior[i] /= counter[i];
            }
        }
        
    } //end if do_map
    printSum(Fweight.wptr(),Fweight.dimzyx,"Fweight after do_map");
    printSum_real(data.wptr(),data.dimzyx,"data_real before normalise");
    printSum_imag(data.wptr(),data.dimzyx,"data_imag brfore normalise");
    // Divide both data and Fweight by normalisation factor to prevent FFT's with very large values....
#ifdef DEBUG_RECONSTRUCT
    std::cerr << " normalise= " << normalise << std::endl;
#endif
    //FOR_ALL_DIRECT_ELEMENTS_IN_MULTIDIMARRAY(Fweight)
    for (int k = 0; k < Fweight.dimz; k++){
        for (int i = 0; i < Fweight.dimy; i++){
#ifdef VECTORIZED_RECONS
            #pragma simd
            for (int j = 0; j < Fweight.dimx; j++)
            {
                ACCESS(Fweight,k,i,j) /= normalise;
            }
            #pragma simd
            for(int j = 0;j < Fweight.dimx; j++){
                ACCESS(data,k,i,j).real /= normalise;
                ACCESS(data,k,i,j).imag /= normalise;
            }
#else
            for (int j = 0; j < Fweight.dimx; j++)
            {
                Fweight(k,i,j) /= normalise;
                data(k,i,j).real /= normalise;
                data(k,i,j).imag /= normalise;
            }
#endif
    }}
        
    printSum_real(data.wptr(),data.dimzyx,"data_real after normalise");
    printSum_imag(data.wptr(),data.dimzyx,"data_imag after normalise");
    printSum(Fweight.wptr(),Fweight.dimzyx,"Fweight after normalise");
    
    if(tmp_folder != "NULL") data.writeToDisk(tmp_folder+"./data_tmp");
    
    // Initialise Fnewweight with 1's and 0's. (also see comments below)
    int weight_origin_z = XMIPP_ORIGIN(weight.dimz);
    int weight_origin_y = XMIPP_ORIGIN(weight.dimy);
    // FOR_ALL_ELEMENTS_IN_ARRAY3D(weight)
    for (int k = 0; k < weight.dimz; k++){
        for (int i = 0; i < weight.dimy; i++){
            int kp = k + weight_origin_z;
            int ip = i + weight_origin_y;
            for (int j = 0; j < weight.dimx; j++)
            {
                int jp = j;
                if (kp * kp + ip * ip + jp * jp < max_r2)
                    weight(k, i, j) = 1.;
                else
                    weight(k, i, j) = 0.;
            }
    }}
    
    decenter(weight, Fnewweight, max_r2);
    printSumDouble(Fnewweight.wptr(),Fnewweight.dimzyx,"Fnewweight after decenter");
    
    //
    if(tmp_folder != "NULL") weight.writeToDisk(tmp_folder+"./weight_tmp");
    if(tmp_folder != "NULL") Fweight.writeToDisk(tmp_folder+"./Fweight_tmp");
    
    TIMEPOINT
    // Iterative algorithm as in  Eq. [14] in Pipe & Menon (1999)
    // or Eq. (4) in Matej (2001)
    for (int iter = 0; iter < max_iter_preweight; iter++)
    {
        TIMEPOINT
        if(tmp_folder != "NULL") Fweight.readFromDisk(tmp_folder+"./Fweight_tmp");
        // Set Fnewweight * Fweight in the transformer
        // In Matej et al (2001), weights w_P^i are convoluted with the kernel,
        // and the initial w_P^0 are 1 at each sampling point
        // Here the initial weights are also 1 (see initialisation Fnewweight above),
        // but each "sampling point" counts "Fweight" times!
        // That is why Fnewweight is multiplied by Fweight prior to the convolution
        // FOR_ALL_DIRECT_ELEMENTS_IN_MULTIDIMARRAY(Fconv)
#ifdef VECTORIZED_RECONS
        #pragma simd
#endif
        for (int n = 0; n < Fconv.dimzyx; n++)
        {
            ACCESS(Fconv,0,0,n).real = ACCESS(Fnewweight,0,0,n) * ACCESS(Fweight,0,0,n);
            ACCESS(Fconv,0,0,n).imag = 0;
        }
		
        if(tmp_folder != "NULL") Fweight.fini();
        
        TIMEPOINT
        // convolute through Fourier-transform (as both grids are rectangular)
        // Note that convoluteRealSpace acts on the complex array inside the transformer
        convoluteBlobRealSpace(Fconv,transformer);
        TIMEPOINT
        FDOUBLE w, corr_min = (std::numeric_limits<FDOUBLE>::max)(), corr_max = (std::numeric_limits<FDOUBLE>::min)(), corr_avg=0., corr_nn=0.;
        // FOR_ALL_ELEMENTS_IN_FFTW_TRANSFORM(Fconv)
#pragma omp parallel for collapse(2) private(w) if(transformer.nr_threads > 1) num_threads(transformer.nr_threads)
        for (int k = 0; k < Fconv.dimz; k++){
            for (int i = 0; i < Fconv.dimy; i++){
                int kp = (k < Fconv.dimx) ? k : k - Fconv.dimz;
                int ip = (i < Fconv.dimx) ? i : i - Fconv.dimy;
                for (int j = 0; j < Fconv.dimx; j++)
                {
                    int jp = j;
                    if (kp * kp + ip * ip + jp * jp < max_r2)
                    {
                        // Make sure no division by zero can occur....
                        w = std::max(1e-6, sqrt(Fconv(k, i, j).real*Fconv(k, i, j).real+Fconv(k, i, j).imag*Fconv(k, i, j).imag));
#ifdef DEBUG_RECONSTRUCT
                        // Monitor min, max and avg conv_weight
                        corr_min = std::min(corr_min, w);
                        corr_max = std::max(corr_max, w);
                        corr_avg += w;
                        corr_nn += 1.;
#endif
                        // Apply division of Eq. [14] in Pipe & Menon (1999)
                        Fnewweight(k, i, j) /= w;
                    }
                }
        }}
        TIMEPOINT
#ifdef DEBUG_RECONSTRUCT
        ERROR_CHECK(transformer.nr_threads > 1, "Remove openmp for loop above for debug.");
        std::cerr << " PREWEIGHTING ITERATION: "<< iter + 1 << " OF " << max_iter_preweight << std::endl;
        // report of maximum and minimum values of current conv_weight
        std::cerr << " corr_avg= " << corr_avg / corr_nn << std::endl;
        std::cerr << " corr_min= " << corr_min << std::endl;
        std::cerr << " corr_max= " << corr_max << std::endl;
#endif
    }
	TIMEPOINT
    printSum_real(Fconv.wptr(),Fconv.dimzyx,"Fconv_real after iter");
    printSum_imag(Fconv.wptr(),Fconv.dimzyx,"Fconv_imag after iter");
    printSumDouble(Fnewweight.wptr(),Fnewweight.dimzyx,"Fnewweight after iter");
#ifdef DEBUG_RECONSTRUCT
    // Image<double> tttt;
    // tttt()=Fnewweight;
    // tttt.write("reconstruct_gridding_weight.spi");
    // FOR_ALL_DIRECT_ELEMENTS_IN_MULTIDIMARRAY(Fconv)
    // {
    //     DIRECT_MULTIDIM_ELEM(ttt(), n) = abs(DIRECT_MULTIDIM_ELEM(Fconv, n));
    // }
    // ttt.write("reconstruct_gridding_correction_term.spi");
#endif
    TIMEPOINT
    // Clear memory
    Fweight.fini();
    
    if(tmp_folder != "NULL") data.readFromDisk(tmp_folder+"./data_tmp");
        
    // Note that Fnewweight now holds the approximation of the inverse of the weights on a regular grid
    
    // Now do the actual reconstruction with the data array
    // Apply the iteratively determined weight
    // to remove any stuff from the input volume
    // Fconv.zero();
    decenter(data, Fconv, max_r2);
    printSum_real(data.wptr(),data.dimzyx,"data_real after decenter");
    printSum_imag(data.wptr(),data.dimzyx,"data_imag after decenter");
    
    // if(tmp_folder != "NULL") data.writeToDisk(tmp_folder+"./data_tmp");
    if(tmp_folder != "NULL") data.fini();
    // FOR_ALL_DIRECT_ELEMENTS_IN_MULTIDIMARRAY(Fconv)
#ifdef VECTORIZED_RECONS
    #pragma simd
#endif
    for (int n = 0; n < Fconv.dimzyx; n++)
    {
#ifdef  FLOAT_PRECISION
        // Prevent numerical instabilities in single-precision reconstruction with very unevenly sampled orientations
        if (Fnewweight(0,0,n) > 1e20)
            Fnewweight(0,0,n) = 1e20;
#endif
        ACCESS(Fconv,0,0,n).real *= ACCESS(Fnewweight,0,0,n);
        ACCESS(Fconv,0,0,n).imag *= ACCESS(Fnewweight,0,0,n);
    }
    printSum_real(Fconv.wptr(),Fconv.dimzyx,"Fconv_real after multiply");
    printSum_imag(Fconv.wptr(),Fconv.dimzyx,"Fconv_imag after multiply");
    // Clear memory
    Fnewweight.fini();
    TIMEPOINT
    // TODO ?????
    // Gridding theory says one now has to interpolate the fine grid onto the coarse one using a blob kernel
    // and then do the inverse transform and divide by the FT of the blob (i.e. do the gridding correction)
    // In practice, this gives all types of artefacts (perhaps I never found the right implementation?!)
    // Therefore, window the Fourier transform and then do the inverse transform
//    #define RECONSTRUCT_CONVOLUTE_BLOB
#ifdef RECONSTRUCT_CONVOLUTE_BLOB
    
    // Apply the same blob-convolution as above to the data array
    // Mask real-space map beyond its original size to prevent aliasing in the downsampling step below
    convoluteBlobRealSpace(Fconv, transformer, true);
    printSum_real(Fconv.wptr(),Fconv.dimzyx,"Fconv_real after Blob");
    printSum_imag(Fconv.wptr(),Fconv.dimzyx,"Fconv_imag after Blob");
    
    // Now just pick every 3rd pixel in Fourier-space (i.e. down-sample)
    // and do a final inverse FT
    // vol_out should be right size
    // if (ref_dim == 2)
    //     vol_out.resize(1, ori_size, ori_size);
    // else
    //     vol_out.resize(ori_size, ori_size, ori_size);
    
    Vol<FDOUBLE> tmp;
    Vol<MKL_Complex> Ftmp;
    if (ref_dim == 3) {
        Ftmp.init(ori_size, ori_size, ori_Fsize);
        tmp.init(ori_size, ori_size, ori_size);
    }
    else{
        Ftmp.init(1, ori_size, ori_Fsize);
        tmp.init(1, ori_size, ori_size);
    }
    transformer.reset_plan(tmp.wptr(), (FourierComplex*)Ftmp.wptr(), tmp.dimx, tmp.dimy, tmp.dimz);
    
    for (int k = 0; k<Ftmp.dimz; k++){
        for (int i = 0; i<Ftmp.dimy; i++){
            int kp = (k < Ftmp.dimx) ? k : k - Ftmp.dimz;
            int ip = (i < Ftmp.dimx) ? i : i - Ftmp.dimy;
            int k2 = (kp*padding_factor) < 0 ? (kp*padding_factor)+Fconv.dimz : (kp*padding_factor);
            int i2 = (ip*padding_factor) < 0 ? (ip*padding_factor)+Fconv.dimy : (ip*padding_factor);
            for (int j = 0; j<Ftmp.dimx; j++)
            {
                int jp = j;
                int j2 = jp*padding_factor;
                if (kp * kp + ip * ip + jp * jp < r_max * r_max)
                {
                    Ftmp(k, i, j).real = Fconv(k2, i2, j2).real;
                    Ftmp(k, i, j).imag = Fconv(k2, i2, j2).imag;
                }
                else
                {
                    Ftmp(k, i, j).real = Ftmp(k, i, j).imag = 0.;
                }
            }
    	}
    }
    
    printSum_real(Ftmp.wptr(),Ftmp.dimzyx,"Ftmp_real after set");
    printSum_imag(Ftmp.wptr(),Ftmp.dimzyx,"Ftmp_imag after set");
    
    // inverse FFT leaves result in vol_out
    transformer.inverseFourierTransform();
    printSum(tmp.wptr(),tmp.dimzyx,"tmp after fft");
    
    // Shift the map back to its origin
    if(ref_dim==3)
        centerFFT3D(tmp.wptr(), vol_out.wptr(), ori_size, false);
    else
        centerFFT2D(tmp.wptr(), vol_out.wptr(), ori_size, false);
    
    printSum(vol_out.wptr(),vol_out.dimzyx,"vol_out after center");
    
    // Un-normalize FFTW (because original FFTs were done with the size of 2D FFTs)
    if (ref_dim==3){
        for (int n = 0; n < vol_out.dimzyx; n++) {
            ACCESS(vol_out,0,0,n) /= ori_size;
        }
    }
    // Mask out corners to prevent aliasing artefacts
    softMaskOutsideMap(vol_out);
    
    printSum(vol_out.wptr(),vol_out.dimzyx,"vol_out after softMaskOutsideMap");
    
    // Gridding correction for the blob
    FDOUBLE normftblob = tab_ftblob(0.);
    // FOR_ALL_ELEMENTS_IN_ARRAY3D(vol_out)
    int vol_out_shift_z = XMIPP_ORIGIN(vol_out.dimz);
    int vol_out_shift_y = XMIPP_ORIGIN(vol_out.dimy);
    int vol_out_shift_x = XMIPP_ORIGIN(vol_out.dimx);
    for (int k = 0; k < vol_out.dimz; k++){
        for (int i = 0; i < vol_out.dimy; i++){
            for (int j = 0; j < vol_out.dimx; j++)
            {
                int kp = k + vol_out_shift_z;
                int ip = i + vol_out_shift_y;
                int jp = j + vol_out_shift_x;
                FDOUBLE r = sqrt((FDOUBLE)(kp*kp+ip*ip+jp*jp));
                FDOUBLE rval = r / (ori_size * padding_factor);
                vol_out(k, i, j) /= tab_ftblob(rval) / normftblob;
                //if (k==0 && i==0)
                //	std::cerr << " j= " << j << " rval= " << rval << " tab_ftblob(rval) / normftblob= " << tab_ftblob(rval) / normftblob << std::endl;
            }
    	}
    }
    tmp.fini();Ftmp.fini();
    printSum(vol_out.wptr(),vol_out.dimzyx,"vol_out after blob.");
    
#else
    // rather than doing the blob-convolution to downsample the data array, do a windowing operation:
    // This is the same as convolution with a SINC. It seems to give better maps.
    // Then just make the blob look as much as a SINC as possible....
    // The "standard" r1.9, m2 and a15 blob looks quite like a sinc until the first zero (perhaps that's why it is standard?)
    //for (FDOUBLE r = 0.1; r < 10.; r+=0.01)
    //{
    //	FDOUBLE sinc = sin(PI * r / padding_factor ) / ( PI * r / padding_factor);
    //	std::cout << " r= " << r << " sinc= " << sinc << " blob= " << blob_val(r, blob) << std::endl;
    //}
    TIMEPOINT
    // Now do inverse FFT and window to original size in real-space
    // Pass the transformer to prevent making and clearing a new one before clearing the one declared above....
    // The latter may give memory problems as detected by electric fence....
    windowToOridimRealSpace(Fconv, vol_out,transformer);
    TIMEPOINT
#endif
    printSum(vol_out.wptr(),vol_out.dimzyx,"vol_out");
#ifdef DEBUG_RECONSTRUCT
    // ttt()=vol_out;
    // ttt.write("reconstruct_before_gridding_correction.spi");
#endif
    
    // Correct for the linear/nearest-neighbour interpolation that led to the data array
    griddingCorrect(vol_out);
    
    printSum(vol_out.wptr(),vol_out.dimzyx,"vol_out after gridding correct");
    // If the tau-values were calculated based on the FSC, then now re-calculate the power spectrum of the actual reconstruction
    if (update_tau2_with_fsc)
    {
        // New tau2 will be the power spectrum of the new map
        std::vector<FDOUBLE> spectrum(vol_out.dimx,0);
        std::vector<FDOUBLE> count(vol_out.dimx,0);
        // Calculate this map's power spectrum
        // Don't call getSpectrum() because we want to use the same transformer object to prevent memory trouble....
        // recycle the same transformer for all images
        Fconv.fini();Fconv.init_nonzero(ori_size, ori_size, ori_Fsize);
        transformer.reset_plan(vol_out.wptr(), (FourierComplex*)Fconv.wptr(), vol_out.dimx, vol_out.dimy, vol_out.dimz);
        transformer.FourierTransform();
        for (int k = 0; k < Fconv.dimz; k++){
            for (int i = 0; i < Fconv.dimy; i++){
                int kp = (k < Fconv.dimx) ? k : k - Fconv.dimz;
                int ip = (i < Fconv.dimx) ? i : i - Fconv.dimy;
                for (int j = 0; j < Fconv.dimx; j++)
                {
                    int jp = j;
                    int idx = round(sqrt(kp*kp + ip*ip + jp*jp));
                    spectrum[idx] += ACCESS(Fconv, k, i, j).real*ACCESS(Fconv, k, i, j).real+ACCESS(Fconv, k, i, j).imag*ACCESS(Fconv, k, i, j).imag;
                    count[idx] += 1.;
                }
        	}
        }
        for(int i = 0;i < ori_Fsize;i++)
            spectrum[i] /= count[i];
        
        // Factor two because of two-dimensionality of the complex plane
        // (just like sigma2_noise estimates, the power spectra should be divided by 2)
        FDOUBLE normfft = (ref_dim == 3) ? (FDOUBLE)(ori_size * ori_size) : 1.;
        for(int i = 0;i < ori_Fsize;i++)
            spectrum[i] *= normfft / 2.;
        
        // New SNR^MAP will be power spectrum divided by the noise in the reconstruction (i.e. sigma2)
        for (int n = 0; n < ori_Fsize; n++)
        {
            tau2[n] =  tau2_fudge * spectrum[n];
        }
    }
    TIMEPOINT
    // Completely empty the transformer object
    transformer.fini();
    
    if(tmp_folder != "NULL") {
        data.readFromDisk(tmp_folder+"./data_tmp");
        weight.readFromDisk(tmp_folder+"./weight_tmp");
    }
    
    TIMEPOINT_FINA
#ifdef DEBUG_RECONSTRUCT
    std::cerr<<"done with reconstruct"<<std::endl;
#endif
    
}

void MyBackProjector::getLowResDataAndWeight(Vol<MKL_Complex > &lowres_data, Vol<FDOUBLE> &lowres_weight,int lowres_r_max)
{
    int lowres_r2_max = padding_factor * padding_factor * lowres_r_max * lowres_r_max;
    int lowres_pad_size = 2 * (padding_factor * lowres_r_max + 1) + 1;
    
    // Check for dimension
    if (ref_dim != 3)
        ERROR_REPORT("MyBackProjector::getLowResDataAndWeight%%ERROR: only implemented for 3D case....");
    
    // Check lowres_r_max is not too big
    if (lowres_r_max > r_max)
        ERROR_REPORT("MyBackProjector::getLowResDataAndWeight%%ERROR: lowres_r_max is bigger than r_max");
    
    // Initialize lowres_data and low_res_weight arrays
    lowres_data.init_zero(lowres_pad_size, lowres_pad_size, lowres_pad_size/2+1);
    lowres_weight.init_zero(lowres_pad_size, lowres_pad_size, lowres_pad_size/2+1);
    
    // fill lowres arrays with relevant values
    int lowres_data_origin_z = XMIPP_ORIGIN(lowres_data.dimz);
    int lowres_data_origin_y = XMIPP_ORIGIN(lowres_data.dimy);
    int data_origin_z = XMIPP_ORIGIN(data.dimz);
    int data_origin_y = XMIPP_ORIGIN(data.dimy);
    for (int k = 0; k < lowres_data.dimz; k++){
        for (int i = 0; i < lowres_data.dimy; i++){
            int kp = k + lowres_data_origin_z;
            int ip = i + lowres_data_origin_y;
            for (int j = 0; j < lowres_data.dimx; j++)
            {
                int jp = j;
                if (kp*kp + ip*ip + jp*jp <= lowres_r2_max)
                {
                    int k2 = kp - data_origin_z;
                    int i2 = ip - data_origin_y;
                    int j2 = jp;
                    ACCESS(lowres_data, k, i, j) = ACCESS(data, k2 , i2, j2);
                    ACCESS(lowres_weight, k, i, j) = ACCESS(weight, k2 , i2, j2);
                }
            }
        }
    }
}

void MyBackProjector::setLowResDataAndWeight(Vol<MKL_Complex > &lowres_data, Vol<FDOUBLE> &lowres_weight,int lowres_r_max)
{
    int lowres_r2_max = padding_factor * padding_factor * lowres_r_max * lowres_r_max;
    int lowres_pad_size = 2 * (padding_factor * lowres_r_max + 1) + 1;
    
    // Check for dimension
    if (ref_dim != 3)
        ERROR_REPORT("MyBackProjector::getLowResDataAndWeight%%ERROR: only implemented for 3D case....");
    
    // Check lowres_r_max is not too big
    if (lowres_r_max > r_max)
        ERROR_REPORT("MyBackProjector::getLowResDataAndWeight%%ERROR: lowres_r_max is bigger than r_max");
    
    // Check sizes of lowres_data and lowres_weight
    if (lowres_data.dimz != lowres_pad_size || lowres_data.dimy != lowres_pad_size || lowres_data.dimx != lowres_pad_size / 2 + 1)
        ERROR_REPORT("MyBackProjector::setLowResDataAndWeight%%ERROR: lowres_data is not of expected size...");
    if (lowres_weight.dimz != lowres_pad_size || lowres_weight.dimy != lowres_pad_size || lowres_weight.dimx != lowres_pad_size / 2 + 1)
        ERROR_REPORT("MyBackProjector::setLowResDataAndWeight%%ERROR: lowres_weight is not of expected size...");
    
    // Re-set origin to the expected place
    
    // Overwrite data and weight with the lowres arrays
    int lowres_data_origin_z = XMIPP_ORIGIN(lowres_data.dimz);
    int lowres_data_origin_y = XMIPP_ORIGIN(lowres_data.dimy);
    int data_origin_z = XMIPP_ORIGIN(data.dimz);
    int data_origin_y = XMIPP_ORIGIN(data.dimy);
    for (int k = 0; k < lowres_data.dimz; k++) {
        for (int i = 0; i < lowres_data.dimy; i++) {
            int kp = k + lowres_data_origin_z;
            int ip = i + lowres_data_origin_y;
            for (int j = 0; j < lowres_data.dimx; j++) {
                int jp = j;
                if (kp*kp + ip*ip + jp*jp <= lowres_r2_max)
                {
                    int k2 = kp - data_origin_z;
                    int i2 = ip - data_origin_y;
                    int j2 = jp;
                    ACCESS(data, k2, i2, j2) = ACCESS(lowres_data, k , i, j);
                    ACCESS(weight, k2, i2, j2) = ACCESS(lowres_weight, k , i, j);
                }
            }
        }
    }
}

void MyBackProjector::getDownsampledAverage(Vol<MKL_Complex > &avg)
{
    Vol<FDOUBLE> down_weight;
    // Pre-set down_data and down_weight sizes
    int down_size = 2 * (r_max + 1) + 1;
    int r2_max = r_max * r_max;
    // Short side of data array
    ERROR_CHECK(ref_dim!=3, "MyBackProjector::getDownsampledAverage%%ERROR: Dimension of the data array should be 3");
    avg.init_zero(down_size,down_size,down_size/2+1);
    // Resize down_weight the same as down_data
    down_weight.init_zero(down_size,down_size,down_size/2+1);
    
    // Now calculate the down-sized sum
    int data_origin_z = XMIPP_ORIGIN(data.dimz);
    int data_origin_y = XMIPP_ORIGIN(data.dimy);
    int avg_origin_z = XMIPP_ORIGIN(avg.dimz);
    int avg_origin_y = XMIPP_ORIGIN(avg.dimy);
    for (int k = 0; k < data.dimz; k++){
        for (int i = 0; i < data.dimy; i++) {
            int k2 = round((FDOUBLE)(k + data_origin_z)/padding_factor)-avg_origin_z;
            int i2 = round((FDOUBLE)(i + data_origin_y)/padding_factor)-avg_origin_y;
            for (int j = 0; j < data.dimx; j++)
            {
                int j2 = round((FDOUBLE)(j+0)/padding_factor)+0;
                ACCESS(avg, k2, i2, j2).real += ACCESS(data, k, i, j).real;
                ACCESS(avg, k2, i2, j2).imag += ACCESS(data, k, i, j).imag;
                ACCESS(down_weight, k2, i2, j2) += ACCESS(weight, k , i, j);
            }
        }
    }
    // Then enforce Hermitian symmetry in the downsampled arrays
    // We already took the average.... so not completely correct, but does not really matter for FSC calculation anyway
    // enforceHermitianSymmetry(avg, down_weight);
    
    // And enforce symmetry in the downsampled arrays
    symmetrise(avg, down_weight, r2_max);
    
    // Calculate the straightforward average in the downsampled arrays
    for (size_t n = 0; n < avg.dimzyx; n++)
    {
        if (ACCESS(down_weight, 0, 0, n) > 0.){
            ACCESS(avg, 0, 0, n).real /= ACCESS(down_weight, 0, 0, n);
            ACCESS(avg, 0, 0, n).imag /= ACCESS(down_weight, 0, 0, n);
        }
        else{
            ACCESS(avg, 0, 0, n).real = ACCESS(avg, 0, 0, n).imag = 0.;
        }
    }
    
    down_weight.fini();
}

void MyBackProjector::calculateDownSampledFourierShellCorrelation(Vol<MKL_Complex > &avg1,
                                                                  Vol<MKL_Complex > &avg2,
                                                                  FDOUBLE* fsc)
{
    if (avg1.dimx!=avg2.dimx || avg1.dimy!=avg2.dimy || avg1.dimz!=avg2.dimz)
        ERROR_REPORT("ERROR MyBackProjector::calculateDownSampledFourierShellCorrelation: two arrays have different sizes");
    
    std::vector<FDOUBLE> num(ori_size/2+1,0);
    std::vector<FDOUBLE> den1(ori_size/2+1,0);
    std::vector<FDOUBLE> den2(ori_size/2+1,0);
    
    int avg1_origin_z = XMIPP_ORIGIN(avg1.dimz);
    int avg1_origin_y = XMIPP_ORIGIN(avg1.dimy);
    for (int k = 0; k < avg1.dimz; k++){
        for (int i = 0; i < avg1.dimy; i++){
            int kp = k + avg1_origin_z;
            int ip = i + avg1_origin_y;
            for (int j = 0;j < avg1.dimx; j++)
            {
                int jp = j;
                FDOUBLE R = sqrt(kp*kp + ip*ip + jp*jp);
                if (R > r_max)
                    continue;
                int idx=round(R);
                auto& z1=ACCESS(avg1, k, i, j);
                auto& z2=ACCESS(avg2, k, i, j);
                FDOUBLE absz1=sqrt(z1.real*z1.real+z1.imag*z1.imag);
                FDOUBLE absz2=sqrt(z2.real*z2.real+z2.imag*z2.imag);
                num[idx]  += z1.real*z2.real+z1.imag*z2.imag;//(conj(z1) * z2).real;
                den1[idx] += absz1*absz1;
                den2[idx] += absz2*absz2;
            }
    	}
    }
    
    for (int i = 0; i < ori_Fsize; i++)
    {
        if (den1[i]*den2[i] > 0.)
            fsc[i] = num[i]/sqrt(den1[i]*den2[i]);
    }
    // Always set zero-resolution shell to FSC=1
    // Raimond Ravelli reported a problem with FSC=1 at res=0 on 13feb2013...
    // (because of a suboptimal normalisation scheme, but anyway)
    fsc[0] = 1.;
    
}