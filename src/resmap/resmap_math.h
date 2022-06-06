/***************************************************************************
 *
 * Authors: "Yongbei(Galow) Ma,Jiayi (Timmy) Wu, Youdong (Jack) Mao"
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

#ifndef _ROME_MATH
#define _ROME_MATH

#include <iostream>
#include <cmath>

#include "./resmap_macros.h"
#include "./resmap_matrix.h"

//
#if !defined(FLT_EPSILON)
#define FLT_EPSILON 1.19209e-07
#endif
//
#define MAP3D_OLD_BASELINE

//
static const double rome_pi = 3.14159265358979323846;
// Radians to degrees
inline double rad2deg(double r) {return r*180./rome_pi;}
//
inline double deg2rad(double d) {return d*rome_pi/180.;}
// ArcCosine in degrees
inline double acosd(double x){return acos(x)*180./rome_pi;}
//
inline double sgn(double x) {return (x >= 0) ? 1 : -1;}
// Wrapping for real numbers between a range
// example : Corrected_angle = realWRAP(angle, 0, 2*PI);
inline double realWrap(double x,double x0,double xF){
    if (x >= x0 && x <= xF)
        return x;
    else if (x < x0)
        return ((x) - (int)(((x) - (x0)) / ((xF) - (x0)) - 1) * ((xF) - (x0)));
    else
        return ((x) - (int)(((x) - (xF)) / ((xF) - (x0)) + 1) * ((xF) - (x0)));
}


/***************************************************************************
 *
 * Authors:     Carlos Oscar S. Sorzano (coss@cnb.csic.es)
 *
 * Unidad de  Bioinformatica of Centro Nacional de Biotecnologia , CSIC
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
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA
 * 02111-1307  USA
 *
 *  All comments concerning this program package may be sent to the
 *  e-mail address 'xmipp@cnb.csic.es'
 ***************************************************************************/

// 1D gaussian value
// This function returns the value of a univariate gaussian function at the
// point x.
FDOUBLE gaussian1D(FDOUBLE x, FDOUBLE sigma, FDOUBLE mu = 0);

// Compute statistics.
//
// The average, standard deviation, minimum and maximum value are
// returned.
void computeStats(const FDOUBLE* data,int size,FDOUBLE& avg, FDOUBLE& stddev, FDOUBLE& minval, FDOUBLE& maxval);

// -------------- euler function -------------- //
// Apply a transformation matrix to Euler angles --------------------------- //
void Euler_apply_transf(const Matrix2D<FDOUBLE> &L,const Matrix2D<FDOUBLE> &R,
                        FDOUBLE rot,FDOUBLE tilt,FDOUBLE psi,
                        FDOUBLE &newrot,FDOUBLE &newtilt,FDOUBLE &newpsi);
// Euler angles --> matrix
void Euler_angles2matrix(FDOUBLE alpha, FDOUBLE beta, FDOUBLE gamma, FDOUBLE A[][3]);
void Euler_angles2matrix(FDOUBLE alpha, FDOUBLE beta, FDOUBLE gamma, Matrix2D<FDOUBLE>& A);

// Matrix --> Euler angles
void Euler_matrix2angles(const FDOUBLE A[][3], FDOUBLE &alpha,FDOUBLE &beta, FDOUBLE &gamma);
void Euler_matrix2angles(const Matrix2D<FDOUBLE>& A, FDOUBLE &alpha,FDOUBLE &beta, FDOUBLE &gamma);
//
void Euler_angles2direction(FDOUBLE alpha, FDOUBLE beta, FDOUBLE v[]);

// Euler direction2angles ------------------------------- //
// gamma is useless but I keep it for simmetry
// with Euler_direction
void Euler_direction2angles(const FDOUBLE v0[],
                            FDOUBLE &alpha, FDOUBLE &beta);

//
class Random_generator {
public:
	Random_generator();
	~Random_generator();

    void init(int seed = -1);
    //
    float rnd_unif(float a = 0., float b = 1.);
    //
    float rnd_gaus(float mu, float sigma);
private:
	class Internals;
	Internals* _internals;
	typedef void rand;
	typedef void srand;
};

extern Random_generator dontShare_Random_generator;

#endif /* defined(_ROME_MATH) */
