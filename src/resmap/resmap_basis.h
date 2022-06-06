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


#ifndef BASIS_H_
#define BASIS_H_

#include <iostream>
#include <cmath>

#include "resmap_macros.h"

namespace GTMBasis
{	
	void set_shift_basis(const double *X,int K,int L,const double *MU,const double *MU_infos,int Mnl,double *PHI,int M);

	void set_shift_basis2(const double *X,int K,int L,const double *MU,const double *MU_infos,int Mnl,double *PHI,int M);

	void set_rotation_basis(const double *X,int K,int L,const double *MU,const double *MU_infos,int Mnl,double *PHI,int M);

	void set_mutilcfm_rotation_basis(const double *X,int K,int L,const double *MU,const double *MU_infos,int Mnl,double *PHI,int M);

	void set_shift_rotation_basis(const double *X,int K,int L,const double *MU,const double *MU_infos,int Mnl,double *PHI,int M);

	void set_oneshift_rotation_basis(const double *X,int K,int L,const double *MU,const double *MU_infos,int Mnl,double *PHI,int M);

	void set_2agnle_healpix_basis(const double *X,int K,int L,const double *MU,double sigma1,double sigma2,int Mnl,double *PHI,int M);

	void set_2angle_basis(const double *X,int K,int L,const double *MU,const double *MU_infos,int Mnl,double *PHI,int M,bool period = false);

	void set_3angle_basis(const double *X,int K,int L,const double *MU,const double *MU_infos,int Mnl,double *PHI,int M);
}

#endif