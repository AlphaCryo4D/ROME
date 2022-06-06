/***************************************************************************
 *
 * Intel® Parallel Computing Center for Structural Biology
 * Principal Investigator : Youdong (Jack) Mao (Youdong_Mao@dfci.harvard.edu)
 * Dana-Farber Cancer Institute, Harvard Medical School and Peking University
 *
 * Authors: "Jian Wang(wj_hust08@hust.edu.cn)"
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

#pragma once

#include "resmap_mpi.h"
#include "resmap_mrcs.h"
#include "resmap_time.h"

#include "res_fft.h"
#include "res_filter.h"
#include "res_spectrum.h"
#include "res_morph.h"
#include "res_utils.h"
#include "res_ndimage.h"
#include "res_linalg.h"
#include "res_optimize.h"
#include "res_chimera.h"


namespace resumap
{
    extern std::string inputFileName;
    extern std::string inputFileName1;
    extern std::string inputFileName2;
    extern double       vxSize;
    extern double       pValue;
    extern double       minRes;
    extern double       maxRes;
    extern double       stepRes;
//    extern double       dataMask;
    extern double       variance;
    extern bool        noiseDiagnostics;

    void calculate();
};

