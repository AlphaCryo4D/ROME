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

#include "res_array.h"
#include "res_fft.h"

namespace cppsci { namespace ndimage {

//    Arrayd zoom(const Arrayd &input, double _zoom, std::string _mode="reflect", int order = 3, double cval=0.0, bool prefilter=true);

//    void spline_filter(Arrayd &input, int order=3);

    Arrayd zoom(const Arrayd &input, double _zoom, std::string _mode="reflect", int order = 3, double cval=0.0, bool prefilter=true);

    void spline_filter(Arrayd &input, int order = 3);

    Arrayd map_coordinates(const Arrayd &input, /*const*/ Arrayd &coordinates, int order=3,
        std::string mode="constant", double cval=0.0, bool prefilter=true);

}} // end namespace cppsci::ndimage

