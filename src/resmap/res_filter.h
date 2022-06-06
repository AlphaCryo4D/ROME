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

#include "res_array.h"
#include "res_fft.h"

namespace resumap {

    /**
     * Same as scipy.ndimage.filters.correlate1d.
     */
    int correlate1d(Arrayd &input, Arrayd &output, const Arrayd &weights,
            int axis = 0, int mode_code = 2, double cval = 0, int origin = 0);

    /**
     * Same as scipy.ndimage.filters.gaussian_filter.
     */
    Arrayd gaussian_filter(const Arrayd &, double sigma, int order=0,
            std::string mode="reflect", double cval=0.0, double truncate=4.0);

    /**
     * Deprecated!!! Filter in frequency domain. 
     */
    inline Arrayd freqFilt3(const Arrayd &in, const Arrayd &filt) {
        return Arrayd::real(Arraycd(in).selfFFT3().selfFFTshift3().selfTimes(filt).selfFFTshift3());
    }

    /**
     * Deprecated!!! Low-pass gaussian filter.
     */
    inline Arrayd lpfGauss(const Shape &shape, double sigma) {
        auto && out = Arrayd::ones(shape);
        int dim = int(shape.size());
        int i = 0;
        ARRAY_EACH(shape, ind) {
            double sum = 0;
            for (int j = 0; j < dim; j++) {
                double n = ind[j] - shape[j]/2.0;
                sum += n*n;
            }
            out[i] = std::exp(-sum/(2*sigma*sigma));
            i++;
        }
        return std::move(out);
    }

}

