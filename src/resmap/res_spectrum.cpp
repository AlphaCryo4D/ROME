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

#include "res_spectrum.h"

namespace resumap {

using spectrum_type = Arrayd;

spectrum_type calculatePowerSpectrum(const Arraycd &m) {
    double epsilon = 1e-10;

    Arrayd dataFabs(m.shape());

    const std::complex<double> *it1 = m.data_rptr();
    double *it2 = dataFabs.data_wptr();

    double min = std::numeric_limits<double>::max();
    double max = std::numeric_limits<double>::min();
    for (int i = 0; i < m.size(); i++) {
        it2[i] = std::abs(it1[i]);
        if (min > it2[i]) min = it2[i];
        if (max < it2[i]) max = it2[i];
    }

    double d = max - min;
    for (int i = 0; i < m.size(); i++) {
        it2[i] = (it2[i] - min) / d;
        it2[i] = (it2[i]) * (it2[i]);
    }

//    auto && dataFabs = Arrayd::abs(m);
//    dataFabs.selfMinus(dataFabs.min());
//    dataFabs.selfDivide(dataFabs.max());
//    dataFabs.selfSquare();
    auto && ls = sphericalAverage(dataFabs);
    return std::move(ls);
}

double isPowerSpectrumLPF(const spectrum_type &m) {
    int l = m.size();

    auto diff = Arrayd::diff(Arrayd::log(m)).selfTimes(-1);

    auto ind = cppsci::signal::find_peaks_cwt(diff, Arrayi::range(1,10), 2);

    if (ind.size() == 0) return 0;
    int max = *std::max_element(ind.begin(), ind.end());
    if (max >= l - 3) return 0;
    auto && sub = diff.sub({{max+2,l}});
    double mean = sub.mean();
    double var = sub.var();
    double thr = 1e-4;
    if (std::abs(mean) < thr && var < thr) return double(max-1)/l;
    else return 0;

}

} // namespace resumap


