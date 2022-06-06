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

namespace cppsci { namespace signal {

    inline Arrayd ricker(int points, double a) {
        double A = 2 / (std::sqrt(3 * a) * (std::pow(PI,0.25)));
        double wsq = a*a;
        Arrayd vec = Arrayd::range(0, points) - (points - 1.0) / 2;
        Arrayd xsq = std::pow(vec,2);
        Arrayd mod = (1 - xsq / wsq);
        Arrayd gauss = std::exp(xsq / (-2 * wsq));
        Arrayd total = A * mod * gauss;
        return std::move(total);
    }

    template<typename T>
    void reverse(Array<T> &in) {
        Shape shape = in.shape();
        Shape v(in.dim(), 0);
        auto symm = [&shape](std::vector<int> v) -> std::vector<int> {
            for (int i = 0; i < v.size(); i++) {
                v[i] = shape[i]-1-v[i];
            }
            return std::move(v);
        };
        ARRAY_EACH(shape, ind) {
            if (ind[0] >= shape[0]/2) break;
            std::swap(in(ind), in(symm(ind)));
        }
    }

    enum {
        CONVOLVE_SAME,
        CONVOLVE_FULL
    };

    inline Arrayd convolve_full(const Arrayd &x, Arrayd y) {
        ERROR_CHECK(x.dim()!=y.dim(), "The dimensions of x and y should be same!");
        int nd = x.dim();
        reverse(y);
        Shape shapex = x.shape();
        Shape shapey = y.shape();
        Shape shape(nd);
        for (int i = 0; i < nd; i++) shape[i] = shapex[i]+shapey[i]-1;
        Arrayd x1(shape, 0), y1(shape, 0);
        for (int i = 0; i < y.size(); i++) y1[i] = y[i];
        int i = 0;
        Shape v(nd);
        ARRAY_EACH(shapex, ind) {
            for (int j = 0; j < nd; j++) v[j] = ind[j]+shapey[j]-1;
            x1(v) = x[i];
            i++;
        }
        return Arrayd::convolve(x1, y1);
    }

    inline Arrayd convolve_same(const Arrayd &x, const Arrayd &y) {
        ERROR_CHECK(x.dim()!=y.dim(), "The dimensions of x and y should be same!");
        int nd = x.dim();
        auto && output = convolve_full(x, y);
        std::vector<std::vector<int>> v(nd);
        for (int i = 0; i < nd; i++) {
            int n = int(std::ceil(y.shape(i)/2.0))-1;
            v[i].push_back(n);
            v[i].push_back(n+x.shape(i));
        }
        return output.sub(v);
    }

    inline Arrayd convolve(const Arrayd &x, const Arrayd &y, int mode=CONVOLVE_FULL) {
        if (mode != CONVOLVE_SAME && mode != CONVOLVE_FULL) {
            ERROR_REPORT("Only 'CONVOLVE_SAME' and 'CONVOLVE_FULL' are supported now!");
        }
        if (mode == CONVOLVE_SAME) return convolve_same(x, y);
        else return convolve_full(x, y);
    }

    template<typename _F = decltype(ricker)>
    Arrayd cwt(const Arrayd &data, const Arrayd &widths, _F && wavelet=ricker) {
        int m = widths.size(), n = data.size();
//        MPI_LOG << "m: " << m << ", n: " << n << std::endl;
        Arrayd output({m, n}, 0);
        for (int i = 0; i < m; i++) {
            auto && wavelet_data = wavelet(std::min(10 * int(widths[i]), n), widths[i]);
//            wavelet_data.identify("wavelet_data");
//            data.identify("data");
            output.sub(i) = convolve(data, wavelet_data, CONVOLVE_SAME);
        }
        return std::move(output);
    }

template<typename _Array, typename _Comp>
Arrayb boolrelextrema(_Array && data, _Comp && comparator, int axis=0, int order=1, int mode=ARRAY_TAKE_MODE_CLIP) {
    ERROR_CHECK(order < 1, "Order must be an int >= 1");

    int datalen = data.shape(axis);
    auto && locs = Arrayi::range(0, datalen);

    Arrayb results(data.shape(), true);
    auto && main = data.takes(locs, axis, mode);
    for (int shift = 1; shift < order+1; shift++) {
        auto && plus = data.takes(locs + shift, axis, mode);
        auto && minus = data.takes(locs - shift, axis, mode);
        results.selfBitAnd(Arrayb::map(comparator, main, plus));
        results.selfBitAnd(Arrayb::map(comparator, main, minus));
        if(!results.any()) {
            return results;
        }
    }
    return results;
}

using Ridge = std::vector<std::array<int, 2>>;

inline Arrayd array_cwt(const Arrayd &ms, const Arrayi &scales) {
    Arrayd psi_xval, psi;
    int nPoints = 1024;
    int n = 8;
    psi_xval = Arrayd::range(-n, n, 0, nPoints+1);
    JN_INFOA(psi_xval);
    auto && psi_xval2 = Arrayd::pow(psi_xval, 2.0);
    JN_INFOA(psi_xval2);
    psi = (2.0/std::sqrt(3.0)*std::pow(PI, -0.25)) * (1 - psi_xval2) * std::exp(-0.5 * psi_xval2);
    JN_INFOA(psi);

    int len = ms.size();
    Arrayd ms2({2*len});
    for (int i = 0; i < len; i++) ms2[i] = ms[i];
    for (int i = 0; i < len; i++) ms2[2*len-1-i] = ms[i];
    JN_INFOA(ms2);
    JN_INFOV(len);
    double scale;
    Arrayd f, wCoefs({scales.size(), len*2}, 0);
    for (int i = 0; i < scales.size(); i++) {
        //MPI_LOG << i << std::endl;
        f = Arrayd({len*2}, 0);
        scale = scales[i];
        int m = scale;
        double d = (double)nPoints / (2*scale);
        for (int j = 0; j <= 2*scale; j++) f[j] = psi[int(j * d)];
        f.selfReverse().selfMinus(f.mean());
        JN_INFOA(f);
        auto ff = Arrayd::convolve(ms2, f).selfTimes(1.0 / std::sqrt(scale));
        auto && sub = wCoefs.sub({{i,i+1},{0,len*2}});
        for (int j = 0; j < len*2-1-m; j++) sub(j) = ff[j+m];
        for (int j = 0; j <= m; j++) sub(len*2-1-m+j) = ff[j];
        //ms.identify("ms");
        //ff.identify("ff");
        //wCoefs.sub({{i, i+1}, {0, len}}) = Arrayd::convolve(ms, f).selfTimes(1.0 / std::sqrt(scale));
    }
    //return std::move(wCoefs);
    return wCoefs.sub({{0,scales.size()},{0,len}});
}

inline Arrayi localMaximum (const Arrayd &x, int winSize) {
    int len = x.size();
    int shift = winSize/2;

    std::vector<int> maxInds;
    double max = x[0], max2 = x[0];
    int ind = 0, ind2 = 0;
    JN_INFOV(winSize);
    for (int i = 0; i < len; i++) {
        double xi = std::abs(x[i]);
        int n = i % winSize;
        if (n == 0)       { max  = xi; ind  = i; }
        if (n == shift+1) { max2 = xi; ind2 = i; }
        if (xi > max+1e-5)   { max  = xi; ind  = i; }
        if (xi > max2+1e-5)  { max2 = xi; ind2 = i; }
        if (n == winSize - 1 || i == len-1) { if (((ind % winSize) != 0) && ((ind % winSize) != (winSize - 1)))   { maxInds.push_back(ind); } }
        if (n == shift       || i == len-1) { if (((ind2 % winSize) != (shift+1)) && ((ind2 % winSize) != shift))       { maxInds.push_back(ind2); } }
        //if (n == winSize - 1 || i == len-1) { { maxInds.push_back(ind); } }
        //if (n == shift       || i == len-1) { { maxInds.push_back(ind2); } }
    }
    std::sort(maxInds.begin(), maxInds.end());
    JN_INFOA(x);
    //for (auto && n : maxInds) MPI_LOG << n << ' '; MPI_LOG << std::endl;
    Arrayi v({len}, 0);
    for (int i = 0; i < maxInds.size(); i++) v[maxInds[i]] = 1;
//    for (int i = 0; i < maxInds.size()-1; i++) {
//        if (maxInds[i+1]-maxInds[i]>shift) v[maxInds[i]] = 1;
//    }
//    if (maxInds.size() >= 1) v[maxInds.back()] = 1;
    return std::move(v);
}

inline Arrayi getLocalMaximumCWT(const Arrayd &wCoefs, int ampTh=0) {
    int len = wCoefs.shape(1);
    Arrayi localMax(wCoefs.shape(), 0);
    for (int i = 0; i < wCoefs.shape(0); i++) {
        localMax.sub({{i,i+1},{0,len}}) = localMaximum(wCoefs.sub({{i,i+1},{0,len}}), 3);
    }
    JN_INFOA(localMax);
    return std::move(localMax);
}

inline std::vector<std::vector<std::array<int, 2>>> getRidge(const Arrayi &localMax, const Arrayd &widths) {
    int gap = int(std::ceil(widths[0]));
//    MPI_LOG << "gap:" << gap << std::endl;
    auto && maxDistances = widths / 4.0;
//    maxDistances.identify("maxDistances");

    std::vector<std::vector<std::array<int, 2>>> v1, v2, v3, v4;
    std::array<int, 2> temp;
    int nRows = localMax.shape(0), nCols = localMax.shape(1);
    for (int i = nRows-1; i >= 0; i--) {
        if (!v2.empty()) v3.resize(v2.size());
        for (int j = 0; j < nCols; j++) {
            if (localMax(i, j) == 1) {
                if (v2.empty()) {
			temp[0] = i;
			temp[1] = j;
			v4.push_back(std::vector<std::array<int, 2> >{temp});
		}
                else {
                    std::vector<int> dists(v2.size());
                    for (int k = 0; k < v2.size(); k++) dists[k] = std::abs(j-v2[k].back()[1]);
                    auto it = std::min_element(dists.begin(), dists.end());
                    if (*it <= maxDistances[i]) {
                        int d = std::distance(dists.begin(), it);
			temp[0] = i;
			temp[1] = j;
                        v3[d].push_back(temp);
                    }
                    else {
			    temp[0] = i;
			    temp[1] = j;
			    v4.push_back(std::vector<std::array<int, 2>>{temp});
                    }
                }
            }
        }
        for (int k = 0; k < v3.size(); k++) {
            if (v3[k].empty()) {
                if (std::abs(v2[k].back()[0]-i) <= gap) v4.push_back(std::move(v2[k]));
                else v1.push_back(std::move(v2[k]));
            }
            else {
                auto it = std::min_element(v3[k].begin(), v3[k].end(), [&v2, &k](const std::array<int, 2> &a, const std::array<int, 2> &b){return std::abs(a[1]-v2[k].back()[1]) < std::abs(b[1]-v2[k].back()[1]); });
                int d = std::distance(v3[k].begin(), it);
                for (int l = 0; l < v3[k].size(); l++) if (l != d) v4.push_back({v3[k][l]});
                v2[k].push_back(v3[k][d]);
                v4.push_back(std::move(v2[k]));
            }
        }
        v2 = std::move(v4);
        v4.clear();
        v3.clear();
    }
    for (auto && v : v2) v1.push_back(std::move(v));
    for (auto && ridge : v1) {
        std::sort(ridge.begin(), ridge.end(), [](const std::array<int, 2> &a1, const std::array<int, 2> &a2){return a1[0]<a2[0];});
    }
    return std::move(v1);
}

inline double quantile(const Arrayd &x, double prob) {
//    auto && x1 = Arrayd::abs(x);
    Arrayd x1 = x;
    std::sort(x1.begin(), x1.end());
    int n = x1.size();
    double h = prob * (n-1);
    int h1 = int(h);
    return x1[h1] + (x1[h1+1]-x1[h1]) * (h - h1);
}

inline double ridgeSignal(const Ridge &ridge, const Arrayd &wCoefs) {
    auto & p = ridge.front();
    return wCoefs(p[0], p[1]);
//    return std::abs(wCoefs(p[0], p[1]));
//    int l = ridge.size();
//    std::vector<double> w(l);
//    for (int i = 0; i < l; i++) w[i] = std::abs(wCoefs(ridge[i][0], ridge[i][1]));
//    return *std::max_element(w.begin(), w.end());
}

inline double ridgeNoise(const Ridge &ridge, const Arrayd &wCoefs/*, const Arrayi &scales*/) {
    auto & p = ridge.front();
    int winSize = std::max(int(std::ceil(wCoefs.shape(1)/20.0)), 2);
//    MPI_LOG << "winSize: " << winSize << std::endl;
    int hfWinSize = winSize / 2;
    int odd = winSize % 2;
//    MPI_LOG << hfWinSize << std::endl;
    int start = std::max(0, p[1]-hfWinSize);
    int end = std::min(p[1]+hfWinSize+odd,wCoefs.shape(1));
//    MPI_LOG << "start: " << start << ", end: " << end << std::endl;
//    rangei.printFull();
//    MPI_LOG << "window: " << start << '-' << end << std::endl;
//    return quantile(wCoefs.sub({{p[0], p[0]+1}, {start, end}}), 0.1);
    return quantile(wCoefs.sub({{0, 1}, {start, end}}), 0.1);
}

inline std::vector<Ridge> identifyMajorPeaks(const std::vector<Ridge> &ridges, const Arrayd &wCoefs, int minSNR) {
    std::vector<Ridge> v;
    for (auto && ridge : ridges) {
        double signal = ridgeSignal(ridge, wCoefs);
        double noise = ridgeNoise(ridge, wCoefs);
        int length = ridge.size();
        double snr = (noise == 0 ? 999 : std::abs(signal / noise));
        //MPI_LOG << signal << ' ' << noise << ' ' << signal / noise << ' ' << length << std::endl;
        if (snr >= minSNR && length >= wCoefs.shape(0)/4.0) v.push_back(ridge);
    }
    return std::move(v);
}

inline Arrayi find_peaks_cwt(const Arrayd &ms, const Arrayd &widths, int SNR = 2) {
    // Perform Continuous Wavelet Transform
    auto && wCoefs = cwt(ms, widths);
    JN_INFOA(wCoefs);
//    wCoefs.printFull();

    // Attach the raw data as the zero level of decomposition
//    wCoefs = cbind(as.vector(ms), wCoefs);
//    colnames(wCoefs) = c(0, scales);

    //-----------------------------------------
    // Identify the local maximum by using a slide window
    // The size of slide window changes over different levels, with the coarse level have bigger window size
    auto && localMax = boolrelextrema(wCoefs, std::greater<double>{}, 1);
//    auto && localMax = getLocalMaximumCWT(wCoefs);
//    localMax.printFull();
//    localMax.identify("localMax");

    //-----------------------------------------
    // Indentify the ridges from coarse level to more detailed levels
    auto && ridgeList = getRidge(localMax, widths);
//    MPI_LOG << "Ridges" << std::endl;
//    for (auto && ridge : ridgeList) {
//        MPI_LOG << "ridge: ";
//        for (auto && p : ridge) MPI_LOG << p[0] << '-' << p[1] << ' ';
//        double signal = ridgeSignal(ridge, wCoefs);
//        double noise = ridgeNoise(ridge, wCoefs);
//        MPI_LOG << signal << ' ' << noise << ' ' << signal / noise << std::endl;
//    }

    //-----------------------------------------
    // Indentify the major peaks and their nearby peaks 
    auto && majorPeaks = identifyMajorPeaks(ridgeList, wCoefs, SNR);
//    MPI_LOG << "majorRidges" << std::endl;
//    for (auto && ridge : majorPeaks) {
//        MPI_LOG << "ridge: ";
//        for (auto && p : ridge) MPI_LOG << p[0] << '-' << p[1] << ' ';
//        MPI_LOG << std::endl;
//    }


    Arrayi v({int(majorPeaks.size())});
    for (int i = 0; i < v.size(); i++) {
        v[i] = majorPeaks[i].front()[1];
    }
    std::sort(v.begin(), v.end());
    return std::move(v);
}

}}
