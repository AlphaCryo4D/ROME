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
#include "res_utils.h"
#include "res_signal.h"
#include "res_poly.h"

namespace resumap
{
    using spectrum_type = Arrayd;

    /**
     * Creates a radius matrix.
     */
    template<typename _Mat = Arrayd>
    _Mat createRmatrix(int n) {
        using value_type = typename _Mat::value_type;
        _Mat r({n, n, n});
        double m = double(n);
        for (int i = 0; i < n; i++) for (int j = 0; j < n; j++) for (int k = 0; k < n; k++) {
            double x = -m/2+i*m/(m-1);
            double y = -m/2+j*m/(m-1);
            double z = -m/2+k*m/(m-1);
            r(i, j, k) = value_type(std::sqrt(x*x+y*y+z*z));
        }
        return std::move(r);
    }

    /**
     *  Calculates the spherically averaged profile of a volume.
     */
    template<typename _Array>
    spectrum_type sphericalAverage(const _Array &image, double binsize = 1.0) {
        int n = image.shape(0);

        auto && r = createRmatrix<Arrayd>(n);

        // int nbins = int(std::ceil(n / 2.0 / binsize));
        int nbins = int(round(((n / 2.0) - 1) / binsize));
        double maxbin = nbins * binsize;
        double new_binsize = maxbin / (nbins - 1);
        //        std::vector<double> bins(nbins);
        //        for (int i = 0; i < bins.size(); i++) bins[i] = i * new_binsize;
        //        for (auto n : bins) MPI_LOG << n << ' '; MPI_LOG << std::endl;

        //MPI_LOG << "nbins:" << nbins << " maxbin:" << maxbin << " new_binsize:" << new_binsize << std::endl;
        spectrum_type sums({nbins}, 0.0);
//        Arrayi which(image.shape(), 0);
        JN_INFOA(sums);
        Arrayi ns({nbins}, 0);
        for (int i = 0; i < n * n * n; i++) {
            int m = int(r(i) / new_binsize);
            if (m < nbins) {
//                which[i] = m;
                sums[m] += image(i);
                ns[m]++;
            }
            else {
//                which[i] = nbins-1;
                sums[nbins-1] += image(i);
                ns[nbins-1]++;
            }
        }
//        Arrayi indices = Arrayi::range(0,n*n*n);
//        std::sort(indices.begin(), indices.end(), [&which](int i, int j){return which[i]<which[j];});
//        Arrayd cumsum({n*n*n},0);
//        cumsum[0] = image[indices[0]];
//        for (int i = 1; i < n*n*n; i++) cumsum[i] = cumsum[i-1]+image[indices[i]];
//        image.identify("image");
//        image.print();
//        MPI_LOG << "cumsum:" << std::endl;
//        cumsum[0] = image[0];
//        for (int i = 1; i < n*n*n; i++) cumsum[i] = cumsum[i-1]+image[i];
//        MPI_LOG << std::setprecision(12);
//:        cumsum.print();
        JN_INFOA(image);
        JN_INFOA(sums);
//        ns.printFull();
//        MPI_LOG << "sum_by_group: ";
        for (int i = 0; i < nbins; i++) {
//            MPI_LOG << sums[i] << ' ';
            if (ns[i] != 0) {
                sums[i] /= double(ns[i]);
            }
            //MPI_LOG << sums[i] << ' ';
        }
//        MPI_LOG << std::endl;
        //MPI_LOG << sums << std::endl;

        JN_INFOA(sums);

        return std::move(sums);
    }

    /**
     * Calculates the radially averaged power spectrum of a volume.
     */
    spectrum_type calculatePowerSpectrum(const Arraycd &m);

    /**
     * Attempts to determine whether there is a low-pass drop in the spectrum.
     */
    double isPowerSpectrumLPF(const spectrum_type &m);

    /**
     * Creates a pre-whitening filter in 3D.
     *
     * Fits a polynomial to the spectrum beyond the "elbowAngstrom" frequency.
     * Returns a whitening filter that can be adjusted using the "rampWeight.
     */
    struct PreWhiteningFilter {
        Arrayd peval;
        Arrayd pcoef;
        Arrayd pWfilter;

        PreWhiteningFilter(const Arrayd &spectrum, double elbowAngstrom, double rampWeight, double vxSize, double n) {
            double epsilon = 1e-10;

            auto && R = createRmatrix(n);

            // Create the x and y variables for the polynomial regression
            Arrayd xpoly = Arrayd::range(1, spectrum.size() + 1);
            Arrayd ypoly = Arrayd::sqrt(spectrum).selfLog();
            ypoly.reshape({ypoly.size(), 1});

            JN_INFOA(spectrum);
            JN_INFOA(xpoly);
            JN_INFOA(ypoly);

            // Create the index of frequencies (depends on vxSize)
            auto Fs     = 1.0 / vxSize;
            Arrayd Findex = Arrayd::linspace(epsilon, 1, xpoly.size()).selfTimes(Fs / 2.0).selfDividedBy(1.0);

            JN_INFOV(Fs);
            JN_INFOA(Findex);

            // Find the points of interest
            int indexElbow = Arrayd::minus(Findex, elbowAngstrom).selfSquare().argmin();
            int indexStart = Arrayd::minus(Findex, (1.05*elbowAngstrom)).selfSquare().argmin();
            int indexNyquist = xpoly[int(xpoly.size()) - 1];

            JN_INFOV(indexElbow);
            JN_INFOV(indexStart);
            JN_INFOV(indexNyquist);
//            MPI_LOG << "indexElbow: " << indexElbow << std::endl;
//            MPI_LOG << "indexStart: " << indexStart << std::endl;
//            MPI_LOG << "indexNyquist: " << indexNyquist << std::endl;

            // Create the weighting function to do a weighted fit
            Arrayd wpoly = Arrayd::bitAnd(Arrayb::gt(xpoly, indexElbow), Arrayb::lt(xpoly, indexNyquist));
            wpoly.selfPlus(Arrayd::bitAnd(Arrayb::gt(xpoly, indexStart), Arrayb::le(xpoly, indexElbow)).selfTimes(0.5));

            JN_INFOA(wpoly);

            // Do the polynomial fit
            pcoef = polyfit(xpoly, ypoly, 2, wpoly).c;
            JN_INFOA(pcoef);

            peval = polyval<double>(xpoly, pcoef);
            JN_INFOA(peval);

            // Don't change any frequencies outside of indexStart to indexNyquist
            for (int i = 0; i < R.size(); i++) if (R[i] < indexStart) R[i] = indexStart;
            for (int i = 0; i < R.size(); i++) if (R[i] > indexNyquist) R[i] = indexNyquist;
            JN_INFOA(R);

            // Create the pre-whitening filter
            pWfilter = polyval<double>(R, Arrayd::times(pcoef, -1.0*rampWeight)).selfExp();
            JN_INFOA(pWfilter);
        }
    };

    /**
     *  Creates a pre-whitening filter in 3D.
     *
     *  Expects a fitted polynomial defined by its "pcoef". Returns a
     *  whitening filter that can be adjusted using the "rampWeight."
     */
    struct PreWhiteningFilterFinal {
        Arrayd pWfilter;

        PreWhiteningFilterFinal(const Arrayd &spectrum, const Arrayd &pcoef, double elbowAngstrom, double rampWeight, double vxSize, int n, int cubeSize) {
            JN_INFOV(n);
            JN_INFOV(cubeSize);
            JN_INFOV(vxSize);

            double epsilon = 1e-10;

            auto && R = createRmatrix(n);

            // Create the x and y variables for the polynomial regression
            auto && xpoly = Arrayd::range(1, spectrum.shape(0) + 1);

            JN_INFOA(xpoly);

            // Create the index of frequencies (depends on vxSize)
            auto Fs     = 1.0 / vxSize;
            auto && Findex = Arrayd::linspace(epsilon, 1, xpoly.shape(0));
            Findex.selfTimes(Fs / 2.0).selfDividedBy(1.0);

            JN_INFOV(vxSize);
            JN_INFOA(Findex);

            // Find the points of interest
            int indexStart    = Arrayd::minus(Findex, 1.05*elbowAngstrom).selfSquare().argmin();
            int indexNyquist  = xpoly[xpoly.shape(0)-1];

            JN_INFOV(indexStart);
            JN_INFOV(indexNyquist);

            // Don't change any frequencies outside of indexStart to indexNyquist
            for (int i = 0; i < R.size(); i++) if (R[i] < indexStart) R[i] = indexStart;
            for (int i = 0; i < R.size(); i++) if (R[i] > indexNyquist) R[i] = indexNyquist;

            // Rescale R such that the polynomial from the cube fit makes sense
            JN_INFOA(R);
            R.selfDivide(double(n) / (cubeSize-1));
            JN_INFOA(R);

            JN_INFOV(rampWeight);
            JN_INFOA(pcoef);

            // Create the pre-whitening filter
            pWfilter = polyval<double>(R, Arrayd::times(pcoef, -1.0*rampWeight)).selfExp();

            JN_INFOA(pWfilter);
        }
    };

    /**
     * Result of prewhiten.
     */
    struct PreWhitenResult {
        Arrayd dataPW;
        Arrayd dataBGPW;
        Arrayd dataPWSpect;
        Arrayd dataPWBGSpect;
        Arrayd peval;
        Arrayd pcoef;
    };

    /**
     * Pre-whitenening using noise estimates from a soft mask of the background. 
     *
     * Returns a the pre-whitened volume and various spectra.
     */
    inline void preWhitenVolumeSoftBG(
            PreWhitenResult &result,
            int n,
            double vxSize,
            double elbowAngstrom,
            double rampWeight,
            Arraycd dataF,
            const Arraycd &dataBGF,
            const Arrayd &dataBGSpect,
            const Arrayb &softBGmask)
    {
        double epsilon = 1e-10;

        PreWhiteningFilter pWfilter(dataBGSpect, elbowAngstrom, rampWeight, vxSize, n);
        result.peval = pWfilter.peval;
        result.pcoef = pWfilter.pcoef;

        // Apply the pre-whitening filter
        dataF = Arraycd::times(pWfilter.pWfilter, dataF);

        auto && dataPWFabs  = Arrayd::abs(dataF);
        dataPWFabs.selfMinus(dataPWFabs.min());
        dataPWFabs.selfDivide(dataPWFabs.max());

        result.dataPWSpect = sphericalAverage(Arrayd::square(dataPWFabs)).selfPlus(epsilon);

        JN_INFOA(result.dataPWSpect);

        result.dataPW = Arrayd::real(Arraycd::fftshift3(dataF).selfIFFT3());

        Arraycd dataPWBG = std::move(Arraycd::times(result.dataPW, softBGmask).selfFFT3().selfFFTshift3());

        JN_INFOA(Arrayd::real(dataPWBG));

        auto && dataPWBGFabs = Arrayd::abs(dataPWBG);

        JN_INFOV(dataPWBGFabs);

        dataPWBGFabs.selfMinus(dataPWBGFabs.min<double>());
        dataPWBGFabs.selfDivide(dataPWBGFabs.max<double>());

        result.dataPWBGSpect = sphericalAverage(Arrayd::square(dataPWBGFabs)).selfPlus(epsilon);
    }


    /**
     * Pre-whitenening using noise estimates from a cube taken from the difference map. 
     *
     * Returns a the pre-whitened volume and various spectra.
     */
    inline void preWhitenCube(
            PreWhitenResult &result,
            int n,
            double vxSize,
            double elbowAngstrom,
            double rampWeight,
            Arraycd dataF,
            Arraycd dataBGF,
            const Arrayd &dataBGSpect)
    {
        double epsilon = 1e-10;
        JN_INFOA(dataBGSpect);
        JN_INFOV(elbowAngstrom);
        JN_INFOV(rampWeight);
        JN_INFOV(vxSize);
        JN_INFOV(n);
        PreWhiteningFilter pWfilter(dataBGSpect, elbowAngstrom, rampWeight, vxSize, n);
        result.peval = pWfilter.peval;
        result.pcoef = pWfilter.pcoef;

        // Apply the pre-whitening filter to the inside cube
//        dataF.selfTimes(pWfilter.pWfilter);
        JN_INFOA(Arrayd::real(dataF));
        dataF = Arraycd::times(pWfilter.pWfilter, dataF);
        JN_INFOA(Arrayd::real(dataF));

        auto && dataPWFabs  = Arrayd::abs(dataF);
        dataPWFabs.selfMinus(dataPWFabs.min<double>());
        dataPWFabs.selfDivide(dataPWFabs.max<double>());
        result.dataPWSpect = sphericalAverage(Arrayd::square(dataPWFabs)).selfPlus(epsilon);

        JN_INFOA(dataPWFabs);
        JN_INFOA(result.dataPWSpect);

//        Arraycd tmp1(dataF.shape()), tmp2(dataF.shape());
//        fftshift3(dataF.data, tmp1.data, dataF.shape(0));
//        ifft3(tmp1.data, tmp2.data, tmp1.shape(0));
//        result.dataPW = Arrayd::real(tmp2);
        result.dataPW = Arrayd::real(Arraycd::fftshift3(dataF).selfIFFT3());
        JN_INFOA(Arrayd::real(dataF));
        JN_INFOA(result.dataPW);

        // Apply the pre-whitening filter to the outside cube
        dataBGF = Arraycd::times(dataBGF, pWfilter.pWfilter);

        auto && dataPWBGFabs  = Arrayd::abs(dataBGF);
        dataPWBGFabs.selfMinus(dataPWBGFabs.min());
        dataPWBGFabs.selfDivide(dataPWBGFabs.max());
        result.dataPWBGSpect = Arrayd::plus(sphericalAverage(Arrayd::square(dataPWBGFabs)), epsilon);
        JN_INFOA(result.dataPWBGSpect);

        result.dataBGPW = Arrayd::real(dataBGF.selfFFTshift3().selfIFFT3());
    }

}


