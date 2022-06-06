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

#include "resmap_mpi.h"
#include "resmap_time.h"

#include "res_fft.h"
#include "res_filter.h"
#include "res_spectrum.h"
#include "res_calculator.h"
#include "res_morph.h"
#include "res_utils.h"


namespace resumap {

    // parameters
    extern double vxSize;
    extern Arrayb dataMask;
    extern double variance;

    // temporary variables
    extern Arrayd data;
    extern Arraycd dataF;
    extern Arrayd dataDiff;
    extern Arrayd dataPowerSpectrum;
    extern Arrayb mask;
    extern int n;
    extern int subVolLPF;
    extern bool splitVolume;
    extern double LPFfactor;
    extern bool estimateVarianceFromBackground;
    extern Arrayb Rinside;
    extern int cubeSize;

    void preWhiten() {
        std::vector<double> times;

        times.push_back(dtime());
        MPI_LOG << "[*] Start to pre-whiten..." << std::endl;

        double oldElbowAngstrom = 0.0;
        double oldRampWeight = 0.0;

        double newElbowAngstrom = std::max(10.0, 2.1*vxSize);
        double newRampWeight = 1.0;

        int widthBox;
        Arrayi dataMaskDistance;
        Arrayi dataOutsideDistance;
        Arrayd cubeOutside;
        Arrayd dataSpect;
        Arrayd dataBGSpect;
        Arrayb dilatedMask;
        Arrayd softBGmask;
        Arraycd dataBGF;
        Arraycd dataDiffF;
        Arrayd dataPowerSpectrumDiff;
        Arrayd dataBG;
        PreWhitenResult preWhitenResult;
        Arrayd dataPW;
        Arrayd dataBGPW;

        times.push_back(dtime());
        MPI_LOG << "[*] Initialize: " << times[times.size()-1]-times[times.size()-2] << " seconds." << std::endl;

        if (variance == 0.0) {
            MPI_LOG << "[*] Variance: 0.0" << std::endl;

            estimateVarianceFromBackground = true;
            if (n > subVolLPF) {

                dataMaskDistance = distance_transform_cdt3(dataMask);

                JN_INFOA(dataMask);
                JN_INFOA(dataMaskDistance);

                times.push_back(dtime());
                MPI_LOG << "[*] Distance transform of mask: " << times[times.size()-1]-times[times.size()-2] << " seconds." << std::endl;

                if (!splitVolume) {
                    auto && dataOutside = Arrayb::logicalAnd(Arrayb::logicalNot(dataMask), Rinside);
                    dataOutsideDistance = distance_transform_cdt3(dataOutside);

                    JN_INFOA(dataOutside);
                    JN_INFOA(dataOutsideDistance);
                }

                times.push_back(dtime());
                MPI_LOG << "[*] Distance transform of ouside data: " << times[times.size()-1]-times[times.size()-2] << " seconds." << std::endl;

                if (splitVolume) {
                    widthBox = dataMaskDistance.max();
                }
                else {
                    widthBox = std::min(dataMaskDistance.max(), dataOutsideDistance.max());
                }
                int halfWidthBox = int(std::floor(widthBox/2));
                cubeSize     = 2*halfWidthBox;

                JN_INFOV(n);
                JN_INFOV(widthBox);
                JN_INFOV(halfWidthBox);
                JN_INFOV(cubeSize);


                auto && insideBox    = unravel_index(dataMaskDistance.argmax(),{n,n,n});
                JN_INFOA(dataMaskDistance);
                JN_INFOV(dataMaskDistance.argmax());
                JN_INFOV(n);
                JN_INFOA(insideBox);

                for (auto && i : insideBox) {
                    i = std::max(i, cubeSize);
                    i = std::min(i, n-cubeSize);
                }

                JN_INFOA(insideBox);

                Arrayd cubeInside   = data.sub({
                        {insideBox[0]-halfWidthBox,insideBox[0]+halfWidthBox},
                        {insideBox[1]-halfWidthBox,insideBox[1]+halfWidthBox},
                        {insideBox[2]-halfWidthBox,insideBox[2]+halfWidthBox}});

                if (splitVolume) {
                    cubeOutside  = dataDiff.sub({
                            {insideBox[0]-halfWidthBox, insideBox[0]+halfWidthBox},
                            {insideBox[1]-halfWidthBox, insideBox[1]+halfWidthBox},
                            {insideBox[2]-halfWidthBox, insideBox[2]+halfWidthBox}});
                }
                else {
                    auto && outsideBox = unravel_index(dataOutsideDistance.argmax(),{n,n,n});
                    JN_INFOA(outsideBox);

                    cubeOutside  = data.sub({
                            {outsideBox[0]-halfWidthBox, outsideBox[0]+halfWidthBox},
                            {outsideBox[1]-halfWidthBox, outsideBox[1]+halfWidthBox},
                            {outsideBox[2]-halfWidthBox, outsideBox[2]+halfWidthBox}});
                }

                auto && hammingWindow1D = hamming(cubeSize);
                auto && hammingWindow2D = array_outer_product(hammingWindow1D,hammingWindow1D);
                auto && hammingWindow3D = array_outer_product2(hammingWindow2D,hammingWindow2D);

                cubeInside.selfTimes(hammingWindow3D);
                cubeOutside.selfTimes(hammingWindow3D);

                JN_INFOA(cubeInside);
                JN_INFOA(cubeOutside);

                times.push_back(dtime());
                MPI_LOG << "[*] Set hamming window: " << times[times.size()-1]-times[times.size()-2] << " seconds." << std::endl;

                dataF = Arraycd(cubeInside).selfFFT3().selfFFTshift3();
                dataSpect = calculatePowerSpectrum(dataF);

                JN_INFOA2(Arrayd::real(dataF), "dataF");
                JN_INFOA2(Arrayd::real(dataSpect), "dataSpect");

                dataBGF = Arraycd(cubeOutside).selfFFT3().selfFFTshift3();
                dataBGSpect = calculatePowerSpectrum(dataBGF);

                JN_INFOA2(Arrayd::real(dataBGF), "dataBGF");
                JN_INFOA2(Arrayd::real(dataBGSpect), "dataBGSpect");

                times.push_back(dtime());
                MPI_LOG << "[*] Calculate power spectrum of the inside and outside cubes: "
                    << times[times.size()-1]-times[times.size()-2] << " seconds." << std::endl;

//                while (newElbowAngstrom != oldElbowAngstrom || oldRampWeight != newRampWeight) {
                    preWhitenCube(preWhitenResult, cubeSize, vxSize, newElbowAngstrom, newRampWeight, dataF, dataBGF, dataBGSpect);                
                    auto && cubeInsidePW = preWhitenResult.dataPW;
                    oldElbowAngstrom = newElbowAngstrom;
                    oldRampWeight    = newRampWeight;
//                }

                times.push_back(dtime());
                MPI_LOG << "[*] Set pre-whiten cube: "
                    << times[times.size()-1]-times[times.size()-2] << " seconds." << std::endl;

                // Apply the pre-whitening filter on the full-sized map
                dataF = std::move(Arraycd(data).selfFFT3().selfFFTshift3());
                dataPowerSpectrum = calculatePowerSpectrum(dataF);
                if (splitVolume) {
                    dataDiffF = std::move(Arraycd(dataDiff).selfFFT3().selfFFTshift3());
                    dataPowerSpectrumDiff = calculatePowerSpectrum(dataDiffF);
                }

                times.push_back(dtime());
                MPI_LOG << "[*] Calculate power spectrum of total data: "
                    << times[times.size()-1]-times[times.size()-2] << " seconds." << std::endl;

                PreWhiteningFilterFinal pwFilterFinal(dataPowerSpectrum, preWhitenResult.pcoef, newElbowAngstrom, newRampWeight, vxSize, n, cubeSize);
                //Arrayd::real(dataF).identify("dataF");
                //pwFilterFinal.pWfilter.identify("pWfilter");
                //Arrayd::real(Arraycd::times(dataF, pwFilterFinal.pWfilter)).identify("dataF2");

                times.push_back(dtime());
                MPI_LOG << "[*] Set final pre-whiten filter: "
                    << times[times.size()-1]-times[times.size()-2] << " seconds." << std::endl;

                dataPW = Arrayd::real(Arraycd::times(dataF, pwFilterFinal.pWfilter).selfFFTshift3().selfIFFT3());
                JN_INFOA(dataPW);
                data = dataPW;
                if (splitVolume) {
                    dataDiff = Arrayd::real(Arraycd::times(pwFilterFinal.pWfilter, dataDiffF).selfFFTshift3().selfIFFT3());
                }

                times.push_back(dtime());
                MPI_LOG << "[*] Apply the pre-whitening filter on the full-sized map: "
                    << times[times.size()-1]-times[times.size()-2] << " seconds." << std::endl;

            }
            else {
                if (!splitVolume) {
                    // Dilate the mask a bit so that we don't seep into the particle when we blur it later
                    Arrayd boxElement({5, 5, 5}, 1);
                    dilatedMask = std::move(binary_dilation(dataMask, boxElement, 3).selfLogicalAnd(Rinside));

                    JN_INFOA(Rinside);
                    JN_INFOA(dataMask);
                    JN_INFOA(boxElement);
                    JN_INFOA(dilatedMask);

                    // Blur the mask
                    softBGmask = gaussian_filter(Arrayb::logicalNot(dilatedMask), double(n)*0.02);

                    JN_INFOV(n);
                    JN_INFOA(softBGmask);

                    // Get the background
                    dataBG = Arrayd::times(data, softBGmask);

                    JN_INFOA(data);
                    JN_INFOA(dataBG);

                }
                else {
                    dataBG = dataDiff;
                }

                // Calculate spectrum of input volume only if downsampled, otherwise use previous computation
                JN_INFOA(Arrayd::real(dataF));
                if (LPFfactor != 0.0 || n > subVolLPF) {
                    dataF = std::move(Arraycd::fft3(Arraycd(data)).selfFFTshift3());
                    dataPowerSpectrum = calculatePowerSpectrum(dataF);
                }
                JN_INFOA(Arrayd::real(dataF));

                // Calculate spectrum of background volume
                dataBGF = std::move(Arraycd::fft3(Arraycd(dataBG)).selfFFTshift3());
                dataBGSpect = calculatePowerSpectrum(dataBGF);

//                while (newElbowAngstrom != oldElbowAngstrom or oldRampWeight != newRampWeight) {

                    if (!splitVolume) {
                        preWhitenVolumeSoftBG(preWhitenResult, n, vxSize, newElbowAngstrom, newRampWeight, dataF, dataBGF, dataBGSpect, softBGmask);
                    }
                    else {
                        preWhitenCube(preWhitenResult, n, vxSize, newElbowAngstrom, newRampWeight, dataF, dataBGF, dataBGSpect);
                    }

                    dataPW   = preWhitenResult.dataPW;
                    if (splitVolume) {
                        dataBGPW = preWhitenResult.dataBGPW;
                    }
                    
                    oldElbowAngstrom = newElbowAngstrom;
                    oldRampWeight    = newRampWeight;
//                }

                JN_INFOA(data);
                JN_INFOA(dataPW);
                data     = dataPW;
                if (splitVolume) {
                    dataDiff = dataBGPW;
                }
            }
        }
        else {
            estimateVarianceFromBackground = false;
//            crudeEstimate = np.var(data[np.logical_and(Rinside,np.logical_not(mask))]);
        }
    }

} // end namespace resumap

