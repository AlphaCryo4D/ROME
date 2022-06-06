/***************************************************************************
 *
 * Intel® Parallel Computing Center for Structural Biology
 * Principal Investigator : Youdong (Jack) Mao (Youdong_Mao@dfci.harvard.edu)
 * Dana-Farber Cancer Institute, Harvard Medical School and Peking University
 *
 * Authors: "Jian Wang(wj_hust08@hust.edu.cn) Yong Bei Ma(galowma@gmail.com)"
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
#include "res_calculator.h"

namespace resumap {

    // parameters
    std::string inputFileName;
    std::string inputFileName1;
    std::string inputFileName2;
    double vxSize;
    double pValue;
    double minRes;
    double maxRes;
    double stepRes;
    Arrayb dataMask;
    double variance;
//    bool  noiseDiagnostics;

    // temporary variables
    Mrcs::MrcVolume volume;
    Mrcs::MrcVolume volume1;
    Mrcs::MrcVolume volume2;
    Arrayd data;
    Arraycd dataF;
    Arrayd dataDiff;
    Arrayd dataPowerSpectrum;
    Arrayb mask;
    int n;
    int subVolLPF;
    double LPFfactor;
    bool splitVolume;
    bool estimateVarianceFromBackground = true;
    int oldSumOfMask;
    Arrayb Rinside;
    int cubeSize;
    double currentRes;

    void preWhiten();

    /**
     * Print parameters.
     */
    void print_pars() {
        MPI_LOG
            << "--------------------- Parameters ------------------------\n"
            << "inputFileName: "  << inputFileName    << "\n"
            << "inputFileName1: " << inputFileName1   << "\n"
            << "inputFileName2: " << inputFileName2   << "\n"
            << "vxSize: "         << vxSize           << "\n"
            << "pValue: "         << pValue           << "\n"
            << "minRes: "         << minRes           << "\n"
            << "maxRes: "         << maxRes           << "\n"
            << "stepRes: "        << stepRes          << "\n"
//            << "dataMask: "       << dataMask         << "\n"
            << "variance: "       << variance         << "\n";
//            << "noiseDiagnostics" << noiseDiagnostics << "\n";
//            << "---------------------------------------------------------\n"
//            << std::endl;
    }

    /**
     * Get the voxel size.
     */
    double getVxSize(const Mrcs::MrcVolume &volume) {
        int mx = volume.head->NX;
        float xlen = volume.head->X_length;
        if (mx > 0 && xlen > 0) return double(xlen) / mx;
        else return 1.0;
    }

    /**
     * Initialization for resolution calculation.
     */
    void init() {
        std::vector<double> times;

        times.push_back(dtime());
        MPI_LOG << "[*] Initialize..." << std::endl;

        splitVolume = false;
        if (!inputFileName.empty()) {
            volume.read(inputFileName);

            times.push_back(dtime());
            MPI_LOG << "[*] Read volume: " << times[times.size()-1]-times[times.size()-2] << " seconds." << std::endl;

            n = volume.size;
//            MPI_LOG << "Volume size: " << n << ", data: " << volume.data << std::endl;
            MapArrayf temp(volume.data, {n, n, n});
            data = Arrayd::minus(temp, temp.mean<double>());

            times.push_back(dtime());
            MPI_LOG << "[*] Set data: " << times[times.size()-1]-times[times.size()-2] << " seconds." << std::endl;

            // Set vxSize
            if (vxSize <= 0) vxSize = getVxSize(volume);

            JN_INFOA(data);
        }
        else if (!inputFileName1.empty() && ! inputFileName2.empty()) {
            splitVolume = true;
            volume1.read(inputFileName1);

            times.push_back(dtime());
            MPI_LOG << "[*] Read volume1: " << times[times.size()-1]-times[times.size()-2] << " seconds." << std::endl;

            volume2.read(inputFileName2);

            times.push_back(dtime());
            MPI_LOG << "[*] Read volume2: " << times[times.size()-1]-times[times.size()-2] << " seconds." << std::endl;

            n = volume1.size;
            data     = Arrayd::plus(MapArrayf(volume1.data, {n, n, n}), MapArrayf(volume2.data, {n, n, n}));
            data.selfTimes(0.5);

            times.push_back(dtime());
            MPI_LOG << "[*] Set data: " << times[times.size()-1]-times[times.size()-2] << " seconds." << std::endl;

            dataDiff = Arrayd::minus(MapArrayf(volume1.data, {n, n, n}), MapArrayf(volume2.data, {n, n, n}));
            dataDiff.selfTimes(0.5);

            times.push_back(dtime());
            MPI_LOG << "[*] Set dataDiff: " << times[times.size()-1]-times[times.size()-2] << " seconds." << std::endl;

            // Set vxSize
            if (vxSize <= 0) vxSize = getVxSize(volume1);

            JN_INFOA(data);
            JN_INFOA(dataDiff);
        }
    }

    /**
     * Check low-passing filtering.
     */
    void checkLPF() {
        std::vector<double> times;

        times.push_back(dtime());
        MPI_LOG << "[*] Check low pass filtration..." << std::endl;

        subVolLPF = 160;
        if (n > subVolLPF) {
            MPI_LOG << "[*] The volume appears to be quite large, the LPF test will" << std::endl;
            MPI_LOG << "    be carried out on a cube of size 160 from the center of the volume." << std::endl;

            int mid  = int(n/2.0);
            int midR = subVolLPF/2;

            Arraycd middleCube = data.sub({{mid-midR, mid+midR}, {mid-midR, mid+midR}, {mid-midR, mid+midR}});

            auto && hammingWindow1D = hamming(subVolLPF);
            auto && hammingWindow2D = array_outer_product(hammingWindow1D, hammingWindow1D);
            auto && hammingWindow3D = array_outer_product2(hammingWindow2D, hammingWindow2D);

            times.push_back(dtime());
            MPI_LOG << "[*] Create hamming window: " << times[times.size()-1]-times[times.size()-2] << " seconds." << std::endl;

            JN_INFOA(data);
            dataF = middleCube.selfTimes(hammingWindow3D).selfFFT3().selfFFTshift3();
            JN_INFOA(data);
            JN_INFOA2(Arrayd::real(dataF), "dataF");

            dataPowerSpectrum = calculatePowerSpectrum(dataF);
            LPFfactor = isPowerSpectrumLPF(dataPowerSpectrum);

            times.push_back(dtime());
            MPI_LOG << "[*] Calculate power spectrum: " << times[times.size()-1]-times[times.size()-2] << " seconds." << std::endl;

            JN_INFOA(dataPowerSpectrum);

        }
        else {
            dataF = std::move(Arraycd::fft3(Arraycd(data)).selfFFTshift3());

            dataPowerSpectrum = calculatePowerSpectrum(dataF);
            LPFfactor = isPowerSpectrumLPF(dataPowerSpectrum);

            JN_INFOA(data);
            JN_INFOA2(Arrayd::real(Arraycd::fft3(Arraycd(data))), "data.fft");
            JN_INFOA2(Arrayd::real(dataF), "dataF");

            times.push_back(dtime());
            MPI_LOG << "[*] Calculate power spectrum: " << times[times.size()-1]-times[times.size()-2] << " seconds." << std::endl;

        }

        if (LPFfactor > 0) {
            MPI_LOG << "The volume appears to be low-pass filtered." << std::endl;

            LPFfactor = round(LPFfactor/0.01)*0.01; // round to the nearest 0.01
            data = cppsci::ndimage::zoom(data, LPFfactor, "reflect");

            times.push_back(dtime());
            MPI_LOG << "[*] Zoom data: " << times[times.size()-1]-times[times.size()-2] << " seconds." << std::endl;

            if (splitVolume) {
                dataDiff = cppsci::ndimage::zoom(dataDiff, LPFfactor, "reflect");

                times.push_back(dtime());
                MPI_LOG << "[*] Zoom dataDiff: " << times[times.size()-1]-times[times.size()-2] << " seconds." << std::endl;

            }
            vxSize = double(vxSize)/LPFfactor;
        }
        else {
            MPI_LOG << "The volume does not appear to be low-pass filtered. Great!" << std::endl;
        }

        print_pars();
    }

    /**
     * Set the mask.
     */
    void setMask() {
        std::vector<double> times;

        times.push_back(dtime());
        MPI_LOG << "[*] Start to set mask..." << std::endl;

        if (minRes <= (2.2*vxSize)) minRes = round((2.2*vxSize)/0.1)*0.1;
        currentRes = minRes;

        if (maxRes == 0.0) maxRes = round((4.0*vxSize)/0.5)*0.5;

        n = data.shape(0);

        auto && R = createRmatrix<Arrayd>(n);
        Rinside = Arrayb::lt(R, n/2 - 1);

        times.push_back(dtime());
        MPI_LOG << "[*] Create Rmatrix: " << times[times.size()-1]-times[times.size()-2] << " seconds." << std::endl;

        JN_INFOA(data);
        if (true) {
            //MPI_LOG << "data:\n" << data << std::endl;
            JN_INFOA(data);
            auto && dataBlurred = gaussian_filter(data, n*0.02);
            JN_INFOA(dataBlurred);
            dataMask = Arrayb::gt(dataBlurred, dataBlurred.max<double>()*5e-2);
            JN_INFOA(dataMask);
//            MPI_LOG << "dataMaskSum: \n" << dataMask.sum<int>() << std::endl;
        }
        else {
            if (LPFfactor == 0.0) {
//                dataMask = dataMask.matrix;
//                dataMask = np.array(dataMask.matrix, dtype=bool);
            }
            else {
//                dataMask = zoom(dataMask.matrix, LPFfactor, "reflect");
//                dataMask = filters.gaussian_filter(dataMask, float(n)*0.02);
//                dataMask = dataMask > np.max(dataMask)*5e-2;
            }
        }

        times.push_back(dtime());
        MPI_LOG << "[*] Gaussian filter: " << times[times.size()-1]-times[times.size()-2] << " seconds." << std::endl;

        mask = Arrayb::bitAnd(dataMask, Arrayb::lt(R, n/2-9));
        JN_INFOA(mask);
        oldSumOfMask = mask.sum<int>();

        times.push_back(dtime());
        MPI_LOG << "[*] Set mask: " << times[times.size()-1]-times[times.size()-2] << " seconds." << std::endl;

    }

    template<typename F>
    void array_slide(Arrayd *pData, Arrayb *pMask, int winSize, F &&f) {
        int x = pData->shape(0);
        int y = pData->shape(1);
        int z = pData->shape(2);
        double *it_d = pData->data_wptr();
        bool *it_m = pMask->data_wptr();
        for (int i = 0; i+winSize-1 < x; i += winSize) {
            double *it_d2 = it_d;
            bool *it_m2 = it_m;
            for (int j = 0; j+winSize-1 < y; j += winSize) {
                double *it_d3 = it_d2;
                bool *it_m3 = it_m2;
                for (int k = 0; k+winSize-1 < z; k += winSize) {
                    f(it_d3, it_m3);
                    it_d3 += winSize;
                    it_m3 += winSize;
                }
                it_d2 += z*winSize;
                it_m2 += z*winSize;
            }
            it_d += y*z*winSize;
            it_m += y*z*winSize;
        }
    }

    template<typename Array_, typename ArrayNum_, typename Func_>
    void array_traverse(Array_ *arr, ArrayNum_ *it, int winSize, Func_ &&f) {
        int x = arr->shape(0);
        int y = arr->shape(1);
        int z = arr->shape(2);
        auto *it1 = it;
        for (int i = 0; i < winSize; i++) {
            auto it2 = it1+i*(y*z);
			//#pragma omp simd collapse(2)
            for (int j = 0; j < winSize; j++) {
                for (int k = 0; k < winSize; k++) {
					auto it3 = it2+j*z+k;
                    f(it3);
                }
            }
        }
    };

    double estimate_variance(Arrayd *pData, Arrayb *pMask, int winSize, Arrayd &LAMBDAd) {
        // Search for all the iters that satisfy specific conditions.
        std::vector<double *> data_iters;
        array_slide(pData, pMask, winSize, [&pMask, &data_iters, &winSize](double *it_d, bool *it_m) {
            bool all = true;
            array_traverse(pMask, it_m, winSize, [&all](bool *it) { all = (all && (*it)); });
            if (all) data_iters.push_back(it_d);
        });

        int sz = int(data_iters.size());
        Arrayd WRSSBG({ sz }, 0);
#ifdef USEMPI
        MPI_Comm mComm;
        int mRank = mpiRank();
        int mSize = mpiSize();
        int mBin = int(std::ceil(sz / (double)mSize));
        int mBeg = mRank*mBin;

#pragma omp parallel for schedule(dynamic)
        for (int i = mBeg; i < std::min(mBeg + mBin, sz); i++) {

#else

#pragma omp parallel for schedule(dynamic)
        for (int i = 0; i < sz; i++) {
            
#endif
            Arrayd temp({ winSize, winSize, winSize });
            double *it_temp = temp.data_wptr();
            array_traverse(pData, data_iters[i], winSize, [&it_temp](double *it_data) {
                *it_temp = *it_data;
                it_temp++;
            });
            temp.selfMinus(temp.mean());
            temp.reshape({ temp.size(),1 });
            WRSSBG[i] = linalg::vdot<double>(temp, linalg::dot(LAMBDAd, temp));
        }
        JN_INFOA(data);
        JN_INFOA(LAMBDAd);

#ifdef USEMPI
        double* local_temp = (double*)aMalloc(sz*sizeof(double), 64);
        for (int index = 0; index < sz; index++) local_temp[index] = WRSSBG[index];
        MPI::COMM_WORLD.Allreduce(local_temp, WRSSBG.data_wptr(), sz, MPI_DOUBLE, MPI::SUM);
        aFree(local_temp);
#endif

        double sum = 0; for (auto && n : WRSSBG) sum += n;
        variance = sum / WRSSBG.size() / LAMBDAd.trace();

        return variance;

    }

    /**
     * Compute the resolution.
     */
    void compute() {
        MPI_LOG << "[*] " << __FUNCTION__ << "..." << std::endl;
        // Initialize the ResMap result volume
        Arrayb maskBG;
        Arrayb maskParticle;
        Arrayd resTOTAL(data.shape(), 0);
        bool moreToProcess;

        // Initialize empty dictionary for histogram plotting
        std::map<int, int> resHisto;

        // Define regions for noise estimation in singleVolume and splitVolume mode
        if (!splitVolume) {
            // Calculate mask of background voxels within Rinside sphere but outside of particle mask
            maskBG = Arrayb::minus(Rinside, mask);
        }
        else {
            maskParticle = mask;
        }

        Arrayd dirs;
        Arrayd kernel, kernelSqrt;
        double kernelSum;
        Arrayd W;
        int i, j, k, l;
        Arrayd A, Ac, Ad, Ack;
        Arrayd H, Hc, Hd;
        Arrayd LAMBDA, LAMBDAc, LAMBDAd, LAMBDAdiff, LAMBDAeig;
        Arrayd LRS;
        double alpha;
        double thrUncorr;
        int kmax;
        int maskSum;
        double maskSumConst;
        int kpoint;
        double thrFDR;
        int newSumOfMask;

        // Continue testing larger and larger scales as long as there is "moreToProcess" (see below)
        std::vector<double> times;

        moreToProcess = true;

        while (moreToProcess) {
            times.push_back(dtime());

            MPI_LOG << "------------------------------------------------------------------" << std::endl;
            MPI_LOG << "[*] Calculating Local Resolution for " << currentRes << " Angstroms" << std::endl;

            // Compute window size and form steerable bases
            int r       = int(std::ceil(0.5*currentRes/vxSize)); // number of pixels around center
            double a       = (2*PI/currentRes) * std::sqrt(2.0/5); // scaling factor so that peak occurs at 1/currentRes Angstroms
            int winSize = 2*r+1;

            times.push_back(dtime());
            MPI_LOG << "[*] Computing window size: " << times[times.size()-1]-times[times.size()-2] << " seconds." << std::endl;

            // Define range of x, y, z for steerable bases
            auto && xyz = mgrid({{-r*vxSize,r*vxSize,0,double(winSize)}, {-r*vxSize,r*vxSize,0,double(winSize)}, {-r*vxSize,r*vxSize,0,double(winSize)}});
            xyz.selfTimes(a);

            auto && x = xyz.sub(0);
            auto && y = xyz.sub(1);
            auto && z = xyz.sub(2);

            dirs = make3DsteerableDirections(x, y, z);

            JN_INFOA(x);
            JN_INFOA(y);
            JN_INFOA(z);
            JN_INFOA(dirs);

            times.push_back(dtime());
            MPI_LOG << "[*] Define range of x, y, z for steerable bases: " << times[times.size()-1]-times[times.size()-2] << " seconds." << std::endl;

            // ResMap has an error here!!! The use -1/2 which is -1, but here -1.0/2 whose value is 0.5 should be used.
            // To be assistent with ResMap, I just use -1 here.
            // Actually, I think -0.5 should be used here.
            
            // Define Gaussian kernel
            kernel = Arrayd::square(x).selfPlus(Arrayd::square(y)).selfPlus(Arrayd::square(z)).selfTimes(-1).selfExp().flatten();
            kernelSqrt = Arrayd::sqrt(kernel);
            kernelSum  = kernel.sum();
            W          = Arrayd::diag(kernel);

            JN_INFOA(kernel);
            JN_INFOA(kernelSqrt);
            JN_INFOA(W);

            times.push_back(dtime());
            MPI_LOG << "[*] Define gaussian kernel: " << times[times.size()-1]-times[times.size()-2] << " seconds." << std::endl;

            // Calculate shape of matrix of bases
            int numberOfPoints = kernel.size();
            int numberOfBases  = dirs.shape(3) + 1;

            JN_INFOV(numberOfPoints);
            JN_INFOV(numberOfBases);

            // Form matrix of Hermite polynomials 
            Arrayd A({numberOfPoints, numberOfBases}, 0);

            JN_INFOA(A);

            for (i = 0, j = 0, l = 0; i < A.shape(0); i++) {
                // Form the G2 (cosine-line terms)
                A[j] = 1; j++;
                for (k = 0; k < 6; k++) { A[j] = 4*square(dirs[l]) - 2; j++; l++; }
                // Form the H2 (sine-like terms)
                for (k = 6; k < 16; k++) { A[j] = std::pow(dirs[l], 3) - 2.254 * dirs[l]; j++; l++; }
            }

            times.push_back(dtime());
            MPI_LOG << "[*] Form matrix of Hermite polynomials: " << times[times.size()-1]-times[times.size()-2] << " seconds." << std::endl;

            JN_INFOA(A);

            // Form matrix of just the constant term
            Arrayd Ac({kernel.size(), 1}, 1);

            JN_INFOA(Ac);

            // Form matrix of all but the constant term
            Ad = A.sub({{0, A.shape(0)}, {1,17}});
            JN_INFOA(Ad);

            //(linalg::dot(A, A.transpose())).identify("testA");
            //(linalg::dot(Ad, Ad.transpose())).identify("testAd");
            //(linalg::dot(Arrayd::diag(kernelSqrt), Ad)).identify("kernelSqrt.diag.dot");
            Hd = linalg::dot(Ad, linalg::dot(linalg::pinv(linalg::dot(Arrayd::diag(kernelSqrt), Ad)), Arrayd::diag(kernelSqrt)));
            JN_INFOA(Hd);

            LAMBDAd = Arrayd::minus(W, linalg::dot(W, Hd));
            JN_INFOA(LAMBDAd);

            // Invert weighted A matrix via SVD
            H = linalg::dot(A, linalg::dot(linalg::pinv(linalg::dot(Arrayd::diag(kernelSqrt), A)), Arrayd::diag(kernelSqrt)));
            JN_INFOA(H);

            // Invert weighted Ac matrix analytically
            Ack = linalg::dot(Arrayd::diag(kernelSqrt), Ac);
            JN_INFOA(Ack);

            Hc = linalg::dot(Ac, linalg::dot(Ack.transpose().selfDivide(square(Ack.norm())), Arrayd::diag(kernelSqrt)));
            JN_INFOA(Hc);

            // Create LAMBDA matrices that correspond to WRSS = Y^T*LAMBDA*Y
            LAMBDA = Arrayd::minus(W, linalg::dot(W, H));
            LAMBDAc = Arrayd::minus(W, linalg::dot(W, Hc));
            LAMBDAdiff = Arrayd::minus(LAMBDAc, LAMBDA);

            JN_INFOA(LAMBDA);
            JN_INFOA(LAMBDAc);
            JN_INFOA(LAMBDAdiff);

            times.push_back(dtime());
            MPI_LOG << "[*] Initialize for variance estimating: " << times[times.size()-1]-times[times.size()-2] << " seconds." << std::endl;
            JN_INFOA(data);

            Arrayd *pData = (splitVolume ? &dataDiff : &data);
            Arrayb *pMask = (splitVolume ? &maskParticle : &maskBG);
            JN_INFOA((*pMask));
            JN_INFOA((*pData));

            // Estimate variance
            variance = estimate_variance(pData, pMask, winSize, LAMBDAd);

            //Arrayd *pData = (splitVolume ? &dataDiff : &data);
            //Arrayb *pMask = (splitVolume ? &maskParticle : &maskBG);
            //std::list<double> WRSSBG;
            //JN_INFOA((*pMask));
            //JN_INFOA((*pData));
            //ArrayWindowRoller roller(pMask->shape(), {winSize, winSize, winSize}, {winSize, winSize, winSize});
            //JN_INFOA(data);
            //do {
            //    auto && winMaskBG = pMask->sub(roller.range);
            //    Arrayd winData = pData->sub(roller.range);
            //    if (winMaskBG.all()) {
            //        winData.selfMinus(winData.mean());
            //        winData.reshape({winData.size(),1});
            //        WRSSBG.push_back(linalg::vdot<double>(winData, linalg::dot(LAMBDAd, winData)));
            //    }
            //} while (roller.roll());
            //JN_INFOA(data);
            //JN_INFOA(LAMBDAd);

            //double sum = 0; for (auto && n : WRSSBG) sum += n;
            //variance = sum / WRSSBG.size() / LAMBDAd.trace();
            // Note: If use the following form, '->double' should not be omitted and 0.0 should not be 0!!!
            // variance = std::accumulate(WRSSBG.begin(), WRSSBG.end(), 0.0, [](double a, double b)->double{return a+b;}) / (WRSSBG.size() * LAMBDAd.trace());

            MPI_LOG << "[*] Variance: " << variance << std::endl;
            times.push_back(dtime());
            MPI_LOG << "[*] Estimate variance: " << times[times.size()-1]-times[times.size()-2] << " seconds." << std::endl;

            //// Compute Likelihood Ratio Statistic

            JN_INFOA(data);
            MPI_LOG << "[*] Begin to compute Likelihood Ratio Statistic" << std::endl;

            // Calculate weighted residual sum of squares difference

            JN_INFOA(data);
            JN_INFOA(LAMBDAdiff);
            JN_INFOA(mask);
            JN_INFOV(r);

            // set up the significant data ii
            // reduce the load imbalance
            //std::vector<int> significant_ii;
            std::vector<std::pair<int, double *>> significant_ii;
            for (int ii = 0; ii < data.size(); ii++) {
                if (mask[ii]) {
                    int a = ii / (n*n);
                    int b = (ii / n) % n;
                    int c = ii % n;
                    if (a-r>=0 && a+r<n && b-r>=0 && b+r<n && c-r>=0 && c+r <n)
                        significant_ii.push_back({ii, data.data_wptr()+ii-r*n*n-r*n-r});
                }
            }
            
            Arrayd WRSSdiff(data.shape(), 0);
            
#ifdef USEMPI

            MPI_Comm mComm;
            int mRank = mpiRank();
            int mSize = mpiSize();
            int mBin = int(std::ceil(significant_ii.size()/(double)mSize));
            int mBeg = mRank*mBin;
            
#pragma omp parallel for schedule(dynamic)
            for (int ii_i = mBeg; ii_i < std::min(mBeg+mBin, (int)significant_ii.size()); ii_i++) {
#else
                
#pragma omp parallel for schedule(dynamic)
            for (int ii_i = 0; ii_i < significant_ii.size(); ii_i++) {
                    
#endif
                int ii;
                double *it_data;
                std::tie(ii, it_data) = significant_ii[ii_i];
                int sz = r*2+1;
                Arrayd winData({sz, sz, sz});
                double *it_temp = winData.data_wptr();
                array_traverse(&data, it_data, sz, [&it_temp](double *it_data){
                    *it_temp = *it_data;
                    it_temp++;
                });
//                int ii = significant_ii[ii_i];
//                assert(mask[ii]);
//                int a = ii / (n*n);
//                int b = (ii / n) % n;
//                int c = ii % n;
//                assert(a-r>=0 && a+r<n && b-r>=0 && b+r<n && c-r>=0 && c+r <n);
//                Arrayd winData = data.sub({{a-r,a+r+1},{b-r,b+r+1},{c-r,c+r+1}});
                winData.reshape({winData.size(), 1});
                WRSSdiff[ii] = linalg::vdot<double>(winData, linalg::dot(LAMBDAdiff, winData));
            }

#ifdef USEMPI
                
            double* local_temp = (double*)aMalloc(data.size()*sizeof(double), 64);
            for (int index = 0; index < data.size(); index++) local_temp[index] = WRSSdiff[index];
            MPI::COMM_WORLD.Allreduce(local_temp,WRSSdiff.data_wptr(),data.size(),MPI_DOUBLE,MPI::SUM);
            aFree(local_temp);
                
            // double *tempDiff;
            // if (MPI_IS_ROOT) tempDiff = new double[mSize*mBin];
            // MPI_Gather(WRSSdiff.data+mBeg, mBin, MPI_DOUBLE, tempDiff, mBin, MPI_DOUBLE, 0, MPI_COMM_WORLD);
            // if (MPI_IS_ROOT) delete [] WRSSdiff.data;
            // if (MPI_IS_ROOT) WRSSdiff.data = tempDiff;
            // MPI_Bcast(WRSSdiff.data, data.size(), MPI_DOUBLE, 0, MPI_COMM_WORLD);

#endif

            times.push_back(dtime());
            MPI_LOG << "[*] Calculate weighted residual sum of squares difference: " << times[times.size()-1]-times[times.size()-2] << " seconds." << std::endl;

            std::list<double> LRSvec;
            LRS = Arrayd({n,n,n}, 0);
            for (i = 0; i < n*n*n; i++) {
                if (mask[i]) {
                    LRSvec.push_back(WRSSdiff[i] / (variance+1e-10));
                    LRS[i] = LRSvec.back();
                }
            }

            //// Numerically Compute Weighted X^2 Statistic

            // Calculate Eigenvalues of LAMBDAdiff
            LAMBDAeig = Arrayd::abs(linalg::eig(LAMBDAdiff).wr);

            JN_INFOA(LAMBDAeig);

            // Remove small values and truncate for numerical stability
            LAMBDAeig = Arrayd::extract(Arrayb::gt(LAMBDAeig, LAMBDAeig.max()*1e-1), LAMBDAeig);

            JN_INFOA(LAMBDAeig);

            // Uncorrected Threshold
            JN_INFOV(pValue);

            alpha = 1-pValue;

            JN_INFOV(alpha);

            auto minResults = minimize_scalar(std::bind(evaluateRuben, std::placeholders::_1, alpha, LAMBDAeig));
            thrUncorr  = minResults.xmin;

            JN_INFOV(thrUncorr);

            // FDR Threshold
            std::vector<double> LRSvecSorted;
            LRSvecSorted.reserve(LRSvec.size());
            for (auto &&n : LRSvec) {
                if (n > thrUncorr) {
                    LRSvecSorted.push_back(n);
                }
            }
            //LRSvecSorted.shrink_to_fit();
            kmax = LRSvecSorted.size();
            std::sort(LRSvecSorted.begin(), LRSvecSorted.end());

            JN_INFOV(kmax);

            times.push_back(dtime());
            MPI_LOG << "[*] Elapsed time: " << times[times.size()-1]-times[times.size()-2] << " seconds." << std::endl;

            if (kmax == 0) {
                thrFDR = 1.7976931348623157e+308;
            }
            else {
                maskSum = mask.sum<int>();
                JN_INFOV(maskSum);
                maskSumConst = Arrayd::range(1, maskSum).selfDividedBy(1.0).sum();
                JN_INFOV(maskSumConst);
                for (k = 1; k < kmax; k += int(std::ceil(kmax / std::min(5e2, double(kmax))))) { // k begins with 1 because it's 1 in python codes
                    auto result = ruben(LAMBDAeig,LRSvecSorted[k]);
                    double tmp    = 1.0-(pValue*((kmax-k)/(maskSum*maskSumConst)));
                    thrFDR = LRSvecSorted[k];
                    if (result.res > tmp) {
                        break;
                    }
                }
            }

            times.push_back(dtime());
            MPI_LOG << "[*] Elapsed time: " << times[times.size()-1]-times[times.size()-2] << " seconds." << std::endl;

            JN_INFOV(thrFDR);

            // Calculate resolution
            auto && res = Arrayi::gt(LRS, thrFDR);
            resTOTAL += res * currentRes;

            JN_INFOA(LRS);
            JN_INFOA(res);
            JN_INFOA(resTOTAL);

            // Update the mask to voxels that failed this level's likelihood test
            for (i = 0; i < mask.size(); i++) {
                if (int(mask[i]) != res[i]) mask[i] = true;
                else mask[i] = false;
            }
            newSumOfMask = mask.sum<int>();

            JN_INFOA(mask);
            JN_INFOV(oldSumOfMask);
            JN_INFOV(newSumOfMask);

            // Update value in histogram
            resHisto[currentRes] = oldSumOfMask - newSumOfMask;
            JN_INFOV(resHisto[currentRes]);

            // Heuristic of telling whether we are likely done
            if (oldSumOfMask-newSumOfMask < n && newSumOfMask < n*n) {
                moreToProcess = false;
            }
            oldSumOfMask = newSumOfMask;

            JN_INFOV(currentRes);
            JN_INFOV(maxRes);
            if (currentRes >= maxRes) {
                moreToProcess = false;
            }

            // Update current resolution
            currentRes += stepRes;

            times.push_back(dtime());
            MPI_LOG << "[*] Calculate resolution: " << times[times.size()-1]-times[times.size()-2] << " seconds." << std::endl;

        }

        // Set all voxels that were outside of the mask or that failed all resolution tests to 100 A
        for (i = 0; i < resTOTAL.size(); i++) if (resTOTAL[i] == 0) resTOTAL[i] = 100;
        std::list<double> resTOTALma;
        double total_mean = 0;
        for (auto && n : resTOTAL) if (n <= currentRes) {
            total_mean += n;
            resTOTALma.push_back(n);
        }
        total_mean /= double(resTOTALma.size());

        // Print results
        MPI_LOG << "MEAN RESOLUTION in MASK: " << total_mean << " Angstrom" << std::endl;
        //MPI_LOG << "MEDIAN RESOLUTION in MASK: " << resTOTALma.median();

        int old_n;
        if (splitVolume) old_n = volume1.size;
        else old_n = volume.size;

        auto && old_coordinates = mgrid({
                {0.0, double(n-1), 0.0, double(old_n)},
                {0.0, double(n-1), 0.0, double(old_n)},
                {0.0, double(n-1), 0.0, double(old_n)}});

        // Up-sample the resulting resolution map if necessary
        if (LPFfactor > 0) {
            resTOTAL = cppsci::ndimage::map_coordinates(resTOTAL, old_coordinates, 1, "nearest");
            for (auto && n : resTOTAL) {
                if (n <= minRes) n = minRes;
                if (n > currentRes) n = 100;
            }
        }

#ifdef USEMPI
    if (MPI_IS_ROOT) {
#endif

        // Write results out as MRC volume
        std::string fname, ext;
        if (splitVolume) {
            std::tie(fname, ext) = os::path::splitext(inputFileName1);
            Arrayf temp(resTOTAL);
            writeMrcData(fname + "_resmap" + ext, temp.data_wptr(), temp.shape(0), *(volume1.head));
        }
        else {
            std::tie(fname, ext) = os::path::splitext(inputFileName);
            Arrayf temp(resTOTAL);
//            temp.identify("resTOTAL");
            writeMrcData(fname + "_resmap" + ext, temp.data_wptr(), temp.shape(0), *(volume.head));
        }

        std::string chimeraScriptFileName ;
        if (splitVolume) {
            chimeraScriptFileName = createChimeraScript(inputFileName1, minRes, maxRes, int(resTOTAL.shape(0)), true);
        }
        else {
            chimeraScriptFileName = createChimeraScript(inputFileName, minRes, maxRes, int(resTOTAL.shape(0)), true);
        }

#ifdef USEMPI
    }
#endif

    }

    /**
     * All the steps for resolution calculation.
     */
    void calculate() {
        std::vector<double> times;

        times.push_back(dtime());

        MPI_LOG << std::endl;
        MPI_LOG << "================== Initializing ================" << std::endl;
        init();
        times.push_back(dtime());
        MPI_LOG << "------------------------------------------------" << std::endl;
        MPI_LOG << "Elapsed time for initializing: " << times[times.size()-1] -  times[times.size()-2] << " seconds." << std::endl;
        MPI_LOG << "================================================\n" << std::endl;

        MPI_LOG << "================== Checking LPF ================" << std::endl;
        checkLPF();
        times.push_back(dtime());
        MPI_LOG << "------------------------------------------------" << std::endl;
        MPI_LOG << "Elapsed time for checking LPF: " << times[times.size()-1] -  times[times.size()-2] << " seconds." << std::endl;
        MPI_LOG << "================================================\n" << std::endl;

        MPI_LOG << "================== Setting Mask ================" << std::endl;
        setMask();
        times.push_back(dtime());
        MPI_LOG << "------------------------------------------------" << std::endl;
        MPI_LOG << "Elapsed time for setting mask: " << times[times.size()-1] -  times[times.size()-2] << " seconds." << std::endl;
        MPI_LOG << "================================================\n" << std::endl;

        MPI_LOG << "================== Pre-Whitening ===============" << std::endl;
        preWhiten();
        times.push_back(dtime());
        MPI_LOG << "------------------------------------------------" << std::endl;
        MPI_LOG << "Elapsed time for pre-whitening: " << times[times.size()-1] -  times[times.size()-2] << " seconds." << std::endl;
        MPI_LOG << "================================================\n" << std::endl;

        MPI_LOG << "============= Resolution Calculation ============" << std::endl;
        compute();
        times.push_back(dtime());
        MPI_LOG << "------------------------------------------------" << std::endl;
        MPI_LOG << "Elapsed time for resolution calculation: " << times[times.size()-1] -  times[times.size()-2] << " seconds." << std::endl;
        MPI_LOG << "================================================\n" << std::endl;
    }

} // end namespace resumap

