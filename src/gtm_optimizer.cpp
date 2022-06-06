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

#include "./gtm_optimizer.h"

// How much memory, total, we are allowed to transfer to Xeon Phi coprocessors
#define OFFLOAD_AMOUNT (offloadlimit)

#define ALLOC alloc_if(1) free_if(0)
#define FREE alloc_if(0) free_if(1)
#define REUSE alloc_if(0) free_if(0)
#define TEMP alloc_if(1) free_if(1)

#define NUM_LOCAL_COMPUTE_ENGINES (numMIC+1)

//TODO - limit work done by all Xeon Phis so that none of them
//ever runs longer than Xeon - we don't want to have to wait on
//sync variables
/// mic_index = 0;  if (iter_index == 1) break; \

#define START_HETERO_OFFLOAD \
    int node_d_start, node_d_end, node_d_chunksize; \
    for(int iter_index = 0;iter_index < NUM_LOCAL_COMPUTE_ENGINES;iter_index++) \
    { \
        int hostamnt; \
        int micstart; \
        \
        int mic_index = NUM_LOCAL_COMPUTE_ENGINES-1-iter_index; \
        \
        if (0==mic_index) \
        { \
            node_d_start = 0; \
            node_d_end = d_interval - estnodesize * numMIC; \
        } \
        else \
        { \
            hostamnt = d_interval - estnodesize * numMIC; \
            micstart = hostamnt;  \
            \
            node_d_start = micstart + (mic_index-1) * estnodesize; \
            node_d_end = ((node_d_start + estnodesize) < d_interval) ? \
            (node_d_start + estnodesize): d_interval; \
        } \
        node_d_chunksize = node_d_end - node_d_start;

#define STOP_HETERO_OFFLOAD  }


namespace GTMoptimizer {
    template<typename T>
    void check(T t,char* what){
        if (std::isnan(t)||std::isinf(t)) {
            std::cout<<"nan ro inf..."<<what<<std::endl<<std::flush;
            exit(1);
        }
    }
    // implicitly initialized
    int node,nodes,numMIC,nodeNameLen,messageTag;
#ifdef USEMPI
    char nodeName[MPI_MAX_PROCESSOR_NAME];
#else
    char nodeName[4096];
#endif
    
    MIC_OFFLOAD_ATTRIBUTES int d_start, d_end, d_interval;

    MIC_OFFLOAD_ATTRIBUTES int maxthreads;
    MIC_OFFLOAD_ATTRIBUTES int numXeonPhithreads;
    
    MIC_OFFLOAD_ATTRIBUTES
        __attribute__((aligned(64))) double* PHI;
    
    MIC_OFFLOAD_ATTRIBUTES int K;
    MIC_OFFLOAD_ATTRIBUTES int M;
    MIC_OFFLOAD_ATTRIBUTES int N;
    MIC_OFFLOAD_ATTRIBUTES int D;

    MIC_OFFLOAD_ATTRIBUTES
        __attribute__((aligned(64))) float* TReal;
    MIC_OFFLOAD_ATTRIBUTES
        __attribute__((aligned(64))) float* TImag;

    // contrast transfer function
    MIC_OFFLOAD_ATTRIBUTES
        __attribute__((aligned(64))) float* CTFvalue;   //DxN

    MIC_OFFLOAD_ATTRIBUTES
        __attribute__((aligned(64))) float* WReal;
    MIC_OFFLOAD_ATTRIBUTES
        __attribute__((aligned(64))) float* WImag;

    MIC_OFFLOAD_ATTRIBUTES
        __attribute__((aligned(64))) double* YReal;
    MIC_OFFLOAD_ATTRIBUTES
        __attribute__((aligned(64))) double* YImag;

    MIC_OFFLOAD_ATTRIBUTES
        __attribute__((aligned(64))) double* R;

    MIC_OFFLOAD_ATTRIBUTES
        __attribute__((aligned(64))) float* Alpha;

    MIC_OFFLOAD_ATTRIBUTES double AverageBeta;
    MIC_OFFLOAD_ATTRIBUTES double AverageAlpha;

    MIC_OFFLOAD_ATTRIBUTES
        __attribute__((aligned(64))) float *avTReal;
    MIC_OFFLOAD_ATTRIBUTES
        __attribute__((aligned(64))) float *avTImag;
    
    MIC_OFFLOAD_ATTRIBUTES
        __attribute__((aligned(64))) int *ipiv;
    
    FineSearch fineSearch;
    
    // (MnlxL)
    int Mnl;
    
    double deviation;
    
    size_t micavail;
    int estnodesize;
    int maxXeonPhi;
    size_t offloadlimit;
    double loadMIC;
    
    int picdim;
    std::string star_fn,result_fn,ref_fn = "NULL";
    double pixel_size;
    bool do_ctf_correction;
    bool only_flip_phases;
    bool ctf_phase_flipped;
    bool intact_ctf_first_peak;
    MetaDataTable metadata;
    int mrcsHead[256];
    char timerstring[16384];
    char buffer[256];
    bool init_by_pca;
    int iter, nr_iter;
    
/// ---------------------------------------------------------------------   
void GTMoptimizer::setupGTMoptimizer(double X_infos[],double Mu_infos[],int set_sampling_dim)
{
    
    double fstart,fstop;
    fstop = 0.;
    
    NODE0ONLY fstart = dtime();
#ifdef USEMPI
    // init MPI outside
    // MPI::Init();
    nodes = MPI::COMM_WORLD.Get_size();
    node = MPI::COMM_WORLD.Get_rank();
    MPI::Get_processor_name(nodeName,nodeNameLen);
#else
    nodes = 1;
    node = 0;
#endif
    
#ifdef __INTEL_OFFLOAD
    // Note - due to a compiler defect, calling this function forces the compiler
    // to assume offload is going to happen - may need to stub out using
    // a build-time define
    ERROR_CHECK(numMIC < 0, "Set the numMic should large than zero.");
    // Note - due to a compiler defect, calling this function forces the compiler
    // to assume offload is going to happen.  Which is why we allow the user to
    // overwrite up-front.
    max_numMIC = _Offload_number_of_devices();
    if (numMIC > max_numMIC)
        numMIC = max_numMIC;
    
#else
    numMIC = 0;
#endif
    
    // How much memory, total, we are allowed to transfer to Xeon Phi coprocessors
    // By default, we are only allowed to use 85% of physical RAM on the card.
    // For 8GB, that's about 6.8GB; for 16GB, about 13.6GB
    // Ideally someday there will be a way to sense the memory capacity of an
    // Intel Xeon Processor.  Until then, we hardcode (or possibly more ideally)
    // make it an input parameter.
    offloadlimit = 6500000000;
    
    maxthreads = omp_get_max_threads();
    numXeonPhithreads = 60;
    
    // result_fn = set_result_fn;
    // do_ctf_correction = _do_ctf_correction;
    // initialize K,M,PHI -----------------------------------
    int L,L2;
    double *X;
    double *Mu;
    switch (set_sampling_dim) {
        case 1:
            X = GTMSampling::sample1L(X_infos,K,L);
            Mu= GTMSampling::sample1L(Mu_infos,Mnl,L2);
            break;
        case 2:
            X = GTMSampling::sample2L(X_infos,K,L);
            Mu= GTMSampling::sample2L(Mu_infos,Mnl,L2);
            break;
        case 3:
            X = GTMSampling::sample3L(X_infos,K,L);
            Mu= GTMSampling::sample3L(Mu_infos,Mnl,L2);
            break;
        default:
            break;
    }
    
    M = Mnl + 1;
    
    PHI = (double*)aMalloc(sizeof(double)*K*M,64);
    
    switch (set_sampling_dim) {
        case 1:
            // Basis rotation_basis;
            GTMBasis::set_rotation_basis(X,K,L,Mu,Mu_infos,Mnl,PHI,M);
            break;
        case 2:
            // Basis shift_basis;
            GTMBasis::set_shift_basis(X,K,L,Mu,Mu_infos,Mnl,PHI,M);
            break;
        case 3:
            // Basis shift_basis;
            GTMBasis::set_shift_basis(X,K,L,Mu,Mu_infos,Mnl,PHI,M);
            break;
        default:
            break;
    }
    
    if (0 == node) std::cout<<"K = "<<K<<",Mnl = "<<Mnl<<std::endl;
    // initialize CTFValue TReal TImag
    
    // initialize picture metadata
    metadata.readFromStar(star_fn);
    
    // std::string mrcs_dir = pathGetDir(star_fn);
    
    if (fineSearch.doFineSearch) {
        metadata.fliterByClass(fineSearch.getCurrentClass());
    }
    
    N = metadata.numberOfParticles();
    
    // TODO - This only works if the initial mrcs file and the output generated by
    // ml2d are in the same directory.  To avoid accidentally over-writing
    // results, recommend that a different directory be allowed for input and
    // ML2D output mrcs files.
    FILE* mrcsFile = fopen((/*mrcs_dir+*/metadata[0].IMAGE.NAME).c_str(),"rb");
    if (NULL == mrcsFile)
        ERROR_REPORT("make sure you have "+/*mrcs_dir+*/metadata[0].IMAGE.NAME);
    
    if(fread((char*)mrcsHead,256*sizeof(float),1,mrcsFile) == NULL)
        ERROR_REPORT("make sure you have "+/*mrcs_dir+*/metadata[0].IMAGE.NAME+".mrcs");
    
    
    fclose(mrcsFile);
    
    // check error
    if(mrcsHead[0] != mrcsHead[1])
        ERROR_REPORT("mrcs image is not square.");
    // checkerror(N > mrcsHead[2], "N > mrcshead[2]");
    if(mrcsHead[3] != 2)
        ERROR_REPORT("mrcs image is not float image.");
    
    picdim = mrcsHead[0];
    
    if(0 == node) std::cout<<"N = "<<N<<",dim = "<<picdim<<std::endl;
    
    D = picdim*(picdim/2+1);
    int D0 = picdim*picdim;
    
    int sub_D = ceil((double)D/nodes);
    d_start = 0 + node*sub_D;
    d_end = (d_start+sub_D) < D?(d_start+sub_D):D;
    d_interval = d_end-d_start;
    
    
    // initialize CTFValue
    CTFvalue = (float*)aMalloc(sizeof(float)*N*d_interval,64);
    // Have the data been CTF phase-flipped?
    // bool ctf_phase_flipped = false;
    // Only perform CTF phase-flipping? (default is full amplitude-correction)
    // bool only_flip_phases = false;
    // Ignore CTFs until their first peak
    // bool intact_ctf_first_peak = false;
    // pixel_size = set_pixel_size;
    
#pragma omp parallel for
    for(int n = 0;n < N;n++){
        if (do_ctf_correction) {
            CTF ctf;
            double *CTFbuffer = (double*)aMalloc(sizeof(double)*D,64);
            
            ctf.setValues(metadata[n].CTF_DEFOCUS_U,metadata[n].CTF_DEFOCUS_V,metadata[n].CTF_DEFOCUS_ANGLE,metadata[n].CTF_VOLTAGE,\
                          metadata[n].CTF_CS,metadata[n].CTF_Q0,metadata[n].CTF_BFAC);
            
            if(n == (N-1)&& 0 == node) std::cout<<metadata[n].CTF_DEFOCUS_U<<" "<<metadata[n].CTF_DEFOCUS_V<<" "<<metadata[n].CTF_DEFOCUS_ANGLE<<" " \
                <<metadata[n].CTF_VOLTAGE<<" "<<metadata[n].CTF_CS<<" "<<metadata[n].CTF_Q0<<" "<<metadata[n].CTF_BFAC<<std::endl;
            
            ctf.getFftwImage(CTFbuffer,picdim,picdim, picdim,pixel_size,ctf_phase_flipped, only_flip_phases, intact_ctf_first_peak, true);
            
            // transposition CTFvalue
            for(int d = 0;d < d_interval;d++)
                CTFvalue[INDEXSZ(INDEXSZ(d)*INDEXSZ(N))+INDEXSZ(n)] = CTFbuffer[d+d_start];
            
            aFree(CTFbuffer);
        }
        else
        {
            for(int d = 0;d < d_interval;d++)
                CTFvalue[INDEXSZ(INDEXSZ(d)*INDEXSZ(N))+INDEXSZ(n)] = 1.;
        }
    }
    
    
    
    // FourierTransformerDFTI fftford_worker(picdim,picdim);
    // MKL_Complex8 *fft_out = (MKL_Complex8*)aMalloc(picdim*(picdim/2+1)*sizeof(MKL_Complex8),64);
    FFTWFTransformer fftford_worker(picdim,picdim);
    SOAComplexArray<float> fft_out(picdim*(picdim/2+1));
    // initialze TReal TImag
    double* T;
    if (init_by_pca) {
        NODE0ONLY std::cout<<"...use PCA to initialize W & averagebelta..."<<std::endl;
        NODE0ONLY std::cout<<"...double precision for Images data,do PCA in one node..."<<std::endl;
        NODE0ONLY std::cout<<"...image data dimension is "<<N<<"x"<<D0<<",size "<<(double(N)*double(D0)*8./1024./1024./1024.)<<" GB..."<<std::endl;
        NODE0ONLY std::cout<<"...need to do for PCA,1)float datatype precision,2)multi-nodes support..."<<std::endl;
        NODE0ONLY T = (double*)aMalloc(sizeof(double)*N*D0,64);
    }
    TReal = (float*)aMalloc(sizeof(float)*d_interval*N,64);
    TImag = (float*)aMalloc(sizeof(float)*d_interval*N,64);
    
    float *buffer = (float*)aMalloc(sizeof(float)*picdim*picdim,64);
    float *buffer2 = (float*)aMalloc(sizeof(float)*picdim*picdim,64);
    // fft data
    int particle_diameter = int(pixel_size*picdim);
    int width_mask_edge = 5;
    
    NODE0ONLY std::cout<<"starting bcast data."<<std::endl;
    double read_time = 0,trans_time = 0,fft_time = 0,bcast_time = 0;
    // open file first
    std::string preFileName = metadata[0].IMAGE.NAME;
    mrcsFile = fopen((/*mrcs_dir+*/metadata[0].IMAGE.NAME).c_str(),"rb");
    
    for(int n = 0; n < N; n++)
    {
        read_time -= dtime();
        
        // NODE0ONLY
        // {
        // if file name is changed,re-open file
        if (metadata[n].IMAGE.NAME != preFileName)
        {
            fclose(mrcsFile);
            mrcsFile = fopen((/*mrcs_dir+*/metadata[n].IMAGE.NAME).c_str(),"rb");
            preFileName = metadata[n].IMAGE.NAME;
        }
        // change type to int may cause bug,avoid offset out of range
        long image_id = metadata[n].IMAGE.INDEX;
        
        long offset = (256+(image_id-1)*picdim*picdim)*sizeof(float);
        
        fseek(mrcsFile,offset,SEEK_SET);
        
        if(fread((char*)buffer,picdim*picdim*sizeof(float),1,mrcsFile) == NULL)
            ERROR_REPORT("read mrcs data failed.");
        
        read_time += dtime();
        
        trans_time -= dtime();
        
        Spider::RT_SF(buffer, buffer2, picdim, picdim, 0, metadata[n].XOFF, metadata[n].YOFF);
        Spider::RT_SF(buffer2, buffer, picdim, picdim, -metadata[n].PSI, 0,0);
        
        if (init_by_pca) for(int d = 0;d < D;d++) T[n*D+d] = buffer[d];
        
        NODE0ONLY {
            if(n == N-1) std::cout<<metadata[n].PSI<<" "<<metadata[n].XOFF<<" "<<metadata[n].YOFF<<" "<<metadata[n].IMAGE<<std::endl;
        }
        
        trans_time += dtime();
        
        fft_time -= dtime();
        // TODO : testing the normalize data case
        fftford_worker.FourierTransform(buffer, fft_out,false);
        // }
        
        fft_time += dtime();
        
        bcast_time -= dtime();
        
// #ifdef USEMPI
//         // MPI::COMM_WORLD.Bcast((float*)fft_out,2*D,MPI::FLOAT,0);
//         MPI::COMM_WORLD.Bcast(fft_out.real,D,MPI::FLOAT,0);
//         MPI::COMM_WORLD.Bcast(fft_out.imag,D,MPI::FLOAT,0);
// #endif
        
        for(int i = d_start;i < d_end;i++){
            TReal[INDEXSZ((INDEXSZ(i)-INDEXSZ(d_start))*INDEXSZ(N))+INDEXSZ(n)] = fft_out.real[i];
            TImag[INDEXSZ((INDEXSZ(i)-INDEXSZ(d_start))*INDEXSZ(N))+INDEXSZ(n)] = fft_out.imag[i];
        }
        
        bcast_time += dtime();
    }
    
    NODE0ONLY std::cout<<" finishing bcast data,read data time : "<<read_time<<", transform image time : "<<trans_time<<std::endl;
    NODE0ONLY std::cout<<" fft time : "<<fft_time<<", bcast data time : "<<bcast_time<<std::endl;
    
    fclose(mrcsFile);
    
    aFree(buffer);
    aFree(buffer2);
    
    
    //  initialize W  -------------------------------------
    WReal = (float*)aMalloc(sizeof(float)*d_interval*M,64);   //should be K*D*M
    WImag = (float*)aMalloc(sizeof(float)*d_interval*M,64);  //////////////////
    
    // ipiv = (int*)aMalloc(sizeof(int)*M,64);
    
    float *W = (float*)aMalloc(sizeof(float)*M*D0,64);
    float *refData = (float*)aMalloc(sizeof(float)*K*D0,64);
    
    // initialize W & belta by PCA if necessary
    // (double *_T,int _N,int _D,double *_X,int _K,int _L,double *_PHI,double *_W,int _M,double* _Belta);
    double beltaTemp;
    double* W_double;
    if (init_by_pca) W_double = (double*)aMalloc(sizeof(double)*M*D0,64);
    Pca pca_initializer(T,N,D0,X,K,L,PHI,W_double,M,&beltaTemp);
    
    if (ref_fn != "NULL")
    {
        // using template to initialize W
        int refHead[256];
        FILE* refFile = fopen((ref_fn).c_str(),"rb");
        
        fread((char*)refHead,256*sizeof(float),1,refFile);
        
        // check the size it should be same as mrcs data
        int ref_N = refHead[2];
        
        if(refHead[0] != refHead[1])
            ERROR_REPORT("reference is not square.");
        if(refHead[0] != picdim)
            ERROR_REPORT("reference dimension is not equal to mrcs image.");
#if 1
        // only need K's reference image
        int ref_N2 = std::min(K,ref_N);
        for (int n = 0; n < ref_N2; n++) {
            if(fread((char*)(refData+n*picdim*picdim),picdim*picdim*sizeof(float),1,refFile) == NULL)
                ERROR_REPORT("read reference data failed.");
        }
        
        if(ref_N < K){// using the first reference image
            if(0 == node) std::cout<<"your import reference data is smaller than K,so throw rest part of data will using only first reference."<<std::endl;
            for(int k = 1;k < K;k++) memcpy(refData+k*D0,refData,sizeof(float)*D0);
        }
        else if(ref_N > K){
            if(0 == node) std::cout<<"your import reference data is larger than K,so throw larger part of data."<<std::endl;
        }
        else{// ref_N == K,do nothing
            // if(node == 0) std::cout<<"reference data right."<<std::endl;
        }
#else
        int select_class = 17;
        long offset = (256+(select_class-1)*picdim*picdim)*sizeof(float);
        fseek(refFile,offset,SEEK_SET);
        float* refData_tmp = (float*)aMalloc(sizeof(float)*picdim*picdim,64);
        fread((char*)refData_tmp,picdim*picdim*sizeof(float),1,refFile);
        for (int k = 0; k < K; k++) {
            memcpy(refData+k*picdim*picdim, refData_tmp, picdim*picdim*sizeof(float));
        }
        aFree(refData_tmp);
        
#endif
        fclose(refFile);
    }
    else
    {
        for (int n = 0; n < K; n++) {
            genGaussCircleRef(refData+n*D0, picdim);
        }
    }
    
    if (init_by_pca){
        NODE0ONLY std::cout<<"...run PCA initializer..."<<std::endl;
        NODE0ONLY pca_initializer.run_pca();
        NODE0ONLY std::cout<<"...complete PCA...copy data..."<<std::endl;
        NODE0ONLY std::cout<<"...belta "<<beltaTemp<<" ..."<<std::endl;
        NODE0ONLY for(int i = 0;i < M*D0;i++) W[i] = W_double[i];
#ifdef USEMPI
        MPI::COMM_WORLD.Bcast(W,M*D0,MPI::FLOAT,0);
        MPI::COMM_WORLD.Bcast(&beltaTemp,1,MPI::DOUBLE,0);
#endif
        // initialize average beta
        AverageBeta = beltaTemp;
        if (node == 1) std::cout<<"...AverageBeta in node 1 "<<AverageBeta<<" ..."<<std::endl;
        NODE0ONLY aFree(T);
        NODE0ONLY aFree(W_double);
        
    }
    else{
        solveW(PHI,K,M,W,refData,D0);
        // initialize average beta
        AverageBeta = 1.;
    }
#ifdef USEMPI
    MPI::COMM_WORLD.Barrier();
#endif
    
    for(int m = 0;m < M;m++){
        
        fftford_worker.FourierTransform(W+m*D0, fft_out,false);
        
        for(int d = 0;d < d_interval;d++){
            WReal[INDEXSZ(INDEXSZ(d)*INDEXSZ(M))+INDEXSZ(m)] = fft_out.real[d+d_start];
            WImag[INDEXSZ(INDEXSZ(d)*INDEXSZ(M))+INDEXSZ(m)] = fft_out.imag[d+d_start];
        }
    }
    
    
    
    aFree(refData);
    aFree(W);
    
    double node0memGB = prepare();
    
    NODE0ONLY fstop=dtime();
    // Prepare timing data string we can easily import into a spreadsheet
    NODE0ONLY std::cout<<" ------- Prepare (data read and allocations sans FFTs): "<<fstop-fstart<<" sec "<<std::endl<<std::flush;
    NODE0ONLY sprintf(timerstring, "XF,%d,%d,%d,%d,%d,%d,%d,%ld,%d,%f,%f", N, picdim,
                      K, M, nodes, maxthreads, numXeonPhithreads, offloadlimit, numMIC, node0memGB, fstop-fstart);
    
}

void GTMoptimizer::destroyGTMoptimizer(){
    
    
    aFree(R);
    
    aFree(YReal);
    aFree(YImag);
    
    aFree(TReal);
    aFree(TImag);
    aFree(WReal);
    aFree(WImag);
    
    aFree(CTFvalue);
    
    aFree(ipiv);
    
    aFree(PHI);
    
}
    
double prepare(){
    
    NODE0ONLY dumpT();
    
    // mkl_solve gets called with size "M"
    ipiv = (int*)aMalloc(sizeof(int)*M,64);
    
	// Delta = (double*)aMalloc(sizeof(double)*K*N,64);
	R = (double*)aMalloc(sizeof(double)*K*N,64);

	YReal = (double*)aMalloc(sizeof(double)*d_interval*K,64);
	YImag = (double*)aMalloc(sizeof(double)*d_interval*K,64);

    NODE0ONLY{
		std::cout<<"space needed for each node: "<<calSpace()<<" GB"<<std::endl;
        std::cout<<"maxthreads = "<<omp_get_max_threads()<<std::endl;
        std::cout<<mrcsHead[2]<<" "<<mrcsHead[1]<<" "<<mrcsHead[0]<<std::endl;
    }
    
    double node0memGB = printDataSize();
    
    // How much space we have for everything except R on the coprocessor
    micavail = (size_t)OFFLOAD_AMOUNT - (size_t)N*(size_t)K*sizeof(double) - (size_t)M*(size_t)K*sizeof(double);
    estnodesize = 0;
    // The maximum fraction of d_interval the coprocessor can traverse
    // given the values of N, M, and K
    estnodesize = (int)(micavail / ((size_t)N*sizeof(float)*(size_t)3 +
                                    (size_t)K*sizeof(double)*(size_t)2 +
                                    (size_t)M*sizeof(float)*(size_t)2));
    
    NODE0ONLY std::cout<<"micavail:  "<<micavail<<" R size: "<<N*K*sizeof(double)<<" PHI size: "<<M*K*sizeof(double)<<std::endl<<std::flush;    NODE0ONLY std::cout<<"estnodesize: "<<estnodesize<<std::endl<<std::flush;
    NODE0ONLY std::cout<<"N*sizeof(float)*3:  "<<N*sizeof(float)*3<<std::endl<<std::flush;
    NODE0ONLY std::cout<<"K*sizeof(double)*2:  "<<K*sizeof(double)*2<<" M*sizeof(float)*2: "<<M*sizeof(float)*2<<std::endl<<std::flush;
    
    // Don't bother with offload if the problem is too small to make
    // it worthwhile.  "Too small" is defined as d_internal is less than
    // estnodesize*(numMIC+1), which would be an even split of the problem
    // between host and coprocessors.
    // double loadMIC = set_loadMIC;
    if (loadMIC != 0)
        estnodesize = d_interval*loadMIC/numMIC;
    else
    {
        if ((d_interval < (estnodesize*(numMIC+1))) && (numMIC>0))
        {
            // If we are here, likely we think we can handle the entire problem on
            // the each coprocessor.  Since we want to spread the work across all compute
            // resources, so let's thread things out
            if ((d_interval/(numMIC+1)) > 1000)
            {
                // Divide the work evenly among compute resources when there are a
                // large number of d iterations on this node
                // TODO - refine this heuristic
                estnodesize = d_interval/(numMIC+1);
            }
            else
            {
                // There really isn't enough work to spread out to coprocessors
                numMIC = 0;   // Turn off coprocessing
                NODE0ONLY std::cout<<"Disabling coprocessing since node's problem size is too small to justify coprocessing"<<std::endl<<std::flush;
            }
        }
    }

    int hostamnt = d_interval - estnodesize * numMIC;
    
    size_t hostsize =   (size_t)hostamnt*(size_t)N*sizeof(float)*(size_t)3 +
                        (size_t)hostamnt*(size_t)K*sizeof(double)*(size_t)2 +
                        (size_t)hostamnt*(size_t)M*sizeof(float)*(size_t)2 +
                        (size_t)K*(size_t)M*sizeof(double) +
                        (size_t)K*(size_t)N*sizeof(double);
    
    NODE0ONLY std::cout<<"Problem on each co-processor: "<<estnodesize<<" bytes"<<std::endl;
    NODE0ONLY std::cout<<"Problem left on host: size "<<hostsize<<" bytes,  local d indices 0 to "<<hostamnt<<std::endl<<std::flush;
    
    

#ifdef USEMPI
	MPI::COMM_WORLD.Barrier();
#endif
    return(node0memGB);
}

/// ---------------------------------------------------------------------

void updateR()
{
  
	static double sumSmall,sumBig;

	shiftMat(R,K,N);
  
#pragma omp parallel for
	for(size_t i = 0;i < K*N;i++)
		R[i] = exp(-R[i]);


#pragma omp parallel for private(sumSmall,sumBig)
	for (int n = 0; n < N; n++){

        sumSmall = 0.;
        sumBig = 0.;
        for (int k = 0; k < K; k++){
            if(R[INDEXSZ(INDEXSZ(k)* INDEXSZ(N))+ INDEXSZ(n)] < 1e-5) 
                sumSmall += R[INDEXSZ(INDEXSZ(k)* INDEXSZ(N))+ INDEXSZ(n)];
            else 
                sumBig += R[INDEXSZ(INDEXSZ(k)* INDEXSZ(N))+ INDEXSZ(n)];   //reduce the round off error
        }

#pragma omp simd   //now suitable
        for (int k = 0; k < K; k++)
            R[INDEXSZ(INDEXSZ(k)* INDEXSZ(N))+ INDEXSZ(n)] = 
                R[INDEXSZ(INDEXSZ(k)* INDEXSZ(N))+ INDEXSZ(n)]/(sumSmall+sumBig);

	}

#ifdef DEBUGGTM
    NODE0ONLY dumpR();
#endif
}

/// ---------------------------------------------------------------------

void iterateUpdateW()
{
    double fstart,fstop;
    
    // Divide up the work between the host and any Xeon Phi cards in the system
    // so that both are working to create the final R at once
    
    START_HETERO_OFFLOAD
    
        fstart = dtime();

        NODE0ONLY std::cout<<"iterateUpdateW section "<<mic_index<<" node_d_start="<<node_d_start<<", node_d_end="<<node_d_end \
                    <<", node_d_chunksize= "<<node_d_chunksize<<std::endl<<std::flush;

#ifndef DONT_USE_OFFLOAD
#ifdef SERIAL_OFFLOAD
    #pragma offload if(mic_index>0) target(mic:(mic_index-1)) \
            nocopy ( TReal[(size_t)node_d_start*(size_t)N:(size_t)node_d_chunksize*(size_t)N] \
            : alloc( TReal[(size_t)node_d_start*(size_t)N:(size_t)node_d_chunksize*(size_t)N]) REUSE) \
            nocopy ( TImag[(size_t)node_d_start*(size_t)N:(size_t)node_d_chunksize*(size_t)N] \
            : alloc( TImag[(size_t)node_d_start*(size_t)N:(size_t)node_d_chunksize*(size_t)N]) REUSE) \
            nocopy ( CTFvalue[(size_t)node_d_start*(size_t)N:(size_t)node_d_chunksize*(size_t)N] \
            : alloc( CTFvalue[(size_t)node_d_start*(size_t)N:(size_t)node_d_chunksize*(size_t)N]) REUSE) \
            nocopy ( PHI:length((size_t)M*(size_t)K) REUSE) \
            nocopy ( WReal[(size_t)node_d_start*(size_t)M:(size_t)node_d_chunksize*(size_t)M] \
            : alloc( WReal[(size_t)node_d_start*(size_t)M:(size_t)node_d_chunksize*(size_t)M]) REUSE) \
            nocopy ( WImag[(size_t)node_d_start*(size_t)M:(size_t)node_d_chunksize*(size_t)M] \
            : alloc( WImag[(size_t)node_d_start*(size_t)M:(size_t)node_d_chunksize*(size_t)M]) REUSE) \
            nocopy ( R:length((size_t)N*(size_t)K) REUSE) \
            nocopy ( ipiv:length((size_t)M) REUSE) \
            in (AverageAlpha, AverageBeta, K, M, N)
#else
    #pragma offload if(mic_index>0) target(mic:(mic_index-1)) signal(&WReal[node_d_start]) \
            nocopy ( TReal[(size_t)node_d_start*(size_t)N:(size_t)node_d_chunksize*(size_t)N] \
            : alloc( TReal[(size_t)node_d_start*(size_t)N:(size_t)node_d_chunksize*(size_t)N]) REUSE) \
            nocopy ( TImag[(size_t)node_d_start*(size_t)N:(size_t)node_d_chunksize*(size_t)N] \
            : alloc( TImag[(size_t)node_d_start*(size_t)N:(size_t)node_d_chunksize*(size_t)N]) REUSE) \
            nocopy ( CTFvalue[(size_t)node_d_start*(size_t)N:(size_t)node_d_chunksize*(size_t)N] \
            : alloc( CTFvalue[(size_t)node_d_start*(size_t)N:(size_t)node_d_chunksize*(size_t)N]) REUSE) \
            nocopy ( PHI:length((size_t)M*(size_t)K) REUSE) \
            nocopy ( WReal[(size_t)node_d_start*(size_t)M:(size_t)node_d_chunksize*(size_t)M] \
            : alloc( WReal[(size_t)node_d_start*(size_t)M:(size_t)node_d_chunksize*(size_t)M]) REUSE) \
            nocopy ( WImag[(size_t)node_d_start*(size_t)M:(size_t)node_d_chunksize*(size_t)M] \
            : alloc( WImag[(size_t)node_d_start*(size_t)M:(size_t)node_d_chunksize*(size_t)M]) REUSE) \
            nocopy ( R:length((size_t)N*(size_t)K) REUSE) \
            nocopy ( ipiv:length((size_t)M) REUSE) \
            in (AverageAlpha, AverageBeta, K, M, N)
#endif
#endif
        {
            int numthreads = maxthreads;
            __attribute__((aligned(64)))  double *AA;
            __attribute__((aligned(64)))  double *bRealArray;
            __attribute__((aligned(64)))  double *bImagArray;
            __attribute__((aligned(64)))  double *A;
            __attribute__((aligned(64)))  double *b;

            if (mic_index > 0)
            {
                // Don't use all threads per core to keep memory subsystem happy
                omp_set_num_threads(numXeonPhithreads);//60
                numthreads = numXeonPhithreads;
            }
            else
            {
            }

            const int totalsizeA = numthreads*sizeof(double)*(M*M + 4096);
            const int totalsizeb = numthreads*sizeof(double)*(M   + 4096);

            AA = (double*)aMalloc(totalsizeA,64);
            bRealArray = (double*)aMalloc(totalsizeb,64);
            bImagArray = (double*)aMalloc(totalsizeb,64);

            A = (double*)aMalloc(sizeof(double)*M*M,64);   //for icpc case
            b = (double*)aMalloc(sizeof(double)*2*M,64);

            // Zero arrays in a way that does proper NUMA first-touch
            // We do not use dynamic scheduling to ensure that the
            // memory is evenly spread among the packages in the system
            // Yes, I know we end up clearing the same memory multiple
            // times.
            // WARN : only do first touch when K is large,otherwise some
            // thread data will not be set to zero....
            if(K >= numthreads)
            {
                #pragma omp parallel for
                for(int k = 0;k < K;k++)   // K is only about 1000
                {
                   int tid = omp_get_thread_num();
                   int offsetA = tid*(M*M + 4096);
                   int offsetb = tid*(M   + 4096);
                   for(int m = 0;m < M;m++)   // M is about 600
                   {  
                       bRealArray[offsetb + m] = 0.;
                       bImagArray[offsetb + m] = 0.;

                        for(int m2 = 0;m2 < M;m2++) { 
                            AA[offsetA + m*M + m2] = 0.;
                        } // m2
                    } // m
                } // k
            }
            else
            {
                memset(AA, 0x00, totalsizeA);
                memset(bRealArray, 0x00, totalsizeb);
                memset(bImagArray, 0x00, totalsizeb);
            }
            // Now update the W array from the latest R
            for (int d = node_d_start; d < node_d_end; d++){
                getWd(d,A,b, numthreads, AA, bRealArray, bImagArray);
                for(int m = 0;m < M;m++){
                    ERROR_CHECK(b[2*m]>std::numeric_limits<float>::max(), "Double to Float conversion error..");
                    ERROR_CHECK(b[2*m+1]>std::numeric_limits<float>::max(), "Double to Float conversion error..");
                    WReal[INDEXSZ(d)* INDEXSZ(M)+ INDEXSZ(m)] = b[2*m];
                    WImag[INDEXSZ(d)* INDEXSZ(M)+ INDEXSZ(m)] = b[2*m+1];
                }
            }

            aFree(A);
            aFree(b);
            aFree(AA);
            aFree(bRealArray);
            aFree(bImagArray);
        } // pragma offload

        fstop = dtime();

        NODE0ONLY std::cout<<"iterateUpdateW body done "<<mic_index<<": "<<fstop-fstart<<" sec"<<std::endl<<std::flush;

    STOP_HETERO_OFFLOAD // iter_index
  
#ifndef DONT_USE_OFFLOAD
#ifndef SERIAL_OFFLOAD
    // Wait for async jobs to complete, if they are not already done
    for(int mic_index = 1;mic_index < NUM_LOCAL_COMPUTE_ENGINES;mic_index++)
    {
        int hostamnt = d_interval - estnodesize * numMIC;
        int micstart = hostamnt;
        
        node_d_start = micstart + (mic_index-1) * estnodesize;
#ifdef __INTEL_OFFLOAD
        if (!_Offload_signaled((mic_index-1), (void *)&WReal[node_d_start]))
        {
            #pragma offload_wait target(mic:(mic_index-1)) wait(&WReal[node_d_start])
            ;
        }
#endif
    }
#endif
#endif

#ifdef DEBUGGTM
        NODE0ONLY dumpW();
#endif
}

/// ---------------------------------------------------------------------

// get vector[M] from matrix W[K][D][M]
MIC_OFFLOAD_ATTRIBUTES
void getWd(const int &d,double *A,double *b, int maxthreads,
           double *AA, double *bRealArray, double *bImagArray)
{
    const int MM = M*M;
    
    __assume_aligned(WReal,64);
    __assume_aligned(WImag,64);
    __assume_aligned(CTFvalue,64);
    __assume_aligned(TReal,64);
    __assume_aligned(TImag,64);
    __assume_aligned(PHI,64);
    __assume_aligned(R,64);    
    __assume_aligned(AA,64);
    __assume_aligned(A,64);
    __assume_aligned(b,64);
    __assume_aligned(bRealArray,64);
    __assume_aligned(bImagArray,64);  
  
    const size_t usedCapacity = size_t(maxthreads);
    bool* used = vNew(bool,usedCapacity);
    for (int i = 0; i < usedCapacity; i++) used[i] = false;

	const int Aoffset = (MM + 4096);	// See const int totalsizeA 
	const int boffset = (M  + 4096);	// See const int totalsizeb
    
    // TODO - blocked implementation
    
#pragma omp parallel for  schedule(dynamic)
    for(int k = 0;k < K;k++)   // K is only about 1000
    {
        const int tid = omp_get_thread_num();
        // assert(tid < maxthreads);
        used[tid] = true;

		double pSA = 0;
        double pSb1 = 0;
        double pSb2 = 0;
        double CTFVdNn;
        double RCTFV;

        // Each thread doing unit stride array accesses
        for(int n = 0; n < N; n++)   // N is massive
        {
            // HOT
            CTFVdNn = CTFvalue[INDEXSZ(INDEXSZ(d) * INDEXSZ(N)) + INDEXSZ(n)];
            RCTFV = R[INDEXSZ(INDEXSZ(k) * INDEXSZ(N)) + INDEXSZ(n)]*CTFVdNn;
            pSb1 += RCTFV*TReal[INDEXSZ(INDEXSZ(d) * INDEXSZ(N)) + INDEXSZ(n)]; //partialSumsb1[k] += RCTFV*TReal[d*N+n];//AverageBeta
            pSA += RCTFV*CTFVdNn;
            pSb2 += RCTFV*TImag[INDEXSZ(INDEXSZ(d) * INDEXSZ(N)) + INDEXSZ(n)];
        }

        const int offsetA = tid*Aoffset;
        const int offsetb = tid*boffset;
        
        double PHIM;
        for(int m = 0;m < M;m++)   // M is about 600
        {  // for each A's rows
            PHIM = PHI[k*M+m];
            
            // Sum this product separately for each thread iterating through
            // all k - we'll boil this down at the end
            bRealArray[offsetb + m] += pSb1*PHIM;
            bImagArray[offsetb + m] += pSb2*PHIM;
            
            for(int m2 = 0;m2 < M;m2++)
            { // for each A's col
                // HOT
                // thread-safe due to use of offset
               AA[offsetA + m*M + m2] += pSA*PHIM*PHI[k*M+m2];//AverageBeta
            } // m2
        } // m
    } // k
    
    // Now we reduce the results from AA into A.   Needs to be in a
    // separate loop for thread-safety
#pragma omp parallel for
    for(int m = 0;m < M;m++)
    {
        // Set A
        for(int m2 = 0;m2 < M;m2++)
        { // for each A's col
            double sumReal = 0;

            // maxthreads is relatively small, and each thread is striding by maxthreads
            // in PHI - hopefully won't hurt the caches too much
            // Iterating with k as the outside loop as above would cause a race
            for(int k = 0;k < maxthreads;k++) {
                // HOT
                sumReal += AA[k*Aoffset + m*M + m2];
                AA[k*Aoffset + m*M + m2] = 0;  // Zero for next time around
            }

            A[m*M + m2] = sumReal;
        }  // m2

        // Set b
        // Also reduce the results for b

		double bReal = 0, bImag = 0;
        for(int t = 0;t < maxthreads;t++)
        {
            if (!used[t]) continue;				// don't waste time adding the zeros
            bReal += bRealArray[t*boffset + m];
            bRealArray[t*boffset + m] = 0; // Zero for next time around
        }
        for(int t = 0;t < maxthreads;t++)
        {
            if (!used[t]) continue;				// don't waste time adding the zeros
            bImag += bImagArray[t*boffset + m];
            bImagArray[t*boffset + m] = 0; // Zero for next time around
        }
        b[2*m] = bReal;   
        b[2*m+1] = bImag;
    } // m
    
    for(int m = 0;m < M-1;m++)   // for each A's rows
        A[m*M + m] += AverageAlpha/AverageBeta;

    mkl_solve(A,M,b,2);

}

/// ---------------------------------------------------------------------

void updateY()
{
    double fstart,fstop;
    
#ifdef DEBUGGTM
    // std::cout<<"node="<<node<<" --------------- Update Y1 -----------------"<<std::endl<<std::flush;
    // dumpW();
    // dumpPHI();
    // std::cout<<"node="<<node<<" --------------- Update Y2 -----------------"<<std::endl<<std::flush;
#endif
    
    // Divide up the work between the host and any Xeon Phi cards in the system
    // so that both are working to create the final R at once
    
    START_HETERO_OFFLOAD
    
        fstart = dtime();

        NODE0ONLY std::cout<<"updateY section "<<mic_index<<" node_d_start="<<node_d_start<<", node_d_end="<<node_d_end \
                    <<", node_d_chunksize= "<<node_d_chunksize<<std::endl<<std::flush;
#ifndef DONT_USE_OFFLOAD
#ifdef SERIAL_OFFLOAD
        #pragma offload if(mic_index>0) target(mic:(mic_index-1))  \
            nocopy ( WReal[(size_t)node_d_start*(size_t)M:(size_t)node_d_chunksize*(size_t)M] \
            : alloc( WReal[(size_t)node_d_start*(size_t)M:(size_t)node_d_chunksize*(size_t)M]) REUSE) \
            nocopy ( WImag[(size_t)node_d_start*(size_t)M:(size_t)node_d_chunksize*(size_t)M] \
            : alloc( WImag[(size_t)node_d_start*(size_t)M:(size_t)node_d_chunksize*(size_t)M]) REUSE) \
            nocopy ( PHI:length((size_t)M*(size_t)K) REUSE) \
            nocopy ( YReal[(size_t)node_d_start*(size_t)K:(size_t)node_d_chunksize*(size_t)K] \
            : alloc( YReal[(size_t)node_d_start*(size_t)K:(size_t)node_d_chunksize*(size_t)K]) REUSE) \
            nocopy ( YImag[(size_t)node_d_start*(size_t)K:(size_t)node_d_chunksize*(size_t)K] \
            : alloc( YImag[(size_t)node_d_start*(size_t)K:(size_t)node_d_chunksize*(size_t)K]) REUSE)
#else
        #pragma offload if(mic_index>0) target(mic:(mic_index-1)) signal(&YReal[node_d_start]) \
            nocopy ( WReal[(size_t)node_d_start*(size_t)M:(size_t)node_d_chunksize*(size_t)M] \
            : alloc( WReal[(size_t)node_d_start*(size_t)M:(size_t)node_d_chunksize*(size_t)M]) REUSE) \
            nocopy ( WImag[(size_t)node_d_start*(size_t)M:(size_t)node_d_chunksize*(size_t)M] \
            : alloc( WImag[(size_t)node_d_start*(size_t)M:(size_t)node_d_chunksize*(size_t)M]) REUSE) \
            nocopy ( PHI:length((size_t)M*(size_t)K) REUSE) \
            nocopy ( YReal[(size_t)node_d_start*(size_t)K:(size_t)node_d_chunksize*(size_t)K] \
            : alloc( YReal[(size_t)node_d_start*(size_t)K:(size_t)node_d_chunksize*(size_t)K]) REUSE) \
            nocopy ( YImag[(size_t)node_d_start*(size_t)K:(size_t)node_d_chunksize*(size_t)K] \
            : alloc( YImag[(size_t)node_d_start*(size_t)K:(size_t)node_d_chunksize*(size_t)K]) REUSE)
#endif
#endif
        {
            if (mic_index > 0)
            {
                // Don't use all threads per core to keep memory subsystem happy
                omp_set_num_threads(numXeonPhithreads);//60
            }
            else
            {
            }

    #pragma omp parallel for
            for(int d = node_d_start;d < node_d_end;d++)
            {
               double PHIWSUMReal,PHIWSUMImag;

               for(int k = 0;k < K;k++)
               {

                  PHIWSUMReal = PHIWSUMImag = 0.;

                  for(int m = 0; m < M; m++)
                  {
                     PHIWSUMReal += PHI[k*M+m]*WReal[INDEXSZ(INDEXSZ(d) * INDEXSZ(M)) + INDEXSZ(m)];
                     PHIWSUMImag += PHI[k*M+m]*WImag[INDEXSZ(INDEXSZ(d) * INDEXSZ(M)) + INDEXSZ(m)];
                  } // m

                  YReal[INDEXSZ(INDEXSZ(d) * INDEXSZ(K)) + INDEXSZ(k)] = PHIWSUMReal;
                  YImag[INDEXSZ(INDEXSZ(d) * INDEXSZ(K)) + INDEXSZ(k)] = PHIWSUMImag;

                } // k
                //        if (d == d_interval-2) std::cout<<"node :"<<node<<" PHIWSUMReal = "<<PHIWSUMReal<<" PHIWSUMImag = "<<PHIWSUMImag<<std::endl<<std::flush;
            } // d

        } // pragma offload

        fstop = dtime();

        NODE0ONLY std::cout<<"updateY body done "<<mic_index<<": "<<fstop-fstart<<" sec"<<std::endl<<std::flush;

    
    STOP_HETERO_OFFLOAD // iter_index
    
#ifndef SERIAL_OFFLOAD
    // Wait for async jobs to complete, if they are not already done
    for(int mic_index = 1;mic_index < NUM_LOCAL_COMPUTE_ENGINES;mic_index++)
    {
        int hostamnt = d_interval - estnodesize * numMIC;
        int micstart = hostamnt;
        
        node_d_start = micstart + (mic_index-1) * estnodesize;
#ifdef __INTEL_OFFLOAD
        if (!_Offload_signaled((mic_index-1), (void *)&YReal[node_d_start]))
        {
#pragma offload_wait target(mic:(mic_index-1)) wait(&YReal[node_d_start])
            ;
        }
#endif
    }
#endif
    
    // Since W and Y are only reference their own subset of d, we don't
    // need to coordinate these results to other nodes  
#ifdef DEBUGGTM
    NODE0ONLY dumpY();
#endif
}

/// ---------------------------------------------------------------------

// iterated update Beta,corresponding to formula (1.5)
double updateBeta()
{
    double tempAverageBeta = 0.0;
    double fstart,fstop;
    double nodeTempAverageBeta[NUM_LOCAL_COMPUTE_ENGINES];
    
    NODE0ONLY std::cout << "updateBeta: N = " <<N<< ", D = " <<D<< ", K = " <<K<< ", M = " <<M<<std::endl<<std::flush;
    
    
    // Divide up the work between the host and any Xeon Phi cards in the system
    // so that both are working to create the final R at once
    
    START_HETERO_OFFLOAD
    
        NODE0ONLY std::cout<<"updateBeta section "<<mic_index<<" node_d_start="<<node_d_start<<", node_d_end="<<node_d_end \
                            <<", node_d_chunksize= "<<node_d_chunksize<<std::endl<<std::flush;

        fstart = dtime();
#ifndef DONT_USE_OFFLOAD
#ifdef SERIAL_OFFLOAD
        #pragma offload if(mic_index>0) target(mic:(mic_index-1)) \
            nocopy ( TReal[(size_t)node_d_start*(size_t)N:(size_t)node_d_chunksize*(size_t)N] \
            : alloc( TReal[(size_t)node_d_start*(size_t)N:(size_t)node_d_chunksize*(size_t)N]) REUSE) \
            nocopy ( TImag[(size_t)node_d_start*(size_t)N:(size_t)node_d_chunksize*(size_t)N] \
            : alloc( TImag[(size_t)node_d_start*(size_t)N:(size_t)node_d_chunksize*(size_t)N]) REUSE) \
            nocopy ( CTFvalue[(size_t)node_d_start*(size_t)N:(size_t)node_d_chunksize*(size_t)N] \
            : alloc( CTFvalue[(size_t)node_d_start*(size_t)N:(size_t)node_d_chunksize*(size_t)N]) REUSE) \
            nocopy ( YReal[(size_t)node_d_start*(size_t)K:(size_t)node_d_chunksize*(size_t)K] \
            : alloc( YReal[(size_t)node_d_start*(size_t)K:(size_t)node_d_chunksize*(size_t)K]) REUSE) \
            nocopy ( YImag[(size_t)node_d_start*(size_t)K:(size_t)node_d_chunksize*(size_t)K] \
            : alloc( YImag[(size_t)node_d_start*(size_t)K:(size_t)node_d_chunksize*(size_t)K]) REUSE) \
            nocopy ( R:length((size_t)N*(size_t)K) REUSE) \
            out (nodeTempAverageBeta[mic_index:1])
#else
        #pragma offload if(mic_index>0) target(mic:(mic_index-1)) signal(&nodeTempAverageBeta[mic_index]) \
            nocopy ( TReal[(size_t)node_d_start*(size_t)N:(size_t)node_d_chunksize*(size_t)N] \
            : alloc( TReal[(size_t)node_d_start*(size_t)N:(size_t)node_d_chunksize*(size_t)N]) REUSE) \
            nocopy ( TImag[(size_t)node_d_start*(size_t)N:(size_t)node_d_chunksize*(size_t)N] \
            : alloc( TImag[(size_t)node_d_start*(size_t)N:(size_t)node_d_chunksize*(size_t)N]) REUSE) \
            nocopy ( CTFvalue[(size_t)node_d_start*(size_t)N:(size_t)node_d_chunksize*(size_t)N] \
            : alloc( CTFvalue[(size_t)node_d_start*(size_t)N:(size_t)node_d_chunksize*(size_t)N]) REUSE) \
            nocopy ( YReal[(size_t)node_d_start*(size_t)K:(size_t)node_d_chunksize*(size_t)K] \
            : alloc( YReal[(size_t)node_d_start*(size_t)K:(size_t)node_d_chunksize*(size_t)K]) REUSE) \
            nocopy ( YImag[(size_t)node_d_start*(size_t)K:(size_t)node_d_chunksize*(size_t)K] \
            : alloc( YImag[(size_t)node_d_start*(size_t)K:(size_t)node_d_chunksize*(size_t)K]) REUSE) \
            nocopy ( R:length((size_t)N*(size_t)K) REUSE) \
            out (nodeTempAverageBeta[mic_index:1])
#endif
#endif
            {
                int D_BLOCK;
                int K_BLOCK;
                int N_BLOCK;

                if (mic_index > 0)
                {
                    // Don't use all threads per core to keep memory subsystem happy
                    omp_set_num_threads(numXeonPhithreads);//60

                    // Experimentally determined best block factors for Xeon Phi 60 threads
                    D_BLOCK = 64;
                    K_BLOCK = 64;
                    N_BLOCK = 1024;
                }
                else
                {
                    // Experimentally determined best block factors for Xeon 24 cores
                    D_BLOCK=128;
                    K_BLOCK=128;
                    N_BLOCK=512;
                }

            double tempNodeAverageBeta = 0.0;

            #pragma omp parallel for collapse(3) schedule(dynamic) reduction(+:tempNodeAverageBeta)
            for (int dd = node_d_start; dd < node_d_end; dd += D_BLOCK)
            {
            for (int kk = 0; kk < K; kk += K_BLOCK)
            {
            for (int nn = 0; nn < N; nn += N_BLOCK)
            {
               int k_end = kk+K_BLOCK;  if (k_end > K) k_end=K; 
               int dd_end = dd+D_BLOCK;  if (dd_end > node_d_end) dd_end=node_d_end; 
               int n_end = nn+N_BLOCK;  if (n_end > N) n_end = N; 

               double dsum = 0.0;
               double delta;
               double YRealk, YImagk;
               double CTFVn;

               for(int d = dd; d < dd_end; d++)
               {
               for(int k = kk; k < k_end; k++)
               {
                  YRealk = YReal[INDEXSZ(d) * INDEXSZ(K) + INDEXSZ(k)];
                  YImagk = YImag[INDEXSZ(d) * INDEXSZ(K) + INDEXSZ(k)];
                  for (int n = nn; n < n_end; n++)
                  {
                     CTFVn = CTFvalue[INDEXSZ(INDEXSZ(d)*INDEXSZ(N))+ INDEXSZ(n)];
                     double a = TReal[INDEXSZ(INDEXSZ(d)*INDEXSZ(N))+ INDEXSZ(n)]-CTFVn*YRealk;
                     double b = TImag[INDEXSZ(INDEXSZ(d)*INDEXSZ(N))+ INDEXSZ(n)]-CTFVn*YImagk;
                     delta = a*a+b*b;
                     dsum += R[INDEXSZ(INDEXSZ(k)*INDEXSZ(N))+ INDEXSZ(n)]*delta; 
                  }  // n
               } } // k & d
               tempNodeAverageBeta += dsum;
            } } }  // pragma omp parallel


                nodeTempAverageBeta[mic_index] = tempNodeAverageBeta; 

            } // pragma offload


        fstop = dtime();
        NODE0ONLY std::cout<<"updateBeta body done "<<mic_index<<": "<<fstop-fstart \
                            <<" sec"<<std::endl<<std::flush;
    
    STOP_HETERO_OFFLOAD // iter_index
    
#ifndef SERIAL_OFFLOAD
    // Wait for async jobs to complete, if they are not already done
    for(int mic_index = 1;mic_index < NUM_LOCAL_COMPUTE_ENGINES;mic_index++)
    {
#ifdef __INTEL_OFFLOAD
       if (!_Offload_signaled((mic_index-1), (void *)&nodeTempAverageBeta[mic_index]))
        {
            #pragma offload_wait target(mic:(mic_index-1)) wait(&nodeTempAverageBeta[mic_index])
            ;
        }
#endif
    }
#endif
    
    // Calculate the final beta for this node
    tempAverageBeta = 0;
    for (int i=0; i < NUM_LOCAL_COMPUTE_ENGINES; i++)
        tempAverageBeta += nodeTempAverageBeta[i];
    
    // AverageBeta  = N*D/tempAverageBeta;
    
#ifdef DEBUGGTM
     NODE0ONLY std::cout<<"iterate update Beta!"<<std::endl<<std::flush;
    // prtMatrix(Beta,N,D);
     NODE0ONLY std::cout<<"tempAverageBeta = "<<tempAverageBeta<<std::endl<<std::flush;
#endif
    
    // return N*D/tempAverageBeta;
    return tempAverageBeta;
}

/// ---------------------------------------------------------------------

double updateAlpha(){

    double tempAverageAlpha = 0.0;
    double fstart,fstop; 
    double nodeTempAverageAlpha[NUM_LOCAL_COMPUTE_ENGINES];

    // Divide up the work between the host and any Xeon Phi cards in the system
    // so that both are working to create the final R at once

   START_HETERO_OFFLOAD
    
        fstart = dtime(); 
        
        NODE0ONLY std::cout<<"updateAlpha section "<<mic_index<<" node_d_start="<<node_d_start<<", node_d_end="<<node_d_end<<", node_d_chunksize= "<<node_d_chunksize<<std::endl<<std::flush;

        // Since this runs really quickly, we will do this processing synchronously so that we can
        // avoid fiddling with synchronization variables
#ifndef DONT_USE_OFFLOAD
		    #pragma offload if(mic_index>0) target(mic:(mic_index-1)) \
            nocopy ( WReal[(size_t)node_d_start*(size_t)M:(size_t)node_d_chunksize*(size_t)M] \
            : alloc( WReal[(size_t)node_d_start*(size_t)M:(size_t)node_d_chunksize*(size_t)M]) REUSE) \
            nocopy ( WImag[(size_t)node_d_start*(size_t)M:(size_t)node_d_chunksize*(size_t)M] \
            : alloc( WImag[(size_t)node_d_start*(size_t)M:(size_t)node_d_chunksize*(size_t)M]) REUSE) \
            out (nodeTempAverageAlpha[mic_index:1])
#endif
        {   
            if (mic_index > 0)
            {
               // Don't use all threads per core to keep memory subsystem happy
               omp_set_num_threads(numXeonPhithreads);
            }
            else
            {
            }

            double tempNodeAverageAlpha = 0;

           #pragma omp parallel for reduction(+:tempNodeAverageAlpha)
               for(int d = node_d_start;d < node_d_end;d++) 
               {
                   double dsum = 0;
                   for(int m = 0;m < M;m++) 
                       dsum += WReal[INDEXSZ(INDEXSZ(d)*INDEXSZ(M)) + INDEXSZ(m)] * 
                               WReal[INDEXSZ(INDEXSZ(d)*INDEXSZ(M)) + INDEXSZ(m)] +
                               WImag[INDEXSZ(INDEXSZ(d)*INDEXSZ(M)) + INDEXSZ(m)] *
                               WImag[INDEXSZ(INDEXSZ(d)*INDEXSZ(M)) + INDEXSZ(m)];   //pow(abs(W[k*M*D+m*D+d]),2);
                   tempNodeAverageAlpha += dsum;
               }

               nodeTempAverageAlpha[mic_index] = tempNodeAverageAlpha;
         } // pragma offload
        
        fstop = dtime();

        NODE0ONLY std::cout<<"updateAlpha body done "<<mic_index<<": "<<fstop-fstart<<" sec"<<std::endl<<std::flush; 

   STOP_HETERO_OFFLOAD // iter_index

   tempAverageAlpha = 0;
    for (int i=0; i < NUM_LOCAL_COMPUTE_ENGINES; i++)
            tempAverageAlpha += nodeTempAverageAlpha[i];

    // return 1./tempAverageAlpha;
    return tempAverageAlpha;
}

    
double update_P_Delta(){
    
    double tempDeviation = 0.0;
    double fstart,fstop;
    double nodeTempDeviation[NUM_LOCAL_COMPUTE_ENGINES];
    
    NODE0ONLY std::cout << "update_P_Delta: N = " <<N<< ", D = " <<D<< ", K = " <<K<< ", M = " <<M<<std::endl<<std::flush;
    // std::cout << "update_P_Delta: d_interval/2 = " <<d_interval/2<<std::endl<<std::flush;
    
    // Divide up the work between the host and any Xeon Phi cards in the system
    // so that both are working to create the final R at once
    
   START_HETERO_OFFLOAD  
        
        fstart = dtime();
        NODE0ONLY std::cout<<"update_P_Delta section "<<mic_index<<" node_d_start="<<node_d_start \
                           <<", node_d_end="<<node_d_end<<", node_d_chunksize= "<<node_d_chunksize<<std::endl<<std::flush;
#ifndef DONT_USE_OFFLOAD       
#ifdef SERIAL_OFFLOAD
        #pragma offload if(mic_index>0) target(mic:(mic_index-1)) \
            nocopy ( TReal[(size_t)node_d_start*(size_t)N:(size_t)node_d_chunksize*(size_t)N] \
            : alloc( TReal[(size_t)node_d_start*(size_t)N:(size_t)node_d_chunksize*(size_t)N]) REUSE) \
            nocopy ( TImag[(size_t)node_d_start*(size_t)N:(size_t)node_d_chunksize*(size_t)N] \
            : alloc( TImag[(size_t)node_d_start*(size_t)N:(size_t)node_d_chunksize*(size_t)N]) REUSE) \
            nocopy ( CTFvalue[(size_t)node_d_start*(size_t)N:(size_t)node_d_chunksize*(size_t)N] \
            : alloc( CTFvalue[(size_t)node_d_start*(size_t)N:(size_t)node_d_chunksize*(size_t)N]) REUSE) \
            nocopy ( YReal[(size_t)node_d_start*(size_t)K:(size_t)node_d_chunksize*(size_t)K] \
            : alloc( YReal[(size_t)node_d_start*(size_t)K:(size_t)node_d_chunksize*(size_t)K]) REUSE) \
            nocopy ( YImag[(size_t)node_d_start*(size_t)K:(size_t)node_d_chunksize*(size_t)K] \
            : alloc( YImag[(size_t)node_d_start*(size_t)K:(size_t)node_d_chunksize*(size_t)K]) REUSE) \
            nocopy ( R:length((size_t)N*(size_t)K) REUSE) \
            out (nodeTempDeviation[mic_index:1])
#else
        #pragma offload if(mic_index>0) target(mic:(mic_index-1)) signal(&nodeTempDeviation[mic_index]) \
            nocopy ( TReal[(size_t)node_d_start*(size_t)N:(size_t)node_d_chunksize*(size_t)N] \
            : alloc( TReal[(size_t)node_d_start*(size_t)N:(size_t)node_d_chunksize*(size_t)N]) REUSE) \
            nocopy ( TImag[(size_t)node_d_start*(size_t)N:(size_t)node_d_chunksize*(size_t)N] \
            : alloc( TImag[(size_t)node_d_start*(size_t)N:(size_t)node_d_chunksize*(size_t)N]) REUSE) \
            nocopy ( CTFvalue[(size_t)node_d_start*(size_t)N:(size_t)node_d_chunksize*(size_t)N] \
            : alloc( CTFvalue[(size_t)node_d_start*(size_t)N:(size_t)node_d_chunksize*(size_t)N]) REUSE) \
            nocopy ( YReal[(size_t)node_d_start*(size_t)K:(size_t)node_d_chunksize*(size_t)K] \
            : alloc( YReal[(size_t)node_d_start*(size_t)K:(size_t)node_d_chunksize*(size_t)K]) REUSE) \
            nocopy ( YImag[(size_t)node_d_start*(size_t)K:(size_t)node_d_chunksize*(size_t)K] \
            : alloc( YImag[(size_t)node_d_start*(size_t)K:(size_t)node_d_chunksize*(size_t)K]) REUSE) \
            nocopy ( R:length((size_t)N*(size_t)K) REUSE) \
            out (nodeTempDeviation[mic_index:1])
#endif
#endif
		{   
            int K_BLOCK;
            int N_BLOCK;
            int D_BLOCK;

            
            if (mic_index > 0)
            {
                // Don't use all threads per core to keep memory subsystem happy
                omp_set_num_threads(120); 
                // TODO - Different number of threads may give different/better
                // performance with different blocking factors

               // Experimentally determined best block factors for Xeon Phi 120 threads
               // May not best for other numbers
               K_BLOCK = 64;
               N_BLOCK = 32; 
            }
            else
            {
               // Experimentally determined best block factors for Xeon 24 cores
               K_BLOCK = 64;
               N_BLOCK = 4096;
               // For KNL 144 threads, K_BLOCK=8, N_BLOCK=512
            }

            // Can we do this in a page-blocked memory clear?
            #pragma omp parallel for
            for(size_t i = 0;i < INDEXSZ(INDEXSZ(K)*INDEXSZ(N));i++)
                R[i] = 0.;       

            double tempNodeDeviation = 0.0;
                
#pragma omp parallel for collapse(2) schedule(dynamic) reduction(+:tempNodeDeviation)
            for (int kk = 0; kk < K; kk += K_BLOCK)
            {
            for (int nn = 0; nn < N; nn += N_BLOCK)
            {
               int k_end = kk+K_BLOCK;  if (k_end > K) k_end=K; 
               int n_end = nn+N_BLOCK;  if (n_end > N) n_end = N; 

               double dsum = 0.0;
               double delta;

               for(int d = node_d_start; d < node_d_end; d++)
               {   
                  for (int k = kk; k < k_end; k++)
                  {                
                     double RealdKk = YReal[INDEXSZ(INDEXSZ(d)* INDEXSZ(K))+ INDEXSZ(k)];
                     double ImagdKk = YImag[INDEXSZ(INDEXSZ(d)* INDEXSZ(K))+ INDEXSZ(k)];

                     for (int n = nn; n < n_end; n++) 
                     {
                         double a = TReal[INDEXSZ(INDEXSZ(d) * INDEXSZ(N)) + INDEXSZ(n)] - 
                                    CTFvalue[INDEXSZ(INDEXSZ(d) * INDEXSZ(N)) + INDEXSZ(n)]*RealdKk;
                         double b = TImag[INDEXSZ(INDEXSZ(d) * INDEXSZ(N)) + INDEXSZ(n)]-
                                    CTFvalue[INDEXSZ(INDEXSZ(d)* INDEXSZ(N)) + INDEXSZ(n)]*ImagdKk;
                         delta = a*a+b*b;

                         R[INDEXSZ(INDEXSZ(k) * INDEXSZ(N)) + INDEXSZ(n)] += delta; 
                                   // Watch out for the race condition here.  Any parallelism
                                   // by d will result in multiple threads updating the same
                                   // R[k*N+n] elements at the same time.
                         dsum += delta; 
                     }  // for n
                  } // for k
               } // for d 
               tempNodeDeviation += dsum;
            } }  // pragma omp parallel
 
           nodeTempDeviation[mic_index] = tempNodeDeviation; 

        } // pragma offload
        fstop = dtime();
        
        NODE0ONLY std::cout<<"update_P_Delta body done "<<mic_index<<": "<<fstop-fstart<<" sec"<<std::endl<<std::flush; 
        
    STOP_HETERO_OFFLOAD // iter_index
    
#ifndef SERIAL_OFFLOAD
    // Wait for async jobs to complete, if they are not already done
    for(int mic_index = 1;mic_index < NUM_LOCAL_COMPUTE_ENGINES;mic_index++)
    {
#ifdef __INTEL_OFFLOAD
       if (!_Offload_signaled((mic_index-1), (void *)&nodeTempDeviation[mic_index]))
        {
            #pragma offload_wait target(mic:(mic_index-1)) wait(&nodeTempDeviation[mic_index])
            ;
        }
#endif
    } 
#endif
    
    // Calculate the final deviation for this node
    tempDeviation = 0;
    for (int i=0; i < NUM_LOCAL_COMPUTE_ENGINES; i++)
            tempDeviation += nodeTempDeviation[i];

    // Combine the R results from all NUM_LOCAL_COMPUTE_ENGINES coprocessors (who now
    // all have their own partial values of R) onto the host into its partial R,
    // thereby building the final R for this node.
    // NOTE - We do this the on host to avoid extra copies of data to the Phi
    if (numMIC > 0)
    {
      double *tempRK = (double*)aMalloc(sizeof(double)*N,64);
      for(int mic_index = 1; mic_index < NUM_LOCAL_COMPUTE_ENGINES; mic_index++)
      {
          size_t chunkStart, chunkSize;
          
          fstart = dtime();
          
          chunkSize = N;
          
          for(int k = 0; k < K; k++)
          {
              // Get a chunk from the MIC card
              chunkStart = k*N;
#ifdef DONT_USE_OFFLOAD
			  std::cerr << "Not using offload" << std::endl;
			  EXIT_ABNORMALLY;
#else
              #pragma offload_transfer target(mic:(mic_index-1)) \
                  out(R[chunkStart:chunkSize] : into(tempRK[0:chunkSize]) REUSE)
#endif
              // Fold it into the host's R
              for (int n = 0; n < N; n++)
                  R[INDEXSZ(INDEXSZ(k) * INDEXSZ(N)) + INDEXSZ(n)] += tempRK[n];
          } // for all K
          
          fstop = dtime();
          NODE0ONLY std::cout<<"Update_P_Delta R gather "<<mic_index<<": "<<fstop-fstart<<" sec"<<std::endl<<std::flush; 

      } // for all cards
      aFree(tempRK);
    }

    // Do the final set up of this node's R
#pragma omp parallel for
    for (int k = 0; k < K; k++)
    {
        for (int n = 0; n < N; n++) 
        {
           R[INDEXSZ(INDEXSZ(k) * INDEXSZ(N)) + INDEXSZ(n)] *= AverageBeta * 0.5;
        }
    } // for k
    
    // NOTE - we do not send R down to the cards here because we need to
    // wait for the reduction between host nodes over MPI first
    
#ifdef DEBUGGTM
    NODE0ONLY std::cout<<"node="<<node<<" -- update_P_Delta R"<<std::endl<<std::flush;
    NODE0ONLY dumpR();
#endif
    
    
#ifdef DEBUGGTM
    // std::cout<<"R before shift"<<std::endl<<std::flush;
    // prtMatrix(R,K,N);
    NODE0ONLY std::cout<<"Delta before shift(10*9)"<<std::endl<<std::flush;
    // prtMatrix(Delta,K,N);
#endif
    
    NODE0ONLY std::cout<<"node="<<node<<" tempDeviation = "<<tempDeviation<<std::endl<<std::flush;
    return tempDeviation;
}
    
/// ---------------------------------------------------------------------
   
// expectation-maximization algorithm
void run(double preci,double probThreshold,double alpha,bool update_beta)
{
	double tstart,tstop,Rtime = 0.,Wtime = 0.,Ytime = 0.,Btime = 0.,Atime = 0.,Dtime = 0.;
	double tempBeta, tempAlpha, tempDeviation;
    double fstart, fstop;
    double sectionStart, sectionEnd;

#ifdef DEBUGGTM
	NODE0ONLY{
        std::cout<<"split info:"<<d_start<<" "<<d_end<<std::endl;
		std::cout<<"###############starting debug!####################"<<std::endl;
		std::cout<<"initialize"<<std::endl;
	}
#endif
	
    iter = 1;
    NODE0ONLY showProgressBar(0, nr_iter);
    sectionStart = dtime();

	// initialize average beta
	// AverageBeta = 1;
    // randRange(0,1);
    AverageAlpha = alpha;//0.01;

    // The key thing to understand here is that each node will will have
    // all arrays sized by d_interval.  So an NxD array is now an Nxd_interval
    // array, whose index starts at 0, not d_start.  
    // Remember, the allocated array has to fit within the node, so it cannot
    // be the size of the overall data set, just a subset.
    // The only exception is R, which looks like it needs to fit within the node

    // --- Transfer persistent data to the Xeon Phi coprocessor(s) if present/enabled ---
    if (numMIC > 0)
    {
        int hostamnt;
        int micstart;
        double ostart, ostop;
        
        hostamnt = d_interval - estnodesize * numMIC;

        ostart = dtime();
        for (int mic_index = 1; mic_index < NUM_LOCAL_COMPUTE_ENGINES; mic_index++)
        {
            int node_d_start, node_d_end, node_d_chunksize;
    
            micstart = hostamnt; 
            
            fstart = dtime();
            
            node_d_start = micstart + (mic_index-1) * estnodesize;
            
            node_d_end = ((node_d_start + estnodesize) < d_interval) ? 
                        (node_d_start + estnodesize): d_interval;
            
            node_d_chunksize = node_d_end - node_d_start;
            
            size_t xfersize = node_d_chunksize*N*sizeof(float)*3 +
                              node_d_chunksize*K*sizeof(double)*2 +
                              node_d_chunksize*M*sizeof(float)*2 +
                              K*M*sizeof(double) +
                              N*K*sizeof(double) +
                              M*sizeof(int);
            
            // We're not going to do async data transfer so that we can get
            // maximum PCI bandwidth for each transfer, rather than sharing
            // the bus with multiple transfers
            NODE0ONLY  std::cout<<"Initial transfer to Phi card "<<mic_index<<": size "<<xfersize<<" bytes,  node_d_start="<<node_d_start
                                <<", node_d_end="<<node_d_end<<", node_d_chunksize= "<<node_d_chunksize<<std::endl<<std::flush;
            
#ifndef DONT_USE_OFFLOAD
			#pragma offload_transfer target(mic:(mic_index-1)) \
            in(      TReal[(size_t)node_d_start*(size_t)N:(size_t)node_d_chunksize*(size_t)N]   \
            : alloc( TReal[(size_t)node_d_start*(size_t)N:(size_t)node_d_chunksize*(size_t)N]) ALLOC) \
            in(      TImag[(size_t)node_d_start*(size_t)N:(size_t)node_d_chunksize*(size_t)N]   \
            : alloc( TImag[(size_t)node_d_start*(size_t)N:(size_t)node_d_chunksize*(size_t)N]) ALLOC) \
            in(      CTFvalue[(size_t)node_d_start*(size_t)N:(size_t)node_d_chunksize*(size_t)N] \
            : alloc( CTFvalue[(size_t)node_d_start*(size_t)N:(size_t)node_d_chunksize*(size_t)N]) ALLOC) \
            nocopy(  YReal[(size_t)node_d_start*(size_t)K:(size_t)node_d_chunksize*(size_t)K] \
            : alloc( YReal[(size_t)node_d_start*(size_t)K:(size_t)node_d_chunksize*(size_t)K]) ALLOC) \
            nocopy(  YImag[(size_t)node_d_start*(size_t)K:(size_t)node_d_chunksize*(size_t)K] \
            : alloc( YImag[(size_t)node_d_start*(size_t)K:(size_t)node_d_chunksize*(size_t)K]) ALLOC) \
            in(      WReal[(size_t)node_d_start*(size_t)M:(size_t)node_d_chunksize*(size_t)M]   \
            : alloc( WReal[(size_t)node_d_start*(size_t)M:(size_t)node_d_chunksize*(size_t)M]) ALLOC) \
            in(      WImag[(size_t)node_d_start*(size_t)M:(size_t)node_d_chunksize*(size_t)M]   \
            : alloc( WImag[(size_t)node_d_start*(size_t)M:(size_t)node_d_chunksize*(size_t)M]) ALLOC) \
            in( PHI:length((size_t)M*(size_t)K) ALLOC) \
            nocopy( R:length((size_t)N*(size_t)K) ALLOC) \
            nocopy( ipiv:length((size_t)M) ALLOC)
#endif

            fstop = dtime();
            NODE0ONLY std::cout<<"Initial transfer to Phi card "<<mic_index<<": "<<fstop-fstart<<" sec"<<std::endl<<std::flush;
        }
        ostop = dtime();
        NODE0ONLY std::cout<<"Total transfer to Phi card: "<<ostop-ostart<<" sec"<<std::endl<<std::flush; 
        NODE0ONLY sprintf(buffer,",%f",ostop-ostart);
        NODE0ONLY strcat(timerstring, buffer); 
        
        
        // At this point our portion of TReal, TImag, CTFvalue, YReal, and YImag are
        // transferred to the card.  Space has been reserved for R on the card, but
        // since one of our first acts is to set it up in update_P_Delta, there
        // is no point in transferring it (since it is uninitialized at this point)
    }
    else
    {
        NODE0ONLY sprintf(buffer,",0.0");
        NODE0ONLY strcat(timerstring, buffer);  
    }

    // --- updateY ----
    
    // Set up YReal and YImag for the first time - requires WReal/WImag and PHI
    NODE0ONLY fstart = dtime();
    updateY();
    NODE0ONLY fstop = dtime();
    NODE0ONLY std::cout<<"Initial updateY: "<<fstop-fstart<<" sec"<<std::endl<<std::flush;
    NODE0ONLY sprintf(buffer,",%f",fstop-fstart);
    NODE0ONLY strcat(timerstring, buffer);
  
    // --- update_P_Delta ----

    NODE0ONLY fstart = dtime();
	update_P_Delta();
    NODE0ONLY fstop = dtime();
    NODE0ONLY std::cout<<"Initial update_P_Delta: "<<fstop-fstart<<" sec"<<std::endl<<std::flush;
    NODE0ONLY sprintf(buffer,",%f",fstop-fstart);
    NODE0ONLY strcat(timerstring, buffer);
    NODE0ONLY fstart = dtime();
    
#ifdef USEMPI
	// reduce R to node 0
	double *tempRK = (double*)aMalloc(sizeof(double)*N,64);
	for(int k = 0;k < K;k++){
        MPI::COMM_WORLD.Reduce(R+INDEXSZ(INDEXSZ(k)*INDEXSZ(N)),tempRK,N,MPI::DOUBLE,MPI::SUM,0);
        MPI::COMM_WORLD.Barrier();//???
        NODE0ONLY memcpy(R+INDEXSZ(INDEXSZ(k)*INDEXSZ(N)),tempRK,sizeof(double)*N); /****node 0****/
	}
	aFree(tempRK);
#endif
    NODE0ONLY fstop = dtime();
    NODE0ONLY std::cout<<"Initial R reduction: "<<fstop-fstart<<" sec"<<std::endl<<std::flush;
    NODE0ONLY sprintf(buffer,",%f",fstop-fstart);
    NODE0ONLY strcat(timerstring, buffer);
    sectionEnd = dtime();
    NODE0ONLY sprintf(buffer,",%f",sectionEnd-sectionStart);
    NODE0ONLY strcat(timerstring, buffer);
    
#ifdef DEBUGGTM
	NODE0ONLY std::cout<<"initialize completed."<<std::endl;
#endif

#ifdef DEBUGGTM
	NODE0ONLY std::cout<<"###############cycle##############################"<<std::endl;
#endif

#define DEBUGGTMTIME

    // ----------- Start Iterations -------------
	do{
        sectionStart = dtime();
#if 1
	NODE0ONLY std::cout<<"#####   the "<<iter<<" round.#######"<<std::endl;
#endif

		///////////    E-step    /////////////////////////////////////////////
#ifdef DEBUGGTM
		tstart = dtime();
#endif

		NODE0ONLY{
            fstart = dtime();
			updateR();
            fstop = dtime();
            std::cout<<"updateR: "<<fstop-fstart<<" sec"<<std::endl<<std::flush;
            sprintf(buffer,",%f",fstop-fstart);
            strcat(timerstring, buffer);  
		}
      
        NODE0ONLY fstart = dtime();
        // Broadcast the current R to all nodes
#ifdef USEMPI
        // reduce R to node 0
        for(int k = 0;k < K;k++){
            MPI::COMM_WORLD.Bcast(R+INDEXSZ(INDEXSZ(k)*INDEXSZ(N)),N,MPI::DOUBLE,0);
        }
        // R > 2GB,crash
		// MPI::COMM_WORLD.Bcast(R,K*N,MPI::DOUBLE,0);
#endif
        NODE0ONLY fstop = dtime();
        NODE0ONLY std::cout<<"R broadcast: "<<fstop-fstart<<" sec"<<std::endl<<std::flush;
        NODE0ONLY sprintf(buffer,",%f",fstop-fstart);
        NODE0ONLY strcat(timerstring, buffer);  
        
        NODE0ONLY fstart = dtime();
        
#ifndef DONT_USE_OFFLOAD
		if (numMIC > 0)
        {
            // Now transfer R to all Xeon Phi cards attached to the host
            for (int mic_index = 1; mic_index < NUM_LOCAL_COMPUTE_ENGINES; mic_index++)
            {
#pragma offload_transfer target(mic:(mic_index-1)) \
                in( R:length((size_t)N*(size_t)K) REUSE)
            }
        } // numMIC
#endif

        NODE0ONLY fstop = dtime();
        NODE0ONLY std::cout<<"R transfer to coprocessors: "<<fstop-fstart<<" sec"<<std::endl<<std::flush;
        NODE0ONLY sprintf(buffer,",%f",fstop-fstart);
        NODE0ONLY strcat(timerstring, buffer);        
       
        // --- iterateUpdateW ----
        
#ifdef DEBUGGTM
        tstop = dtime();
        Rtime += (tstop - tstart);
        tstart = dtime();
#endif
        // Updates WReal/WImag from changed R
        // Does not need YReal/YImag
        NODE0ONLY fstart = dtime();
		iterateUpdateW();
        NODE0ONLY fstop = dtime();
        NODE0ONLY std::cout<<" ** iterateUpdateW: "<<fstop-fstart<<" sec"<<std::endl<<std::flush;
        NODE0ONLY sprintf(buffer,",%f",fstop-fstart);
        NODE0ONLY strcat(timerstring, buffer);

#ifdef DEBUGGTM
        tstop = dtime();
        Wtime += (tstop - tstart);
        tstart = dtime();
#endif        
        // --- updateY ----

        // Update YReal and YImag now that WReal and WImag have changed        
        NODE0ONLY fstart = dtime();

		updateY();
        NODE0ONLY fstop = dtime();
        NODE0ONLY std::cout<<"updateY: "<<fstop-fstart<<" sec"<<std::endl<<std::flush;
        NODE0ONLY sprintf(buffer,",%f",fstop-fstart);
        NODE0ONLY strcat(timerstring, buffer);    
        
        std::string iter_str = num2str(iter);
        writeClassAverageCTF(result_fn+"_iter_"+iter_str,result_fn+"_iter_"+iter_str,probThreshold);
        
#ifdef USEMPI
		MPI::COMM_WORLD.Barrier();
#endif
        
#ifdef DEBUGGTM
		tstop = dtime();
		Ytime += (tstop - tstart);
		tstart = dtime();
#endif
        
		if(++iter > nr_iter)
			break;

        // --- updateBeta ----
        
        NODE0ONLY fstart = dtime();
        // Use the updated YReal and YImag to get a new value for Beta
		if(update_beta) tempBeta = updateBeta();  
        NODE0ONLY fstop = dtime();
        NODE0ONLY std::cout<<"updateBeta: "<<fstop-fstart<<" sec"<<std::endl<<std::flush;
        NODE0ONLY sprintf(buffer,",%f",fstop-fstart);
        NODE0ONLY strcat(timerstring, buffer);
        
        NODE0ONLY fstart = dtime();
     
        // TODO - fix - if !update_beta then tempBeta is not initialized
#ifdef USEMPI
		MPI::COMM_WORLD.Allreduce(&tempBeta,&AverageBeta,1,MPI::DOUBLE,MPI::SUM);
#else
        AverageBeta = tempBeta;
#endif
        
		AverageBeta = INDEXSZ(INDEXSZ(N)*INDEXSZ(D))/AverageBeta;
        
        if(!update_beta) AverageBeta = 0.01;
        NODE0ONLY fstop = dtime();
        NODE0ONLY std::cout<<"AverageBeta(="<<AverageBeta<<") reduction: "<<fstop-fstart<<" sec"<<std::endl<<std::flush;
        NODE0ONLY sprintf(buffer,",%f",fstop-fstart);
        NODE0ONLY strcat(timerstring, buffer);
        
#ifdef DEBUGGTM
        tstop = dtime();
        Btime += (tstop - tstart);
        tstart = dtime();
#endif
        
        // --- updateAlpha ----
        
        // if(update_alpha) tempAlpha = updateAlpha();
        
#ifdef USEMPI
		MPI::COMM_WORLD.Allreduce(&tempAlpha,&AverageAlpha,1,MPI::DOUBLE,MPI::SUM);
#endif
        
		AverageAlpha = K/AverageAlpha;
        
		if(alpha > 0) AverageAlpha = alpha;
        
        NODE0ONLY fstop = dtime();
        NODE0ONLY std::cout<<"AverageAlpha reduction: "<<fstop-fstart<<" sec"<<std::endl<<std::flush;
        NODE0ONLY sprintf(buffer,",%f",fstop-fstart);
        NODE0ONLY strcat(timerstring, buffer);
        NODE0ONLY fstart = dtime();
        
#ifdef DEBUGGTM
        tstop = dtime();
        Atime += (tstop - tstart);
        tstart = dtime();
#endif
        // --- update_P_Delta ----
        
        // rebuild R now that YReal/YImag have changed
        // Requires new YReal/YImag
		tempDeviation = update_P_Delta();  //update P and Delta
        NODE0ONLY fstop = dtime();
        NODE0ONLY std::cout<<"loop update_P_Delta: "<<fstop-fstart<<" sec"<<std::endl<<std::flush;
        NODE0ONLY sprintf(buffer,",%f",fstop-fstart);
        NODE0ONLY strcat(timerstring, buffer);
        NODE0ONLY fstart = dtime();        
#ifdef USEMPI
		MPI::COMM_WORLD.Allreduce(&tempDeviation,&deviation,1,MPI::DOUBLE,MPI::SUM);
#else
        deviation = tempDeviation;
#endif

        NODE0ONLY fstop = dtime();
        NODE0ONLY std::cout<<"deviation reduction: "<<fstop-fstart<<" sec"<<std::endl<<std::flush;
        NODE0ONLY sprintf(buffer,",%f",fstop-fstart);
        NODE0ONLY strcat(timerstring, buffer);
        NODE0ONLY fstart = dtime();
        
#ifdef USEMPI
		//reduce R to node 0
		double *tempRK = (double*)aMalloc(sizeof(double)*N,64);
		for(int k = 0;k < K;k++){
			MPI::COMM_WORLD.Reduce(R+INDEXSZ(INDEXSZ(k)*INDEXSZ(N)),tempRK,N,MPI::DOUBLE,MPI::SUM,0);
			NODE0ONLY memcpy(R+INDEXSZ(INDEXSZ(k)*INDEXSZ(N)),tempRK,sizeof(double)*N); /****node 0****/
		}
		aFree(tempRK);
#endif
        NODE0ONLY fstop = dtime();
        NODE0ONLY std::cout<<"R reduction: "<<fstop-fstart<<" sec"<<std::endl<<std::flush;
        NODE0ONLY sprintf(buffer,",%f",fstop-fstart);
        NODE0ONLY strcat(timerstring, buffer);
        NODE0ONLY fstart = dtime();
        
#ifdef DEBUGGTM
        tstop = dtime();
        Dtime += (tstop - tstart);
#endif
        
#ifdef DEBUGGTM
		NODE0ONLY{
			std::cout<<"AverageAlpha = "<<AverageAlpha<<std::endl;
			std::cout<<"averageBeta = "<<AverageBeta<<std::endl;
			std::cout<<"AverageAlpha/AverageBeta = "<<AverageAlpha/AverageBeta<<std::endl;
			std::cout<<"Deviation = "<<deviation<<std::endl;
		}
#endif

        sectionEnd = dtime();
        NODE0ONLY std::cout<<" ------- One iteration completed - "<<sectionEnd-sectionStart<<" sec"<<std::endl<<std::flush;
        NODE0ONLY sprintf(buffer,",%f",sectionEnd-sectionStart);
        NODE0ONLY strcat(timerstring, buffer);
	NODE0ONLY showProgressBar(iter-1, nr_iter);

	}while(!checkConvergeBeta(preci) || !checkConvergeAlpha(preci));
    // ----------- End Iterations -------------

	if(iter > nr_iter)
		NODE0ONLY std::cout<<"#### ^_^ not reach precision. ~_~ ######"<<std::endl;
	else
		NODE0ONLY std::cout<<"#### ^_^ reach presion. ~_~ ######"<<std::endl;

#ifdef DEBUGGTM
	NODE0ONLY{
		std::cout<<"after "<<(iter-1)<<" iterations,time cost : "<<std::endl;
		std::cout<<"update R       : "<<Rtime<<std::endl;
		std::cout<<"update W       : "<<Wtime<<std::endl;
        std::cout<<"update Y       : "<<Ytime<<std::endl;
		std::cout<<"update Beta    : "<<Btime<<std::endl;
		std::cout<<"update Alpha   : "<<Atime<<std::endl;
        std::cout<<"update Delta&R : "<<Dtime<<std::endl;
		std::cout<<"total          : "<<(Rtime+Wtime+Ytime+Btime+Atime+Dtime)<<std::endl;
	}
#endif

        NODE0ONLY showProgressBar(nr_iter, nr_iter);
    
}



MIC_OFFLOAD_ATTRIBUTES
void mkl_solve(double *A,int N,double *b,int M){
	LAPACKE_dgetrf(LAPACK_ROW_MAJOR,N,N,A,N,ipiv);
	LAPACKE_dgetrs(LAPACK_ROW_MAJOR,'N',N,M,A,N,ipiv,b,M);
}


// shift Matrix by column mean
void shiftMat(double *Mat,int M,int N){
	double max,min,distCorr;

	#pragma omp parallel for private(max,min,distCorr)
	for(int n = 0;n < N;n++){   //for each column
		max = (std::numeric_limits<double>::min)();
		min = (std::numeric_limits<double>::max)();

		for(int m = 0;m < M;m++){
			max = Mat[INDEXSZ(INDEXSZ(m) * INDEXSZ(N)) + INDEXSZ(n)] > 
                max ? Mat[INDEXSZ(INDEXSZ(m) * INDEXSZ(N)) + INDEXSZ(n)] : max;
			min = Mat[INDEXSZ(INDEXSZ(m) * INDEXSZ(N)) + INDEXSZ(n)] < 
                min ? Mat[INDEXSZ(INDEXSZ(m) * INDEXSZ(N)) + INDEXSZ(n)] : min;
	}
	distCorr = (max+min)/2.;
	// exp(709) < realmax < exp(710), plus a few digits margin to avoid
	distCorr = distCorr < (min+700)?distCorr:(min+700);

	#pragma simd
	for(int m = 0;m < M;m++)
		Mat[INDEXSZ(INDEXSZ(m) * INDEXSZ(N)) + INDEXSZ(n)] -= distCorr;
	}
}

/// ---------------------------------------------------------------------

void analysisR(double threshold/* = 0.98*/){

	NODE0ONLY std::cout<<"each k corresponding to n with probability big than"<<threshold<<std::endl;
	for(int k = 0;k < K;k++){

		NODE0ONLY std::cout<<"k = "<<std::setw(3)<<(k+1)<<"|";
		for(int n = 0;n < N;n++)
			if(R[INDEXSZ(INDEXSZ(k)* INDEXSZ(N)) + INDEXSZ(n)] > threshold) 
                NODE0ONLY std::cout<<" "<<std::setw(4)<<(n+1)<<" ";
		NODE0ONLY std::cout<<std::endl;
	}

	NODE0ONLY std::cout<<"each n corresponding to k with biggest probability."<<std::endl;
	NODE0ONLY std::cout<<"n"<<" K "<<"probability"<<std::endl;
	double big_prob = 0;
	double big_prob_k = 0;
	for(int n = 0;n < N;n++){
	big_prob = 0.;
    
	for(int k = 0;k < K;k++)
		if(R[INDEXSZ(INDEXSZ(k)* INDEXSZ(N)) + INDEXSZ(n)] > big_prob){
			big_prob = R[INDEXSZ(INDEXSZ(k)* INDEXSZ(N)) + INDEXSZ(n)];
			big_prob_k = k;
		}

		NODE0ONLY std::cout<<(n+1)<<" "<<(big_prob_k+1)<<" "<<big_prob<<std::endl;
	}

}

/// ---------------------------------------------------------------------

void writeClassAverageCTF(std::string filename_mrcs,std::string fn_metadata,double probThreshold){

	// all class average picture
	// ------- get each average picture corresponding picture number ----------
	NODE0ONLY{

        int classes_count[K];
        for (int k = 0; k < K; k++) {
            classes_count[k] = 0;
        }

        for (int n = 0; n < N; n++) {
            double max_R = (std::numeric_limits<double>::min)();
            for (int k = 0; k < K; k++) {
                if (R[INDEXSZ(INDEXSZ(k) * INDEXSZ(N)) + INDEXSZ(n)] > max_R) {
                    metadata[n].CLASS = k+1;
                    max_R = R[INDEXSZ(INDEXSZ(k) * INDEXSZ(N)) + INDEXSZ(n)];
                }
            }
            classes_count[(int)metadata[n].CLASS-1] ++;
            // metadata[n].PMAX = max_R;
        }

        if (fineSearch.doFineSearch) {
            fineSearch.append(metadata);
        }
        else{
            metadata.writeToStar(fn_metadata);
            std::string fn_metadata_info = fn_metadata+".info";
            std::ofstream infofile(fn_metadata_info.c_str());
            infofile<<"iclass"<<" images number"<<std::endl;
            for (int k = 0; k < K; k++) {
                infofile<<k+1<<" "<<classes_count[k]<<std::endl;
            }
            infofile.close();
        }
	}
 
	// ---------  get the average picture   --------------
	avTReal = (float*)aMalloc(sizeof(float)*d_interval*K,64);
	avTImag = (float*)aMalloc(sizeof(float)*d_interval*K,64);

	float *data;
    
    FFTWFTransformer fftback_worker(mrcsHead[1],mrcsHead[0]);
    
    NODE0ONLY {
        data = (float*)aMalloc(sizeof(float)*K*mrcsHead[1]*mrcsHead[0],64);
    }
	
    ERROR_CHECK(probThreshold < 0, "weightedsum should be in 0~1.");
    if (probThreshold >= 1.0)
        writeClassAverageCTFsub_maxPro();
    else
        writeClassAverageCTFsub(probThreshold);// write the weighted sum classaverage
    
    // Now, one output image at a time, reconstruct it from the data
    // scattered across the cluster, do the inverse FFT, and write the
    // result to the output file
    for(int k = 0;k < K;k++)
    {
        int bufsize = sizeof(float) * D;
        float RFFTbuffer[D], IFFTbuffer[D];
        float mergedRFFTbuffer[D], mergedIFFTbuffer[D];
        
        // Clear buffers so that the reduce works
        memset(&RFFTbuffer, 0x00, bufsize);
        memset(&IFFTbuffer, 0x00, bufsize);
        
        // Now copy our slice of the average image into the buffer
        for (int d = d_start; d < d_end; d++)
        {
            RFFTbuffer[d] = avTReal[INDEXSZ(INDEXSZ(d-d_start) * INDEXSZ(K)) + INDEXSZ(k)];
            IFFTbuffer[d] = avTImag[INDEXSZ(INDEXSZ(d-d_start) * INDEXSZ(K)) + INDEXSZ(k)];
        }
        
#ifdef USEMPI
        // Merge all this together on the head node using a reduce
        // Since each node will have only filled in its slide of RFFTbuffer and
        // IFFTbuffer, this is a whole lot of adding elements containing zero to
        // a single element containing the data from a given node
        MPI::COMM_WORLD.Reduce(RFFTbuffer,mergedRFFTbuffer,D,MPI::FLOAT,MPI::SUM,0);
        MPI::COMM_WORLD.Reduce(IFFTbuffer,mergedIFFTbuffer,D,MPI::FLOAT,MPI::SUM,0);
        
        MPI::COMM_WORLD.Barrier();
#else
        memcpy(mergedRFFTbuffer, RFFTbuffer, sizeof(float)*D);
        memcpy(mergedIFFTbuffer, IFFTbuffer, sizeof(float)*D);
#endif
        // Now Node 0 transforms the result from FFT space back to an image
        // and writes that image out
        NODE0ONLY
        {

            // data are coming in transposed, so
            fftback_worker.inverseFourierTransform(mergedRFFTbuffer, mergedIFFTbuffer, data+k*mrcsHead[1]*mrcsHead[0]);
            // Now write data back to the output MRCS file - note that we are writing
            // a single image here, so we need to append to what has been written
            // so far
			//write a single image

        }  // node 0 processing

	}

	NODE0ONLY{
        Mrcs::MrcsImages listOfImages(data,mrcsHead[1], K);
        
        for (int k = 0; k < K; k++)
            normalizeData(listOfImages.image_ptr(k), mrcsHead[1]*mrcsHead[0]);
        if (fineSearch.doFineSearch) {
            fineSearch.append(data);
        }
        else{
            Mrcs::writeMrcsData(filename_mrcs,listOfImages);
        }
	}

	NODE0ONLY aFree(data);
	aFree(avTReal);
	aFree(avTImag);

#ifdef USEMPI
	MPI::COMM_WORLD.Barrier();
#endif

}

/// ---------------------------------------------------------------------

void writeClassAverageCTFsub(double probThreshold){

	//------------ all class average picture -----------------
   
	// class average picture
	for(int d = 0;d < d_interval;d++)
		for(int k = 0;k < K;k++)
            avTReal[INDEXSZ(INDEXSZ(d) * INDEXSZ(K)) + INDEXSZ(k)] = 
                avTImag[INDEXSZ(INDEXSZ(d) * INDEXSZ(K)) + INDEXSZ(k)] = 0;

    // TODO : too large probThreshold will make more classaverage black.
    // so need use largest probability image instead.
  
	#pragma omp parallel for
	for(int d = 0;d < d_interval;d++){
        double temp;
		for(int k = 0;k < K;k++){
            
			for(int n = 0;n < N;n++){
                if (R[INDEXSZ(INDEXSZ(k) * INDEXSZ(N)) + INDEXSZ(n)] > probThreshold) {
                    avTReal[INDEXSZ(INDEXSZ(d) * INDEXSZ(K)) + INDEXSZ(k)] += 
                        R[INDEXSZ(INDEXSZ(k) * INDEXSZ(N)) + INDEXSZ(n)]*AverageBeta*
                        CTFvalue[INDEXSZ(INDEXSZ(d)*INDEXSZ(N))+ INDEXSZ(n)]*
                        TReal[INDEXSZ(INDEXSZ(d)*INDEXSZ(N))+ INDEXSZ(n)];
                    avTImag[INDEXSZ(INDEXSZ(d) * INDEXSZ(K)) + INDEXSZ(k)] += 
                        R[INDEXSZ(INDEXSZ(k) * INDEXSZ(N)) + INDEXSZ(n)]*AverageBeta*
                        CTFvalue[INDEXSZ(INDEXSZ(d)*INDEXSZ(N))+ INDEXSZ(n)]*
                        TImag[INDEXSZ(INDEXSZ(d)*INDEXSZ(N))+ INDEXSZ(n)];
                }
			}

			temp = 0;
			for(int n = 0;n < N;n++)
                if (R[INDEXSZ(INDEXSZ(k) * INDEXSZ(N)) + INDEXSZ(n)] > probThreshold) {
                    temp += R[INDEXSZ(INDEXSZ(k) * INDEXSZ(N)) + INDEXSZ(n)]*AverageBeta*
                        CTFvalue[INDEXSZ(INDEXSZ(d)*INDEXSZ(N))+ INDEXSZ(n)]*
                        CTFvalue[INDEXSZ(INDEXSZ(d)*INDEXSZ(N))+ INDEXSZ(n)];
                }

			temp += AverageAlpha;
			avTReal[INDEXSZ(INDEXSZ(d) * INDEXSZ(K)) + INDEXSZ(k)] = 
                avTReal[INDEXSZ(INDEXSZ(d) * INDEXSZ(K)) + INDEXSZ(k)]/temp;
			avTImag[INDEXSZ(INDEXSZ(d) * INDEXSZ(K)) + INDEXSZ(k)] = 
                avTImag[INDEXSZ(INDEXSZ(d) * INDEXSZ(K)) + INDEXSZ(k)]/temp;
		}
    }

}

/// ---------------------------------------------------------------------

void writeClassAverageCTFsub_maxPro(){
    
    // ---------------- all class average picture -----------------
    
    // class average picture
    for(int d = 0;d < d_interval;d++)
        for(int k = 0;k < K;k++)
            avTReal[INDEXSZ(INDEXSZ(d) * INDEXSZ(K)) + INDEXSZ(k)] = 
            avTImag[INDEXSZ(INDEXSZ(d) * INDEXSZ(K)) + INDEXSZ(k)] = 0;
    
    // double temp;
    double *max_R = (double*)aMalloc(sizeof(double)*N,64);
    int *max_R_k = (int*)aMalloc(sizeof(int)*N,64);
    
#pragma omp parallel for
    for (int n = 0; n < N; n++) {
        max_R[n] = (std::numeric_limits<double>::min)();
        max_R_k[n] = 0;
        for (int k = 0; k < K; k++) {
            if (R[INDEXSZ(INDEXSZ(k) * INDEXSZ(N)) + INDEXSZ(n)] > max_R[n]) {
                max_R[n] = R[INDEXSZ(INDEXSZ(k) * INDEXSZ(N)) + INDEXSZ(n)];
                max_R_k[n] = k;
            }
        }
        // std::cout<<max_R_k[n]<<" "<<max_R[n]<<std::endl;
    }
    
#pragma omp parallel for
    for(int d = 0;d < d_interval;d++)
        for(int k = 0;k < K;k++){
            
            for(int n = 0;n < N;n++){
                if (max_R_k[n] == k) {
                    avTReal[INDEXSZ(INDEXSZ(d) * INDEXSZ(K)) + INDEXSZ(k)] += 
                        max_R[n]*AverageBeta*
                        CTFvalue[INDEXSZ(INDEXSZ(d)*INDEXSZ(N))+ INDEXSZ(n)]*
                        TReal[INDEXSZ(INDEXSZ(d)*INDEXSZ(N))+ INDEXSZ(n)];
                    avTImag[INDEXSZ(INDEXSZ(d) * INDEXSZ(K)) + INDEXSZ(k)] += 
                        max_R[n]*AverageBeta*
                        CTFvalue[INDEXSZ(INDEXSZ(d)*INDEXSZ(N))+ INDEXSZ(n)]*
                        TImag[INDEXSZ(INDEXSZ(d)*INDEXSZ(N))+ INDEXSZ(n)];
                }
            }
            
            double temp = 0;
            for(int n = 0;n < N;n++)
                if (max_R_k[n] == k) {
                    temp += max_R[n]*AverageBeta*
                        CTFvalue[INDEXSZ(INDEXSZ(d)*INDEXSZ(N))+ INDEXSZ(n)]*
                        CTFvalue[INDEXSZ(INDEXSZ(d)*INDEXSZ(N))+ INDEXSZ(n)];
                }
            
            temp += AverageAlpha;
            avTReal[INDEXSZ(INDEXSZ(d) * INDEXSZ(K)) + INDEXSZ(k)] = 
                avTReal[INDEXSZ(INDEXSZ(d) * INDEXSZ(K)) + INDEXSZ(k)]/temp;
            avTImag[INDEXSZ(INDEXSZ(d) * INDEXSZ(K)) + INDEXSZ(k)] = 
                avTImag[INDEXSZ(INDEXSZ(d) * INDEXSZ(K)) + INDEXSZ(k)]/temp;
            
        }
    aFree(max_R);
    aFree(max_R_k);
}

/// ---------------------------------------------------------------------

// if iterated 3 times and each change smaller than precision,then convergence
bool checkConvergeBeta(double precision){
    static double a = -1.,b = 0.,c = 1.;
    // static double precision = 10e-10;
    
    // if AverageBeta is not change
    if(c == AverageBeta) return false;
    a = b;
    b = c;
    c = AverageBeta;
    // std::cout<<a<<" "<<b<<" "<<c<<" "<<std::endl;
    return (fabs(a-b)/a < precision) && (fabs(c-b)/b < precision);
}

/// ---------------------------------------------------------------------

bool checkConvergeAlpha(double precision){
    static double a = -1.,b = 0.,c = 1.;
    // static double precision = 10e-10;
    // if AverageBeta is not change
    if(c == AverageAlpha) return false;
    a = b;
    b = c;
    c = AverageAlpha;
    // std::cout<<a<<" "<<b<<" "<<c<<" "<<std::endl;
    return (fabs(a-b)/a < precision) && (fabs(c-b)/b < precision);
}

/// ---------------------------------------------------------------------

bool checkConvergeDeviation(double precision/* = 10e-20*/){
    static double a = -1.,b = 0.,c = 1.;
    // static double precision = 10e-10;
    //if deviation is not change
    if(c == deviation) return false;
    a = b;
    b = c;
    c = deviation;
    // std::cout<<a<<" "<<b<<" "<<c<<" "<<std::endl;
    return (fabs(a-b)/a < precision) && (fabs(c-b)/b < precision);
}

/// ---------------------------------------------------------------------
//  Debugging Functions
/// ---------------------------------------------------------------------

void dumpT(){
	std::cout<<"node = "<<node<<std::endl<<std::flush;
	std::cout<<"fft data and fft complete."<<std::endl<<std::flush;
	std::cout<<"here is some info:"<<std::endl<<std::flush;
	std::cout<<"TReal[0] = "<<TReal[0]<<",TImag[0] = "<<TImag[0]<<std::endl<<std::flush;
	std::cout<<"TReal[1] = "<<TReal[1]<<",TImag[1] = "<<TImag[1]<<std::endl<<std::flush;
	std::cout<<"TReal[2] = "<<TReal[2]<<",TImag[2] = "<<TImag[2]<<std::endl<<std::flush;
	std::cout<<"TReal[3] = "<<TReal[3]<<",TImag[3] = "<<TImag[3]<<std::endl<<std::flush;
	std::cout<<"TReal[4] = "<<TReal[4]<<",TImag[4] = "<<TImag[4]<<std::endl<<std::flush;
	std::cout<<"TReal[5] = "<<TReal[5]<<",TImag[5] = "<<TImag[5]<<std::endl<<std::flush;
    
	std::cout<<"N = "<<N<<",D = "<<D<<std::endl<<std::flush;
	std::cout<<"WReal[0] = "<<WReal[0]<<",WImag[0] = "<<WImag[0]<<std::endl<<std::flush;
	std::cout<<"WReal[1] = "<<WReal[1]<<",WImag[1] = "<<WImag[1]<<std::endl<<std::flush;
	std::cout<<"WReal[2] = "<<WReal[2]<<",WImag[2] = "<<WImag[2]<<std::endl<<std::flush;
	std::cout<<"WReal[3] = "<<WReal[3]<<",WImag[3] = "<<WImag[3]<<std::endl<<std::flush;
	std::cout<<"WReal[4] = "<<WReal[4]<<",WImag[4] = "<<WImag[4]<<std::endl<<std::flush;
	std::cout<<"WReal[5] = "<<WReal[5]<<",WImag[5] = "<<WImag[5]<<std::endl<<std::flush;
}

void dumpR(){
    std::cout<<"node = "<<node<<",update R!"<<std::endl<<std::flush;
    // prtMatrix(R,1,10);
    int rows = 1;
    int cols = 10;
    std::cout.flags(std::ios::left);
    std::cout.precision(8);
    for(int i = 0;i < rows;i++){
        for(int j = 0;j < cols;j++)
            std::cout<<std::setw(15)<<R[i*cols+j]<<" ";
        std::cout<<std::endl;
    }
    
}

void dumpW(){
	std::cout<<"node = "<<node<<std::endl<<std::flush;
	std::cout<<"update WReal."<<std::endl<<std::flush;
	for(int m = 0;m < 2;m++){
		for(int d = 0;d < 10;d++)
			std::cout<<WReal[d*M+m]<<" ";
			std::cout<<std::endl<<std::flush;
	}
	std::cout<<"update WImag."<<std::endl<<std::flush;
	for(int m = 0;m < 2;m++){
		for(int d = 0;d < 10;d++)
			std::cout<<WImag[d*M+m]<<" ";
			std::cout<<std::endl<<std::flush;
	}
}

void dumpY(){
	std::cout<<"node = "<<node<<std::endl<<std::flush;
	std::cout<<"update YReal."<<std::endl<<std::flush;
    for(int k = 0;k < 2;k++){
		for(int d = 0;d < 10;d++)
			std::cout<<YReal[d*K+k]<<" ";
		std::cout<<std::endl<<std::flush;
    }
    std::cout<<"update YImag."<<std::endl<<std::flush;
    for(int k = 0;k < 2;k++){
		for(int d = 0;d < 10;d++)
			std::cout<<YImag[d*K+k]<<" ";
		std::cout<<std::endl<<std::flush;
    }
}

void dumpBelta(){
	// std::cout<<"node = "<<node<<std::endl<<std::flush;
	// std::cout<<"iterate update Beta!"<<std::endl<<std::flush;
	// prtMatrix(Beta,N,D);
}

void TestOffloadMemory(int mic_index)
{
    double fstart, fstop;
    
    START_HETERO_OFFLOAD
    
    fstart = dtime();
 
#ifndef DONT_USE_OFFLOAD
#ifdef SERIAL_OFFLOAD
#pragma offload if(mic_index>0) target(mic:(mic_index-1))  \
    nocopy( TReal[node_d_start*N:node_d_chunksize*N]    : alloc( TReal[node_d_start*N:node_d_chunksize*N]) REUSE) \
    nocopy( TImag[node_d_start*N:node_d_chunksize*N]    : alloc( TImag[node_d_start*N:node_d_chunksize*N]) REUSE) \
    nocopy( CTFvalue[node_d_start*N:node_d_chunksize*N] : alloc( CTFvalue[node_d_start*N:node_d_chunksize*N]) REUSE) \
    nocopy( YReal[node_d_start*K:node_d_chunksize*K]: alloc( YReal[node_d_start*K:node_d_chunksize*K]) REUSE) \
    nocopy( YImag[node_d_start*K:node_d_chunksize*K]: alloc( YImag[node_d_start*K:node_d_chunksize*K]) REUSE) \
    nocopy( WReal[node_d_start*M:node_d_chunksize*M]    : alloc( WReal[node_d_start*M:node_d_chunksize*M]) REUSE) \
    nocopy( WImag[node_d_start*M:node_d_chunksize*M]    : alloc( WImag[node_d_start*M:node_d_chunksize*M]) REUSE) \
    nocopy( PHI:length(M*K) REUSE) \
    nocopy( R:length(N*K) REUSE)
#else
#pragma offload if(mic_index>0) target(mic:(mic_index-1)) signal(&YReal[node_d_start]) \
    nocopy( TReal[node_d_start*N:node_d_chunksize*N]    : alloc( TReal[node_d_start*N:node_d_chunksize*N]) REUSE) \
    nocopy( TImag[node_d_start*N:node_d_chunksize*N]    : alloc( TImag[node_d_start*N:node_d_chunksize*N]) REUSE) \
    nocopy( CTFvalue[node_d_start*N:node_d_chunksize*N] : alloc( CTFvalue[node_d_start*N:node_d_chunksize*N]) REUSE) \
    nocopy( YReal[node_d_start*K:node_d_chunksize*K]: alloc( YReal[node_d_start*K:node_d_chunksize*K]) REUSE) \
    nocopy( YImag[node_d_start*K:node_d_chunksize*K]: alloc( YImag[node_d_start*K:node_d_chunksize*K]) REUSE) \
    nocopy( WReal[node_d_start*M:node_d_chunksize*M]    : alloc( WReal[node_d_start*M:node_d_chunksize*M]) REUSE) \
    nocopy( WImag[node_d_start*M:node_d_chunksize*M]    : alloc( WImag[node_d_start*M:node_d_chunksize*M]) REUSE) \
    nocopy( PHI:length(M*K) REUSE) \
    nocopy( R:length(N*K) REUSE)
#endif
#endif
    {
        double testval = 0;
        
        for (int d = node_d_start; d < node_d_end; d++)
        {
            for (int n = 0; n < N; n++)
            {
                testval += TReal[d*N+n];
                testval += TImag[d*N+n];
                testval += CTFvalue[d*N+n];
            }  // N
            for (int m = 0; m < M; m++)
            {
                testval += WReal[d*M+m];
                testval += WImag[d*M+m];
            }  // M
            for (int k = 0; k < K; k++)
            {
                testval += YReal[d*K+k];
                testval += YImag[d*K+k];
            }  // K
        }  // d
        for (int k = 0; k < K; k++)
        {
            for (int n = 0; n < N; n++)
            {
                testval += R[k*N+n];
            }
            for (int m = 0; m < M; m++)
            {
                testval += PHI[k*M+m];
            }
        }
    }  // pragma offload
    fstop = dtime();
    
    std::cout<<"TestOffloadMemory body done "<<mic_index<<": "<<fstop-fstart<<" sec"<<std::endl<<std::flush;
    
    STOP_HETERO_OFFLOAD
}
    
    
double calSpace(){
    double spaceNeed = 0;
    // some space keep constant
    spaceNeed += sizeof(*PHI)*K*M;//for PHI
    spaceNeed += sizeof(*TReal)*d_interval*N*2;//for TReal and TImag
    spaceNeed += sizeof(*CTFvalue)*d_interval*N;//for CTFvalue
    spaceNeed += sizeof(*WReal)*d_interval*M*2;//for WReal and WImag
    spaceNeed += sizeof(*YReal)*d_interval*K*2;//for YReal and YImag
    spaceNeed += sizeof(*R)*K*N;//for R
    spaceNeed += sizeof(*Alpha)*K;//for Alpha
    // some space for write data
    spaceNeed += sizeof(*avTReal)*K*d_interval*2;//for avTReal and avTImag
    spaceNeed += sizeof(float)*K*D;//for write mrcs data
    // some space for iterateUpdateW()
    spaceNeed += maxthreads*sizeof(double)*M*M;//for AA
    spaceNeed += maxthreads*sizeof(double)*M*2;//for bRealArray and bImagArray
    spaceNeed += sizeof(double)*M*M;//for A
    spaceNeed += sizeof(double)*2*M;//for b
    
    spaceNeed /= (1024*1024*1024);
    return spaceNeed;
}
    
double printDataSize(){
    // Print out resulting memory requirements
    NODE0ONLY
    {
        // Memory needed on all nodes
        size_t PHImem = sizeof(double)*K*M;
        size_t Rmem = sizeof(double)*K*N;
        size_t TRImem = 2 * sizeof(float)*d_interval*N;
        size_t YRImem = 2 * sizeof(double)*d_interval*K;
        size_t WRImem = 2 * sizeof(float)*d_interval*M;
        size_t FFTbuffermem = 2 * sizeof(float)*D;
        size_t CTFmem = sizeof(float)*N*d_interval;
        size_t CTFresmem = sizeof(double)*D;
        size_t AAmem = maxthreads*sizeof(double)*M*M;
        size_t bRImem = 2 * maxthreads*sizeof(double)*M;
        size_t Amem = sizeof(double)*M*M;
        size_t bmem = sizeof(double)*2*M;
        size_t tRKmem = sizeof(double)*N;
        size_t ipivmem = M*sizeof(int);
        // Ignoring memory use in writeY2File, currently unused
        // TODO - add that to the calculation if used
        size_t oriensmem = 0;
        size_t Ymem = 0;
        size_t Ydatamem = 0;
        size_t avTmem = sizeof(double)*K*2;
        size_t avTRImem = 2 * sizeof(float)*d_interval*K;
        size_t writedatamem = sizeof(float)*picdim*picdim;
        
        // Memory needed only on node 0
        size_t Wmem = sizeof(float)*M*picdim*picdim;
        
        size_t computemem = PHImem + Rmem + YRImem + TRImem + WRImem +
            FFTbuffermem + CTFmem + CTFresmem + AAmem + bRImem +
            Amem + bmem + tRKmem + ipivmem + oriensmem + Ymem +
            Ydatamem + avTmem + avTRImem + writedatamem;
        double computememGB = (double)computemem / 1073741824.0;
        size_t node0mem = computemem + Wmem;
        double node0memGB = (double)node0mem / 1073741824.0;
        
        std::cout << "-------------------------------------------------"<<std::endl<<std::flush;
        std::cout << "N = " <<N<< ", D = " <<D<< ", K = " <<K<< ", M = " <<M<<std::endl<<std::flush;
        std::cout << "Number of nodes="<<nodes<<" max threads="<<maxthreads<<std::endl<<std::flush;
        std::cout << "Memory use (bytes):"<<std::endl<<std::flush;
        std::cout << "Total on compute nodes: " << computemem << " (" << computememGB << " GB)" <<std::endl<<std::flush;
        std::cout << "Total on Node 0: " << node0mem << " (" << node0memGB << " GB)" <<std::endl<<std::flush;
        std::cout << "   PHI: " << PHImem << " (" << ((double)PHImem / 1073741824.0) << " GB)" <<std::endl<<std::flush;
        std::cout << "   R: " << Rmem << " (" << ((double)Rmem / 1073741824.0) << " GB)" <<std::endl<<std::flush;
        std::cout << "   TReal & TImag: " << TRImem << " (" << ((double)TRImem / 1073741824.0) << " GB)" <<std::endl<<std::flush;
        std::cout << "   YReal & YImag: " << YRImem << " (" << ((double)YRImem / 1073741824.0) << " GB)" <<std::endl<<std::flush;
        std::cout << "   WReal & WImag: " << WRImem << " (" << ((double)WRImem / 1073741824.0) << " GB)" <<std::endl<<std::flush;
        std::cout << "   FFTbuffers: " << FFTbuffermem << " (" << ((double)FFTbuffermem / 1073741824.0) << " GB)" <<std::endl<<std::flush;
        std::cout << "   CTFvalue: " << CTFmem << " (" << ((double)CTFmem / 1073741824.0) << " GB)" <<std::endl<<std::flush;
        std::cout << "   CTFvalue read memory: " << CTFresmem << " (" << ((double)CTFresmem / 1073741824.0) << " GB)" <<std::endl<<std::flush;
        std::cout << "   AA: " << AAmem << " (" << ((double)AAmem / 1073741824.0) << " GB)" <<std::endl<<std::flush;
        std::cout << "   bReal and bImag: " << bRImem << " (" << ((double)bRImem / 1073741824.0) << " GB)" <<std::endl<<std::flush;
        std::cout << "   A: " << Amem << " (" << ((double)Amem / 1073741824.0) << " GB)" <<std::endl<<std::flush;
        std::cout << "   b: " << bmem << " (" << ((double)bmem / 1073741824.0) << " GB)" <<std::endl<<std::flush;
        std::cout << "   Partial R buffer: " << tRKmem << " (" << ((double)tRKmem / 1073741824.0) << " GB)" <<std::endl<<std::flush;
        std::cout << "   MKL ipiv: " << ipivmem << " (" << ((double)ipivmem / 1073741824.0) << " GB)" <<std::endl<<std::flush;
        std::cout << "   orientations: " << oriensmem << " (" << ((double)oriensmem / 1073741824.0) << " GB)" <<std::endl<<std::flush;
        std::cout << "   complex Y: " << Ymem << " (" << ((double)Ymem / 1073741824.0) << " GB)" <<std::endl<<std::flush;
        std::cout << "   Y fft buffer: " << Ydatamem << " (" << ((double)Ydatamem / 1073741824.0) << " GB)" <<std::endl<<std::flush;
        std::cout << "   avTnumbers: " << avTmem << " (" << ((double)avTmem / 1073741824.0) << " GB)" <<std::endl<<std::flush;
        std::cout << "   avTReal and avTImag: " << avTRImem << " (" << ((double)avTRImem / 1073741824.0) << " GB)" <<std::endl<<std::flush;
        std::cout << "   picture write buffer: " << writedatamem << " (" << ((double)writedatamem / 1073741824.0) << " GB)" <<std::endl<<std::flush;
        std::cout << "   W buffer on node 0: " << Wmem << " (" << ((double)Wmem / 1073741824.0) << " GB)" <<std::endl<<std::flush;
        std::cout << "-------------------------------------------------"<<std::endl<<std::flush;
        return(node0memGB);
    }
    return(0.0);
}

    
    
} // namespace GTMoptimizer
