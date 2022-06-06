/***************************************************************************
 *
 * Authors: "Yongbei(Glow) Ma,Jiayi (Timmy) Wu, Youdong (Jack) Mao"
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
#ifndef MAP3D_OPTIMIZER_OLD_H_
#define MAP3D_OPTIMIZER_OLD_H_

// #include "tbb/parallel_sort.h"
#include "./map_optimizer_base.h"

namespace Map3dOptimizer_old
{
    using namespace MapOptimizerBase_old;
    
    // model
    extern HealpixSampler sampler3d;
    extern MAPModel mapModel;
    extern MLModel mlModel;
    extern ParticleModel particleModel;
    extern HiddenVariableMonitor hiddenVarMonitor;
    
    // image data
    extern Images images_data;
    extern MetaDataTable metadata;
    
    // sampling
    extern double offset_step;
    extern double offset_range;
    extern int sampler3d_healpix_order;
    extern std::string sampler3d_fn_sym;
    extern SamplingGrid samplingGrid;
    
    // ------------ link  some variable ------------- //
    static auto& particle_diameter = particleModel.particle_diameter;
    static auto& only_flip_phases = particleModel.only_flip_phases;
    static auto& ctf_phase_flipped = particleModel.ctf_phase_flipped;
    static auto& intact_ctf_first_peak = particleModel.intact_ctf_first_peak;
        
    // ------------ variable for expectation step  ------------- //
#define SEP
#define ELTONE(T,N,S1,S2) extern T N;
#define ELTVE1(T,N,S1,S2) extern VectorOfArray1d<T> N;
#define ELTVE2(T,N,S1,S2) extern VectorOfArray2d<T> N;
#define ELTVE3(T,N,S1,S2) extern Aligned3dArray <T> N;
    MAPOPTIMIZER_OLD_EXP_VARS
    MAPOPTIMIZER_OLD_THREAD_VARS
#undef SEP
#undef ELTONE
#undef ELTVE1
#undef ELTVE2
#undef ELTVE3
    //
    extern VectorOfStruct<MetaDataElem> exp_metadata;
    //
    //
    // ---------------   thread variable   ------------ //
    extern Aligned3dArray<FDOUBLE> thread_Frefctf_real;
    extern Aligned3dArray<FDOUBLE> thread_Frefctf_imag;
    extern Aligned3dArray<FDOUBLE> thread_Fimg_real;
    extern Aligned3dArray<FDOUBLE> thread_Fimg_imag;
    extern Aligned3dArray<FDOUBLE> thread_Fweight;
    extern Aligned3dArray<FDOUBLE> thread_wsum_pdf_direction;
    extern Aligned2dArray< char > threadfake_do_scale_norm_class;
    extern std::vector<std::vector<GridIndex>> thread_exp_max_weight_index;
    
    //  ---------------   setup      ---------------------- //
    void setupMLoptimizer();
    
    void destroyMLoptimizer();
    
    // Interpret command line for the initial start of a run
    void prepare();
    
    // Perform expectation-maximization iterations
    void iterate();
    
    // ------------   EM-Iteration     ----------------- //
    
    // (de)allocate memory space for each expectation step
    void prepareExpData();
    void endExpData();
    
    void expectation();
    
    // expectation nr_pool image each time
    void expectationSomeParticles();
    
    // ------------ some function in expectationsomeparticles functon    --------- //
    
    // get all rotated reference  and the significant rotation
    void getReferenceAllOrientations();
    
    // get all reference and images 's squared differences
    void getAllSquaredDifferences(bool do_coarse_search);
    
    // calculates exp_sum_weight and, for adaptive approach, also exp_significant_weight
    void findAllSignificantPoints(Exp_Mweight_old& exp_Mweight);
    
    // convert all squared difference to weight(P = exp(-x))
    void convertSquaredDifferencesToWeights(Exp_Mweight_old& exp_Mweight);
    
    // calculate norm_correction,dLL and Pmax
    void storeWeightedSums();
    
    // update mlModel
    void updateOtherParams();
    
    // add all shifted and rotated images back
    void backProjection();
    
    // ------------  Maximization step    ------------ //
    
    void maximization();
    
    // Updates all other model parameters (besides the reconstructions)
    void maximizationOtherParameters();
    
    // ------------  read and Write files     ------------- //
    
    void readResult();
    void writeAndCheckResult();
    void printMem(int set_nr_pool);
    
    // ------------  some help function     ------------- //
    
    inline bool isSignificantAnyParticleAnyTranslation(int iorient);
    
    //
    // use StatusTracer to write data to disk
    // to compare the different.....
    
};


#endif
