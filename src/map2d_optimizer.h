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
#ifndef MAP2D_OPTIMIZER_OLD_H_
#define MAP2D_OPTIMIZER_OLD_H_

#include "./map_optimizer_base.h"

namespace Map2dOptimizer_old
{
    using namespace MapOptimizerBase_old;
    
    // model
    extern HealpixSampler sampler2d;
    extern MAPModel mapModel;
    extern MLModel mlModel;
    extern ParticleModel particleModel;
    extern HiddenVariableMonitor hiddenVarMonitor;
    
    // image data
    extern Images images;
    extern MetaDataTable metadata;
	
    // sampling
    extern double offset_step;
    extern double offset_range;
    extern double psi_step;
    extern SamplingGrid samplingGrid;
    
	// ------------ link  some variable ------------- //
	static auto& particle_diameter = particleModel.particle_diameter;
    static auto& only_flip_phases = particleModel.only_flip_phases;
    static auto& ctf_phase_flipped = particleModel.ctf_phase_flipped;
    static auto& intact_ctf_first_peak = particleModel.intact_ctf_first_peak;
    
    // ------------ variable for expectation step -------------- //
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
    extern Aligned3dArray<FDOUBLE> exp_Fimgs_shifted_real;
    extern Aligned3dArray<FDOUBLE> exp_Fimgs_shifted_imag;
    extern Aligned3dArray<FDOUBLE> exp_Frefs_Rot_real;
    extern Aligned3dArray<FDOUBLE> exp_Frefs_Rot_imag;
    extern Aligned3dArray<FDOUBLE> exp_Fweight_Rot;
    
    // metadata of the weight,structure is iclass->rot->iover_rot->trans->iover_trans
    extern double* exp_Mweight;
    extern int exp_Mweight_xsize;
    
    // ------------ thread variable ----------------- //
    extern Aligned2dArray<FDOUBLE> thread_Frefctf_real;
    extern Aligned2dArray<FDOUBLE> thread_Frefctf_imag;
    extern Aligned2dArray<FDOUBLE> thread_wsum_pdf_direction;
    extern std::vector<std::vector<GridIndex>> thread_exp_max_weight_index;
    
    //
    void readImages();
    
    void setupMap2dOptimizer();
    
    void destroyMap2dOptimizer();

    // Interpret command line for the initial start of a run
    void prepare();
    
    // Perform expectation-maximization iterations
    void iterate();
    
	// ------------   EM-Iteration     ----------------- //
    
	void expectation();
    
    // expectation nr_pool image each time
    void expectationSomeParticles();
    
    // (de)allocate memory space for each expectation step
    void prepareExpData();
    void endExpData();
    
    // ------------ some function in expectationsomeparticles functon    --------- //
    
    // get all rotated reference  and the significant rotation
    void getReferenceAllOrientations();

	// get all reference and images 's squared differences
	void getAllSquaredDifferences();
    
    // calculates exp_sum_weight and, for adaptive approach, also exp_significant_weight
    void findAllSignificantPoints();
    
	// convert all squared difference to weight(P = exp(-x))
	void convertSquaredDifferencesToWeights();

    // calculate norm_correction,dLL and Pmax
    void storeWeightedSums();
    
	// update other parameters for refine the model:norm_correction,sigma2_noise,pdf_class,pdf_direction,prior_offsetx(y)_class
	void updateOtherParams();
    
    // add all shifted and rotated images back
    void backProjection();
    
    // ------------  Maximization step    ------------ //
    
    void maximization();
    
    // Perform the actual reconstructions
    void maximizationReconstructClass(int iclass);
    
    // Updates all other model parameters (besides the reconstructions)
    void maximizationOtherParameters();
    
    // ------------  Write files to disc     ------------- //

    void writeClassesAndMetadata(std::string fn_class,std::string fn_metadata);
    
    // ------------  some help function     ------------- //
    inline bool isSignificantAnyParticleAnyTranslation(int iorient);
    
    // --------------------------------------------------- //
    void calspace(int ori_size,int current_size,int set_nr_pool);
    
};


#endif
