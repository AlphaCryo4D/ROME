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

#ifndef MAP_OPTIMIZERBASE_H_
#define MAP_OPTIMIZERBASE_H_

#include "mkl.h"
#include <omp.h>
//
#include "./resmap_util.h"
#include "./resmap_error.h"
#include "./resmap_string.h"
#include "./resmap_mpi.h"
#include "./resmap_time.h"
#include "./resmap_mrcs.h"
#include "./resmap_image.h"
#include "./resmap_initialize.h"
#include "./resmap_ctf.h"
#include "./resmap_fft_fftw3.h"
#include "./resmap_sampling.h"
#include "./resmap_metadata.h"
#include "./resmap_exp_model.h"
#include "./resmap_map_model.h"
#include "./resmap_ml_model.h"
//
#include "./resmap_progressbar.h"
#include "./resmap_exp_mweight.h"

#include "../../src_tools/statusTracer.h"

namespace MapOptimizerBase
{
#ifdef USEMPI
    extern int nodeNameLen;
    extern char nodeName[MPI_MAX_PROCESSOR_NAME];
#endif
	//
    // ------------      base variables         	 -------------- //
#define MAPOPTIMIZER_BASE_VARS \
    /* (TYPE		, NAME							,VALUE ) */	\
    ELT(int			, maxthreads					, 0		)	SEP \
    /* how to divide image across the nodes */						\
    ELT(int			, node							, 0		)	SEP \
    ELT(int			, nodes							, 0		)	SEP \
    ELT(bool		, debug_flag					, 0		)	SEP \
    ELT(bool		, tile_size_select				, 0		)	SEP \
	ELT(int			, nr_local_images				, 0		)	SEP \
    ELT(int			, nr_global_images				, 0		)	SEP \
    ELT(int			, nr_pool						, 0		)	SEP \
	ELT(int			, first_local_image				, 0		)	SEP \
    ELT(int			, last_local_image				, 0		)	SEP \
	ELT(std::string	, mrcs_dir						, "NULL")	SEP \
	ELT(std::string	, star_fn						, "NULL")	SEP \
    ELT(std::string	, write_path					, "NULL")	SEP \
    ELT(std::string	, write_fn						, "NULL")	SEP \
    ELT(std::string	, ref_fn						, "NULL")	SEP \
    ELT(std::string	, mask_fn						, "NULL")	SEP \
    ELT(std::string	, tile_size_fn					, "NULL")	SEP \
    ELT(std::string	, compare_matrix_fn				, "NULL")	SEP \
    ELT(std::string	, guidance_fn					, "NULL")	SEP \
	ELT(std::string	, continue_fn					, "NULL")	SEP \
    ELT(std::string	, data_stream_in_fn				, "NULL")	SEP \
    ELT(std::string	, data_stream_out_fn			, "NULL")	SEP \
    /* Seed for random number generator */							\
    ELT(int			, random_seed					, 0		)	SEP \
    ELT(int			, adaptive_oversampling			, 0		)	SEP \
	/* Use images only up to a certain resolution in the expectation step */ \
    ELT(int			, coarse_size					, 0		)	SEP \
    ELT(int			, ori_size						, 0		)	SEP \
    ELT(int			, current_size					, 0		)	SEP \
    ELT(int			, nr_classes					, 0		)	SEP \
	/* Pixel size (in Angstrom) */									\
	ELT(double		, pixel_size					, 0		)	SEP \
    /* Total number iterations and current iteration */				\
	ELT(int			, iter							, 0		)	SEP \
    ELT(int			, nr_iter						, 20	)	SEP \
	/* Flag whether to do CTF correction */							\
	ELT(bool		, do_ctf_correction				, false	)	SEP \
    /* Flag whether current references are ctf corrected */			\
	ELT(bool		, refs_are_ctf_corrected		, false	)	SEP \
    /* Flag whether to do image-wise intensity-scale correction */	\
    ELT(bool		, do_norm_correction			, false	)	SEP \
	/* Flag to indicate the refinement has converged */				\
    ELT(bool		, has_converged					, false	)	SEP \
    /* Flag whether to use maximum a posteriori (MAP) estimation */	\
    ELT(bool		, do_map						, true	)	SEP \
	/* Do zero-filled soft-masks instead of noisy masks on experimental particles */ \
	/* This will increase SNR but introduce correlations that are not modelled... */ \
	/* Until now the best refinements have used the noisy mask, not the soft mask.... */ \
	ELT(bool		, do_zero_mask					, false	)	SEP \
	/* Flag to flatten solvent */									\
    ELT(bool		, do_solvent					, false	)	SEP \
	/* get shifted images on-fly */									\
    ELT(bool		, do_shifts_onthefly			, false	)	SEP \
	/* first iteration do cross-correlation,use to get 3D volume seed*/\
    ELT(bool		, do_firstiter_cc				, false	)	SEP \
	/* always do cross-correlation */								\
    ELT(bool		, do_always_cc					, false	)	SEP \
	/* do scale correction ,I have not add this to map2d */			\
    ELT(bool		, do_scale_correction			, false	)	SEP \
    /* Flag whether to split data from the beginning into two random halves */ \
    ELT(bool		, do_split_random_halves		, false	)	SEP \
	/* Initial low-pass filter for all references (in digital frequency) */ \
    ELT(double		, ini_high						, 0		)	SEP \
	/* Number of iterations for gridding preweighting reconstruction */	\
	ELT(int			, gridding_nr_iter				, 10	)	SEP \
	/* Multiplicative fdge factor for the sigma estimates */		\
	ELT(double		, sigma2_fudge					, 1		)	SEP \
	/* Variance angle for the orientational pdf */					\
    ELT(double		, sigma2_angle					, 0		)	SEP \
    /* Fudge factor to adjust estimated tau2_class spectra */		\
    ELT(double		, tau2_fudge_factor				, 0		)	SEP \
	/* Minimum resolution to perform Bayesian estimates of the model */ \
	ELT(int			, minres_map					, 5		)	SEP \
    /* Width of the soft-edges of the circular masks */ 			\
	ELT(int			, width_mask_edge				, 5		)	SEP \
    ELT(int			, width_fmask_edge				, 2		)	SEP \
	/* Fraction of the weights for which to consider the finer sampling */ \
	/* The closer to one, the more orientations will be oversampled */ \
	/* The default is 0.999. */ 									\
	ELT(double		, adaptive_fraction				, 0.999f)	// end of MAPOPTIMIZER_BASE_VARS
    
    
#define AUTOREFINE_BASE_VARS \
	/* Flag whether to use the auto-refine procedure */	\
	ELT(bool		, do_auto_refine				,false	)	SEP \
    /* In auto-sampling mode, use local searches from this sampling rate on */ \
    ELT(int			, autosampling_hporder_local_searches, 4)	SEP	\
    /* resolution (in Angstrom) to join the two random halves */	\
    ELT(double		, low_resol_join_halves			, -1	)	SEP	\
	/* Flag to join random halves again */							\
    ELT(bool		, do_join_random_halves			, false	)	SEP \
	/* Flag to indicate to use all data out to Nyquist */			\
    ELT(bool		, do_use_all_data				, false	)	SEP \
	/* Number of iterations without a resolution increase */		\
    ELT(int			, nr_iter_wo_resol_gain			, 0		)	SEP \
    /* Best resolution obtained thus far */							\
    ELT(double		, best_resol_thus_far			, 1./999.)	SEP	\
	/* Is the FSC still high at the resolution limit? */			\
    ELT(bool 		, has_high_fsc_at_limit			, false	)	SEP \
	/* Flag to indicate that angular sampling in auto-sampling has reached its limit */	\
    ELT(bool		, has_fine_enough_angular_sampling, false)	SEP	\
	/* Flag to always join random halves, this is a developmental option for testing of sub-optimal FSC-usage only! */ \
	ELT(bool		, do_always_join_random_halves	, false	)	// end of AUTOREFINE_BASE_VARS
    
#define SEP
#define ELT(T,N,V) extern T N;
    MAPOPTIMIZER_BASE_VARS
#undef SEP
#undef ELT
    
    // ------------ variable for expectation step -------------- //
#define MAPOPTIMIZER_EXP_VARS \
	ELTONE(int		, exp_first_image					, 1				, 1				) SEP \
	ELTONE(int		, exp_last_image					, 1				, 1				) SEP \
    ELTONE(int		, exp_nr_images						, 1				, 1				) SEP \
	/*Loop from iclass_min to iclass_max to deal with seed generation in the first iteration*/\
    ELTONE(int		, exp_iclass_min					, 1				, 1				) SEP \
	ELTONE(int		, exp_iclass_max					, 1				, 1				) SEP \
    ELTONE(int		, exp_ipass							, 1				, 1				) SEP \
    ELTONE(int		, exp_current_size					, 1				, 1				) SEP \
    ELTONE(int		, exp_current_oversampling			, 1				, 1				) SEP \
    ELTONE(int		, exp_nr_over_trans					, 1				, 1				) SEP \
    ELTONE(int		, exp_nr_over_rot					, 1				, 1				) SEP \
    ELTONE(int		, exp_nr_dir						, 1				, 1				) SEP \
    ELTONE(int		, exp_nr_psi						, 1				, 1				) SEP \
    ELTONE(int		, exp_nr_trans						, 1				, 1				) SEP \
	ELTVE1(FDOUBLE	, exp_highres_Xi2_imgs				, 1				, exp_nr_images	) SEP \
    ELTVE1(double	, exp_min_diff2						, 1				, exp_nr_images ) SEP \
    ELTVE1(double	, exp_sum_weight					, 1				, exp_nr_images ) SEP \
    ELTVE1(double	, exp_significant_weight			, 1				, exp_nr_images ) SEP \
    ELTVE1(FDOUBLE	, exp_old_offsetx					, 1				, exp_nr_images	) SEP \
    ELTVE1(FDOUBLE	, exp_old_offsety					, 1				, exp_nr_images	) SEP \
    ELTVE1(FDOUBLE	, exp_wsum_norm_correction			, 1				, exp_nr_images	) SEP \
    ELTVE1(FDOUBLE	, exp_local_sqrtXi2					, 1				, exp_nr_images ) SEP \
    ELTVE2(FDOUBLE	, exp_local_Fctfs					, exp_nr_images	, current_Fsize2) SEP \
    ELTVE2(FDOUBLE	, exp_local_Minvsigma2s				, exp_nr_images	, current_Fsize2) SEP \
    ELTVE2(FDOUBLE	, exp_wsum_scale_correction_XA		, exp_nr_images	, ori_Fsize		) SEP \
    ELTVE2(FDOUBLE	, exp_wsum_scale_correction_AA		, exp_nr_images	, ori_Fsize		) SEP \
    ELTVE2(FDOUBLE	, exp_imgs							, exp_nr_images , ori_size2		) SEP \
    ELTVE2(FDOUBLE	, exp_power_imgs					, exp_nr_images	, ori_Fsize		) SEP \
	ELTVE2(char		, exp_Mcoarse_significant			, exp_nr_images	, nr_classes*sampler3d.NrPoints(0)) SEP \
    ELTVE2(char		, exp_Rot_significant				, nr_classes	, exp_nr_dir*exp_nr_psi) SEP \
    ELTVE3(FDOUBLE	, exp_Fimgs_shifted_real			, exp_nr_images*shifts, current_Fsize2) SEP \
    ELTVE3(FDOUBLE	, exp_Fimgs_shifted_imag			, exp_nr_images*shifts, current_Fsize2) SEP \
	/*Array with pointers to the resolution of each point in a Fourier-space FFTW-like array*/\
    ELTVE1(int		, Npix_per_shell					, 1				, ori_Fsize 	) SEP \
    ELTVE1(int		, Mresol_coarse						, 1				, current_Fsize2) SEP \
    ELTVE1(int		, Mresol_fine						, 1				, coarse_Fsize2 ) // end of MAPOPTIMIZER_EXP_VARS
    //
	//
    //
    //
    // ------------ 					thread variable 		-------------- //
#define MAPOPTIMIZER_THREAD_VARS \
	ELTVE2(double	, thread_exp_sum_weight				, maxthreads	, exp_nr_images	) SEP \
    ELTVE2(double	, thread_exp_max_weight				, maxthreads	, exp_nr_images	) SEP \
    ELTVE2(double	, thread_exp_min_diff2				, maxthreads	, exp_nr_images	) SEP \
    ELTVE2(FDOUBLE	, thread_sumw_group					, maxthreads	, exp_nr_images	) SEP \
    ELTVE2(FDOUBLE	, thread_wsum_norm_correction		, maxthreads	, exp_nr_images	) SEP \
    ELTVE2(FDOUBLE	, thread_wsum_sigma2_offset			, maxthreads	, exp_nr_images	) SEP \
    ELTVE2(FDOUBLE	, thread_wsum_pdf_class				, maxthreads	, nr_classes	) SEP \
    ELTVE2(FDOUBLE	, thread_wsum_prior_offsetx_class	, maxthreads	, nr_classes	) SEP \
    ELTVE2(FDOUBLE	, thread_wsum_prior_offsety_class	, maxthreads	, nr_classes	) SEP \
    ELTVE3(FDOUBLE	, thread_wsum_sigma2_noise			, maxthreads*exp_nr_images, current_Fsize2) SEP \
    ELTVE3(FDOUBLE	, thread_wsum_scale_correction_XA	, maxthreads*exp_nr_images, current_Fsize2) SEP \
    ELTVE3(FDOUBLE	, thread_wsum_scale_correction_AA	, maxthreads*exp_nr_images, current_Fsize2) // end of MAPOPTIMIZER_THREAD_VARS
    //
    //
    //
    //   ------------ track the change of hidden variable ---------------   //
    class HiddenVariableMonitor{
#define MAPOPTIMIZER_HIDDEN_VAR_ELT \
        /* Smallest changes thus far in the optimal translational offsets, orientations and classes */ \
        ELT(double	, smallest_changes_optimal_offsets			, 999		)			SEP \
        ELT(double	, smallest_changes_optimal_orientations		, 999		)			SEP \
        ELT(int		, smallest_changes_optimal_classes			, 9999999	)			SEP \
        /* Changes from one iteration to the next in the angles,translations,class assignments */ \
        ELT(double	, current_changes_optimal_offsets			, 999		)			SEP \
        ELT(double	, current_changes_optimal_orientations		, 999		)			SEP \
        ELT(double	, current_changes_optimal_classes			, 9999999	)			SEP \
        ELT(double	, sum_changes_optimal_offsets				, 0			)			SEP \
        ELT(double	, sum_changes_optimal_orientations			, 0			)			SEP \
        ELT(double	, sum_changes_optimal_classes				, 0			)			SEP \
        /* Just count how often the optimal changes are summed */							\
        ELT(double	, sum_changes_count							, 0			)			SEP \
        /* Number of iterations without a decrease in OffsetChanges */						\
        ELT(int		, nr_iter_wo_large_hidden_variable_changes	, 0			)		// end of HIDDEN_VAR_ELT
    public:
#define SEP
#define ELT(T,N,V) T N;
        MAPOPTIMIZER_HIDDEN_VAR_ELT
#undef ELT
#undef SEP
        HiddenVariableMonitor(){
#define SEP
#define ELT(T,N,V) N = V;
            MAPOPTIMIZER_HIDDEN_VAR_ELT
#undef ELT
#undef SEP
        }
        ~HiddenVariableMonitor(){}
        // Monitor the changes in the optimal translations, orientations and class assignments for some particles
        void monitorHiddenVariableChanges(HealpixSampler& sampling,const MetaDataTable& metadata_old,int first_image,
                                          const MetaDataElem* metadata_new,int nr_images);
        //
        void reduceData();
        //
        void bcastData(int bcastNode);
        // Updates the overall changes in the hidden variables and keeps track of nr_iter_wo_large_changes_in_hidden_variables
        void updateOverallChangesInHiddenVariables();
    };
    //
    static inline int getPixelFromResolution(double resol)	{ return (int)(resol * pixel_size * ori_size); }
    static inline double getResolution(int ipix)	{ return (double)ipix/(pixel_size * ori_size); }
    static inline double getResolutionAngstrom(int ipix)	{ return (ipix==0) ? 999. : (pixel_size * ori_size)/(double)ipix; }
    //
    //  ---------------   setup   TODO   ---------------------- //
    //
    void setupMLoptimizer();
    
    //
    void destroyMLoptimizer();
    
    // Interpret command line for the initial start of a run
    void prepare();
    
    // Perform expectation-maximization iterations
    void iterate();
    
    //
    void writeOutOptimizer(std::string fn_optimiser);
    //
    std::pair<std::string,std::string> readFromOptimizer(std::string fn_optimiser,bool debug_flag = false);
    //
    // ---------------------------- below is all useful data structure    ----------------------------- //
    //
    struct GridIndex{size_t iimage;size_t iclass;size_t idir;size_t ipsi;size_t iover_rot;size_t itrans;size_t iover_trans;};

    class IndexCoder{
    private:
        static const int max_index_num = 10;
        int index_num;
        int index_size[max_index_num];
    public:
        IndexCoder():index_num(0){
            for (int i = 0; i < max_index_num; i++) index_size[i] = 0;
        }
        ~IndexCoder(){}
        void init(std::initializer_list<int> _index_size);
        
        void fini();
        template<typename ... Types>// &index,&level,..index(inner to outer)...
        void decode(int index,int& first, Types& ... rest){
            int level = index_num;
            decodeLevel(index/index_size[level-1],level-1,rest...);
            first = index%index_size[level-1];
        }
        int code(std::initializer_list<int> _index);
        void unitTest();
    private:
        // NOTE : notice the peel order,from inner to outer
        template<typename ... Types>
        void decodeLevel(int index,int level,int &rest){
            assert(level==1);
            assert(index<index_size[0]);
            rest = index;
        }
        template<typename ... Types>
        void decodeLevel(int index,int level,int& first, Types& ... rest){
            decodeLevel(index/index_size[level-1],level-1,rest...);
            first = index%index_size[level-1];
        }
    };
    class IntDoublePair{
    private:
        int* first;
        double* second;
        int current_size;
        int total_size;
    public:
        IntDoublePair():first(NULL),second(NULL),current_size(0),total_size(0){}
        ~IntDoublePair(){}
        void init(int _size){
            current_size = 0;total_size = _size;
            assert(!first);first 	= (int*)	aMalloc(sizeof(int)		*total_size,64);
            assert(!second);second 	= (double*)	aMalloc(sizeof(double)	*total_size,64);
        }
        void fini(){
            total_size = current_size = 0;
            if (first) aFree(first);if (second) aFree(second);
        }
        inline void push_back(int _first,double _second){
            assert(current_size < total_size);
            first[current_size] = _first;second[current_size++] = _second;
        }
        inline void set_second(int _i,int _first,double _second){
            assert(_i < current_size);assert(first[_i]==_first);
            second[_i] = _second;
        }
        inline void add_second(int _i,int _first,double _second){
            assert(_i < current_size);assert(first[_i]==_first);
            second[_i] += _second;
        }
        inline int get_first(int _i){
            assert(_i < current_size);
            return first[_i];
        }
        inline double get_second(int _i){
            assert(_i < current_size);
            return second[_i];
        }
        void clear(){current_size = 0;}
        int size(){return current_size;}
    };

    //
    inline void setImagesDistribution(int& nr_global_images,int& nr_local_images,int nodes,int node,
                                      int& first_local_image,int& last_local_image,
                                      const MetaDataTable& metadata,bool do_split_random_halves)
    {
        if (do_split_random_halves)
        {
            if(node < nodes/2)// half one
            {
                nr_global_images = metadata.numberOfParticles(1);
                nr_local_images = divide_equally(nr_global_images, nodes/2, node, first_local_image, last_local_image);
            }
            else// half two
            {
                int nr_global_images_half_one = metadata.numberOfParticles(1);
                nr_global_images = metadata.numberOfParticles(2);
                nr_local_images = divide_equally(nr_global_images, nodes/2, node-nodes/2, first_local_image, last_local_image);
            }
        }
        else
        {
            nr_global_images = metadata.numberOfParticles();
            // NOTE : SAME_IMAGE_DIVISION uses to produce the same result for different nodes
#define SAME_IMAGE_DIVISION
#ifdef SAME_IMAGE_DIVISION
            //divided like not mpi version,just to easily debug
            nr_local_images = divide_equally_pool(nr_global_images, nr_pool, nodes, node, first_local_image, last_local_image);
#else
            nr_local_images = divide_equally(nr_global_images, nodes, node, first_local_image, last_local_image);
#endif
        }
    }
#ifdef USEMPI
    /*******  MPI function  *********/
    
    // reduce the metadata (ROT,TILT,PSI,XOFF,YOFF,ZOFF,CLASS,DLL,PMAX,NR_SIGN,NORM) to node0
    // or maybe donot reduce this,directly write metadata to disk
    void gatherMetaDataToMaster(MetaDataTable& metadata,std::vector<int> masterNodes = {0},bool do_split_random_halves = false);
    
#endif
    
}
#endif /* defined(MAP_OPTIMIZER_BASE_H_) */
