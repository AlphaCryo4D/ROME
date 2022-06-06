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

#ifndef MAP_MODEL_H_
#define MAP_MODEL_H_

#include "./resmap_reconstruct.h"
#include "./resmap_mpi.h"
#include "./resmap_macros.h"
#include "./resmap_time.h"
#include "./resmap_string.h"
#include "./resmap_array_vector.h"

class MAPModel{
    
#define DEBUG_MAP_MODEL

public:
     //(Datatype,  Name                          ,Value      ,Size                       ,Index)
#define MAPMODEL_ELTS \
    ELTONE(int      , ori_size                      , 0         , 1                         , 0   ) SEP \
    ELTONE(int      , ori_Fsize                     , 0         , 1                         , 1   ) SEP \
    ELTONE(int      , ref_dim                     	, 0         , 1                         , 2   ) SEP \
    ELTONE(int      , nr_classes                    , 0         , 1                         , 3   ) SEP \
    ELTONE(double   , pixel_size                    , 0         , 1                         , 4   ) SEP \
    ELTONE(double   , particle_diameter             , 0         , 1                         , 5   ) SEP \
    ELTONE(double   , ini_high                    	, 0         , 1                         , 6   ) SEP \
    ELTONE(double   , current_resolution            , 0         , 1                         , 7   ) SEP \
    ELTONE(double   , width_fmask_edge            	, 0         , 1                         , 8   ) SEP \
    /* Width of the soft-edges of the circular masks */													\
    ELTONE(double   , width_mask_edge	            , 0         , 1                         , 9   ) SEP \
    /* Minimum resolution to perform Bayesian estimates of the model */									\
    ELTONE(int	    , minres_map            		, 0         , 1                         , 10  ) SEP \
    /* How many fourier shells should be included beyond the highest shell where evidenceVsPriorRatio < 1? */ \
    ELTONE(int 		, incr_size						, 10		, 1							, 11  )	SEP \
    /* Flag whether to use maximum a posteriori (MAP) estimation */										\
    ELTONE(bool	    , do_map            			, 0         , 1                         , 12  ) SEP \
    ELTONE(bool	    , do_gridding            		, 0         , 1                         , 13  ) SEP \
    /* Original reference */																			\
    ELTVEC(Vol<FDOUBLE>, Irefs            			, 0         , 1                         , 14  ) SEP \
    /* Symmetry group */																				\
    ELTONE(std::string, fn_sym            			, "C1"      , 1                         , 15  ) SEP \
    /* projector and back-projector for references(2D or 3D) */                                         \
    ELTVEC(MyProjector, projector            		, 0         , 1                         , 16  ) SEP \
    ELTVEC(MyBackProjector, backprojector           , 0         , 1                         , 17  ) SEP \
	ELTONE(bool	    , reset_nrthreads_save_memory	, 0         , 1                         , 18  ) SEP \
    ELTONE(int	    , nr_threads            		, 0         , 1                         , 19  ) SEP \
    ELTVEC(Vol<MKL_Complex>, thread_data          	, 0         , 1                         , 20  ) SEP \
    ELTVEC(Vol<FDOUBLE>, thread_weight            	, 0         , 1                         , 21  ) // end of macro

#define SEP
#define ELTONE(T,N,V,S,I) T N;
#define ELTVEC(T,N,V,S,I) std::vector<T> N;
    MAPMODEL_ELTS
#undef ELTONE
#undef ELTVEC
#undef SEP
    //
    Mrcs::MrcsHead refHead;
    
    MAPModel(){
#define SEP
#define ELTONE(T,N,V,S,I) N = V;
#define ELTVEC(T,N,V,S,I)
        MAPMODEL_ELTS
#undef ELTONE
#undef ELTVEC
#undef SEP
    }
    ~MAPModel(){
        finalize();
    }
    
    void finalize(){
#define SEP
#define ELTONE(T,N,V,S,I) N = V;
#define ELTVEC(T,N,V,S,I) N.resize(0);
        MAPMODEL_ELTS
#undef ELTONE
#undef ELTVEC
#undef SEP
    }
    
    static inline double getResolution(int ipix,double pixel_size,double ori_size)	{ return (double)ipix/(pixel_size * ori_size); }
    // maybe cause some round error,like resol * pixel_size * ori_size=4.99999999,be careful
    static inline int getPixelFromResolution(double resol,double pixel_size,double ori_size)	{ return (int)(resol * pixel_size * ori_size+0.0001); }
    
    //
    void initialize(int _nr_classes,int _ori_size,double _particle_diameter,
                    double _pixel_size,double _ini_high,int _nr_threads,
                    double _width_mask_edge = 5,double _width_fmask_edge = 2,
                    int _minres_map = 5,bool _do_gridding = true,bool _do_map = true,
                    int _ref_dim = 3,std::string _fn_sym = "C1");
    
    // initialize reference from file
    void initializeRef(std::string ref_fn);
    void initializeRef(const std::vector<std::string>& ref_fns);
    
    // initialize reference from average image
    void initializeRef(Image& averageImage);
    
    //
    void setFourierTransformMaps(std::vector<VectorOfFDOUBLE>& power_spectrum,int current_size,bool reset_nrthreads_save_memory = false);
    
    // get 2D slice form 3D map
    void get2DFourierTransform(int iclass,FDOUBLE* img_out_real, FDOUBLE* img_out_imag,int img_out_size,FDOUBLE A[][3], bool inv);
    
    // get 2D slice form 3D map
    void get2DFourierTransformOneTile(int iclass,FDOUBLE* img_out_real, FDOUBLE* img_out_imag,
                                      int n_start,int n_end,int img_out_size,const FDOUBLE A[][3],
                                      const int* nIndex = nullptr);
    
    // set 2D slice to 3D map
    void set2DFourierTransform(int thread,int iclass,const FDOUBLE* img_in_real,const FDOUBLE* img_in_imag,int img_in_size,
                               const FDOUBLE A[][3], bool inv,const FDOUBLE* Mweight = NULL);
    
    // set 2D slice to 3D map
    void set2DFourierTransformOneTile(int thread,const FDOUBLE* img_in_real,const FDOUBLE* img_in_imag,int n_start,
                                      int n_end,int img_in_size,const FDOUBLE A[][3],const FDOUBLE* Mweight,const int* nIndex);
    
    //
    void reduceThreadMap(int iclass);
    
    //
    void reconstruction(int iclass,int max_iter_preweight,bool do_map,double tau2_fudge,FDOUBLE* tau2,FDOUBLE* sigma2,FDOUBLE* data_vs_prior,
                        const FDOUBLE* fsc /* only input*/, double normalise = 1., bool update_tau2_with_fsc = false,
                        bool is_whole_instead_of_half = false,int _nr_threads = 1);
    
    //
    void writeResult(std::string filename);
    
    // apply low-pass filter to all map
    void applyLowPassFilter();
    
    // Apply a solvent flattening to all map
    void applySolventFlatten(std::string fn_mask);
    
    // Updates the current resolution (from data_vs_prior array) and keeps track of best resolution thus far
    // and precalculates a 2D Fourier-space array with pointers to the resolution of each point in a FFTW-centered array
    void updateCurrentResolution(std::vector<VectorOfFDOUBLE>& data_vs_prior_class,bool set_by_ini_high,int& nr_iter_wo_resol_gain,FDOUBLE& best_resol_thus_far);
    
    //
    void printResolution(std::vector<VectorOfFDOUBLE>& data_vs_prior_class,bool set_by_ini_high,bool debug_flag = false);
    
    // Update the current and coarse image size
    // and precalculates a 2D Fourier-space array with pointers to the resolution of each point in a FFTW-centered array
    void updateImageSizeAndResolutionPointers(VectorOfInt& Npix_per_shell,VectorOfInt& Mresol_coarse,VectorOfInt& Mresol_fine,
                                              int& coarse_size,int& current_size,
                                              int adaptive_oversampling,double angularSampling,FDOUBLE model_ave_Pmax,
                                              bool has_high_fsc_at_limit = false,bool do_use_all_data = false,bool debug_flag = false);
    
    // reduce Fourier-space sampling data to each node
    void reduceData(const MPI::Intracomm& currentWorld,bool use_allreduce = true);
    
    // broadcast Volume data to each nodes
    void broadcastData();
    
    //
    void printSpaceInfo(std::ostream &os);
};

#endif /* defined(MAP_MODEL_H_) */
