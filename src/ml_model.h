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

#ifndef _ROME_ML_MODEL
#define _ROME_ML_MODEL

#include "./sampling.h"
#include "./image.h"
#include "./map_model.h"
#include "./mpi.h"
#include "./string.h"
#include "./metadata.h"
#include "./array_vector.h"

class MLModel{
public:
    //(Datatype,  Name                          ,Value      ,Size                       ,Index)
#define MLMODEL_ELTS \
ELTONE(int      , ori_size                      , 0         , 1                         , 0   ) SEP \
ELTONE(int      , ori_Fsize                     , 0         , 1                         , 1   ) SEP \
ELTONE(int      , nr_classes                    , 0         , 1                         , 2   ) SEP \
ELTONE(int      , nr_groups                     , 0         , 1                         , 3   ) SEP \
ELTONE(int      , nr_directions                 , 1         , 1                         , 4   ) SEP \
ELTONE(FDOUBLE  , sigma2_offset                 , 3*3       , 1                         , 5   ) SEP \
ELTONE(FDOUBLE  , avg_norm_correction           , 1.0       , 1                         , 6   ) SEP \
ELTONE(FDOUBLE  , ave_Pmax                      , 0         , 1                         , 7   ) SEP \
ELTONE(FDOUBLE  , LL                            , 0         , 1                         , 8   ) SEP \
ELTONE(int      , orientational_prior_mode      , NOPRIOR   , 1                         , 9   ) SEP \
/* Variance angle for the orientational pdf */														\
ELTONE(FDOUBLE	, sigma2_rot					, 0			, 1							, 10  ) SEP \
ELTONE(FDOUBLE	, sigma2_tilt					, 0			, 1							, 11  ) SEP \
ELTONE(FDOUBLE	, sigma2_psi					, 0			, 1							, 12  ) SEP \
ELTVE2(FDOUBLE  , tau2_class                    , 0         , nr_classes*ori_Fsize      , 13  ) SEP \
ELTVE2(FDOUBLE  , sigma2_noise                  , 0         , nr_groups*ori_Fsize       , 14  ) SEP \
ELTVE2(FDOUBLE  , sigma2_class                  , 0         , nr_classes*ori_Fsize      , 15  ) SEP \
ELTVE2(FDOUBLE  , data_vs_prior_class           , 0         , nr_classes*ori_Fsize      , 16  ) SEP \
ELTVE2(FDOUBLE  , fsc_halves_class              , 0         , nr_classes*ori_Fsize      , 17  ) SEP \
ELTVE1(FDOUBLE  , scale_correction              , 0         , nr_groups                 , 18  ) SEP \
ELTVE1(FDOUBLE  , prior_offsetx_class           , 0         , nr_classes                , 19  ) SEP \
ELTVE1(FDOUBLE  , prior_offsety_class           , 0         , nr_classes                , 20  ) SEP \
ELTVE1(FDOUBLE  , pdf_class                     , 0         , nr_classes                , 21  ) SEP \
ELTVE1(FDOUBLE  , acc_rot                     	, 0         , nr_classes                , 22  ) SEP \
ELTVE1(FDOUBLE  , acc_trans                     , 0         , nr_classes                , 23  ) SEP \
ELTVE2(FDOUBLE  , pdf_direction                 , 0         , nr_classes*nr_directions  , 24  ) SEP \
ELTVE2(FDOUBLE	, orientability_contrib			, 0			, nr_classes*ori_Fsize		, 25  ) SEP \
ELTVE1(int      , nr_particles_group            , 0         , nr_groups                 , 26  ) SEP \
/*  For the refinement of group intensity scales and bfactors */                                    \
/*  each group store weighted sums of experimental image times reference image */                   \
/*    as a function of resolution */                                                                \
ELTVE2(FDOUBLE  , wsum_signal_product_spectra   , 0         , nr_groups*ori_Fsize       , 27  ) SEP \
/* each group store weighted sums of squared reference as a function of resolution */               \
ELTVE2(FDOUBLE  , wsum_reference_power_spectra  , 0         , nr_groups*ori_Fsize       , 28  ) SEP \
ELTVE1(FDOUBLE  , wsum_sumw_group               , 0         , nr_groups                 , 29  ) SEP \
ELTONE(FDOUBLE  , wsum_LL                       , 0         , 1                         , 30  ) SEP \
ELTONE(FDOUBLE  , wsum_ave_Pmax                 , 0         , 1                         , 31  ) SEP \
ELTONE(FDOUBLE  , wsum_avg_norm_correction      , 0         , 1                         , 32  ) SEP \
ELTONE(FDOUBLE  , wsum_sigma2_offset            , 0         , 1                         , 33  ) SEP \
ELTVE2(FDOUBLE  , wsum_sigma2_noise             , 0         , nr_groups*ori_Fsize       , 34  ) SEP \
ELTVE1(FDOUBLE  , wsum_pdf_class                , 0         , nr_classes                , 35  ) SEP \
ELTVE1(FDOUBLE  , wsum_prior_offsetx_class      , 0         , nr_classes                , 36  ) SEP \
ELTVE1(FDOUBLE  , wsum_prior_offsety_class      , 0         , nr_classes                , 37  ) SEP \
ELTVE2(FDOUBLE  , wsum_pdf_direction            , 0         , nr_classes*nr_directions  , 38  ) // end of macro
    
#define SEP
#define ELTONE(T,N,V,S,I) T N;
#define ELTVE1(T,N,V,S,I) VectorOfArray1d<T> N;
#define ELTVE2(T,N,V,S,I) VectorOfArray2d<T> N;
    MLMODEL_ELTS
#undef ELTONE
#undef ELTVE1
#undef ELTVE2
#undef SEP
    
    MLModel()
    {
#define SEP
#define ELTONE(T,N,V,S,I) N = V;
#define ELTVE1(T,N,V,S,I)
#define ELTVE2(T,N,V,S,I)
        MLMODEL_ELTS
#undef ELTONE
#undef ELTVE1
#undef ELTVE2
#undef SEP
    };
    ~MLModel(){};
    
    //
    void initialize(int _ori_size,int _nr_classes,int _nr_groups,int _nr_directions,double sigma2_angle = 0);
    
    //
    void finalize();
    
    //
    void resetZero();
    
    //
    template<typename T>
    void calculateSumOfPowerSpectraAndAverageImage(T& Mavg,std::vector<T>& images,bool do_zero_mask,
                                                   const MetaDataTable& metadata,int first_local_image,const MAPModel& mapModel)
    {
        Mavg.zero();
        if (images.size() == 0) return;
        
        int nr_local_images = images.size();
        int ori_size2 = ori_size*ori_size;
        int ori_Fsize2 = ori_size*ori_Fsize;
        
        // For spectrum calculation: recycle the transformer (so do not call getSpectrum all the time)
        FourierTransformer transformer(ori_size,ori_size);
        
        auto ind_spectrum	= Heap::allocScalars<FDOUBLE>(ori_size, __FILE__, __LINE__);
        auto count			= Heap::allocScalars<FDOUBLE>(ori_size, __FILE__, __LINE__);
        
        auto Faux_real = Heap::allocScalars<FDOUBLE>(ori_Fsize2, __FILE__, __LINE__);
        auto Faux_imag = Heap::allocScalars<FDOUBLE>(ori_Fsize2, __FILE__, __LINE__);
        
        Image img_temp(ori_size*ori_size);
        
        auto wsum_sigma2_noise_group0 = wsum_sigma2_noise[0].wptr(ori_Fsize);
        
        for (int iimage = 0; iimage < nr_local_images; iimage++)
        {
            int igroup = metadata[iimage+first_local_image].GROUP_NO-1;
            // Copy the image (only needed if do_zero_mask, but since that is true, not optimized)
            //
            img_temp = images[iimage];
            
            // Apply a similar softMask as below (assume zero translations)
            if (do_zero_mask){
                softMaskOutsideMap(img_temp.wptr(img_temp.size()), ori_size,
                                   mapModel.particle_diameter / (2. * mapModel.pixel_size),
                                   (double)mapModel.width_mask_edge,(FDOUBLE*)nullptr);
            }
            
            // Calculate the sum of the masked images
            //
            Mavg += img_temp;
            
            // Calculate the masked image's power spectrum and noise
            //
            transformer.FourierTransform(img_temp.wptr(img_temp.size()), Faux_real, Faux_imag);
            
            memset(ind_spectrum, 0, sizeof(FDOUBLE)*ori_size);
            memset(count, 0, sizeof(FDOUBLE)*ori_size);
            for (int i = 0; i<ori_size; i++)
                for (int j = 0; j<ori_Fsize; j++)
                {
                    int ip = (i < ori_Fsize) ? i : i - ori_size;
                    int jp = j;
                    
                    long int idx = round(sqrt(ip*ip + jp*jp));
                    ind_spectrum[idx] += Faux_real[i*ori_Fsize+j]*Faux_real[i*ori_Fsize+j] + Faux_imag[i*ori_Fsize+j]*Faux_imag[i*ori_Fsize+j];
                    count[idx] += 1.;
                }
            
            assert(ori_Fsize <= ori_size);
            for (int i = 0; i < ori_Fsize; i++) {
                CHECK_NOT_NAN_OR_INF(ind_spectrum[i] /= count[i]);
                wsum_sigma2_noise_group0[i] += CHECK_NOT_NAN_OR_INF(ind_spectrum[i]);
            }
            
            for (int i = ori_Fsize; i < ori_size; i++) {
                CHECK_NOT_NAN_OR_INF(ind_spectrum[i] /= count[i]);
            }
            
            wsum_sumw_group.wptrAll()[0] += 1.;
            nr_particles_group.wptrAll()[igroup] += 1;
            
        } // end loop iimage
        
        Heap::freeScalars(ind_spectrum);
        Heap::freeScalars(count);
        Heap::freeScalars(Faux_real);
        Heap::freeScalars(Faux_imag);
    }
    
    //
    void setSigmaNoiseEstimates(Image& Iref_avg,bool do_calculate_initial_sigma_noise = true);
    
    //
    void initialiseDataVersusPrior(MAPModel& myMapModel,double tau2_fudge_factor);
    
    //
    void initialisePdfDirection(int newsize);
    
    //
    void writeOutModel(std::string fn_model,const MAPModel& mapModel,const MetaDataTable& metadata,int random_subset = -1);
    
    //
    void writeOutBild(std::string fn_class, MAPModel& mapModel, HealpixSampler& sampling);
    
    //
    void readFromModel(std::string fn_model,MAPModel& mapModel,bool debug_flag = false);
    
    // print needed space for MLModel data
    void printSpaceInfo(std::ostream &os);
    
    // write MLModel data to disk
    // working with -continue
    void writeToDisk(std::string fn_ml_model);
    
    // read MLModel data from disk
    // working with -continue
    void readFromDisk(std::string fn_ml_model);
    
    // reduce MLModel data for all nodes
    void reduceData(const MPI::Intracomm& currentWorld);
    
    //
    void broadcastData();
    
    // it may return all zero if using float precision
    inline double calculatePdfOffset(FDOUBLE offset_x,FDOUBLE offset_y, FDOUBLE prior_x,FDOUBLE prior_y)
    {
        double tdiff2 = (offset_x-prior_x)*(offset_x-prior_x)+(offset_y-prior_y)*(offset_y-prior_y);
        if (sigma2_offset < 0.0001)
        {
            return (tdiff2 > 0.) ? 0. : 1.;
        }
        else
        {
            return exp ( tdiff2 / (-2. * sigma2_offset) ) / ( 2. * rome_pi * sigma2_offset );
        }
    }
    
};

#endif /* defined(_ROME_ML_MODEL) */
