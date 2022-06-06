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

#ifndef EXP_MODEL_H_
#define EXP_MODEL_H_


#include "./ctf.h"
#include "./debug.h"
#include "./fft_fftw3.h"
#include "./image.h"
#include "./map_model.h"
#include "./metadata.h"
#include "./mrcs.h"
#include "./ml_model.h"
#include "./sampling.h"

// 4~10% 10e-6 error if turn on
// so default use single translation for comparison
//#define DOUBLE_TRANSLATION
//#define TRIPLE_TRANSLATION

// rearrange the image array to resolution increase array....
class FourierShellTranslation : public NoCopy {
	int* nIndexFine;// old index for each new pixel point
	int* iResFine;// resolution for each new pixel points
	std::vector<double*> TempDataThread;
	int oldFsize2Fine, newFsize2Fine;
	// coase image
	int* nIndexCoarse;
	int oldFsize2Coarse, newFsize2Coarse;// assert(coarseSize <= fineSize);
										 // norm correction
	int min_iRes, max_iRes;
	int normCorrLo, normCorrHi;
	// scale correction
	std::vector<int> scaleCorrThresholdRes;
	std::vector<int> scaleCorrLo;
	std::vector<int> scaleCorrHi;
    std::vector<int*> scaleCorrExcludingIndex;// If data_va_prior is not monotone decreasing
public:
	FourierShellTranslation() :nIndexFine(nullptr), iResFine(nullptr), oldFsize2Fine(0), newFsize2Fine(0),
		nIndexCoarse(nullptr), oldFsize2Coarse(0), newFsize2Coarse(0),
		min_iRes(0), max_iRes(0), normCorrLo(0), normCorrHi(0),
		scaleCorrThresholdRes(0), scaleCorrLo(0), scaleCorrHi(0), scaleCorrExcludingIndex(0) {}
	~FourierShellTranslation() { clear(); }
	//
	void checksum(Checksum<int> & checksumCoarse, Checksum<int> & checksumFine) const {
		checksumCoarse.init(nIndexCoarse, newFsize2Coarse);
		checksumFine.init(nIndexFine, newFsize2Fine);
	}
	//
	void setCoarseImageBoundary(const int* _Mresol_coarse, int _size);
	//
	void setNormCorrectionBoundary(const int* _Mresol_fine, int _size);
	inline int getNormCorrLo() { return normCorrLo; }
	inline int getNormCorrHi() { return normCorrHi; }
	//
	void setScaleCorrectionBoundary(const VectorOfArray2d<FDOUBLE>& data_vs_prior_class);
	inline int getScaleCorrLo(int iclass) { return scaleCorrLo[iclass]; }
	inline int getScaleCorrHi(int iclass) { return scaleCorrHi[iclass]; }
	//
	void setTransformFineAndCoarse(int nr_threads) {
		TempDataThread.resize(nr_threads);
		for (int thread = 0; thread < nr_threads; thread++) { TempDataThread[thread] = (double*)aMalloc(sizeof(double)*newFsize2Fine, 64); }
	}

	static IntPerformanceCounter transformFinePerformanceCounter;
	static IntPerformanceCounter transformCoarsePerformanceCounter;

	void transformFine(int tid, double* data, int size) const {
		transformFinePerformanceCounter.count.v++;
		assert(oldFsize2Fine == size);
		for (int i = 0; i < newFsize2Fine; i++) { (TempDataThread[tid])[i] = data[nIndexFine[i]]; }
		for (int i = 0; i < newFsize2Fine; i++) { data[i] = (TempDataThread[tid])[i]; }
	}
	void transformBackFine(int tid, double* data, int size) const {
		assert(oldFsize2Fine == size);
		for (int i = 0; i < newFsize2Fine; i++) { (TempDataThread[tid])[i] = data[i]; }
		for (int i = 0; i < newFsize2Fine; i++) { data[nIndexFine[i]] = (TempDataThread[tid])[i]; }
	}

	void transformCoarse(int tid, double* data, int size) const {
		transformCoarsePerformanceCounter.count.v++;
		assert(oldFsize2Coarse == size);
		for (int i = 0; i < newFsize2Coarse; i++) { (TempDataThread[tid])[i] = data[nIndexCoarse[i]]; }
		for (int i = 0; i < newFsize2Coarse; i++) { data[i] = (TempDataThread[tid])[i]; }
	}
	//
	void clear();
	inline int* rptr_nIndexFine() const { return nIndexFine; }
	inline int* rptr_nIndexCoarse() const { return nIndexCoarse; }
	inline int fineFsize2() { return newFsize2Fine; }
	inline int coarseFsize2() { return newFsize2Coarse; }
};


class ParticleModel
{
public:
     	//(Datatype,  Name                          ,Value      ,Size                       )
#define SIZE_ELTS \
    ELTONE(int      , ori_size                      , 0         , 1                         ) SEP \
    ELTONE(int      , ori_Fsize                     , 0         , 1                   		) SEP \
    ELTONE(int      , ori_Fsize2                    , 0         , 1                         ) SEP \
    ELTONE(int      , ori_size2                     , 0         , 1                         ) SEP \
    ELTONE(int      , coarse_size                   , 0         , 1                         ) SEP \
    ELTONE(int      , coarse_Fsize                  , 0         , 1                         ) SEP \
    ELTONE(int      , coarse_Fsize2                 , 0         , 1                         ) SEP \
    ELTONE(int      , coarse_size2                  , 0         , 1                         ) SEP \
    ELTONE(int      , fine_size                  	, 0         , 1                         ) SEP \
    ELTONE(int      , fine_Fsize                  	, 0         , 1                         ) SEP \
    ELTONE(int      , fine_Fsize2                  	, 0         , 1                         ) SEP \
    ELTONE(int      , fine_size2                  	, 0         , 1                         ) SEP \
    ELTONE(int      , exp_nr_images                 , 0         , 1                         ) // end of SIZE_ELTS

#define CONFIG_ELTS \
    ELTONE(double   , pixel_size            		, 0         , 1                         ) SEP \
    /* Particle diameter (in Ang) */													   		  \
    ELTONE(double   , particle_diameter            	, 0         , 1                         ) SEP \
    ELTONE(int   	, width_mask_edge            	, 0         , 1                         ) SEP \
    ELTONE(double   , sigma2_fudge            		, 0         , 1                         ) SEP \
    ELTONE(int   	, random_seed            		, 0         , 1                         ) SEP \
    ELTONE(bool     , do_norm_correction            , 0         , 1                         ) SEP \
    ELTONE(bool     , do_zero_mask            		, 0         , 1                         ) SEP \
    /* Pnly perform phase-flipping CTF correction */											  \
    ELTONE(bool     , only_flip_phases            	, 0         , 1                         ) SEP \
	/* Images have been CTF phase flipped al;ready */											  \
    ELTONE(bool     , ctf_phase_flipped            	, 0         , 1                         ) SEP \
    ELTONE(bool     , intact_ctf_first_peak         , 0         , 1                         ) // end of CONFIG_ELTS

#define IMAGE_ELTS \
    ELTONE(int   	, exp_first_image            	, 0         , 1                         ) SEP \
    ELTVEC(MetaDataElem, exp_metadata            	, 0         , exp_nr_images             ) SEP \
    ELTVE2(FDOUBLE  , exp_images            		, 0         , exp_nr_images*ori_size2   ) SEP \
    ELTVE2(FDOUBLE  , Fctfs            				, 0         , exp_nr_images*fine_Fsize2 ) SEP \
    ELTVE2(FDOUBLE  , FctfsOne            			, 0         , exp_nr_images*fine_Fsize2 ) SEP \
    /* masked image,fine search size (current_size),for calculate weight */ 					  \
    ELTVE2(FDOUBLE  , Fimages_mask_fine_real        , 0         , exp_nr_images*fine_Fsize2 ) SEP \
    ELTVE2(FDOUBLE  , Fimages_mask_fine_imag        , 0         , exp_nr_images*fine_Fsize2 ) SEP \
    /* masked image,coarse search size (coarse_size),for calculate weight*/ 					  \
    ELTVE2(FDOUBLE  , Fimages_mask_coarse_real      , 0         , exp_nr_images*coarse_Fsize2) SEP \
    ELTVE2(FDOUBLE  , Fimages_mask_coarse_imag      , 0         , exp_nr_images*coarse_Fsize2) SEP \
    /* un-mask image,for backprojection */ 														  \
    ELTVE2(FDOUBLE  , Fimages_nomask_real      		, 0         , exp_nr_images*fine_Fsize2 ) SEP \
    ELTVE2(FDOUBLE  , Fimages_nomask_imag      		, 0         , exp_nr_images*fine_Fsize2 ) SEP \
    ELTONE(int   	, nr_threads            		, 0         , 1                         ) SEP \
    ELTVE2(FDOUBLE  , threadFimages_real            , 0         , nr_threads*ori_size2      ) SEP \
    ELTVE2(FDOUBLE  , threadFimages_imag            , 0         , nr_threads*ori_size2      ) SEP \
    ELTVE2(FDOUBLE  , threadImages            		, 0         , nr_threads*ori_size2      ) // end of IMAGE_ELTS
    
    //
    DataStream* global_data_stream;
#if defined(FLOAT_PRECISION)
    std::vector< FFTWFTransformer* > threadFFTtransformer;
#else
    std::vector< FFTWTransformer* > threadFFTtransformer;
#endif
    std::vector< CTF*			  > threadCTFer;
    bool do_shifts_onthefly;

//#define DEBUG_SHIFTER_POOL
	struct State {
		enum Status {
			uninited,
			initTableDone,
			initDone,
			transformDone
		} status;
		State() { init(); }
		void init() { 
			// Set all so State can be compared
			status = uninited; 
#if defined(DEBUG_SHIFTER_POOL)
			debugTransformed = false;
#endif
		}
		// initTableDone
		int fine_Fsize2;
#if defined(DEBUG_SHIFTER_POOL)
		bool debugTransformed;
#endif
		// initDone
		int current_size; double shiftx; double shifty; int ori_size;
		// transformDone
		FourierShellTranslation const * fourierShellTrans;
		int fourierShellTrans_uid;
		bool transformToCoarse;
		int exp_current_Fsize2;
#if defined(DEBUG_SHIFTER_POOL)
		Checksum<double> aTable_checksum_before, aTable_checksum_after;
		Checksum<int>    fourierShellTrans_checksum_coarse, fourierShellTrans_checksum_fine;
#endif

		bool operator==(State const & rhs) const {
#define EQL(N) (N == rhs.N)
			if (!EQL(status)) return false;
			if (status < initTableDone) return true;
			if (!EQL(fine_Fsize2)) return false;
			if (status < initDone) return true;
			if (!( EQL(current_size)
				&& EQL(shiftx)
				&& EQL(shifty)
				&& EQL(ori_size))) return false;
			if (status < transformDone) return true;
			if (!( EQL(fourierShellTrans)
				&& (fourierShellTrans_uid == fourierShellTrans->uid())
				&& EQL(transformToCoarse)
				&& EQL(exp_current_Fsize2)
				)) return false;
			return true;
#undef EQL
		}
	};

	class ShifterPool;
	class ShiftImageAssistor;

	class Shifter : public ShiftImageInFourierTransformNew<FDOUBLE, FDOUBLE> {
	public:
		Shifter() : index(-1), acquiredCount(0), freeNext(nullptr), freePrev(nullptr) {}
		virtual ~Shifter() {}
	protected:
		Shifter*	freeNext;
		Shifter*	freePrev;
		int			index;
		int			acquiredCount;
		State		currentState;
		friend class ParticleModel::ShifterPool;
		friend class ParticleModel::ShiftImageAssistor;
		friend void doTransformNow(
			const ParticleModel::State & desiredState,
			ParticleModel::Shifter* shifter);
		bool acquire(State const & state, int index);
		void release();
	};

	class ShiftImageAssistor : public NoCopy {
		//
		// We are working towards having this table 
		// be a subset of just the tables that are needed and fit in the L2 cache
		// recomputing the ones that are needed rather than polluting the cache with unneeded ones
		//
		int				         statesCapacity;
		int				         statesSize;
		State*                   desiredStates;
		ShifterPool*	  const  shifterPool;

		int max_dim0, max_dim1, dim0, dim1;
	public:
		ShiftImageAssistor(ParticleModel & parent);
		~ShiftImageAssistor();
		void setMaxDims(int max_exp_nr_trans, int max_exp_nr_over_trans);
		void initTable (int fine_Fsize2);
		void setCurrDims(int exp_nr_trans, int exp_nr_over_trans);
		void setInitArgs(int itrans, int iover_trans, int itrans_over_trans, int current_size, double shiftx, double shifty, int ori_size);
		void transform(int itrans, int iover_trans, int itrans_over_trans, const int exp_current_Fsize2, FourierShellTranslation const & fourierShellTrans, bool useCoarse);
		Shifter const * acquireShifter (int itrans, int iover_trans, int itrans_over_trans);
		void releaseShifter(Shifter const * & shifter);
	private:
		Shifter       * acquireShifterW(int itrans, int iover_trans, int itrans_over_trans);
		void releaseShifter(Shifter       * & shifter);

		int index(int itrans, int iover_trans, int itrans_over_trans) const {
			if (itrans < 0 || dim0 <= itrans || iover_trans < 0 || dim1 <= iover_trans)
				EXIT_ABNORMALLY;
			auto result = itrans*dim1 + iover_trans;
			if (result != itrans_over_trans) {
				std::cerr << "ShiftImageAssistor::index problem "
					<< " itrans:" << itrans
					<< " dim1:" << dim1
					<< " iover_trans:" << iover_trans
					<< " itrans_over_trans:" << itrans_over_trans
					<< std::endl;
				EXIT_ABNORMALLY;
			}
			return result;
		}
	} shiftImageAssistor;

    ShiftImageAssistor largeShiftedABTable;
    ShiftImageAssistor smallShiftedABTable;
    ShiftImageAssistor largeXShiftedABTable;
    ShiftImageAssistor largeYShiftedABTable;
    
	void shiftAssistorsSetMaxDims(int max_exp_nr_trans, int max_exp_nr_over_trans) {
		shiftImageAssistor       .setMaxDims(max_exp_nr_trans, max_exp_nr_over_trans);
#if defined(DOUBLE_TRANSLATION) || defined(TRIPLE_TRANSLATION)
		largeShiftedABTable      .setMaxDims(max_exp_nr_trans,                     1);
		smallShiftedABTable      .setMaxDims(1,                max_exp_nr_over_trans);
		largeXShiftedABTable     .setMaxDims(max_exp_nr_trans,                     1);
		largeYShiftedABTable     .setMaxDims(max_exp_nr_trans,                     1);
#endif
	}

	void shiftAssistorsSetCurrDims(int exp_nr_trans, int exp_nr_over_trans) {
		shiftImageAssistor		 .setCurrDims(exp_nr_trans, exp_nr_over_trans);
#if defined(DOUBLE_TRANSLATION) || defined(TRIPLE_TRANSLATION)
		largeShiftedABTable		 .setCurrDims(exp_nr_trans,					1);
		smallShiftedABTable		 .setCurrDims(1,			exp_nr_over_trans);
		largeXShiftedABTable	 .setCurrDims(exp_nr_trans,					1);
		largeYShiftedABTable	 .setCurrDims(exp_nr_trans,					1);
#endif
	}

	void shiftAssistorsInitTable(int fine_Fsize2) {
		shiftImageAssistor.initTable(fine_Fsize2);
#if defined(DOUBLE_TRANSLATION) || defined(TRIPLE_TRANSLATION)
		largeShiftedABTable .initTable(fine_Fsize2);
		largeXShiftedABTable.initTable(fine_Fsize2);
		largeYShiftedABTable.initTable(fine_Fsize2);
		smallShiftedABTable .initTable(fine_Fsize2);
#endif
	}

#define SEP
#define ELTONE(T,N,V,S) T N;
#define ELTVEC(T,N,V,S) VectorOfStruct<T> N;
#define ELTVE1(T,N,V,S) VectorOfArray1d<T> N;
#define ELTVE2(T,N,V,S) VectorOfArray2d<T> N;
    SIZE_ELTS
    CONFIG_ELTS
    IMAGE_ELTS
#undef ELTVE2
#undef ELTVE1
#undef ELTVEC
#undef ELTONE
#undef SEP
    
    ParticleModel() : 
		shiftImageAssistor(*this),
		largeShiftedABTable(*this),
		smallShiftedABTable(*this),
		largeXShiftedABTable(*this),
		largeYShiftedABTable(*this)
	{
#define ELTONE(T,N,V,S) N = V;
#define ELTVEC(T,N,V,S)
#define ELTVE1(T,N,V,S)
#define ELTVE2(T,N,V,S)
#define SEP
        SIZE_ELTS
        CONFIG_ELTS
        IMAGE_ELTS
#undef ELTVE2
#undef ELTVE1
#undef ELTVEC
#undef ELTONE
#undef SEP
    }
    ~ParticleModel(){}
    //
    void initialize(int _ori_size,double _pixel_size,double _particle_diameter,int _width_mask_edge,
                    double _sigma2_fudge,int _random_seed,bool _do_norm_correction,bool _do_zero_mask,
                    bool _do_shifts_onthefly,int _nr_threads,DataStream* global_data_stream);
    //
    void finalize();
    
    // Set up the Images Transformer
    void setup(int _nr_pool,int _current_size,int _coarse_size,int exp_nr_trans = 0,int exp_nr_over_trans = 0);
    
    // destroy the Images Transformer
    void destroy();
    
    // set up the expectation step image data
    // NOTE : the exp_nr_images maybe small than nr_images in final search pool
    template<typename T>
    void prepare(std::vector<T>& _exp_images,VectorOfStruct<MetaDataElem>& _exp_metadata,int _exp_first_image,int _exp_nr_images){
        exp_nr_images = _exp_nr_images;
        exp_first_image = _exp_first_image;
#pragma omp parallel for
        for (int iimage = 0; iimage < exp_nr_images; iimage++) {
            ::copy(exp_images[iimage].wptrAll(), ori_size2, _exp_images[iimage].rptrAll(), ori_size2);
            exp_metadata[iimage] = _exp_metadata[iimage];
        }
    }
    // get Image(with mask and nomask) and CTF in Fourier Space
    // it will downsize the image from ori_size to current_size
    // Set the mask and nomask image
    template<typename T1,typename T2>
    void setFourierTransforms(std::vector<T1>& exp_power_imgs,VectorOfScalar<T2>& exp_highres_Xi2_imgs,
                              VectorOfScalar<T2>& exp_old_offsetx,VectorOfScalar<T2>& exp_old_offsety,
                              MLModel& mlModel,bool do_ctf_correction)
    {
#pragma omp parallel for
        for (int iimage = 0; iimage < exp_nr_images; iimage++)
        {
#ifdef DATA_STREAM
            global_data_stream->foutInt(iimage, "setFourierTransforms()_start------", __FILE__, __LINE__);
#endif
            int tid = omp_get_thread_num();
            
            int igroup = exp_metadata[iimage].GROUP_NO-1;
            // thread temporary data
            auto tempImage = threadImages[tid].wptr(ori_size2);
            auto tempFimage_real = threadFimages_real[tid].wptr(ori_size2);
            auto tempFimage_imag = threadFimages_imag[tid].wptr(ori_size2);
            
            // We never need any resolutions higher than current_size So resize the Fourier transforms
            // to current_size
            auto exp_Fimage_nomask_real = Fimages_nomask_real[iimage].wptr(fine_Fsize2);
            auto exp_Fimage_nomask_imag = Fimages_nomask_imag[iimage].wptr(fine_Fsize2);
            auto exp_Fimage_mask_fine_real = Fimages_mask_fine_real[iimage].wptr(fine_Fsize2);
            auto exp_Fimage_mask_fine_imag = Fimages_mask_fine_imag[iimage].wptr(fine_Fsize2);
            auto exp_Fimage_mask_coarse_real = Fimages_mask_coarse_real[iimage].wptr(coarse_Fsize2);
            auto exp_Fimage_mask_coarse_imag = Fimages_mask_coarse_imag[iimage].wptr(coarse_Fsize2);
            auto exp_image = exp_images[iimage].wptr(ori_size2);
            
            for (int i = 0; i < ori_size2; i++)
                tempImage[i] = exp_image[i];
            
#ifdef DATA_STREAM
            global_data_stream->foutDouble(tempImage, ori_size2, "getFourierTransformsAndCtfs()_original", __FILE__, __LINE__);
            global_data_stream->check();global_data_stream->flush();
#endif
            
            // Apply the norm_correction term
            if (do_norm_correction)
            {
                // Get the norm_correction
                double normcorr = exp_metadata[iimage].NORM;
                
                for (int i = 0; i < ori_size2; i++)
                    tempImage[i] = tempImage[i] * mlModel.avg_norm_correction / normcorr;
            }
            
            CHECK_NOT_IND(tempImage[0]);
            if (tempImage[0] * tempImage[0] < 0) {          // not true for negative indefinites
                std::cerr << "exp_metadata[iimage].NORM:" << exp_metadata[iimage].NORM << " model_avg_norm_correction:" << mlModel.avg_norm_correction << std::endl;
                assert(tempImage[0] * tempImage[0] >= 0);   // not true for negative indefinites
            }
            
            // Apply (rounded) old offsets first then translate the image
            exp_old_offsetx[iimage] = round(exp_metadata[iimage].XOFF);
            exp_old_offsety[iimage] = round(exp_metadata[iimage].YOFF);
            translate(tempImage, exp_image, ori_size, exp_old_offsetx[iimage],exp_old_offsety[iimage],DONT_WRAP);
            CHECK_NOT_IND(exp_image[0]);
            
#ifdef DATA_STREAM
            global_data_stream->foutDouble(tempImage, ori_size2, "getFourierTransformsAndCtfs()_norm_img", __FILE__, __LINE__);
            global_data_stream->foutDouble(exp_metadata[iimage].XOFF, "getFourierTransformsAndCtfs()_my_old_offsetx", __FILE__, __LINE__);
            global_data_stream->foutDouble(exp_metadata[iimage].YOFF, "getFourierTransformsAndCtfs()_my_old_offsety", __FILE__, __LINE__);
            global_data_stream->foutDouble(exp_old_offsetx[iimage], "getFourierTransformsAndCtfs()_my_old_offsetx_round", __FILE__, __LINE__);
            global_data_stream->foutDouble(exp_old_offsety[iimage], "getFourierTransformsAndCtfs()_my_old_offsety_round", __FILE__, __LINE__);
            global_data_stream->foutDouble(exp_image, ori_size2, "getFourierTransformsAndCtfs()_translate", __FILE__, __LINE__);
            global_data_stream->check();global_data_stream->flush();
#endif
            
            // Inside Projector and Backprojector the origin of the Fourier Transform is centered!
            CenterFFT(exp_image,tempImage,ori_size,true);
            CHECK_NOT_IND(tempImage[0]);
            
            threadFFTtransformer[tid]->FourierTransform(tempImage, tempFimage_real,tempFimage_imag);
            CHECK_NOT_IND(tempFimage_real[0]);CHECK_NOT_IND(tempFimage_imag[0]);
            
            // get the Fourier Transform of the image Fimg nomask
            // We never need any resolutions higher than current_size
            // so resize the Fourier transforms to current_size
            //
            // BEVIN This is often an unnecessary copy
            windowFourierTransform(tempFimage_real, tempFimage_imag, ori_size,
                                   exp_Fimage_nomask_real, exp_Fimage_nomask_imag, fine_size);
            
#ifdef DATA_STREAM
            global_data_stream->foutDouble(tempImage, ori_size2, "getFourierTransformsAndCtfs()_centerFFT", __FILE__, __LINE__);
            global_data_stream->foutDouble(tempFimage_real, ori_Fsize2, "getFourierTransformsAndCtfs()_FFT_real", __FILE__, __LINE__);
            global_data_stream->foutDouble(tempFimage_imag, ori_Fsize2, "getFourierTransformsAndCtfs()_FFT_imag", __FILE__, __LINE__);
            global_data_stream->foutDouble(exp_Fimage_nomask_real, fine_Fsize2, "getFourierTransformsAndCtfs()_expFFT_real", __FILE__, __LINE__);
            global_data_stream->foutDouble(exp_Fimage_nomask_imag, fine_Fsize2, "getFourierTransformsAndCtfs()_expFFT_imag", __FILE__, __LINE__);
            global_data_stream->check();global_data_stream->flush();
#endif
            
            if (!do_zero_mask)
            {
                // Create noisy image for outside the mask
                auto Mnoise = tempImage;
                auto Fnoise_real = tempFimage_real;
                auto Fnoise_imag = tempFimage_imag;
                
                auto sigma2_noise_igroup = mlModel.sigma2_noise[igroup].wptr(ori_Fsize);
                // Different MPI-distributed subsets may otherwise have different instances of the random noise below,
                // because work is on an on-demand basis and therefore variable with the timing of distinct nodes...
                // Have the seed based on the ipart, so that each particle has a different instant of the noise
                // Do this all inside a mutex for the threads, because they all use the same static variables inside ran1...
                // (So the mutex only goal is to make things exactly reproducible with the same random_seed.)
#pragma omp critical
                {
                    Random_generator seeded_Random_generator;
                    seeded_Random_generator.init(random_seed + exp_first_image + iimage);
                    
                    // Fill Fnoise with random numbers, use power spectrum of the noise for its variance
                    for (int i = 0; i<ori_size; i++){
                        
                        int ip = (i < ori_Fsize) ? i : i - ori_size;
                        
                        for (int j = 0; j<ori_Fsize; j++)
                        {
                            int jp = j;
                            int ires = round( sqrt( (double)(ip * ip + jp * jp) ) );
                            if (ires >= 0 && ires < ori_Fsize)
                            {
                                // If we're doing running averages, then the sigma2_noise was already adjusted for the running averages.
                                // Make a noisy background image with the same spectrum as the sigma2_noise
                                double power_noise = sigma2_fudge*sigma2_noise_igroup[ires];
                                double sigma = sqrt(power_noise);
                                Fnoise_real[i*ori_Fsize+j] = seeded_Random_generator.rnd_gaus(0., sigma);
                                Fnoise_imag[i*ori_Fsize+j] = seeded_Random_generator.rnd_gaus(0., sigma);
                            }
                            else
                            {
                                Fnoise_real[i*ori_Fsize+j] = 0;
                                Fnoise_imag[i*ori_Fsize+j] = 0;
                            }
                        }
                    }
                }
                
                // Back to real space Mnoise
                threadFFTtransformer[tid]->inverseFourierTransform(Fnoise_real, Fnoise_imag, Mnoise);
                
#ifdef DATA_STREAM
                global_data_stream->foutDouble(Mnoise, ori_size2, "getFourierTransformsAndCtfs()_Mnoise", __FILE__, __LINE__);
                global_data_stream->check();global_data_stream->flush();
#endif
                softMaskOutsideMap(exp_image,ori_size,particle_diameter / (2. * pixel_size), (double)width_mask_edge, Mnoise);
                
            }
            else
            {
                softMaskOutsideMap(exp_image,ori_size,particle_diameter / (2. * pixel_size), (double)width_mask_edge, (FDOUBLE*)nullptr);
            }
            // Inside Projector and Backprojector the origin of the Fourier Transform is centered!
            CenterFFT(exp_image,tempImage,ori_size,true);
            
            threadFFTtransformer[tid]->FourierTransform(tempImage, tempFimage_real, tempFimage_imag);
            
            // get the Fourier Transform of the image Fimg
            // BEVIN This is a waste of time when ori_size == current_size
            windowFourierTransform(tempFimage_real, tempFimage_imag, ori_size,
                                   exp_Fimage_mask_fine_real, exp_Fimage_mask_fine_imag, fine_size);
            
            // BEVIN This is a waste of time when ori_size == current_size
            windowFourierTransform(exp_Fimage_mask_fine_real, exp_Fimage_mask_fine_imag, fine_size,
                                   exp_Fimage_mask_coarse_real, exp_Fimage_mask_coarse_imag, coarse_size);
            
#ifdef DATA_STREAM
            global_data_stream->foutDouble(exp_image, ori_size2, "getFourierTransformsAndCtfs()_maskimg", __FILE__, __LINE__);
            global_data_stream->foutDouble(tempImage, ori_size2, "getFourierTransformsAndCtfs()_maskimg_centerFFT", __FILE__, __LINE__);
            global_data_stream->foutDouble(tempFimage_real, ori_Fsize2, "getFourierTransformsAndCtfs()_maksFFT_real", __FILE__, __LINE__);
            global_data_stream->foutDouble(tempFimage_imag, ori_Fsize2, "getFourierTransformsAndCtfs()_maskFFT_imag", __FILE__, __LINE__);
#endif
            // Store the power_class spectrum of the whole image (to fill sigma2_noise between current_size and ori_size
            if (fine_size < ori_size)
            {
                exp_power_imgs[iimage].zero();
                auto spectrum = exp_power_imgs[iimage].wptr(ori_Fsize);
                
                double highres_Xi2 = 0.;
                
                for (int i = 0; i<ori_size; i++){
                    
                    int ip = (i < ori_Fsize) ? i : i - ori_size;
                    
                    for (int j = 0; j<ori_Fsize; j++)
                    {
                        int jp = j;
                        int ires = round( sqrt( (double)(ip*ip + jp*jp) ) );
                        // Skip Hermitian pairs in the x==0 column
                        if (ires > 0 && ires < ori_Fsize && !(jp==0 && ip < 0) )
                        {
                            double normFaux = tempFimage_real[i*ori_Fsize+j]*tempFimage_real[i*ori_Fsize+j] + \
                            tempFimage_imag[i*ori_Fsize+j]*tempFimage_imag[i*ori_Fsize+j];
                            spectrum[ires] += normFaux;
                            // Store sumXi2 from current_size until ori_size
                            if (ires >= fine_Fsize)
                                highres_Xi2 += normFaux;
                        }
                    }
                }
                // Let's use .at() here instead of [] to check whether we go outside the vectors bounds
                exp_highres_Xi2_imgs[iimage] = highres_Xi2;
#ifdef DATA_STREAM
                global_data_stream->foutDouble(exp_power_imgs[iimage].wptr(ori_Fsize), ori_Fsize, "getFourierTransformsAndCtfs()_exp_power_imgs", __FILE__, __LINE__);
#endif
            }
            else
            {
                exp_highres_Xi2_imgs[iimage] = 0.;
            }
            
#ifdef DATA_STREAM
            global_data_stream->foutDouble(exp_highres_Xi2_imgs[iimage], "getFourierTransformsAndCtfs()_exp_highres_Xi2_imgs", __FILE__, __LINE__);
            global_data_stream->foutDouble(exp_Fimage_mask_fine_real, fine_Fsize2, "getFourierTransformsAndCtfs()_expmaksFFT_real", __FILE__, __LINE__);
            global_data_stream->foutDouble(exp_Fimage_mask_fine_imag, fine_Fsize2, "getFourierTransformsAndCtfs()_expmaskFFT_imag", __FILE__, __LINE__);
            global_data_stream->check();global_data_stream->flush();
#endif
            
            // Now calculate the actual CTF
            if (do_ctf_correction)
            {
                threadCTFer[tid]->setValues(exp_metadata[iimage].CTF_DEFOCUS_U,
                                            exp_metadata[iimage].CTF_DEFOCUS_V,
                                            exp_metadata[iimage].CTF_DEFOCUS_ANGLE,
                                            exp_metadata[iimage].CTF_VOLTAGE,
                                            exp_metadata[iimage].CTF_CS,
                                            exp_metadata[iimage].CTF_Q0,
                                            exp_metadata[iimage].CTF_BFAC);
                
                threadCTFer[tid]->getFftwImage(Fctfs[iimage].wptr(fine_Fsize2), fine_size, ori_size, ori_size, pixel_size,
                                               ctf_phase_flipped, only_flip_phases, intact_ctf_first_peak, true);
            }
            else
            {
                // TODO : check whetehr can be removed
                Fctfs[iimage].fill(1.);
            }
#ifdef DATA_STREAM
            global_data_stream->foutDouble(Fctfs[iimage].wptr(fine_Fsize2), fine_Fsize2, "getFourierTransformsAndCtfs()_CTF", __FILE__, __LINE__);
            global_data_stream->foutInt(iimage, "setFourierTransforms()_end------", __FILE__, __LINE__);
            global_data_stream->check();global_data_stream->flush();
#endif
        }
    }
    
    //
    void getLargeShiftedMaskImageOneTile(int iimage,FDOUBLE* Fimgs_shifted_mask_real,FDOUBLE* Fimgs_shifted_mask_imag,
                                  		 int n_start,int n_end,int itrans);
    //
    void getLargeShiftedNomaskImageOneTile(int iimage,FDOUBLE* Fimgs_shifted_nomask_real,FDOUBLE* Fimgs_shifted_nomask_imag,
                                           int n_start,int n_end,int itrans);
    //
    void getLargeShiftedMaskImageDecompOneTile(int iimage,FDOUBLE* Fimgs_shifted_mask_real,FDOUBLE* Fimgs_shifted_mask_imag,
                                               int n_start,int n_end,int itrans,SamplingGrid& samplingGrid);
    //
    void getLargeShiftedNomaskImageDecompOneTile(int iimage,FDOUBLE* Fimgs_shifted_nomask_real,FDOUBLE* Fimgs_shifted_nomask_imag,
                                                 int n_start,int n_end,int itrans,SamplingGrid& samplingGrid);
    //
    void getShiftedMaskImageOneTileCoarse(int iimage,FDOUBLE* Fimgs_shifted_mask_real,FDOUBLE* Fimgs_shifted_mask_imag,
                                          int n_start,int n_end, int itrans, int iover_trans, int itrans_over_trans);
    //
    void getShiftedMaskImageOneTileFine(int iimage,FDOUBLE* Fimgs_shifted_mask_real,FDOUBLE* Fimgs_shifted_mask_imag,
                                        int n_start,int n_end, int itrans, int iover_trans, int itrans_over_trans);
    //
    void getShiftedNomaskImageOneTile(int iimage,FDOUBLE* Fimgs_shifted_nomask_real,FDOUBLE* Fimgs_shifted_nomask_imag,
                                      int n_start,int n_end, int itrans, int iover_trans, int itrans_over_trans);
    //
    void testDoubleOrTripleTranslation(SamplingGrid& samplingGrid,int exp_current_size);
    
    // get all shifted Fimg,Fimg_nomask and CTF,also get inversed sigma^2
    void preShiftedImagesCtfsAndInvSigma2s(Aligned3dArray<FDOUBLE>& exp_Fimgs_shifted_real,
                                           Aligned3dArray<FDOUBLE>& exp_Fimgs_shifted_imag,
                                           VectorOfArray2d<FDOUBLE>& exp_local_Minvsigma2s,
                                           VectorOfArray1d<FDOUBLE>& exp_local_sqrtXi2,
                                           VectorOfArray2d<FDOUBLE>& exp_local_Fctfs,
                                           int exp_current_size,bool do_cc,bool do_coarse_search,
                                           SamplingGrid& samplingGrid,MLModel& mlModel,
                                           VectorOfInt& Mresol_coarse,VectorOfInt& Mresol_fine);
    
    // TODO
    void pregetShiftedImagesNomask(Aligned3dArray<FDOUBLE>& exp_nomask_Fimgs_shifted_real,
                                   Aligned3dArray<FDOUBLE>& exp_nomask_Fimgs_shifted_imag,
                                   int exp_current_size,SamplingGrid& samplingGrid);
    
    //
    void unitTest();
private:
    void inline shiftOneFimageByXYabTable(const FDOUBLE* Fimgs_in_real,const FDOUBLE* Fimgs_in_imag,
                                          FDOUBLE* Fimgs_out_real,FDOUBLE* Fimgs_out_imag,
                                          int n_start,int n_end,int itrans,SamplingGrid& samplingGrid);
};


#endif /* defined(EXP_MODEL_H_) */
