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

#include "util.h"		// used for building precompiled headers on Windows

#include "ml_model.h"

void MLModel::initialize(int _ori_size,int _nr_classes,int _nr_groups,int _nr_directions,double sigma2_angle)
{
    //
    ori_size = _ori_size;ori_Fsize = ori_size/2+1;
    nr_classes = _nr_classes;
    nr_groups = _nr_groups;
    nr_directions = _nr_directions;
    if (sigma2_angle != 0) {
        orientational_prior_mode = PRIOR_ROTTILT_PSI;
        sigma2_rot = sigma2_tilt = sigma2_psi = sigma2_angle * sigma2_angle;
    }
    //
    tau2_class			.init(nr_classes, ori_Fsize);
    sigma2_noise		.init(nr_groups,  ori_Fsize);
    sigma2_class		.init(nr_classes, ori_Fsize);	sigma2_class.fill(0.);
    data_vs_prior_class	.init(nr_classes, ori_Fsize);
    fsc_halves_class	.init(nr_classes, ori_Fsize);	fsc_halves_class.fill(0.);
    
    orientability_contrib.init(nr_classes, ori_Fsize);
    acc_rot				.init(nr_classes);
    acc_trans			.init(nr_classes);
    
    scale_correction	.init(nr_groups);				scale_correction.fill(1.);
    
    sigma2_offset = 3*3;
    avg_norm_correction = 1.;
    ave_Pmax = 0;
    LL = 0;
    
    prior_offsetx_class	.init(nr_classes);				prior_offsetx_class.fill(0.);
    prior_offsety_class	.init(nr_classes);				prior_offsety_class.fill(0.);
    
    pdf_class			.init(nr_classes);				pdf_class.fill(1./(double)nr_classes);
    pdf_direction		.init(nr_classes, nr_directions);pdf_direction.fill(1./((double)nr_classes*nr_directions));
    
    nr_particles_group	.init(nr_groups);				nr_particles_group.fill(0.);
    
    wsum_sigma2_noise	.init(nr_groups, ori_Fsize);	wsum_sigma2_noise.fill(0.);
    wsum_pdf_class		.init(nr_classes);				wsum_pdf_class.fill(0.);
    wsum_sumw_group		.init(nr_groups);				wsum_sumw_group.fill(0.);
    wsum_signal_product_spectra.init(nr_groups, ori_Fsize);wsum_signal_product_spectra.fill(0.);
    wsum_reference_power_spectra.init(nr_groups, ori_Fsize);wsum_reference_power_spectra.fill(0.);
    
    wsum_LL = 0.;
    wsum_ave_Pmax = 0.;
    wsum_avg_norm_correction = 0.;
    
    //
    wsum_prior_offsetx_class.init(nr_classes);
    wsum_prior_offsety_class.init(nr_classes);
    wsum_pdf_direction.init(nr_classes, nr_directions);
}

void MLModel::finalize()
{
#define SEP
#define ELTONE(T,N,V,S,I)
#define ELTVE1(T,N,V,S,I) N.fini();
#define ELTVE2(T,N,V,S,I) N.fini();
    MLMODEL_ELTS
#undef ELTONE
#undef ELTVE1
#undef ELTVE2
#undef SEP
}

void MLModel::resetZero()
{
    wsum_pdf_class				.fill(0.);
    wsum_sigma2_noise			.fill(0.);
    wsum_signal_product_spectra	.fill(0.);
    wsum_reference_power_spectra.fill(0.);
    wsum_sumw_group				.fill(0.);
    
    wsum_sigma2_offset = 0.;
    wsum_avg_norm_correction = 0.;
    wsum_LL = 0.;
    wsum_ave_Pmax = 0.;
    
    wsum_prior_offsetx_class	.fill(0.);
    wsum_prior_offsety_class	.fill(0.);
    wsum_pdf_direction			.fill(0.);
    
    orientability_contrib		.fill(0.);
    acc_rot						.fill(0.);
    acc_trans					.fill(0.);
}

//
void MLModel::setSigmaNoiseEstimates(Image& Iref_avg,bool do_calculate_initial_sigma_noise /*= true*/)
{
    // Calculate sigma2_noise estimates as average of power class spectra, and subtract power spectrum of the average image from that
    if (do_calculate_initial_sigma_noise)
    {
        auto wsum_sigma2_noise_group0 = wsum_sigma2_noise[0].wptr(ori_Fsize);
        auto sigma2_noise_group0 = sigma2_noise[0].wptr(ori_Fsize);
        // Factor 2 because of 2-dimensionality of the complex plane
        for (int i = 0; i < ori_Fsize; i++) {
            CHECK_NOT_NAN_OR_INF(sigma2_noise_group0[i] = wsum_sigma2_noise_group0[i] / ( 2. * wsum_sumw_group.rptrAll()[0] ));
        }
        
        // Calculate power spectrum of the average image
        //
        auto spect = Heap::allocScalars<FDOUBLE>(ori_size, __FILE__, __LINE__);
        
        getSpectrum(Iref_avg.rptrAll(), ori_size, spect, POWER_SPECTRUM);
        
        for (int i = 0; i < ori_Fsize; i++) {
            //because of 2-dimensionality of the complex plane
            spect[i] /= 2;
            // Now subtract power spectrum of the average image from the average power spectrum of the individual images
            CHECK_NOT_NAN_OR_INF(sigma2_noise_group0[i]  -= spect[i]);
        }
#if defined(FLOAT_PRECISION)
        Heap::freeFloats(spect);
#else
        Heap::freeDoubles(spect);
#endif
        // set same specturm for all groups
        for (int igroup = 0; igroup < nr_groups; igroup++){
            auto sigma2_noise_group = sigma2_noise[igroup].wptr(ori_Fsize);
            for (int i = 0; i < ori_Fsize; i++)
                sigma2_noise_group[i] = sigma2_noise_group0[i];
        }
    }
}

//
void MLModel::initialiseDataVersusPrior(MAPModel& mapModel,double tau2_fudge_factor)
{
    bool fix_tau = false;
    int nr_particles = 0;
    // Get total number of particles
    for (int igroup = 0; igroup < nr_groups; igroup++)
        nr_particles += nr_particles_group.rptrAll()[igroup];
    
    // Calculate average sigma2_noise over all image groups,now it just equal to sigma2_noise
    auto avg_sigma2_noise = Heap::allocScalars<FDOUBLE>(ori_Fsize, __FILE__, __LINE__);
    for (int i = 0; i < ori_Fsize; i++) avg_sigma2_noise[i] = 0.;
    // each group has its own sigma,the average sigma will be weighted sum
    for (int igroup = 0; igroup < nr_groups; igroup++) {
        auto sigma2_noise_igroup = sigma2_noise[igroup].wptr(ori_Fsize);
        for (int i = 0; i < ori_Fsize; i++)
            avg_sigma2_noise[i] += nr_particles_group.rptrAll()[igroup]*sigma2_noise_igroup[i];
    }
    for (int i = 0; i < ori_Fsize; i++)
        avg_sigma2_noise[i] /= nr_particles;
    
    // Get the FT of all reference structures
    // The Fourier Transforms are all "normalised" for 2D transforms of size = ori_size x ori_size
    // And spectrum is squared, so ori_size*ori_size in the 3D case!
    FDOUBLE normfft = (mapModel.ref_dim == 3) ? (FDOUBLE)(ori_size * ori_size) : 1.;
    
    // Get the power spectrum of the reference,ignore the ori_Fsize element
    auto spectrum = Heap::allocScalars<FDOUBLE>(ori_size, __FILE__, __LINE__);
    
    for (int iclass = 0; iclass < nr_classes; iclass++)
    {
        auto tau2_class_iclass = tau2_class[iclass].wptr(ori_Fsize);
        auto data_vs_prior_class_iclass = data_vs_prior_class[iclass].wptr(ori_Fsize);
        // Update the tau2_class spectrum for this reference
        // This is only for writing out in the it000000_model.star file
        if (!fix_tau)
        {
            // Initialise output arrays to correct size
            getSpectrum(mapModel.Irefs[iclass], spectrum, POWER_SPECTRUM);
            
            // TODO a more powerful assignment operator
            
            // Factor two because of two-dimensionality of the complex plane
            // (just like sigma2_noise estimates, the power spectra should be divided by 2)
            for (int i = 0; i < ori_size; i++) {
                spectrum[i] *= normfft / 2.;
            }
            
            for (int i = 0; i < ori_Fsize; i++)
                tau2_class_iclass[i] = tau2_fudge_factor * spectrum[i];
        }
        
        // Calculate data_vs_prior_class as spectral_nr_observations_per_class/sigma2_noise vs 1/tau2_class
        for (int i = 0; i < ori_Fsize; i++)
        {
            double evidence = CHECK_NOT_NAN_OR_INF(nr_particles * pdf_class.rptrAll()[iclass] / avg_sigma2_noise[i]);
            // empirical accounting for ratio of pixels in 3D shells compared to 2D shells
            if (mapModel.ref_dim == 3 && i > 0) evidence /= (2. * (FDOUBLE)i);// ?????Yongbei
            FDOUBLE prior = 1. /  tau2_class_iclass[i];
            FDOUBLE myssnr = evidence / prior;
            data_vs_prior_class_iclass[i] = myssnr;
            // Also initialise FSC-halves here (...)
            // DIRECT_A1D_ELEM(fsc_halves_class[iclass], i ) = myssnr / (myssnr + 1);
        }
    } // end loop iclass
#if defined(FLOAT_PRECISION)
    Heap::freeFloats(avg_sigma2_noise);
    Heap::freeFloats(spectrum);
#else
    Heap::freeDoubles(avg_sigma2_noise);
    Heap::freeDoubles(spectrum);
#endif
}

void MLModel::initialisePdfDirection(int newsize)
{
    // If the pdf_direction were already filled (size!=0), and newsize=oldsize then leave them as they were
    assert(nr_directions!=newsize);
    // If they were still empty, or if the size changes, then initialise them with an even distribution
    nr_directions = newsize;
    pdf_direction		.fini();
    pdf_direction		.init(nr_classes, nr_directions);pdf_direction.fill(1./((double)nr_classes*nr_directions));
    wsum_pdf_direction	.fini();
    wsum_pdf_direction	.init(nr_classes, nr_directions);wsum_pdf_direction.fill(0.);
}

//
void MLModel::writeOutModel(std::string fn_model,const MAPModel& mapModel,const MetaDataTable& metadata,int random_subset)
{
    std::ofstream modelFile;
    modelFile.open((fn_model+".star").c_str(), std::ios::out);
    //
    {
        modelFile << std::endl;
        modelFile << "data_model_general" <<std::endl;
        modelFile << std::endl;
#define COUTMETADATA(v1,v2) modelFile << "_rln" << std::left<<std::setw(30) << v1 << std::right<<std::setw(18) << v2 <<std::endl;
        //
        COUTMETADATA("ReferenceDimensionality"	, mapModel.ref_dim				)
        // COUTMETADATA("DataDimensionality"		, 2								) // support > 1.4
        COUTMETADATA("OriginalImageSize"		, ori_size						)// equal particle model fine size
        COUTMETADATA("CurrentResolution"		, mapModel.current_resolution	)
        COUTMETADATA("CurrentImageSize"			, ori_size						)
        COUTMETADATA("PaddingFactor"			, 2								)
        COUTMETADATA("FourierSpaceInterpolator"	, 1								)
        COUTMETADATA("MinRadiusNnInterpolation" , 10							)
        COUTMETADATA("PixelSize"				, mapModel.pixel_size			)
        COUTMETADATA("NrClasses"				, nr_classes					)
        COUTMETADATA("NrGroups"					, nr_groups						)
        COUTMETADATA("Tau2FudgeFactor"			, 4								)
        COUTMETADATA("NormCorrectionAverage"	, avg_norm_correction			)
        COUTMETADATA("SigmaOffsets"				, sqrt(sigma2_offset)			)
        COUTMETADATA("OrientationalPriorMode"	, 0								)
        COUTMETADATA("SigmaPriorRotAngle"		, 0								)
        COUTMETADATA("SigmaPriorTiltAngle"		, 0								)
        COUTMETADATA("SigmaPriorPsiAngle"		, 0								)
        COUTMETADATA("LogLikelihood"			, LL							)
        COUTMETADATA("AveragePmax"				, ave_Pmax						)
        //
#undef COUTMETADATA
        modelFile << std::endl << std::endl;
    }
    //
    std::vector<int> SpectralIndex(ori_Fsize,0);
    std::vector<double> Resolution(ori_Fsize,0);
    std::vector<double> AngstromResolution(ori_Fsize,0);
    for (int ipix = 0; ipix < ori_Fsize; ipix++){
        SpectralIndex[ipix] = ipix;
        Resolution[ipix] = ipix/(mapModel.pixel_size * ori_size);
        AngstromResolution[ipix] = ( (ipix==0) ? 999. : (mapModel.pixel_size * ori_size)/(double)ipix );
    }
    // data_model_classes
    {
        SmallMetataDataTable data_model_classes("data_model_classes");
        data_model_classes.appendName({"ReferenceImage","ClassDistribution","AccuracyRotations","AccuracyTranslations"});
        data_model_classes.appendType({ElemTypeChar,ElemTypeDouble,ElemTypeDouble,ElemTypeDouble});
        std::vector<std::string*> ReferenceImage(nr_classes);
        std::vector<double> AccuracyRotations(nr_classes);
        std::vector<double> AccuracyTranslations(nr_classes);
        for (int iclass = 0; iclass < nr_classes; iclass++) {
            ReferenceImage[iclass] = 
#include "./util_heap_undefs.h"
				born(new std::string(fn_model.substr(0,fn_model.find("_model"))+"_class" + num2str(iclass+1) + ".mrc"));
#include "./util_heap_defs.h"
			AccuracyRotations[iclass] = 999.;
            AccuracyTranslations[iclass] = 999.;
        }
        data_model_classes.appendElem({ReferenceImage.data(),pdf_class.wptrAll(),AccuracyRotations.data(),AccuracyTranslations.data()}, nr_classes);
        data_model_classes.print(modelFile);
    }
    // data_model_class_*
    for (int iclass = 0; iclass < nr_classes; iclass++)
    {
        SmallMetataDataTable data_model_classes_one("data_model_class_"+std::to_string((long long)(iclass+1)));
        data_model_classes_one.appendName({ "SpectralIndex","Resolution","AngstromResolution","SsnrMap",
            "GoldStandardFsc","ReferenceSigma2","ReferenceTau2","SpectralOrientabilityContribution"});
        data_model_classes_one.appendType({	ElemTypeInt,ElemTypeDouble,ElemTypeDouble,ElemTypeDouble,
            ElemTypeDouble,ElemTypeDouble,ElemTypeDouble,ElemTypeDouble});
        std::vector<double> orientability_contrib(ori_Fsize,0);
        orientability_contrib[0] = 1.;
        data_model_classes_one.appendElem({	SpectralIndex.data(),Resolution.data(),AngstromResolution.data(),data_vs_prior_class[iclass].wptrAll(),
            								fsc_halves_class[iclass].wptrAll(),sigma2_class[iclass].wptrAll(),
            								tau2_class[iclass].wptrAll(),orientability_contrib.data()}, ori_Fsize);
        data_model_classes_one.print(modelFile);
    }
    // data_model_groups
    {
        SmallMetataDataTable data_model_groups("data_model_groups");
        data_model_groups.appendName({"GroupNumber","GroupName","GroupNrParticles","GroupScaleCorrection"});
        data_model_groups.appendType({ElemTypeInt,ElemTypeChar,ElemTypeInt,ElemTypeDouble});
        std::vector<int> GroupNumber(nr_groups);
        for (int i = 0; i < nr_groups; i++)
            GroupNumber[i] = i+1;
        std::vector<std::string*> GroupName(nr_groups);
        std::vector<int> GroupNrParticles(nr_groups,0);
        // write all group_name and group_number
        int nr_particles_all = metadata.numberOfParticles(-1);
        for (int i = 0; i < nr_particles_all; i++){
            int group_no = metadata.accessAll(i).GROUP_NO-1;
            GroupName[group_no] =
#include "./util_heap_undefs.h"
            born(new std::string(metadata.accessAll(i).GROUP_NAME));
#include "./util_heap_defs.h"
        }
        // only write this half particles number
        int nr_particles_half = metadata.numberOfParticles(random_subset);
        for (int i = 0; i < nr_particles_half; i++){
            int group_no = metadata[i].GROUP_NO-1;
            GroupNrParticles[group_no]++;
		}
        data_model_groups.appendElem({GroupNumber.data(),GroupName.data(),GroupNrParticles.data(),scale_correction.wptrAll()}, nr_groups);
        data_model_groups.print(modelFile);
    }
    // data_model_groups_*
    for (int igroup = 0; igroup < nr_groups; igroup++)
    {
        SmallMetataDataTable data_model_groups_one("data_model_group_"+std::to_string((long long)(igroup+1)));
        data_model_groups_one.appendName({"SpectralIndex","Resolution","Sigma2Noise"});
        data_model_groups_one.appendType({ElemTypeInt,ElemTypeDouble,ElemTypeDouble});
        data_model_groups_one.appendElem({SpectralIndex.data(),Resolution.data(),sigma2_noise[igroup].wptrAll()}, ori_Fsize);
        data_model_groups_one.print(modelFile);
    }
    // data_model_pdf_orient_class_*
    for (int iclass = 0; iclass < nr_classes; iclass++) {
        SmallMetataDataTable data_model_classes_one("data_model_pdf_orient_class_"+std::to_string((long long)(iclass+1)) );
        data_model_classes_one.appendName({"OrientationDistribution"});
        data_model_classes_one.appendType({ElemTypeDouble});
        data_model_classes_one.appendElem({pdf_direction[iclass].wptrAll()}, nr_directions);
        data_model_classes_one.print(modelFile);
    }
    modelFile.close();
}

void MLModel::writeOutBild(std::string fn_class, MAPModel& mapModel, HealpixSampler& sampling)
{
    for (int iclass = 0; iclass < nr_classes; iclass++) {
        std::string fn_bild = fn_class + num2str(iclass+1) + "_angdist.bild";
        FDOUBLE offset = ori_size * mapModel.pixel_size / 2.;
        sampling.writeBildFileOrientationalDistribution(pdf_direction[iclass], fn_bild, offset, offset);
    }
}

void MLModel::readFromModel(std::string fn_model,MAPModel& mapModel,bool debug_flag /*= false*/)
{
    ifstreamCheckingExistence modelFile(fn_model.c_str());
    {
        bool startingRead = false;
        std::string line;
        while (true) {
            if (startingRead) {
                double doubleTmp;
#define CINMETADATA(V) modelFile >> line;modelFile >> doubleTmp;if(debug_flag) {MASTERNODE std::cout<<std::setw(30)<<line<<" "<<doubleTmp<<std::endl;}
                //
                CINMETADATA("ReferenceDimensionality"	)
                CINMETADATA("DataDimensionality"		)
                if (line.find("DataDimensionality") != std::string::npos) { // RELION 1.4,support 3D data
                    CINMETADATA("OriginalImageSize"		)
                }
                else{} // RELION 1.3
                CINMETADATA("CurrentResolution"			);mapModel.current_resolution = doubleTmp;
                CINMETADATA("CurrentImageSize"			)
                CINMETADATA("PaddingFactor"				)
                CINMETADATA("FourierSpaceInterpolator"	)
                CINMETADATA("MinRadiusNnInterpolation"	)
                CINMETADATA("pixel_size"				);
                if(fabs(doubleTmp-mapModel.pixel_size)>1e-6){
                    std::cerr<<doubleTmp<<" "<<mapModel.pixel_size<<std::endl;
                    ERROR_REPORT("Set pixel_size.");
                }
                CINMETADATA("NrClasses"					);ERROR_CHECK(nr_classes!=(int)doubleTmp, "Set mlModel.nr_classes.");
                CINMETADATA("NrGroups"					);ERROR_CHECK(nr_groups!=(int)doubleTmp, "Set mlModel.nr_groups.");
                CINMETADATA("Tau2FudgeFactor"			)
                CINMETADATA("NormCorrectionAverage"		);avg_norm_correction = doubleTmp;
                CINMETADATA("SigmaOffsets"				);sigma2_offset = doubleTmp*doubleTmp;
                CINMETADATA("OrientationalPriorMode"	)
                CINMETADATA("SigmaPriorRotAngle"		)
                CINMETADATA("SigmaPriorTiltAngle"		)
                CINMETADATA("SigmaPriorPsiAngle"		)
                CINMETADATA("LogLikelihood"				)
                CINMETADATA("AveragePmax"				);ave_Pmax = doubleTmp;
                //
#undef CINMETADATA
                break;
            }
            else{
                getline(modelFile,line);
                if(debug_flag) MASTERNODE std::cout<<line<<std::endl;
                if ((line.find("data_model_general") !=std::string::npos) ){
                    startingRead = true;
                    getline(modelFile,line);assert(line=="");// escape a empty line
                }
                ERROR_CHECK(modelFile.eof(), "end of model file,can not find data_model_general.");
            }
        }
    }
    //
    std::vector<int> SpectralIndex(ori_Fsize,0);
    std::vector<double> Resolution(ori_Fsize,0);
    std::vector<double> AngstromResolution(ori_Fsize,0);
    // data_model_classes
    {
        SmallMetataDataTable data_model_classes("data_model_classes");
        data_model_classes.appendName({"ReferenceImage","ClassDistribution","AccuracyRotations","AccuracyTranslations"});
        data_model_classes.appendType({ElemTypeChar,ElemTypeDouble,ElemTypeDouble,ElemTypeDouble});
        std::vector<std::string*> ReferenceImage(nr_classes);
        std::vector<double> AccuracyRotations(nr_classes);
        std::vector<double> AccuracyTranslations(nr_classes);
        data_model_classes.appendElem({ReferenceImage.data(),pdf_class.wptrAll(),AccuracyRotations.data(),AccuracyTranslations.data()}, nr_classes);
        bool metadata_fit = data_model_classes.read(modelFile);
        ERROR_CHECK(!metadata_fit, "ml_model class number is not correct for pdf_class");
        // read the reference
        std::vector<std::string> ref_fns(nr_classes);
        for (int iclass = 0; iclass < nr_classes; iclass++) ref_fns[iclass] = *ReferenceImage[iclass];
        mapModel.initializeRef(ref_fns);
    }
    // data_model_class_*
    {
        SmallMetataDataTable data_model_classes_one("data_model_class_"+std::to_string((long long)(0+1)));
        if (data_model_classes_one.contain(modelFile, "SpectralOrientabilityContribution"))
        {
            for (int iclass = 0; iclass < nr_classes; iclass++)
            {
                SmallMetataDataTable data_model_classes_one("data_model_class_"+std::to_string((long long)(iclass+1)));
                data_model_classes_one.appendName({ "SpectralIndex","Resolution","AngstromResolution","SsnrMap",
                    "GoldStandardFsc","ReferenceSigma2","ReferenceTau2","SpectralOrientabilityContribution"});
                data_model_classes_one.appendType({	ElemTypeInt,ElemTypeDouble,ElemTypeDouble,ElemTypeDouble,
                    ElemTypeDouble,ElemTypeDouble,ElemTypeDouble,ElemTypeDouble});
                std::vector<double> orientability_contrib(ori_Fsize,0);
                data_model_classes_one.appendElem({	SpectralIndex.data(),Resolution.data(),AngstromResolution.data(),data_vs_prior_class[iclass].wptrAll(),
                    fsc_halves_class[iclass].wptrAll(),sigma2_class[iclass].wptrAll(),
                    tau2_class[iclass].wptrAll(),orientability_contrib.data()}, ori_Fsize);
                bool metadata_fit = data_model_classes_one.read(modelFile);
                ERROR_CHECK(!metadata_fit, "ml_model class number is not correct for fsc_halves_class[iclass],sigma2_class,tau2_class");
            }
        }
        else
        {
            for (int iclass = 0; iclass < nr_classes; iclass++)
            {
                SmallMetataDataTable data_model_classes_one("data_model_class_"+std::to_string((long long)(iclass+1)));
                data_model_classes_one.appendName({ "SpectralIndex","Resolution","AngstromResolution","SsnrMap",
                    "GoldStandardFsc","ReferenceSigma2","ReferenceTau2"/*,"SpectralOrientabilityContribution"*/});
                data_model_classes_one.appendType({	ElemTypeInt,ElemTypeDouble,ElemTypeDouble,ElemTypeDouble,
                    ElemTypeDouble,ElemTypeDouble,ElemTypeDouble/*,ElemTypeDouble*/});
                std::vector<double> orientability_contrib(ori_Fsize,0);
                data_model_classes_one.appendElem({	SpectralIndex.data(),Resolution.data(),AngstromResolution.data(),data_vs_prior_class[iclass].wptrAll(),
            								fsc_halves_class[iclass].wptrAll(),sigma2_class[iclass].wptrAll(),
            								tau2_class[iclass].wptrAll()/*,orientability_contrib.data()*/}, ori_Fsize);
                bool metadata_fit = data_model_classes_one.read(modelFile);
                ERROR_CHECK(!metadata_fit, "ml_model class number is not correct for fsc_halves_class[iclass],sigma2_class,tau2_class");
            }
        }
    }
    // data_model_groups
    {
        SmallMetataDataTable data_model_groups("data_model_groups");
        data_model_groups.appendName({"GroupNumber","GroupName","GroupNrParticles","GroupScaleCorrection"});
        data_model_groups.appendType({ElemTypeInt,ElemTypeChar,ElemTypeInt,ElemTypeDouble});
        std::vector<int> GroupNumber(nr_groups);
        std::vector<std::string*> GroupName(nr_groups);
        std::vector<int> GroupNrParticles(nr_groups,0);
        data_model_groups.appendElem({GroupNumber.data(),GroupName.data(),GroupNrParticles.data(),scale_correction.wptrAll()}, nr_groups);
        bool metadata_fit = data_model_groups.read(modelFile);
        ERROR_CHECK(!metadata_fit, "ml_model group number is not correct for scale_correction");
    }
    // data_model_groups_*
    for (int igroup = 0; igroup < nr_groups; igroup++)
    {
        SmallMetataDataTable data_model_groups_one("data_model_group_"+std::to_string((long long)(igroup+1)));
        data_model_groups_one.appendName({"SpectralIndex","Resolution","Sigma2Noise"});
        data_model_groups_one.appendType({ElemTypeInt,ElemTypeDouble,ElemTypeDouble});
        data_model_groups_one.appendElem({SpectralIndex.data(),Resolution.data(),sigma2_noise[igroup].wptrAll()}, ori_Fsize);
        bool metadata_fit = data_model_groups_one.read(modelFile);
        ERROR_CHECK(!metadata_fit, "ml_model half-size is not correct for sigma2_noise[igroup]");
    }
    // data_model_pdf_orient_class_*
    for (int iclass = 0; iclass < nr_classes; iclass++) {
        SmallMetataDataTable data_model_classes_one("data_model_pdf_orient_class_"+std::to_string((long long)(iclass+1)) );
        data_model_classes_one.appendName({"OrientationDistribution"});
        data_model_classes_one.appendType({ElemTypeDouble});
        data_model_classes_one.appendElem({pdf_direction[iclass].wptrAll()}, nr_directions);
        bool metadata_fit = data_model_classes_one.read(modelFile);
        // use increase or decrease the heal_pix_order,reset the pdf_direction
        // directions_have_changed
        if (!metadata_fit) {
            if(iclass==0) MASTERNODE std::cout<<"##### reset pdf_direction because increase the healpix order!"<<std::endl;
            if(debug_flag) MASTERNODE std::cout<<"????? why not reset pdf_direction.... TODO..."<<std::endl;
            for (int idir = 0; idir < nr_directions; idir++)
                pdf_direction[iclass].wptrAll()[idir] = 1./(nr_classes * nr_directions);
        }
    }
    modelFile.close();
}

//
void MLModel::writeToDisk(std::string fn_ml_model)
{
    fn_ml_model = pathRemoveSuffix(fn_ml_model);
    std::ofstream mLModel_file;
    mLModel_file.open((fn_ml_model+".mlmodel").c_str(),std::ios::out | std::ios::binary | std::ios::trunc);
    
#define SEP
#define ELTONE(T,N,V,S,I) mLModel_file.write((char*)&N,(S)*sizeof(T));
#define ELTVE1(T,N,V,S,I) mLModel_file.write((char*)N.wptr(N.size()),(S)*sizeof(T));
#define ELTVE2(T,N,V,S,I) for(int i=0;i<N.size();i++) mLModel_file.write((char*)N[i].wptr(N[i].size()),(S)*sizeof(T));
    MLMODEL_ELTS
#undef ELTONE
#undef ELTVE1
#undef ELTVE2
#undef SEP
    
    mLModel_file.close();
}

//
void MLModel::readFromDisk(std::string fn_ml_model)
{
    fn_ml_model = pathRemoveSuffix(fn_ml_model);
    
    std::ifstream mLModel_file;
    mLModel_file.open((fn_ml_model+".mlmodel").c_str(),std::ios::in | std::ios::binary);
    
#define SEP
#define ELTONE(T,N,V,S,I) mLModel_file.read((char*)&N,(S)*sizeof(T));
#define ELTVE1(T,N,V,S,I) mLModel_file.read((char*)N.wptr(N.size()),(S)*sizeof(T));
#define ELTVE2(T,N,V,S,I) for(int i=0;i<N.size();i++) mLModel_file.read((char*)N[i].wptr(N[i].size()),(S)*sizeof(T));
    MLMODEL_ELTS
#undef ELTONE
#undef ELTVE1
#undef ELTVE2
#undef SEP
    
    mLModel_file.close();
}

//
void MLModel::printSpaceInfo(std::ostream &os)
{
    double perGB = 1024.*1024.*1024.;
    double totalSpace = 0.;
    os<<"Memory space for MLModel : "<<std::endl;
#define SEP
#define ELTONE(T,N,V,S,I)
#define ELTVE1(T,N,V,S,I) os<<std::setw(40)<<#N<<": "<<sizeof(T)*S/perGB<<" GB."<<std::endl;totalSpace += sizeof(T)*S/perGB;
#define ELTVE2(T,N,V,S,I) os<<std::setw(40)<<#N<<": "<<sizeof(T)*S/perGB<<" GB."<<std::endl;totalSpace += sizeof(T)*S/perGB;
    MLMODEL_ELTS
#undef ELTONE
#undef ELTVE1
#undef ELTVE2
#undef SEP
    os<<std::setw(40)<<"MLMODEL TotalSpace"<<": "<<totalSpace<<" GB."<<std::endl;
    std::cout<<"----------------------------------------------"<<std::endl;
}

// TODO
void MLModel::reduceData(const MPI::Intracomm& currentWorld)
{
#ifdef USEMPI
    int local_temp_size = std::max(std::max(nr_groups,ori_Fsize), std::max(nr_classes,nr_directions));
    FDOUBLE* local_temp = (FDOUBLE*)aMalloc(sizeof(FDOUBLE)*local_temp_size,64);
    
    for (int igroup = 0; igroup < nr_groups; igroup++)
    {
        //
        memcpy(local_temp, wsum_sigma2_noise[igroup].rptrAll(), sizeof(FDOUBLE)*ori_Fsize);
        currentWorld.Allreduce(local_temp,wsum_sigma2_noise[igroup].wptrAll(),ori_Fsize,MPI_FDOUBLE,MPI::SUM);
        //
        memcpy(local_temp, wsum_signal_product_spectra[igroup].wptrAll(), sizeof(FDOUBLE)*ori_Fsize);
        currentWorld.Allreduce(local_temp,wsum_signal_product_spectra[igroup].wptrAll(),ori_Fsize,MPI_FDOUBLE,MPI::SUM);
        //
        memcpy(local_temp, wsum_reference_power_spectra[igroup].rptrAll(), sizeof(FDOUBLE)*ori_Fsize);
        currentWorld.Allreduce(local_temp,wsum_reference_power_spectra[igroup].wptrAll(),ori_Fsize,MPI_FDOUBLE,MPI::SUM);
    }
    //
    memcpy(local_temp, wsum_pdf_class.rptrAll(), sizeof(FDOUBLE)*nr_classes);
    currentWorld.Allreduce(local_temp,wsum_pdf_class.wptrAll(),nr_classes,MPI_FDOUBLE,MPI::SUM);
    //
    for (int iclass = 0; iclass < nr_classes; iclass++) {
        memcpy(local_temp, wsum_pdf_direction[iclass].rptrAll(), sizeof(FDOUBLE)*nr_directions);
        currentWorld.Allreduce(local_temp,wsum_pdf_direction[iclass].wptrAll(),nr_directions,MPI_FDOUBLE,MPI::SUM);
    }
    // wsum_prior_offsetx_class and wsum_prior_offsety_class is not needed for 3D
    memcpy(local_temp, wsum_prior_offsetx_class.rptrAll(), sizeof(FDOUBLE)*nr_classes);
    currentWorld.Allreduce(local_temp,wsum_prior_offsetx_class.wptrAll(),nr_classes,MPI_FDOUBLE,MPI::SUM);
    //
    memcpy(local_temp, wsum_prior_offsety_class.rptrAll(), sizeof(FDOUBLE)*nr_classes);
    currentWorld.Allreduce(local_temp,wsum_prior_offsety_class.wptrAll(),nr_classes,MPI_FDOUBLE,MPI::SUM);
    //
    memcpy(local_temp, wsum_sumw_group.rptrAll(), sizeof(FDOUBLE)*nr_groups);
    currentWorld.Allreduce(local_temp,wsum_sumw_group.wptrAll(),nr_groups,MPI_FDOUBLE,MPI::SUM);
    //
    local_temp[0] = wsum_sigma2_offset;
    currentWorld.Allreduce(local_temp,&wsum_sigma2_offset,1,MPI_FDOUBLE,MPI::SUM);
    
    local_temp[0] = wsum_avg_norm_correction;
    currentWorld.Allreduce(local_temp,&wsum_avg_norm_correction,1,MPI_FDOUBLE,MPI::SUM);
    
    local_temp[0] = wsum_LL;
    currentWorld.Allreduce(local_temp,&wsum_LL,1,MPI_FDOUBLE,MPI::SUM);
    
    local_temp[0] = wsum_ave_Pmax;
    currentWorld.Allreduce(local_temp,&wsum_ave_Pmax,1,MPI_FDOUBLE,MPI::SUM);
    
    aFree(local_temp);
#endif
}

// TODO
void MLModel::broadcastData()
{
//#ifdef USEMPI
//    int nodes = MPI::COMM_WORLD.Get_size();
//    int node = MPI::COMM_WORLD.Get_rank();
//    int ranks = nodes;
//    int first_local_class,last_local_class;
//    int nr_local_classes = divide_equally(nr_classes, nodes, node, first_local_class, last_local_class);
//    // some MPI datatype
//    struct DataIclassAtom {
//        double sigma2_iclass;
//        double tau2_iclass;
//        double data_vs_prior_iclass;
//        double fsc_halves_iclass;
//    };
//    int num_members = 4;
//    int lengths[] = { 1, 1, 1, 1 };
//    MPI::Aint offsets[] = { offsetof(DataIclassAtom, sigma2_iclass),
//                            offsetof(DataIclassAtom, tau2_iclass),
//        					offsetof(DataIclassAtom, data_vs_prior_iclass),
//        					offsetof(DataIclassAtom, fsc_halves_iclass)};
//    MPI::Datatype type[] = {MPI::DOUBLE,MPI::DOUBLE,MPI::DOUBLE,MPI::DOUBLE};
//    MPI::Datatype MPIDataIclassAtom = MPI::Datatype::Create_struct(num_members, lengths, offsets, type);
//    MPIDataIclassAtom.Commit();
//    //
//    MPI::Datatype MPIDataIclass = MPIDataIclassAtom.Create_contiguous(ori_Fsize);
//    MPIDataIclass.Commit();
//    // manually pack the data
//    DataIclassAtom* sendBuf = (DataIclassAtom*)aMalloc(sizeof(DataIclassAtom)*nr_local_classes*ori_Fsize,64);
//    for (int i = 0; i < nr_local_classes*ori_Fsize; i++) {
//        sendBuf[i].sigma2_iclass = sigma2_class[first_local_class*ori_Fsize+i];
//        sendBuf[i].tau2_iclass = tau2_class[first_local_class*ori_Fsize+i];
//        sendBuf[i].data_vs_prior_iclass = data_vs_prior_class[first_local_class*ori_Fsize+i];
//        sendBuf[i].fsc_halves_iclass = fsc_halves_class[first_local_class*ori_Fsize+i];
//    }
//    // recv data
//    DataIclassAtom* recvBuf = (DataIclassAtom*)aMalloc(sizeof(DataIclassAtom)*nr_classes*ori_Fsize,64);
//    int* rcounts = (int*)aMalloc(sizeof(int)*ranks,64);
//    int* displs = (int*)aMalloc(sizeof(int)*ranks,64);
//    for	(int rank = 0; rank < ranks; rank++){
//        int first_local_class,last_local_class;
//        int nr_local_classes = divide_equally(nr_classes, ranks, rank, first_local_class, last_local_class);
//        if (rank == 0){
//            rcounts[rank] = nr_local_classes;
//            displs[rank] = 0;
//        }else{
//            rcounts[rank] = nr_local_classes;
//            displs[rank] = displs[rank-1]+rcounts[rank-1];
//        }
//    }
//    //
//    MPI::COMM_WORLD.Allgatherv(sendBuf,nr_local_classes,MPIDataIclass,recvBuf,rcounts,displs,MPIDataIclass);
//    // unpack the data
//    for (int i = 0; i < nr_classes*ori_Fsize; i++) {
//        sigma2_class[i] = recvBuf[i].sigma2_iclass;
//        tau2_class[i] = recvBuf[i].tau2_iclass;
//        data_vs_prior_class[i] = recvBuf[i].data_vs_prior_iclass;
//        fsc_halves_class[i] = recvBuf[i].fsc_halves_iclass;
//    }
//    //
//    aFree(sendBuf);aFree(recvBuf);aFree(rcounts);aFree(displs);
//    MPIDataIclass.Free();
//    MPIDataIclassAtom.Free();
//#endif
//    // TODO,directly memory copy and allgather
//    // which is faster...
}