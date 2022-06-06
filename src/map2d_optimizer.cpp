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
#include "./map2d_optimizer.h"

namespace Map2dOptimizer_old {
    
    // use datastream to write out and compare the Intermediate result
    DataStream global_data_stream;
    
    // model
    HealpixSampler sampler2d;
    MAPModel mapModel;
    MLModel mlModel;
    ParticleModel particleModel;
    HiddenVariableMonitor hiddenVarMonitor;
    
    // image data
    Images images;
    MetaDataTable metadata;
    
    // sampling
    double offset_step;
    double offset_range;
    double psi_step;
    SamplingGrid samplingGrid;
    
    // ------------ variable for expectation step -------------- //
#define SEP
#define ELTONE(T,N,S1,S2) T N;
#define ELTVE1(T,N,S1,S2) VectorOfArray1d<T> N;
#define ELTVE2(T,N,S1,S2) VectorOfArray2d<T> N;
#define ELTVE3(T,N,S1,S2) Aligned3dArray <T> N;
    MAPOPTIMIZER_OLD_EXP_VARS
    MAPOPTIMIZER_OLD_THREAD_VARS
#undef SEP
#undef ELTONE
#undef ELTVE1
#undef ELTVE2
#undef ELTVE3
	//
    VectorOfStruct<MetaDataElem> exp_metadata;
	//
    Aligned3dArray<FDOUBLE> exp_Frefs_Rot_real;
    Aligned3dArray<FDOUBLE> exp_Frefs_Rot_imag;
    Aligned3dArray<FDOUBLE> exp_Fweight_Rot;
    
    // metadata of the weight,structure is iclass->rot->iover_rot->trans->iover_trans
    double* exp_Mweight;
    int exp_Mweight_xsize;
    
    // ------------ thread variable ----------------- //
    Aligned2dArray<FDOUBLE> thread_Frefctf_real;
    Aligned2dArray<FDOUBLE> thread_Frefctf_imag;
    Aligned2dArray<FDOUBLE> thread_wsum_pdf_direction;
    std::vector<std::vector<GridIndex>> thread_exp_max_weight_index;
    //

void setupMap2dOptimizer() {
#ifdef DATA_STREAM
    global_data_stream.init(data_stream_out_fn, data_stream_in_fn);
    global_data_stream.foutInt(&nr_iter, 1, "nr_iter", __FILE__, __LINE__);
    global_data_stream.foutInt(&nr_classes, 1, "nr_classes", __FILE__, __LINE__);
    global_data_stream.foutDouble(&pixel_size, 1, "pixel_size", __FILE__, __LINE__);
    global_data_stream.foutInt(&random_seed, 1, "random_seed", __FILE__, __LINE__);
    global_data_stream.foutDouble(&ini_high, 1, "ini_high", __FILE__, __LINE__);
    global_data_stream.foutDouble(&tau2_fudge_factor, 1, "tau2_fudge_factor", __FILE__, __LINE__);
    global_data_stream.foutDouble(&particle_diameter, 1, "particle_diameter", __FILE__, __LINE__);
    global_data_stream.foutInt(&adaptive_oversampling, 1, "adaptive_oversampling", __FILE__, __LINE__);
    global_data_stream.foutInt(&sampler2d.healpix_order, 1, "sampling.healpix_order", __FILE__, __LINE__);
    global_data_stream.foutDouble(&sampler2d.psi_step, 1, "sampling.psi_step", __FILE__, __LINE__);
    global_data_stream.foutDouble(-91, "sampling.limit_tilt", __FILE__, __LINE__);
    global_data_stream.foutDouble(&offset_range, 1, "sampling.offset_range", __FILE__, __LINE__);
    global_data_stream.foutDouble(&offset_step, 1, "sampling.offset_step", __FILE__, __LINE__);
    global_data_stream.foutDouble(0.5, "sampling.perturbation_factor", __FILE__, __LINE__);
    global_data_stream.check();global_data_stream.flush();
#endif
}

void readImages() {
    
#ifdef USEMPI
    // MPI::Init();//init mpi outside
    nodes = MPI::COMM_WORLD.Get_size();
    node = MPI::COMM_WORLD.Get_rank();
    MPI::Get_processor_name(nodeName,nodeNameLen);
#else
    nodes = 1;
    node = 0;
#endif
    
    // omp_set_num_threads(8);
    if(0 == node) std::cout<<"number of threads = "<<omp_get_max_threads()<<std::endl;
    
    metadata.readFromStar(star_fn);
    nr_global_images = metadata.numberOfParticles();
    
    Mrcs::MrcsHead mrcsHead;
    Mrcs::readMrcsHead(metadata[0].IMAGE.NAME, mrcsHead);
    NODE0ONLY std::cout<<"x = "<<mrcsHead.NC<<",y = "<<mrcsHead.NC<<",z = "<<mrcsHead.NS<<",mode = "<<mrcsHead.MODE<<std::endl;
    
    ori_size = mrcsHead.NC;
    
    metadata.shuffle();
    
    // NOTE : SAME_IMAGE_DIVISION uses to produce the same result for different nodes
#ifdef SAME_IMAGE_DIVISION
    //divided like not mpi version,just to easily debug
    nr_local_images = divide_equally_pool(nr_global_images,nr_pool,nodes,node, first_local_image, last_local_image);
#else
    nr_local_images = divide_equally(nr_global_images,nodes,node, first_local_image, last_local_image);
#endif

    images.resize(nr_local_images);
    for (auto& image : images) image.init(ori_size*ori_size);
    
    // All nodes read the file at once rather than broadcasting
    float *buffer = (float*)aMalloc(sizeof(float)*ori_size*ori_size,64);
    std::string preFileName = metadata[first_local_image].IMAGE.NAME;
    FILE * mrcsFile = fopen((metadata[first_local_image].IMAGE.NAME).c_str(),"rb");
    
    for (int iimage = 0; iimage < nr_local_images; iimage++) {
        
        if (metadata[iimage+first_local_image].IMAGE.NAME != preFileName)
        {
            fclose(mrcsFile);
            mrcsFile = fopen((metadata[iimage+first_local_image].IMAGE.NAME).c_str(),"rb");
            preFileName = metadata[iimage+first_local_image].IMAGE.NAME;
        }
        
        //change type to int may cause bug,avoid offset out of range
        long image_id = metadata[iimage+first_local_image].IMAGE.INDEX;
        
        long offset = (256+(image_id-1)*ori_size*ori_size)*sizeof(float);
        
        //std::cout<<image_id<<" "<<offset<<std::endl;
        
        fseek(mrcsFile,offset,SEEK_SET);
        
        if(fread((char*)buffer,ori_size*ori_size*sizeof(float),1,mrcsFile) == NULL){
            std::cerr<<"read file failed."<<std::endl;
			EXIT_ABNORMALLY;
        }
        
        // copy
        auto image_ptr = images[iimage].wptr(ori_size*ori_size);
        for (int i = 0; i < ori_size*ori_size; i++) {
            image_ptr[i] = buffer[i];
        }
        
    }
    
    aFree(buffer);
    
    fclose(mrcsFile);
    
    // ------------ initialize sampling ---------- //
    
    // Initialise the sampling object (sets prior mode and fills translations and rotations inside sampling object)
    sampler2d.initialize(offset_step,offset_range,psi_step);
    
    // Also randomize random-number-generator for perturbations on the angles
    dontShare_Random_generator.init(random_seed);
    
}
    
void prepare()
{
    maxthreads = omp_get_max_threads();
    NODE0ONLY std::cout<<"maxthreads = "<<maxthreads<<std::endl;
    
	// --------------- initialize Particle Model ----------------- //
    if (particle_diameter < 0.)
        particle_diameter = (ori_size - width_mask_edge) * pixel_size;
    
    particleModel.initialize(ori_size, pixel_size, particle_diameter, width_mask_edge,
                             sigma2_fudge, random_seed, do_norm_correction, do_zero_mask,
                             false,maxthreads,&global_data_stream);
    
    // --------------- initialize MAP Model ---------------- //
    // For do_average_unaligned, always use initial low_pass filter
    // By default, use 0.07 dig.freq. low-pass filter
    // See S.H.W. Scheres (2010) Meth Enzym.
    ini_high = 1./getResolution(round(0.07 * ori_size));
    
    // ini_high ,user set
    mapModel.initialize(nr_classes,ori_size,particle_diameter,pixel_size,ini_high,maxthreads,
                        width_mask_edge,width_fmask_edge,5,true,do_map,2);
    // ------------ initialize model --------------- //
    // do_scale_correction = false;
    // NOTE : not do scale correction for 2D
    assert(sampler2d.NrDir()==1);
    mlModel.initialize(ori_size, nr_classes, metadata.numberOfGroups(), sampler2d.NrDir());
    
    // Calculate initial sigma noise model from power_class spectra of the individual images
    Image Mavg;
    Mavg.init(ori_size*ori_size); Mavg.zero();
    
#ifdef DATA_STREAM
    //TODO : check why different with relion.
    global_data_stream.foutInt(metadata.numberOfGroups(), "numberOfGroups", __FILE__, __LINE__);
    global_data_stream.foutInt(metadata.numberOfMicrographs(), "numberOfMicrographs", __FILE__, __LINE__);
    global_data_stream.foutInt(metadata.numberOfParticles(), "numberOfOriginalParticles", __FILE__, __LINE__);
    global_data_stream.foutDouble(Mavg.wptrAll(), ori_size*ori_size, "ref_1", __FILE__, __LINE__);// zero
    global_data_stream.foutDouble(Mavg.wptrAll(), ori_size*ori_size, "ref_1", __FILE__, __LINE__);// zero
    global_data_stream.foutInt(0, "mymodel.orientational_prior_mode", __FILE__, __LINE__);
    global_data_stream.foutInt(mapModel.ref_dim, "mymodel.ref_dim", __FILE__, __LINE__);
    global_data_stream.foutInt(sampler2d.NrDir(), "sampling.NrDirections()", __FILE__, __LINE__);
    global_data_stream.foutDouble(mlModel.pdf_direction[0].wptr(sampler2d.NrDir()), sampler2d.NrDir(), "mymodel.pdf_direction[0]", __FILE__, __LINE__);
    global_data_stream.foutDouble(mlModel.pdf_direction[nr_classes-1].wptr(sampler2d.NrDir()), sampler2d.NrDir(), "mymodel.pdf_direction_nr_class", __FILE__, __LINE__);
    global_data_stream.check();global_data_stream.flush();
#endif
    
    // check the data whether normalized
    if (true/* check_norm */) {
        checkNormalize(images, ori_size, particle_diameter, pixel_size);
    }
    
    // initialize Mavg,wsum_sigma2_noise,wsum_sumw_group,
    mlModel.calculateSumOfPowerSpectraAndAverageImage(Mavg, images, do_zero_mask, metadata, first_local_image, mapModel);
    
#ifdef USEMPI
    // reduce data
    auto local_temp = Heap::allocScalars<FDOUBLE>(std::max(mlModel.nr_groups,ori_size*ori_size), __FILE__, __LINE__);
    memcpy(local_temp, Mavg.rptrAll(), sizeof(FDOUBLE)*ori_size*ori_size);
    MPI::COMM_WORLD.Reduce(local_temp,Mavg.wptrAll(),ori_size*ori_size,MPI_FDOUBLE,MPI::SUM,0);
    
    // NOTE : because we only get average sigma2_noise,do not need to reduce all nr_groups's sigma2 noise
    memcpy(local_temp, mlModel.wsum_sigma2_noise[0].rptrAll(), sizeof(FDOUBLE)*(ori_size/2+1));
    MPI::COMM_WORLD.Reduce(local_temp,mlModel.wsum_sigma2_noise[0].wptrAll(),ori_size/2+1,MPI_FDOUBLE,MPI::SUM,0);
    
    memcpy(local_temp, mlModel.wsum_sumw_group.rptrAll(), sizeof(FDOUBLE)*mlModel.nr_groups);
    MPI::COMM_WORLD.Reduce(local_temp,mlModel.wsum_sumw_group.wptrAll(),mlModel.nr_groups,MPI_FDOUBLE,MPI::SUM,0);
    
    Heap::freeScalars(local_temp);
    
    auto local_temp2 = Heap::allocScalars<int>(mlModel.nr_groups, __FILE__, __LINE__);
    
    memcpy(local_temp2, mlModel.nr_particles_group.rptrAll(), sizeof(int)*mlModel.nr_groups);
    MPI::COMM_WORLD.Allreduce(local_temp2,mlModel.nr_particles_group.wptrAll(),mlModel.nr_groups,MPI::INT,MPI::SUM);
    
    Heap::freeScalars(local_temp2);

    if (false/*node == 1*/) {
        for (int i = 0; i < mlModel.nr_groups; i++) {
            std::cerr<<mlModel.nr_particles_group.rptrAll()[i]<<" ";
        }
        std::cerr<<std::endl;
    }
#endif

#ifdef DATA_STREAM
    global_data_stream.foutDouble(Mavg.wptrAll(), ori_size*ori_size, "Mavg", __FILE__, __LINE__);
    global_data_stream.foutDouble(mlModel.wsum_sigma2_noise[0].wptr(ori_size/2+1), ori_size/2+1, "wsum_model_sigma2_noise", __FILE__, __LINE__);
    global_data_stream.foutDouble(mlModel.wsum_sumw_group.wptrAll()[0], "wsum_model_sumw_group1", __FILE__, __LINE__);
    global_data_stream.foutDouble(mlModel.wsum_sumw_group.wptrAll()[mlModel.nr_groups-1], "wsum_model_sumw_groupN", __FILE__, __LINE__);
    global_data_stream.foutInt((int)mlModel.nr_particles_group.wptrAll()[0], "mymodel_sumw_nr_particles_group1", __FILE__, __LINE__);
    global_data_stream.foutInt((int)mlModel.nr_particles_group.wptrAll()[mlModel.nr_groups-1], "mymodel_sumw_nr_particles_groupN", __FILE__, __LINE__);
    global_data_stream.check();global_data_stream.flush();
#endif
    
    NODE0ONLY{
        // First calculate average image
        for (int i = 0; i < ori_size*ori_size; i++)
            Mavg[i] /= mlModel.wsum_sumw_group.rptrAll()[0];
        mapModel.initializeRef(Mavg);
    }
    // Set model_sigma2_noise and model_Iref from averaged poser spectra and Mavg
    NODE0ONLY mlModel.setSigmaNoiseEstimates(Mavg);
    
#ifdef DATA_STREAM
    global_data_stream.foutDouble(mlModel.sigma2_noise[0].wptr(ori_size/2+1), ori_size/2+1, "sigma2_noise1", __FILE__, __LINE__);
    global_data_stream.foutDouble(mlModel.sigma2_noise[mlModel.nr_groups-1].wptr(ori_size/2+1), ori_size/2+1, "sigma2_noiseN", __FILE__, __LINE__);
    global_data_stream.foutInt(ori_size/2+1, "sigma2_size", __FILE__, __LINE__);
    global_data_stream.check();global_data_stream.flush();
#endif
    
    // First low-pass filter the initial references
    NODE0ONLY mapModel.applyLowPassFilter();
    
#ifdef DATA_STREAM
    global_data_stream.foutDouble(mapModel.Irefs[0].wptr(), ori_size*ori_size, "ref_1", __FILE__, __LINE__);
    global_data_stream.foutDouble(mapModel.Irefs[nr_classes-1].wptr(), ori_size*ori_size, "ref_N", __FILE__, __LINE__);
    global_data_stream.check();global_data_stream.flush();
#endif
    
    // Initialise the model_data_versus_prior ratio to get the initial current_size right
    // model_tau2_class,model_data_vs_prior_class
    NODE0ONLY mlModel.initialiseDataVersusPrior(mapModel,tau2_fudge_factor);
    
#ifdef DATA_STREAM
    global_data_stream.foutDouble(mlModel.data_vs_prior_class[0].wptr(ori_size/2+1), ori_size/2+1, "data_vs_prior_class_1", __FILE__, __LINE__);
    global_data_stream.foutDouble(mlModel.data_vs_prior_class[nr_classes-1].wptr(ori_size/2+1), ori_size/2+1, "data_vs_prior_class_N", __FILE__, __LINE__);
    global_data_stream.check();global_data_stream.flush();
#endif
    
#ifdef USEMPI
    for (int igroup = 0; igroup < mlModel.nr_groups; igroup++)
        MPI::COMM_WORLD.Bcast(mlModel.sigma2_noise[igroup].wptrAll(),mlModel.ori_Fsize,MPI_FDOUBLE,0);
    for (int iclass = 0; iclass < nr_classes; iclass++){
        MPI::COMM_WORLD.Bcast(mapModel.Irefs[iclass].wptr(),mapModel.Irefs[iclass].dimzyx,MPI_FDOUBLE,0);
        MPI::COMM_WORLD.Bcast(mlModel.tau2_class[iclass].wptrAll(),mlModel.ori_Fsize,MPI_FDOUBLE,0);
        MPI::COMM_WORLD.Bcast(mlModel.data_vs_prior_class[iclass].wptrAll(),mlModel.ori_Fsize,MPI_FDOUBLE,0);
    }
#endif
    Mavg.fini();
}

void destroyMap2dOptimizer(){
    images			.resize(0);
    particleModel	.finalize();
    mapModel		.finalize();
    mlModel			.finalize();
}

StatusTracer statusTracer;

void setupStatusTracer()
{
    statusTracer.clear();
    for (int i = 0; i < mapModel.Irefs.size(); i++) {
#if defined(FLOAT_PRECISION)
        statusTracer.appendFloatPtr(mapModel.Irefs[i].wptr(), mapModel.Irefs[i].dimzyx, "mapmodel_Irefs", true);
#else
        statusTracer.appendDoublePtr(mapModel.Irefs[i].wptr(), mapModel.Irefs[i].dimzyx, "mapmodel_Irefs");
#endif
    }
    statusTracer.appendDoublePtr(&mapModel.current_resolution, 1, "mapModel_current_resolution");
    //
    for (int igroup = 0; igroup < mlModel.nr_groups; igroup++) {
#if defined(FLOAT_PRECISION)
        statusTracer.appendFloatPtr(mlModel.sigma2_noise[igroup].wptr(mlModel.ori_Fsize), mlModel.ori_Fsize, "mlmodel_sigma2_noise", true);
        statusTracer.appendFloatPtr(mlModel.wsum_signal_product_spectra[igroup].wptr(mlModel.ori_Fsize), mlModel.ori_Fsize, "mlmodel_wsum_signal_product_spectra", true);
        statusTracer.appendFloatPtr(mlModel.wsum_reference_power_spectra[igroup].wptr(mlModel.ori_Fsize), mlModel.ori_Fsize, "mlmodel_wsum_reference_power_spectra", true);
#else
        statusTracer.appendDoublePtr(mlModel.sigma2_noise[igroup].wptr(mlModel.ori_Fsize),
                                     mlModel.ori_Fsize, "mlmodel_sigma2_noise_"+num2str(igroup,6));
        statusTracer.appendDoublePtr(mlModel.wsum_signal_product_spectra[igroup].wptr(mlModel.ori_Fsize),
                                     mlModel.ori_Fsize, "mlmodel_wsum_signal_product_spectra_"+num2str(igroup,6));
        statusTracer.appendDoublePtr(mlModel.wsum_reference_power_spectra[igroup].wptr(mlModel.ori_Fsize),
                                     mlModel.ori_Fsize, "mlmodel_wsum_reference_power_spectra_"+num2str(igroup,6));
#endif
    }
    for (int iclass = 0; iclass < mlModel.nr_classes; iclass++) {
#if defined(FLOAT_PRECISION)
        statusTracer.appendFloatPtr(mlModel.tau2_class[iclass].wptr(mlModel.ori_Fsize), mlModel.ori_Fsize, "mlmodel_tau2_class", true);
        statusTracer.appendFloatPtr(mlModel.sigma2_class[iclass].wptr(mlModel.ori_Fsize), mlModel.ori_Fsize, "mlmodel_sigma2_class", true);
        statusTracer.appendFloatPtr(mlModel.data_vs_prior_class[iclass].wptr(mlModel.ori_Fsize), mlModel.ori_Fsize, "mlmodel_data_vs_prior_class", true);
        statusTracer.appendFloatPtr(mlModel.fsc_halves_class[iclass].wptr(mlModel.ori_Fsize), mlModel.ori_Fsize, "mlmodel_fsc_halves_class", true);
        statusTracer.appendFloatPtr(mlModel.pdf_direction[iclass].wptr(mlModel.nr_directions), mlModel.nr_directions, "mlmodel_pdf_direction", true);
#else
        statusTracer.appendDoublePtr(mlModel.tau2_class[iclass].wptr(mlModel.ori_Fsize),
                                     mlModel.ori_Fsize, "mlmodel_tau2_class_"+num2str(iclass,6));
        statusTracer.appendDoublePtr(mlModel.sigma2_class[iclass].wptr(mlModel.ori_Fsize),
                                     mlModel.ori_Fsize, "mlmodel_sigma2_class_"+num2str(iclass,6));
        statusTracer.appendDoublePtr(mlModel.data_vs_prior_class[iclass].wptr(mlModel.ori_Fsize),
                                     mlModel.ori_Fsize, "mlmodel_data_vs_prior_class_"+num2str(iclass,6));
        statusTracer.appendDoublePtr(mlModel.fsc_halves_class[iclass].wptr(mlModel.ori_Fsize),
                                     mlModel.ori_Fsize, "mlmodel_fsc_halves_class_"+num2str(iclass,6));
        statusTracer.appendDoublePtr(mlModel.pdf_direction[iclass].wptr(mlModel.nr_directions),
                                     mlModel.nr_directions, "mlmodel_pdf_direction_"+num2str(iclass,6));
#endif
    }
#if defined(FLOAT_PRECISION)
    statusTracer.appendFloatPtr(mlModel.scale_correction.wptr(mlModel.nr_groups), mlModel.nr_groups, "mlmodel_scale_correction", true);
    statusTracer.appendFloatPtr(mlModel.prior_offsetx_class.wptr(mlModel.nr_classes), mlModel.nr_classes, "mlmodel_prior_offsetx_class", true);
    statusTracer.appendFloatPtr(mlModel.prior_offsety_class.wptr(mlModel.nr_classes), mlModel.nr_classes, "mlmodel_prior_offsety_class", true);
    statusTracer.appendFloatPtr(mlModel.pdf_class.wptr(mlModel.nr_classes), mlModel.nr_classes, "mlmodel_pdf_class", true);
    statusTracer.appendFloatPtr(&mlModel.avg_norm_correction, 1, "mlmodel_avg_norm_correction", true);
    statusTracer.appendFloatPtr(&mlModel.ave_Pmax, 1, "mlmodel_ave_Pmax", true);
    statusTracer.appendFloatPtr(&mlModel.sigma2_offset, 1, "mlmodel_sigma2_offset", true);
#else
    statusTracer.appendDoublePtr(mlModel.scale_correction.wptr(mlModel.nr_groups),
                                 mlModel.nr_groups, "mlmodel_scale_correction");
    statusTracer.appendDoublePtr(mlModel.prior_offsetx_class.wptr(mlModel.nr_classes),
                                 mlModel.nr_classes, "mlmodel_prior_offsetx_class");
    statusTracer.appendDoublePtr(mlModel.prior_offsety_class.wptr(mlModel.nr_classes),
                                 mlModel.nr_classes, "mlmodel_prior_offsety_class");
    statusTracer.appendDoublePtr(mlModel.pdf_class.wptr(mlModel.nr_classes),
                                 mlModel.nr_classes, "mlmodel_pdf_class");
    statusTracer.appendDoublePtr(&mlModel.avg_norm_correction, 1, "mlmodel_avg_norm_correction");
    statusTracer.appendDoublePtr(&mlModel.ave_Pmax, 1, "mlmodel_ave_Pmax");
    statusTracer.appendDoublePtr(&mlModel.sigma2_offset, 1, "mlmodel_sigma2_offset");
#endif
    for (int igroup = 0; igroup < mlModel.nr_groups; igroup++) {
#if defined(FLOAT_PRECISION)
        statusTracer.appendFloatPtr(mlModel.wsum_sigma2_noise[igroup].wptr(mlModel.ori_Fsize), mlModel.ori_Fsize, "mlmodel_wsum_sigma2_noise", true);
#else
        statusTracer.appendDoublePtr(mlModel.wsum_sigma2_noise[igroup].wptr(mlModel.ori_Fsize),
                                     mlModel.ori_Fsize, "mlmodel_wsum_sigma2_noise_"+num2str(igroup,6));
#endif
    }
#if defined(FLOAT_PRECISION)
    statusTracer.appendFloatPtr(mlModel.wsum_sumw_group.wptr(mlModel.nr_groups), mlModel.nr_groups, "mlmodel_wsum_sumw_group", true);
    statusTracer.appendFloatPtr(mlModel.wsum_pdf_class.wptr(mlModel.nr_classes), mlModel.nr_classes, "mlmodel_wsum_pdf_class", true);
    statusTracer.appendFloatPtr(mlModel.wsum_prior_offsetx_class.wptr(mlModel.nr_classes), mlModel.nr_classes, "mlmodel_wsum_prior_offsetx_class", true);
    statusTracer.appendFloatPtr(mlModel.wsum_prior_offsety_class.wptr(mlModel.nr_classes), mlModel.nr_classes, "mlmodel_wsum_prior_offsety_class", true);
#else
    statusTracer.appendDoublePtr(mlModel.wsum_sumw_group.wptr(mlModel.nr_groups),
                                 mlModel.nr_groups, "mlmodel_wsum_sumw_group");
    statusTracer.appendDoublePtr(mlModel.wsum_pdf_class.wptr(mlModel.nr_classes),
                                 mlModel.nr_classes, "mlmodel_wsum_pdf_class");
    statusTracer.appendDoublePtr(mlModel.wsum_prior_offsetx_class.wptr(mlModel.nr_classes),
                                 mlModel.nr_classes, "mlmodel_wsum_prior_offsetx_class");
    statusTracer.appendDoublePtr(mlModel.wsum_prior_offsety_class.wptr(mlModel.nr_classes),
                                 mlModel.nr_classes, "mlmodel_wsum_prior_offsety_class");
#endif
    for (int iclass = 0; iclass < mlModel.nr_classes; iclass++) {
#if defined(FLOAT_PRECISION)
        statusTracer.appendFloatPtr(mlModel.wsum_pdf_direction[iclass].wptr(mlModel.nr_directions), mlModel.nr_directions, "mlmodel_wsum_pdf_direction", true);
#else
        statusTracer.appendDoublePtr(mlModel.wsum_pdf_direction[iclass].wptr(mlModel.nr_directions),
                                     mlModel.nr_directions, "mlmodel_wsum_pdf_direction"+num2str(iclass,6));
#endif
    }
    //
#if defined(FLOAT_PRECISION)
    statusTracer.appendFloatPtr(&sampler2d.random_perturbation,1,"sampling3d.random_perturbation",true);
#else
    statusTracer.appendDoublePtr(&sampler2d.random_perturbation,1,"sampling3d.random_perturbation");
#endif
    //
    for (int iimage = 0; iimage < nr_global_images; iimage++) {
#if defined(FLOAT_PRECISION)
        statusTracer.appendFloatPtr(&metadata[iimage].NORM, 1, "metadata["+num2str(iimage,6)+"].NORM", true);
        statusTracer.appendFloatPtr(&metadata[iimage].XOFF, 1, "metadata["+num2str(iimage,6)+"].XOFF", true);
        statusTracer.appendFloatPtr(&metadata[iimage].YOFF, 1, "metadata["+num2str(iimage,6)+"].YOFF", true);
        statusTracer.appendFloatPtr(&metadata[iimage].PSI, 1, "metadata["+num2str(iimage,6)+"].PSI", true);
#else
        statusTracer.appendDoublePtr(&metadata[iimage].NORM, 1, "metadata["+num2str(iimage,6)+"].NORM");
        statusTracer.appendDoublePtr(&metadata[iimage].XOFF, 1, "metadata["+num2str(iimage,6)+"].XOFF");
        statusTracer.appendDoublePtr(&metadata[iimage].YOFF, 1, "metadata["+num2str(iimage,6)+"].YOFF");
        statusTracer.appendDoublePtr(&metadata[iimage].PSI, 1, "metadata["+num2str(iimage,6)+"].PSI");
#endif
        statusTracer.appendIntPtr(&metadata[iimage].CLASS, 1, "metadata["+num2str(iimage,6)+"].CLASS");
    }
}
    
// ------------------------- EM-Iteration  ------------------------- //
void iterate()
{
	// Update the current resolution and image sizes, and precalculate resolution pointers
	// The rest of the time this will be done after maximization and before writing output files,
	// so that current resolution is in the output files of the current iteration
    bool set_by_ini_high = ini_high > 0. && (iter == 0);
    static int nr_iter_wo_resol_gain = 0;
    mapModel.updateCurrentResolution(mlModel.data_vs_prior_class,set_by_ini_high,nr_iter_wo_resol_gain);
    
#ifdef DATA_STREAM
    global_data_stream.foutDouble(mapModel.current_resolution, "updateCurrentResolution()_current_resolution", __FILE__, __LINE__);
#endif
    
    NODE0ONLY std::cout<<"current_resolution = "<<mapModel.current_resolution<<std::endl;
	bool has_already_reached_convergence = false;
    //
    setupStatusTracer();
    
    // continue
    if (iter != 0) {
        NODE0ONLY std::cout<<"------------------ from "<<iter<<"  --------------------"<<std::endl;
        std::string statusFn = pathRemoveSuffix(star_fn);
        statusTracer.recoveryStatus(statusFn);
        // After the first iteration the references are always CTF-corrected
        if (do_ctf_correction)
            refs_are_ctf_corrected = true;
    }
	for (iter = iter + 1; iter <= nr_iter; iter++)
    {
#ifdef DATA_STREAM
        global_data_stream.foutInt(iter, "iterate()_iter", __FILE__, __LINE__);
#endif
        const double starttime = dtime();
        // update coarse_size,current_size,Npix_per_shell,Mresol_coarse,Mresol_fine
        double angularSampler = sampler2d.getAngularSampling();
        mapModel.updateImageSizeAndResolutionPointers(Npix_per_shell, Mresol_coarse, Mresol_fine,
                                                      coarse_size, current_size,
                                                      adaptive_oversampling, angularSampler,
                                                      mlModel.ave_Pmax);
#ifdef DATA_STREAM
        global_data_stream.foutInt(current_size, "updateImageSizeAndResolutionPointers()_current_size", __FILE__, __LINE__);
        global_data_stream.foutInt(coarse_size, "updateImageSizeAndResolutionPointers()_coarse_size", __FILE__, __LINE__);
        global_data_stream.foutInt(Npix_per_shell.wptrAll(), (ori_size/2+1), "updateImageSizeAndResolutionPointers()_Npix_per_shell", __FILE__, __LINE__);
        global_data_stream.foutInt(Mresol_fine.wptrAll(), current_size*(current_size/2+1), "updateImageSizeAndResolutionPointers()_Mresol_fine", __FILE__, __LINE__);
        global_data_stream.foutInt(Mresol_coarse.wptrAll(), coarse_size*(coarse_size/2+1), "updateImageSizeAndResolutionPointers()_Mresol_coarse", __FILE__, __LINE__);
        global_data_stream.check();global_data_stream.flush();
#endif
        NODE0ONLY std::cout<<"current_size = "<<current_size<<" coarse_size = "<<coarse_size<<std::endl;
        
        // set the exp_metadata and exp_image_data and some other data
        // allocate the maximum data
        prepareExpData();

		expectation();
        
#ifdef USEMPI
        
        mlModel.reduceData(MPI::COMM_WORLD);
        
        mapModel.reduceData(MPI::COMM_WORLD);
        
        gatherMetaDataToMaster(metadata);
        
        MPI::COMM_WORLD.Barrier();
        
#endif
        
		maximization();

		// Apply masks to the reference images
		// At the last iteration, do not mask the map for validation purposes
		if (do_solvent && !has_converged)
			mapModel.applySolventFlatten("NULL");
        
#ifdef DATA_STREAM
        global_data_stream.foutDouble(mapModel.Irefs[0].wptr(), mapModel.Irefs[0].dimzyx, "iterate()_Iref[0]", __FILE__, __LINE__);
        global_data_stream.foutDouble(mapModel.Irefs[nr_classes-1].wptr(), mapModel.Irefs[nr_classes-1].dimzyx, "iterate()_Iref[nr_classes-1]", __FILE__, __LINE__);
#endif
        
		// Re-calculate the current resolution, do this before writing to get the correct values in the output files
        bool set_by_ini_high = false;
		mapModel.updateCurrentResolution(mlModel.data_vs_prior_class,set_by_ini_high,nr_iter_wo_resol_gain);
        
#ifdef DATA_STREAM
        global_data_stream.foutDouble(mapModel.current_resolution, "updateCurrentResolution()_current_resolution", __FILE__, __LINE__);
#endif
        
        NODE0ONLY{
        	if(guidance_fn!= "NULL" && guidance_fn != "BACK")
                statusTracer.checkStatus(guidance_fn+"_iter"+num2str(iter));
        }
		// Write output files
        if(write_path != "NULL"){
            std::string iterStr = num2str(iter);
            NODE0ONLY {
                std::string fn_class    = write_path+write_fn+"_iter"+iterStr;
                std::string fn_metadata = write_path+write_fn+"_iter"+iterStr;
                writeClassesAndMetadata(fn_class,fn_metadata);
                std::string statusFn = write_path+write_fn+"_iter"+iterStr;
				if (guidance_fn == "BACK")
					statusTracer.backupStatus(statusFn);
            }
        }
        // update the metadata,free exp_metadata and exp_image_data
        endExpData();
        const double endtime = dtime();
        NODE0ONLY std::cout<<"**** iteration "<<iter<<" completed in "<<endtime-starttime<<" seconds ****"<<std::endl<<std::flush;
    }
    
}

void prepareExpData()
{
    
#ifdef DATA_STREAM
    // ------------ initialize model and wsum ---------------- //
    global_data_stream.foutDouble(mlModel.tau2_class[0].wptr(ori_size/2+1), ori_size/2+1, "expectationSetup()_tau2_class", __FILE__, __LINE__);
    global_data_stream.foutDouble(mlModel.tau2_class[nr_classes-1].wptr(ori_size/2+1), ori_size/2+1, "expectationSetup()_tau2_class", __FILE__, __LINE__);
    global_data_stream.check();global_data_stream.flush();
#endif
    
    // ------------ initialize model and wsum ---------------- //
    mapModel.setFourierTransformMaps(mlModel.tau2_class, current_size);
    
#ifdef DATA_STREAM
    global_data_stream.foutDouble(mlModel.tau2_class[0].wptr(ori_size/2+1), ori_size/2+1, "expectationSetup()_tau2_class", __FILE__, __LINE__);
    global_data_stream.foutDouble(mlModel.tau2_class[nr_classes-1].wptr(ori_size/2+1), ori_size/2+1, "expectationSetup()_tau2_class", __FILE__, __LINE__);
    global_data_stream.check();global_data_stream.flush();
#endif
    
    // Initialise all weighted sums to zero
    mlModel.resetZero();
    
    dontShare_Random_generator.init(random_seed + iter);
    // Reset the random perturbation for this sampling
    sampler2d.resetRandomlyPerturbedSampling();
    
    NODE0ONLY {
        // check error
        int n_trials_acc = 10;
        n_trials_acc = std::min(n_trials_acc, (int)nr_local_images);
        // calculateExpectedAngularErrors(0, n_trials_acc-1);
    }
    // the first time we use small nr_pool to initialize model_Iref,
    // in order to peak much performance,set it larger after iter 1.
    // if(iter == 2) nr_pool = std::min(300,nr_local_images);
    NODE0ONLY calspace(ori_size,current_size,nr_pool);
    
    int current_Fsize2 = current_size*(current_size/2 + 1);
	int nr_images = nr_pool;
    exp_highres_Xi2_imgs		.init(nr_images);
    exp_min_diff2				.init(nr_images);
    exp_sum_weight				.init(nr_images);
    exp_significant_weight		.init(nr_images);
    exp_old_offsetx				.init(nr_images);
    exp_old_offsety				.init(nr_images);
    exp_wsum_norm_correction	.init(nr_images);
    exp_local_sqrtXi2			.init(nr_images);
    exp_local_Fctfs				.init(nr_images, current_Fsize2);
    exp_local_Minvsigma2s		.init(nr_images, current_Fsize2);
    exp_imgs					.init(nr_images, ori_size*ori_size);
    exp_power_imgs				.init(nr_images, ori_size/2+1);
    //
    exp_metadata				.init(nr_images);
    
    // Number of rotational and translational sampling points
    samplingGrid.initialize(sampler2d, adaptive_oversampling);
    exp_nr_trans = samplingGrid.exp_nr_trans();
    exp_nr_psi = samplingGrid.exp_nr_psi();
    exp_nr_dir = samplingGrid.exp_nr_dir();
    
    //intialzie some shift images and its ctf
    int exp_nr_all_trans    	= sampler2d.NrTrans(adaptive_oversampling);
    exp_Fimgs_shifted_real		.init(nr_images, exp_nr_all_trans, current_Fsize2);
    exp_Fimgs_shifted_imag		.init(nr_images, exp_nr_all_trans, current_Fsize2);
    
    int exp_nr_all_rot      	= sampler2d.NrPsi(adaptive_oversampling);
    exp_Frefs_Rot_real			.init(nr_classes, exp_nr_all_rot, current_Fsize2);
    exp_Frefs_Rot_imag			.init(nr_classes, exp_nr_all_rot, current_Fsize2);
    exp_Fweight_Rot				.init(nr_classes, exp_nr_all_rot, current_Fsize2);
    
    exp_Rot_significant			.init(nr_classes, exp_nr_psi);
    exp_Mcoarse_significant		.init(nr_images, nr_classes*sampler2d.NrPoints(0));
    exp_Mweight             	= (double*)aMalloc(sizeof(double)*nr_images*nr_classes*sampler2d.NrPoints(adaptive_oversampling),64);
    
    // some thread variable
    thread_Frefctf_real				.init(maxthreads, current_Fsize2);
    thread_Frefctf_imag				.init(maxthreads, current_Fsize2);
    
    // ori_sizexori_size for image
    thread_exp_max_weight_index		.resize(maxthreads);
    for (auto& index : thread_exp_max_weight_index) index.resize(nr_images);
    thread_sumw_group				.init(maxthreads, nr_images);
    thread_exp_min_diff2			.init(maxthreads, nr_images);
    thread_exp_sum_weight			.init(maxthreads, nr_images);
    thread_exp_max_weight			.init(maxthreads, nr_images);
    thread_wsum_norm_correction		.init(maxthreads, nr_images);
    thread_wsum_sigma2_noise		.init(maxthreads, nr_images, current_Fsize2);
    thread_wsum_pdf_direction		.init(maxthreads, nr_classes);
    thread_wsum_pdf_class			.init(maxthreads, nr_classes);
    thread_wsum_prior_offsetx_class	.init(maxthreads, nr_classes);
    thread_wsum_prior_offsety_class	.init(maxthreads, nr_classes);
    thread_wsum_sigma2_offset		.init(maxthreads, 1);
    //
    particleModel.setup(nr_images, current_size, coarse_size);
}

void endExpData()
{
    //
    Npix_per_shell					.fini();
    Mresol_coarse					.fini();
    Mresol_fine						.fini();
    //
    exp_highres_Xi2_imgs			.fini();
    exp_min_diff2					.fini();
    exp_sum_weight					.fini();
    exp_significant_weight			.fini();
    exp_old_offsetx					.fini();
    exp_old_offsety					.fini();
    exp_wsum_norm_correction		.fini();
    exp_local_sqrtXi2				.fini();
    exp_local_Fctfs					.fini();
    exp_local_Minvsigma2s			.fini();
    exp_imgs						.fini();
    exp_power_imgs					.fini();
    //
    exp_metadata					.fini();
    
    samplingGrid					.finalize();
    
    exp_Fimgs_shifted_real			.fini();
    exp_Fimgs_shifted_imag			.fini();
    
	//
    exp_Frefs_Rot_real				.fini();
    exp_Frefs_Rot_imag				.fini();
    exp_Fweight_Rot					.fini();
    //
    exp_Rot_significant				.fini();
    exp_Mcoarse_significant			.fini();
    aFree(exp_Mweight);
    
    // some thread variable
    thread_exp_max_weight_index		.resize(0);
    thread_Frefctf_real				.fini();
    thread_Frefctf_imag				.fini();
    thread_exp_min_diff2			.fini();
    thread_sumw_group				.fini();
    thread_exp_sum_weight			.fini();
    thread_exp_max_weight			.fini();
    thread_wsum_norm_correction		.fini();
    thread_wsum_sigma2_noise		.fini();
    thread_wsum_pdf_direction		.fini();
    thread_wsum_pdf_class			.fini();
    thread_wsum_prior_offsetx_class	.fini();
    thread_wsum_prior_offsety_class	.fini();
    thread_wsum_sigma2_offset		.fini();
    //
    particleModel					.destroy();
}

void expectation()
{
    
    NODE0ONLY std::cout << " Expectation iteration " << iter<< " of " << nr_iter<<std::endl;
    
    int my_first_image,my_last_image,nr_images_done = 0;

    NODE0ONLY showProgressBar(0, nr_local_images);
    
    while (nr_images_done < nr_local_images)
    {
        my_first_image = first_local_image + nr_images_done;
        
        if (iter == 1){
            // divied the images that can be equally split by nr_classes
            // only use in first iteration to generate reference seed
            int iclass = divide_equally_which_group(nr_global_images, nr_classes, my_first_image);
            int first,last;
            divide_equally(nr_global_images, nr_classes, iclass, first, last);
            int suitable_pool = std::min(nr_pool, last-my_first_image+1);
            my_last_image = std::min(first_local_image+nr_local_images - 1, my_first_image + suitable_pool - 1);
        }
        else{
            my_last_image = std::min(first_local_image+nr_local_images - 1, my_first_image + nr_pool - 1);
        }
        exp_first_image = my_first_image;
        exp_last_image = my_last_image;
        exp_nr_images = my_last_image - my_first_image + 1;
        
#ifdef DATA_STREAM
        global_data_stream.foutInt(exp_first_image, "expectation()_my_first_ori_particle", __FILE__, __LINE__);
        global_data_stream.foutInt(exp_last_image, "expectation()_my_last_ori_particle", __FILE__, __LINE__);
        global_data_stream.foutInt(exp_nr_images, "expectation()_nr_pool", __FILE__, __LINE__);
#endif
        
        // first time divided the images by nr_classes peiece,make the initialized reference different
        exp_iclass_min = 0;
        exp_iclass_max = nr_classes - 1;
        if (iter == 1){
            exp_iclass_min = exp_iclass_max = divide_equally_which_group(nr_global_images, nr_classes, exp_first_image);
        }
        
#ifdef DATA_STREAM
        global_data_stream.foutInt(0, "expectationOneParticle()_do_firstiter_cc", __FILE__, __LINE__);
        global_data_stream.foutInt(exp_iclass_min, "expectationOneParticle()_exp_iclass_min", __FILE__, __LINE__);
        global_data_stream.foutInt(exp_iclass_max, "expectationOneParticle()_exp_iclass_max", __FILE__, __LINE__);
        global_data_stream.check();global_data_stream.flush();
#endif
        
        // prepare exp_image_data and exp_metadata
        // each node keep nr_local_images's images data and all(nr_global_images) metadata
        for (int iimage = 0; iimage < exp_nr_images; iimage++)
        {
            auto exp_img_ptr = exp_imgs[iimage].wptr(ori_size*ori_size);
            auto image_ptr = images[(my_first_image-first_local_image)+iimage].wptrAll();
            for (int i = 0; i < ori_size*ori_size; i++)
                exp_img_ptr[i] = image_ptr[i];
            //
            exp_metadata[iimage] = metadata[iimage+my_first_image];
        }
        // get Image(with mask and nomask) and CTF in Fourier Space
        // it will downsize the image from ori_size to current_size
        particleModel.prepare(exp_imgs, exp_metadata, exp_first_image, exp_nr_images);
        particleModel.setFourierTransforms(exp_power_imgs, exp_highres_Xi2_imgs,
                                           exp_old_offsetx, exp_old_offsety,
                                           mlModel,do_ctf_correction);
        
        // perform the actual expectation step on several particles
        expectationSomeParticles();
        
        // Also monitor the changes in the optimal orientations and classes
        hiddenVarMonitor.monitorHiddenVariableChanges(sampler2d, metadata, my_first_image, &exp_metadata[0],exp_nr_images);
        
        // update the metadata
        for (int iimage = 0; iimage < exp_nr_images; iimage++) metadata[iimage+my_first_image] = exp_metadata[iimage];
        //
        nr_images_done += my_last_image - my_first_image + 1;
        
        NODE0ONLY showProgressBar(nr_images_done, nr_local_images);
    }
    NODE0ONLY showProgressBar(nr_local_images, nr_local_images);
}


void expectationSomeParticles()
{

	// Only perform a second pass when using adaptive oversampling
	int nr_sampling_passes = (adaptive_oversampling > 0) ? 2 : 1;
    
	// Pass twice through the sampling of the entire space of rot, tilt and psi
	// The first pass uses a coarser angular sampling and possibly smaller FFTs than the second pass.
	// Only those sampling points that contribute to the highest x% of the weights in the first pass are oversampled in the second pass
	// Only those sampling points will contribute to the weighted sums in the third loop below
	for (exp_ipass = 0; exp_ipass < nr_sampling_passes; exp_ipass++)
	{
        // Use smaller images in the first pass, larger ones in the second pass
        exp_current_size = (exp_ipass == 0 && adaptive_oversampling > 0) ? coarse_size : current_size;
        
		// Use coarse sampling in the first pass, oversampled one the second pass
		exp_current_oversampling = (exp_ipass == 0) ? 0 : adaptive_oversampling;
        samplingGrid.computeGrid2D(sampler2d, exp_current_oversampling);
        exp_nr_over_rot = samplingGrid.exp_nr_over_rot();
        exp_nr_over_trans = samplingGrid.exp_nr_over_trans();
        
        // get all shifted Fimg,Fimg_nomask and CTF,also get inversed sigma^2
        bool do_cross_correlation = false;
        particleModel.preShiftedImagesCtfsAndInvSigma2s(exp_Fimgs_shifted_real, exp_Fimgs_shifted_imag,
                                                        exp_local_Minvsigma2s, exp_local_sqrtXi2, exp_local_Fctfs,
                                                        exp_current_size, do_cross_correlation, exp_ipass == 0,
                                                        samplingGrid, mlModel, Mresol_coarse, Mresol_fine);
        
        // get all rotated reference  and the significant rotation
        getReferenceAllOrientations();
        
        // get all reference and images 's squared differences
		getAllSquaredDifferences();

		// convert all squared differences to weight,and find significant(maximum) weight for each image
        convertSquaredDifferencesToWeights();
        
		findAllSignificantPoints();

	}// end loop over 2 exp_ipass iterations
    
    //
    // last ipass iteration,update some parameter for maximization
    // For the reconstruction step use model current_size!
    exp_current_size = current_size;
    
    updateOtherParams();
    
    backProjection();
    
    storeWeightedSums();

}
    
void getReferenceAllOrientations()
{
    int exp_current_Fsize2 = exp_current_size*(exp_current_size/2+1);

#pragma omp parallel for collapse(2)
    for (int iclass = exp_iclass_min; iclass <= exp_iclass_max; iclass++)
    {
        for (int ipsi = 0; ipsi < exp_nr_psi; ipsi++)
        {
            if (mlModel.pdf_class.rptrAll()[iclass] <= 0.)
                continue;
            
            int iorientclass_offset = iclass * exp_nr_psi;
            int iorientclass = iorientclass_offset + ipsi;
            // In the first pass, always proceed
            // In the second pass, check whether one of the translations for this orientation of any of
            // the particles had a significant weight in the first pass
            // if so, proceed with projecting the reference in that direction
            if (exp_ipass == 0)
                exp_Rot_significant[iclass].wptrAll()[ipsi] = true;
            else
                exp_Rot_significant[iclass].wptrAll()[ipsi] = isSignificantAnyParticleAnyTranslation(iorientclass);
        }// end loop of ipsi
    }// end loop of iclass
    
    
#pragma omp parallel for collapse(3)
    for (int iclass = exp_iclass_min; iclass <= exp_iclass_max; iclass++)
    {
        for (int ipsi = 0; ipsi < exp_nr_psi; ipsi++)
        {
            for (int iover_rot = 0; iover_rot < exp_nr_over_rot; iover_rot++)
            {
                auto pdf_orientation = mlModel.pdf_direction[iclass][0];
                
                if (!exp_Rot_significant[iclass].rptrAll()[ipsi] || pdf_orientation <= 0. || mlModel.pdf_class.rptrAll()[iclass] <= 0.)
                    continue;
                auto exp_Frefs_Rot_real_aux = exp_Frefs_Rot_real.wptr(iclass, (ipsi*exp_nr_over_rot+iover_rot), exp_current_Fsize2);
                auto exp_Frefs_Rot_imag_aux = exp_Frefs_Rot_imag.wptr(iclass, (ipsi*exp_nr_over_rot+iover_rot), exp_current_Fsize2);
                // transposition matrix
                FDOUBLE A[3][3];
                // Take tilt-series into account
                Euler_angles2matrix(0, 0, samplingGrid.exp_over_psi[ipsi][iover_rot], A);
                mapModel.get2DFourierTransform(iclass, exp_Frefs_Rot_real_aux, exp_Frefs_Rot_imag_aux, exp_current_size, A, false);
            }// end loop of iover_rot
        }//end loop of ipsi
    }// end loop of iclass
    
}


void getAllSquaredDifferences()
{
    int exp_current_Fsize2 = exp_current_size*(exp_current_size/2+1);
    
    // trick!the real needed coarse maybe smaller than exp_Mweight_xsize depend on two thing:
    // 1)in first EM-interation,not all class used in big while loop
    // 2)some not significant rot,direction,class,trans will be escapsed
    // so we need some minimum integer number to flag this
    exp_Mweight_xsize = nr_classes * sampler2d.NrPoints(exp_current_oversampling);
    assert(exp_Mweight_xsize > 0);
    for (int i = 0; i < exp_nr_images*exp_Mweight_xsize; i++) {
        exp_Mweight[i] = (std::numeric_limits<int>::min)();
    }
    thread_exp_min_diff2.fill_with_first_touch((std::numeric_limits<double>::max)());
    
#ifdef DATA_STREAM
    // turn off data stream inside openmp for loop
    // global_data_stream.turnOff();
    for (int iimage = 0; iimage < exp_nr_images; iimage++){
        for (int iclass = exp_iclass_min; iclass <= exp_iclass_max; iclass++){
            for (int ipsi = 0; ipsi < exp_nr_psi; ipsi++){
                for (int iover_rot = 0; iover_rot < exp_nr_over_rot; iover_rot++)
                {
                    if (!exp_Rot_significant[iclass].rptrAll()[ipsi]
                        || mlModel.pdf_direction[iclass].rptrAll()[0] <= 0.
                        || mlModel.pdf_class.rptrAll()[iclass] <= 0.) continue;
                    
                    auto exp_Frefs_Rot_real_aux = exp_Frefs_Rot_real.wptr(iclass, ipsi*exp_nr_over_rot+iover_rot, exp_current_Fsize2);
                    auto exp_Frefs_Rot_imag_aux = exp_Frefs_Rot_imag.wptr(iclass, ipsi*exp_nr_over_rot+iover_rot, exp_current_Fsize2);
                    auto Minvsigma2_aux 		= exp_local_Minvsigma2s[iimage].wptr(exp_current_Fsize2);
                    // write out data stream
                    global_data_stream.foutInt(iclass, "getAllSquaredDifferences()_iclass", __FILE__, __LINE__);
                    global_data_stream.foutInt(iimage, "getAllSquaredDifferences()_iimage", __FILE__, __LINE__);
                    global_data_stream.foutInt(0/*idir*/, "getAllSquaredDifferences()_idir", __FILE__, __LINE__);
                    global_data_stream.foutInt(ipsi, "getAllSquaredDifferences()_ipsi", __FILE__, __LINE__);
                    global_data_stream.foutInt(iover_rot, "getAllSquaredDifferences()_iover_rot", __FILE__, __LINE__);
                    global_data_stream.foutDouble(mlModel.pdf_class.wptrAll()[iclass], "getAllSquaredDifferences()_pdf_class[iclass]", __FILE__, __LINE__);
                    global_data_stream.foutDouble(mlModel.pdf_direction[iclass][0/*idir*/], "getAllSquaredDifferences()_pdf_direction[iclass][idir]", __FILE__, __LINE__);
                    global_data_stream.foutDouble(0/*samplingGrid.exp_over_rot[ipsi][iover_rot]*/, "getAllSquaredDifferences()_exp_over_rot[iorient][iover_rot]", __FILE__, __LINE__);
                    global_data_stream.foutDouble(0/*samplingGrid.exp_over_tilt[ipsi][iover_rot]*/, "getAllSquaredDifferences()_exp_over_tilt[iorient][iover_rot]", __FILE__, __LINE__);
                    global_data_stream.foutDouble(samplingGrid.exp_over_psi[ipsi][iover_rot], "getAllSquaredDifferences()_exp_over_psi[iorient][iover_rot]", __FILE__, __LINE__);
                    global_data_stream.foutInt(0/*do_firstiter_cc*/, "getAllSquaredDifferences()_do_firstiter_cc", __FILE__, __LINE__);
                    global_data_stream.foutInt(0/*do_always_cc*/, "getAllSquaredDifferences()_do_always_cc", __FILE__, __LINE__);
                    global_data_stream.foutInt(0/*do_scale_correction*/, "getAllSquaredDifferences()_do_scale_correction", __FILE__, __LINE__);
                    global_data_stream.foutInt(do_ctf_correction, "getAllSquaredDifferences()_do_ctf_correction", __FILE__, __LINE__);
                    global_data_stream.foutInt(refs_are_ctf_corrected, "getAllSquaredDifferences()_refs_are_ctf_corrected", __FILE__, __LINE__);
                    global_data_stream.foutInt(exp_metadata[iimage].GROUP_NO-1, "getAllSquaredDifferences()_GROUP_NO", __FILE__, __LINE__);
                    global_data_stream.foutDouble(mlModel.scale_correction.wptrAll()[exp_metadata[iimage].GROUP_NO-1], "getAllSquaredDifferences()_scale_correction", __FILE__, __LINE__);
                    global_data_stream.foutDouble(exp_Frefs_Rot_real_aux, exp_current_Fsize2, "getAllSquaredDifferences()_Fref_real", __FILE__, __LINE__);
                    global_data_stream.foutDouble(exp_Frefs_Rot_imag_aux, exp_current_Fsize2, "getAllSquaredDifferences_Fref_imag", __FILE__, __LINE__);
                    global_data_stream.foutDouble(Minvsigma2_aux, exp_current_Fsize2, "getAllSquaredDifferences_Fref_imag()_Minvsigma2", __FILE__, __LINE__);
                    global_data_stream.check();global_data_stream.flush();
                }
            }
        }
    }
#endif
#pragma omp parallel for collapse(4) schedule(dynamic)
    for (int iimage = 0; iimage < exp_nr_images; iimage++)
    {
        for (int iclass = exp_iclass_min; iclass <= exp_iclass_max; iclass++)
        {
            for (int ipsi = 0; ipsi < exp_nr_psi; ipsi++)
            {
                for (int iover_rot = 0; iover_rot < exp_nr_over_rot; iover_rot++)
                {
                    if (!exp_Rot_significant[iclass].rptrAll()[ipsi]
                        || mlModel.pdf_direction[iclass].rptrAll()[0] <= 0.
                        || mlModel.pdf_class.rptrAll()[iclass] <= 0.) continue;
                    
                    // exp_iclass loop does not always go from 0 to nr_classes!
                    int iclass_rot = iclass * exp_nr_psi + ipsi;
                    int iclass_rot_over = iclass_rot * exp_nr_over_rot + iover_rot;
                    int tid = omp_get_thread_num();
                    auto thread_exp_min_diff2_tid = thread_exp_min_diff2[tid].wptr(exp_nr_images);
                    auto Frefctf_real 			= thread_Frefctf_real.wptr(tid, exp_current_Fsize2);
                    auto Frefctf_imag 			= thread_Frefctf_imag.wptr(tid, exp_current_Fsize2);
                    auto Minvsigma2_aux 		= exp_local_Minvsigma2s[iimage].wptr(exp_current_Fsize2);
                    auto exp_Frefs_Rot_real_aux = exp_Frefs_Rot_real.wptr(iclass, ipsi*exp_nr_over_rot+iover_rot, exp_current_Fsize2);
                    auto exp_Frefs_Rot_imag_aux = exp_Frefs_Rot_imag.wptr(iclass, ipsi*exp_nr_over_rot+iover_rot, exp_current_Fsize2);

                    // Apply CTF to reference projection
                    // after first iteration refs_are_ctf_corrected = true
                    if (do_ctf_correction && refs_are_ctf_corrected)
                    {
                        auto exp_local_Fctfs_aux = exp_local_Fctfs[iimage].wptr(exp_current_Fsize2);
                        
                        for (int n = 0; n < exp_current_Fsize2; n++) {
                            Frefctf_real[n] = exp_Frefs_Rot_real_aux[n]*exp_local_Fctfs_aux[n];
                            Frefctf_imag[n] = exp_Frefs_Rot_imag_aux[n]*exp_local_Fctfs_aux[n];
                        }
                    }
                    else{
                        for (int n = 0; n < exp_current_Fsize2; n++) {
                            Frefctf_real[n] = exp_Frefs_Rot_real_aux[n];
                            Frefctf_imag[n] = exp_Frefs_Rot_imag_aux[n];
                        }
                    }
                    
                    for (int itrans = 0; itrans < exp_nr_trans; itrans++)//100 or 400
                    {
                        for (int iover_trans = 0; iover_trans < exp_nr_over_trans; iover_trans++)//1,4
                        {
                            int iclass_rot_trans = iclass_rot * exp_nr_trans + itrans;
                            
                            // In the first pass, always proceed
                            // In the second pass, check whether this translations (&orientation) had a significant weight in the first pass
                            if (exp_ipass != 0 && exp_Mcoarse_significant[iimage].rptrAll()[iclass_rot_trans] == false) {
                                continue;
                            }
                            
                            // Get the shifted image
                            int itrans_over_trans = itrans*exp_nr_over_trans+iover_trans;
                            auto Fimg_shift_real_aux = exp_Fimgs_shifted_real.wptr(iimage, itrans_over_trans, exp_current_Fsize2);
                            auto Fimg_shift_imag_aux = exp_Fimgs_shifted_imag.wptr(iimage, itrans_over_trans, exp_current_Fsize2);
                            
                            double diff2 = 0;
                            
                            for(int n = 0;n < exp_current_Fsize2;n++)
                            {
                                auto diff_real = Frefctf_real[n] - Fimg_shift_real_aux[n];
                                auto diff_imag = Frefctf_imag[n] - Fimg_shift_imag_aux[n];
                                diff2 += (diff_real * diff_real + diff_imag * diff_imag) * Minvsigma2_aux[n];
                            }
                            
                            // Calculate the actual squared difference term of the Gaussian probability function
                            // If current_size < model ori_size diff2 is initialised to the sum of
                            // all |Xij|2 terms that lie between current_size and ori_size
                            // Factor two because of factor 2 in division below, NOT because of 2-dimensionality of the complex plane!
                            diff2 = (diff2 + exp_highres_Xi2_imgs.rptrAll()[iimage]) / 2;
                            // Store all diff2 in exp_Mweight
                            int iclass_rot_trans_over = iclass_rot_over*exp_nr_trans*exp_nr_over_trans + itrans*exp_nr_over_trans+iover_trans;
                            
                            //remember the structure of exp_Mweight class->rot->over_rot->trans->over_trans
                            exp_Mweight[iimage*exp_Mweight_xsize+iclass_rot_trans_over] = diff2;
                            //
                            if (diff2 < thread_exp_min_diff2_tid[iimage]) {
                                thread_exp_min_diff2_tid[iimage] = diff2;
                            }
                        } // end loop iover_trans
                    } // end loop itrans
                }// end loop iover_rot
            } // end loop ipsi
        }//end loop iclass
    } // end loop iimage
    
    //
    exp_min_diff2.fill( (std::numeric_limits<double>::max)() );
    for (int thread = 0; thread < maxthreads; thread++){
        auto thread_exp_min_diff2_tid = thread_exp_min_diff2[thread].wptr(exp_nr_images);
        for (int iimage = 0; iimage < exp_nr_images; iimage++){
            if (thread_exp_min_diff2_tid[iimage] < exp_min_diff2.rptrAll()[iimage]) {
                exp_min_diff2.wptrAll()[iimage] = thread_exp_min_diff2_tid[iimage];
            }
        }
    }
    
    
#ifdef DATA_STREAM
    //
    global_data_stream.turnOn();
    //
    for (int iimage = 0; iimage < exp_nr_images; iimage++)
        global_data_stream.foutDouble(exp_highres_Xi2_imgs.wptrAll()[iimage], "getAllSquaredDifferences()_exp_highres_Xi2_imgs", __FILE__, __LINE__);
    global_data_stream.foutDouble(exp_min_diff2.wptrAll()[0], "getAllSquaredDifferences()_exp_min_diff2[0]", __FILE__, __LINE__);
    global_data_stream.check();global_data_stream.flush();
#endif
}


void convertSquaredDifferencesToWeights()
{
    thread_exp_max_weight.fill_with_first_touch((std::numeric_limits<double>::min)());
    thread_exp_sum_weight.fill_with_first_touch(0.);
#pragma omp parallel for
    for (int thread = 0; thread < maxthreads; thread++)
    {
        auto thread_exp_max_weight_index_tid = thread_exp_max_weight_index[thread].data();
        for (int iimage = 0; iimage < exp_nr_images; iimage++)
            thread_exp_max_weight_index_tid[iimage] = {0,0,0,0,0,0,0};
    }
#ifdef DATA_STREAM
    // turn off data stream inside openmp for loop
    // global_data_stream.turnOff();
#endif
#pragma omp parallel for collapse(6)
    for (int iimage = 0; iimage < exp_nr_images;iimage++)
    {
        for (int iclass = exp_iclass_min; iclass <= exp_iclass_max; iclass++)
        {
            for (int ipsi = 0; ipsi < exp_nr_psi; ipsi++)
            {
                for (int itrans = 0; itrans < exp_nr_trans; itrans++)
                {
                	for (int iover_rot = 0; iover_rot < exp_nr_over_rot; iover_rot++)
                	{
                        for (int iover_trans = 0; iover_trans < exp_nr_over_trans; iover_trans++)
                        {
                            // exp_Mweight order iclass->irot->iover_rot->itrans->iover_trans
                            int iclass_rot = iclass * exp_nr_psi + ipsi;
                            int iclass_rot_over = iclass_rot * exp_nr_over_rot + iover_rot;
                            int iclass_rot_trans_over = iclass_rot_over*exp_nr_trans*exp_nr_over_trans + itrans*exp_nr_over_trans+iover_trans;
                            
                            // Only exponentiate for determined values of exp_Mweight
                            // (this is always true in the first pass, but not so in the second pass)
                            // Only deal with this sampling point if its weight was significant
                            if (exp_Mweight[iimage*exp_Mweight_xsize+iclass_rot_trans_over] < 0.)
                            {
                                exp_Mweight[iimage*exp_Mweight_xsize+iclass_rot_trans_over] = 0.;
                            }
                            else
                            {
                                FDOUBLE offset_x,offset_y;
                                sampler2d.getTranslation(itrans, offset_x, offset_y);
                                offset_x += exp_old_offsetx.rptrAll()[iimage];
                                offset_y += exp_old_offsety.rptrAll()[iimage];
                                // Convert offsets back to Angstroms to calculate PDF!
                                double pdf_offset = mlModel.calculatePdfOffset(offset_x,offset_y,mlModel.prior_offsetx_class.rptrAll()[iclass],
                                                                               mlModel.prior_offsety_class.rptrAll()[iclass]);
                                
                                auto pdf_orientation = mlModel.pdf_direction[iclass][0];
                                
                                auto weight = pdf_orientation * pdf_offset;
                                
                                auto diff2 = exp_Mweight[iimage*exp_Mweight_xsize+iclass_rot_trans_over] - exp_min_diff2.rptrAll()[iimage];
                                // next line because of numerical precision of exp-function
                                if (diff2 > 700.) weight = 0.;
                                // TODO: use tabulated exp function?
                                else weight *= exp(-diff2);
                                
                                // Store the weight
                                exp_Mweight[iimage*exp_Mweight_xsize+iclass_rot_trans_over] = weight;
                                //
                                int tid = omp_get_thread_num();
                                auto thread_exp_sum_weight_tid = thread_exp_sum_weight[tid].wptr(exp_nr_images);
                                auto thread_exp_max_weight_tid = thread_exp_max_weight[tid].wptr(exp_nr_images);
                                auto thread_exp_max_weight_index_tid = thread_exp_max_weight_index[tid].data();
                                thread_exp_sum_weight_tid[iimage] += weight;
                                if (weight > thread_exp_max_weight_tid[iimage])
                                {
                                    thread_exp_max_weight_tid[iimage] = weight;
                                    thread_exp_max_weight_index_tid[iimage] = {iimage,iclass,0/*idir*/,ipsi,iover_rot,itrans,iover_trans};
                                }
#ifdef DATA_STREAM
                                global_data_stream.foutInt(iimage, "convertSquaredDifferencesToWeights()_iimage", __FILE__, __LINE__);
                                global_data_stream.foutInt(iclass, "convertSquaredDifferencesToWeights()_iclass", __FILE__, __LINE__);
                                global_data_stream.foutInt(0/*idir*/, "convertSquaredDifferencesToWeights()_idir", __FILE__, __LINE__);
                                global_data_stream.foutInt(ipsi, "convertSquaredDifferencesToWeights()_ipsi", __FILE__, __LINE__);
                                global_data_stream.foutInt(itrans, "convertSquaredDifferencesToWeights()_itrans", __FILE__, __LINE__);
                                global_data_stream.foutInt(iover_rot, "convertSquaredDifferencesToWeights()_iover_rot", __FILE__, __LINE__);
                                global_data_stream.foutInt(iover_trans, "convertSquaredDifferencesToWeights()_iover_trans", __FILE__, __LINE__);
                                global_data_stream.foutDouble(offset_x, "convertSquaredDifferencesToWeights()_offset_x", __FILE__, __LINE__);
                                global_data_stream.foutDouble(offset_y, "convertSquaredDifferencesToWeights()_offset_y", __FILE__, __LINE__);
                                global_data_stream.foutDouble(pdf_offset, "convertSquaredDifferencesToWeights()_pdf_offset", __FILE__, __LINE__);
                                global_data_stream.foutDouble(pdf_orientation, "convertSquaredDifferencesToWeights()_pdf_orientation", __FILE__, __LINE__);
                                global_data_stream.foutDouble(diff2, "convertSquaredDifferencesToWeights()_diff2", __FILE__, __LINE__);
                                global_data_stream.foutDouble(weight, "convertSquaredDifferencesToWeights()_weight", __FILE__, __LINE__);
                                global_data_stream.foutDouble(mlModel.prior_offsetx_class.wptrAll()[iclass], "convertSquaredDifferencesToWeights()_prior_offsetx_class[iclass]", __FILE__, __LINE__);
                                global_data_stream.foutDouble(mlModel.prior_offsety_class.wptrAll()[iclass], "convertSquaredDifferencesToWeights()_prior_offsety_class[iclass]", __FILE__, __LINE__);
                                global_data_stream.foutDouble(mlModel.sigma2_offset, "convertSquaredDifferencesToWeights()_sigma2_offset", __FILE__, __LINE__);
                                global_data_stream.foutDouble(exp_old_offsetx.wptrAll()[iimage], "convertSquaredDifferencesToWeights()_exp_old_offsetx[iimage]", __FILE__, __LINE__);
                                global_data_stream.foutDouble(exp_old_offsety.wptrAll()[iimage], "convertSquaredDifferencesToWeights()_exp_old_offsety[iimage]", __FILE__, __LINE__);
                                global_data_stream.check();global_data_stream.flush();
#endif
                            } // end if/else exp_Mweight < 0.
                        } // end loop iover_trans
                    }// end loop iover_rot
                } // end loop itrans
            } // end loop iorient
        } //end loop iclass
    }//end loop iimage
#ifdef DATA_STREAM
    global_data_stream.turnOn();
#endif
    exp_sum_weight.fill(0.);
#pragma omp parallel for
    for (int iimage = 0; iimage < exp_nr_images; iimage++)
    {
        double exp_max_weight = 0;
        GridIndex exp_max_weight_index;
        
        for (int thread = 0; thread < maxthreads; thread++)
        {
            const auto thread_exp_sum_weight_tid = thread_exp_sum_weight[thread].wptr(exp_nr_images);
            exp_sum_weight.wptrAll()[iimage] += thread_exp_sum_weight_tid[iimage];
            const auto thread_exp_max_weight_tid = thread_exp_max_weight[thread].wptr(exp_nr_images);
            const auto thread_exp_max_weight_index_tid = thread_exp_max_weight_index[thread].data();
            if(thread_exp_max_weight_tid[iimage] > exp_max_weight){
                exp_max_weight = thread_exp_max_weight_tid[iimage];
                exp_max_weight_index = thread_exp_max_weight_index_tid[iimage];
            }
        }
        //
        int iimage = exp_max_weight_index.iimage;
        int iclass = exp_max_weight_index.iclass;
        int idir = exp_max_weight_index.idir;
        int ipsi = exp_max_weight_index.ipsi;
        int iover_rot = exp_max_weight_index.iover_rot;
        int itrans = exp_max_weight_index.itrans;
        int iover_trans = exp_max_weight_index.iover_trans;
        int iorient = idir * exp_nr_psi + ipsi;
        //
        exp_metadata[iimage].PMAX = exp_max_weight/exp_sum_weight.rptrAll()[iimage];
        exp_metadata[iimage].CLASS = iclass+1;
        exp_metadata[iimage].ROT = 0;//samplingGrid.exp_over_rot[iorient][iover_rot];
        exp_metadata[iimage].TILT = 0;//samplingGrid.exp_over_tilt[iorient][iover_rot];
        exp_metadata[iimage].PSI = samplingGrid.exp_over_psi[iorient][iover_rot];
        if (exp_metadata[iimage].PSI > 180.) exp_metadata[iimage].PSI -= 360.;
        if (exp_metadata[iimage].PSI < -180.) exp_metadata[iimage].PSI += 360.;
        exp_metadata[iimage].XOFF = exp_old_offsetx.rptrAll()[iimage] + samplingGrid.exp_over_trans_x[itrans][iover_trans];
        exp_metadata[iimage].YOFF = exp_old_offsety.rptrAll()[iimage] + samplingGrid.exp_over_trans_y[itrans][iover_trans];
    }
    
#ifdef DATA_STREAM
    global_data_stream.foutDouble(exp_sum_weight.wptrAll()[0], "convertSquaredDifferencesToWeights()_exp_sum_weight[0]", __FILE__, __LINE__);
    global_data_stream.check();global_data_stream.flush();
#endif
}

void findAllSignificantPoints()
{
#ifdef DATA_STREAM
    // turn off data stream inside openmp for loop
    // global_data_stream.turnOff();
#endif
    // Now, for each image,find the exp_significant_weight that encompasses adaptive_fraction of exp_sum_weight
#pragma omp parallel for
    for (int iimage = 0;iimage < exp_nr_images; iimage++)
    {
        double* sorted_weight = (double*)aMalloc(sizeof(double)*exp_Mweight_xsize,64);
        int sorted_weight_count = 0;
        
        double frac_weight = 0.;
        double my_significant_weight;
        int my_nr_significant_coarse_samples = 0;
        
        
        for (int ihidden = 0; ihidden < exp_Mweight_xsize; ihidden++)
        {
            if (exp_Mweight[iimage*exp_Mweight_xsize+ihidden] > 0.)
            {
                sorted_weight[sorted_weight_count] = exp_Mweight[iimage*exp_Mweight_xsize+ihidden];
                sorted_weight_count++;
            }
        }
        
        // sort from small to large
        std::sort(sorted_weight, sorted_weight + sorted_weight_count);
        
        for (int i = sorted_weight_count - 1; i >= 0; i--)
        {
            if (exp_ipass==0) my_nr_significant_coarse_samples++;
            my_significant_weight = sorted_weight[i];
            frac_weight += my_significant_weight;
            if (frac_weight > adaptive_fraction * exp_sum_weight.rptrAll()[iimage])
                break;
        }
        
        if (exp_ipass==0 && my_nr_significant_coarse_samples == 0)
        {
            std::cerr << " ipart= " << iimage << " adaptive_fraction= " << adaptive_fraction << std::endl;
            std::cerr << " frac-weight= " << frac_weight << std::endl;
            std::cerr << " exp_sum_weight[iimage] = " << exp_sum_weight.rptrAll()[iimage] << std::endl;
            std::cerr << " sorted_weight_count = " << sorted_weight_count << std::endl;
            if (sorted_weight_count > 0)
            {
                std::cout<<"sum of abs(sorted_weight) = "<<sumvec(sorted_weight, sorted_weight_count)<<std::endl;
                std::cerr << "written sorted_weight.spi" << std::endl;
            }
            ERROR_REPORT("my_nr_significant_coarse_samples == 0");
        }
        
        if (exp_ipass==0)
        {
            exp_metadata[iimage].NR_SIGN = (double)my_nr_significant_coarse_samples;
            assert(exp_Mweight_xsize==exp_Mcoarse_significant[0].size());
            // Keep track of which coarse samplings were significant were significant for this particle
            // in ipass = 0,exp_Mcoarse_xsize eqaul to exp_Mweight_xsize
            for (int ihidden = 0; ihidden < exp_Mweight_xsize; ihidden++)
            {
                if (exp_Mweight[iimage*exp_Mweight_xsize+ihidden] >= my_significant_weight)
                    exp_Mcoarse_significant[iimage].wptrAll()[ihidden] = true;
                else
                    exp_Mcoarse_significant[iimage].wptrAll()[ihidden] = false;
            }
            
        }
        
        exp_significant_weight.wptrAll()[iimage] = my_significant_weight;
        
#ifdef DATA_STREAM
        global_data_stream.foutDouble(sorted_weight_count, "findAllSignificantPoints()_np", __FILE__, __LINE__);
        global_data_stream.foutDouble(&sorted_weight[0], sorted_weight_count, "findAllSignificantPoints()_sorted_weight", __FILE__, __LINE__);
        global_data_stream.foutDouble(my_nr_significant_coarse_samples, "findAllSignificantPoints()_my_nr_significant_coarse_samples", __FILE__, __LINE__);
        global_data_stream.foutDouble(my_significant_weight, "findAllSignificantPoints()_my_significant_weight", __FILE__, __LINE__);
        global_data_stream.foutDouble(frac_weight, "findAllSignificantPoints()_frac_weight", __FILE__, __LINE__);
        global_data_stream.foutDouble(adaptive_fraction, "findAllSignificantPoints()_adaptive_fraction", __FILE__, __LINE__);
        global_data_stream.check();global_data_stream.flush();
#endif
        aFree(sorted_weight);
    } // end loop iimage
#ifdef DATA_STREAM
    global_data_stream.turnOn();
#endif
}
    

void storeWeightedSums()
{
    int current_Fsize = current_size/2+1;
    int current_Fsize2 = current_size*(current_size/2+1);
    int ori_Fsize = ori_size/2+1;
#ifdef DATA_STREAM
    // turn off data stream inside openmp for loop
    global_data_stream.turnOff();
#endif
    // Extend norm_correction and sigma2_noise estimation to higher resolutions
    // NOTE : some memory conflict for group-data,like wsum_sigma2_noise,wsum_signal_product_spectra...
    // so removing openmp
    for (int iimage = 0; iimage < exp_nr_images; iimage++)
    {
        int igroup = exp_metadata[iimage].GROUP_NO-1;
        auto exp_power_imgs_iimage = exp_power_imgs[iimage].wptr(ori_Fsize);
        auto wsum_sigma2_noise_igroup = mlModel.wsum_sigma2_noise[igroup].wptr(mlModel.ori_Fsize);
        // If the current images were smaller than the original size, fill the rest of wsum_model.sigma2_noise with the power_class spectrum of the images
        for (int ires = current_Fsize; ires < ori_Fsize; ires++)
        {
            wsum_sigma2_noise_igroup[ires] += exp_power_imgs_iimage[ires];
            // Also extend the weighted sum of the norm_correction
            exp_wsum_norm_correction.wptrAll()[iimage] += exp_power_imgs_iimage[ires];
        }
        
        // Store norm_correction
        // Multiply by old value because the old norm_correction term was already applied to the image
        if (do_norm_correction)
        {
            auto old_norm_correction = exp_metadata[iimage].NORM / mlModel.avg_norm_correction;
            // Now set the new norm_correction in the relevant position of exp_metadata22
            // The factor two below is because exp_wsum_norm_correctiom is similar to sigma2_noise, which is the variance for the real/imag components
            // The variance of the total image (on which one normalizes) is twice this value!
            exp_metadata[iimage].NORM = old_norm_correction * sqrt(exp_wsum_norm_correction.rptrAll()[iimage] * 2.);
            if (exp_metadata[iimage].NORM > 10.)
            {
#pragma omp critical
                std::cerr<<"Warning in storeWeightedSums(),please debug this function."<<std::endl;
            }
#pragma omp atomic
            mlModel.wsum_avg_norm_correction += old_norm_correction * sqrt(exp_wsum_norm_correction.rptrAll()[iimage] * 2.);
        }
        
        // Some analytics...
        // Calculate normalization constant for dLL
        // loop over all particles inside this ori_particle
        double logsigma2 = 0.;
        auto sigma2_noise_igroup = mlModel.sigma2_noise[igroup].wptr(mlModel.ori_Fsize);
        for (int n = 0; n < current_Fsize2; n++)
        {
            int ires = Mresol_fine.rptrAll()[n];
            // Note there is no sqrt in the normalisation term because of the 2-dimensionality of the complex-plane
            // Also exclude origin from logsigma2, as this will not be considered in the P-calculations
            if (ires > 0)
                logsigma2 += log( 2. * PI * sigma2_noise_igroup[ires]);
        }
        
        
        if (exp_sum_weight.rptrAll()[iimage]==0)
        {
            ERROR_REPORT("ERROR: exp_sum_weight[ipart]==0");
        }
        
        double dLL = log(exp_sum_weight.rptrAll()[iimage]) - exp_min_diff2.rptrAll()[iimage] - logsigma2;
        
        // Also store dLL of each image in the output array
        exp_metadata[iimage].DLL = dLL;
        
        mlModel.wsum_LL += dLL;

        mlModel.wsum_ave_Pmax += exp_metadata[iimage].PMAX;
        
#ifdef DATA_STREAM
        global_data_stream.foutInt(iimage, "storeWeightedSums()_iimage", __FILE__, __LINE__);
        global_data_stream.foutDouble(exp_wsum_norm_correction.wptrAll()[iimage], "storeWeightedSums()_exp_wsum_norm_correction", __FILE__, __LINE__);
        global_data_stream.foutDouble(exp_metadata[iimage].NORM, "storeWeightedSums()_exp_metadata[iimage].NORM", __FILE__, __LINE__);
        global_data_stream.foutDouble(exp_metadata[iimage].DLL, "storeWeightedSums()_exp_metadata[iimage].DLL", __FILE__, __LINE__);
        global_data_stream.check();global_data_stream.flush();
#endif
    } // end loop iimage
#ifdef DATA_STREAM
    global_data_stream.turnOn();
    global_data_stream.foutDouble(mlModel.wsum_LL, "storeWeightedSums()_wsum_LL", __FILE__, __LINE__);
    global_data_stream.foutDouble(mlModel.wsum_ave_Pmax, "storeWeightedSums()_wsum_ave_Pmax", __FILE__, __LINE__);
    global_data_stream.foutDouble(mlModel.wsum_avg_norm_correction, "storeWeightedSums()_wsum_avg_norm_correction", __FILE__, __LINE__);
    global_data_stream.foutDouble(mlModel.wsum_signal_product_spectra[0].wptr(ori_Fsize), ori_Fsize, "storeWeightedSums()_wsum_signal_product_spectra1", __FILE__, __LINE__);
    global_data_stream.foutDouble(mlModel.wsum_signal_product_spectra[mlModel.nr_groups-1].wptr(ori_Fsize), ori_Fsize, "toreWeightedSums()_wsum_signal_product_spectraN", __FILE__, __LINE__);
    global_data_stream.foutDouble(mlModel.wsum_reference_power_spectra[0].wptr(ori_Fsize),ori_Fsize , "storeWeightedSums()_wsum_reference_power_spectra1", __FILE__, __LINE__);
    global_data_stream.foutDouble(mlModel.wsum_reference_power_spectra[mlModel.nr_groups-1].wptr(ori_Fsize), ori_Fsize, "storeWeightedSums()_wsum_reference_power_spectraN", __FILE__, __LINE__);
    global_data_stream.foutDouble(mlModel.wsum_sigma2_noise[0].wptr(ori_Fsize), ori_Fsize, "storeWeightedSums()_wsum_sigma2_noise1", __FILE__, __LINE__);
    global_data_stream.foutDouble(mlModel.wsum_sigma2_noise[mlModel.nr_groups-1].wptr(ori_Fsize), ori_Fsize, "storeWeightedSums()_wsum_sigma2_noiseN", __FILE__, __LINE__);
    global_data_stream.foutDouble(mlModel.wsum_sumw_group.wptrAll()[0], "storeWeightedSums()_wsum_sumw_group1", __FILE__, __LINE__);
    global_data_stream.foutDouble(mlModel.wsum_sumw_group.wptrAll()[mlModel.nr_groups-1], "storeWeightedSums()_wsum_sumw_groupN", __FILE__, __LINE__);
    global_data_stream.foutDouble(mlModel.wsum_sigma2_offset, "storeWeightedSums()_wsum_sigma2_offset", __FILE__, __LINE__);
    global_data_stream.foutDouble(mlModel.wsum_pdf_class.wptrAll()[0], "storeWeightedSums()_wsum_pdf_class0", __FILE__, __LINE__);
    global_data_stream.foutDouble(mlModel.wsum_pdf_class.wptrAll()[nr_classes-1], "storeWeightedSums()_wsum_pdf_classN", __FILE__, __LINE__);
    global_data_stream.foutDouble(mlModel.wsum_pdf_direction[0].wptr(exp_nr_dir), exp_nr_dir, "storeWeightedSums()_wsum_pdf_direction0", __FILE__, __LINE__);
    global_data_stream.foutDouble(mlModel.wsum_pdf_direction[nr_classes-1].wptr(exp_nr_dir), exp_nr_dir, "storeWeightedSums()_wsum_pdf_directionN", __FILE__, __LINE__);
    global_data_stream.check();global_data_stream.flush();
#endif
}


void updateOtherParams()
{
    int exp_current_Fsize2 = exp_current_size*(exp_current_size/2+1);
    
    thread_wsum_norm_correction	.fill_with_first_touch(0.);
    thread_wsum_sigma2_noise	.fill_with_first_touch(0.);
    
#pragma omp parallel for collapse(4) schedule(dynamic)
    for (int iimage = 0; iimage < exp_nr_images; iimage++)
    {
        for (int iclass = exp_iclass_min; iclass <= exp_iclass_max; iclass++)
        {
            for (int ipsi = 0; ipsi < exp_nr_psi; ipsi++)
            {
                for (int iover_rot = 0; iover_rot < exp_nr_over_rot; iover_rot++)
                {
                    if (!exp_Rot_significant[iclass].rptrAll()[ipsi])
                        continue;
                    
                    int iclass_rot = iclass * exp_nr_psi + ipsi;
                    int iclass_rot_over = iclass_rot * exp_nr_over_rot + iover_rot;
                    
					int tid = omp_get_thread_num();
                    auto Frefctf_real				=	thread_Frefctf_real.wptr(tid, exp_current_Fsize2);
                    auto Frefctf_imag 				=	thread_Frefctf_imag.wptr(tid, exp_current_Fsize2);
                    auto wsum_sigma2_noise_tid 		=	thread_wsum_sigma2_noise.wptr(tid, iimage, exp_current_Fsize2);
                    auto wsum_norm_correction_tid 	=	thread_wsum_norm_correction[tid].wptr(exp_nr_images);
                    auto Frefs_Rot_real		 		=	exp_Frefs_Rot_real.wptr(iclass, ipsi*exp_nr_over_rot+iover_rot, exp_current_Fsize2);
                    auto Frefs_Rot_imag		 		=	exp_Frefs_Rot_imag.wptr(iclass, ipsi*exp_nr_over_rot+iover_rot, exp_current_Fsize2);
                    
                    // Apply CTF to reference
                    if (do_ctf_correction && refs_are_ctf_corrected)
                    {
                        auto exp_local_Fctfs_aux = exp_local_Fctfs[iimage].wptr(exp_current_Fsize2);
                        #pragma ivdep
                        for (int n = 0; n < exp_current_Fsize2; n++)
                        {
                            Frefctf_real[n] = Frefs_Rot_real[n] * exp_local_Fctfs_aux[n];
                            Frefctf_imag[n] = Frefs_Rot_imag[n] * exp_local_Fctfs_aux[n];
                        }
                    }
                    else
                    {
                        #pragma ivdep
                        for (int n = 0; n < exp_current_Fsize2; n++)
                        {
                            Frefctf_real[n] = Frefs_Rot_real[n];
                            Frefctf_imag[n] = Frefs_Rot_imag[n];
                        }
                    }
                    
                    
                    for (int itrans = 0; itrans < exp_nr_trans; itrans++)//100 or 400
                    {
                        for (int iover_trans = 0; iover_trans < exp_nr_over_trans; iover_trans++)//4
                        {
                            
                            int iclass_rot_trans_over = iclass_rot_over*exp_nr_trans*exp_nr_over_trans + itrans*exp_nr_over_trans+iover_trans;
                            
                            // Only deal with this sampling point if its weight was significant
                            double weight = exp_Mweight[iimage*exp_Mweight_xsize+iclass_rot_trans_over];
                            // Only sum weights for non-zero weights
                            if (weight >= exp_significant_weight.rptrAll()[iimage])
                            {
                                // Normalise the weight (do this after the comparison with exp_significant_weight!)
                                weight /= exp_sum_weight.rptrAll()[iimage];
                                
                                // Get the shifted image
                                int itrans_over_trans = itrans*exp_nr_over_trans+iover_trans;
                                auto Fimg_shift_real = exp_Fimgs_shifted_real.wptr(iimage, itrans_over_trans, exp_current_Fsize2);
                                auto Fimg_shift_imag = exp_Fimgs_shifted_imag.wptr(iimage, itrans_over_trans, exp_current_Fsize2);

                                FDOUBLE sum_wdiff2 = 0.;
                                // Store weighted sum of squared differences for sigma2_noise estimation
                                #pragma ivdep
                                for (int n = 0; n < exp_current_Fsize2; n++)
                                {
                                    int ires = Mresol_fine.rptrAll()[n];
                                    if (ires > -1)
                                    {
                                        // Use FT of masked image for noise estimation!
                                        auto diff_real = Frefctf_real[n] - Fimg_shift_real[n];
                                        auto diff_imag = Frefctf_imag[n] - Fimg_shift_imag[n];
                                        auto wdiff2 = weight * (diff_real*diff_real + diff_imag*diff_imag);
                                        // group-wise sigma2_noise
                                        wsum_sigma2_noise_tid[n] += wdiff2;
                                        // For norm_correction
                                        sum_wdiff2 += wdiff2;
                                    }
                                }
                                
                                wsum_norm_correction_tid[iimage] += sum_wdiff2;
                                
                            } // end if weight >= exp_significant_weight
                        } // end loop iover_trans
                    } // end loop itrans
                }// end if iover_rot
            } // end loop iorient
        }// end of iclass
    } // end loop iimage
    
    exp_wsum_norm_correction	.fill(0.);

    for (int thread = 0; thread < maxthreads; thread++)
    {
        auto thread_wsum_norm_correction_tid = thread_wsum_norm_correction[thread].wptr(exp_nr_images);
        for (int iimage = 0; iimage < exp_nr_images; iimage++)
            exp_wsum_norm_correction.wptrAll()[iimage] += thread_wsum_norm_correction_tid[iimage];
        for (int iimage = 0; iimage < exp_nr_images; iimage++)
        {
            auto thread_wsum_sigma2_noise_iimage = thread_wsum_sigma2_noise.wptr(thread, iimage, exp_current_Fsize2);
            int igroup = exp_metadata[iimage].GROUP_NO-1;
            auto wsum_sigma2_noise_igroup = mlModel.wsum_sigma2_noise[igroup].wptr(mlModel.ori_Fsize);
            for (int n = 0; n < exp_current_Fsize2; n++)
            {
                int ires = Mresol_fine.rptrAll()[n];
                if (ires > -1)
                    wsum_sigma2_noise_igroup[ires] += thread_wsum_sigma2_noise_iimage[n];
            }
        }
    }
    
#ifdef DATA_STREAM
    global_data_stream.foutDouble(exp_wsum_norm_correction.wptrAll()[0], "updateOtherParams()_exp_wsum_norm_correction", __FILE__, __LINE__);
    // global_data_stream.foutDouble(mlModel.wsum_sigma2_noise, mlModel.ori_Fsize, "updateOtherParams()_wsum_sigma2_noise1", __FILE__, __LINE__);
    // global_data_stream.foutDouble(mlModel.wsum_sigma2_noise+(mlModel.nr_groups-1)*mlModel.ori_Fsize, mlModel.ori_Fsize, "updateOtherParams()_wsum_sigma2_noiseN", __FILE__, __LINE__);
    // todo,print different group
    if (/*do_scale_correction*/false) {
        // global_data_stream.foutDouble(exp_wsum_scale_correction_XA, mlModel.ori_Fsize, "updateOtherParams()_exp_wsum_scale_correction_XA", __FILE__, __LINE__);
        // global_data_stream.foutDouble(exp_wsum_scale_correction_AA, mlModel.ori_Fsize, "updateOtherParams()_exp_wsum_scale_correction_AA", __FILE__, __LINE__);
    }
    global_data_stream.check();global_data_stream.flush();
#endif
    
    thread_wsum_sigma2_offset		.fill_with_first_touch(0.);
    thread_sumw_group				.fill_with_first_touch(0.);
    thread_wsum_pdf_direction		.fill_with_first_touch(0.);
    thread_wsum_pdf_class			.fill_with_first_touch(0.);
    thread_wsum_prior_offsetx_class	.fill_with_first_touch(0.);
    thread_wsum_prior_offsety_class	.fill_with_first_touch(0.);
    
#pragma omp parallel for collapse(4) schedule(dynamic)
    for (int iimage = 0; iimage < exp_nr_images; iimage++)
    {
        for (int iclass = exp_iclass_min; iclass <= exp_iclass_max; iclass++)
        {
            for (int ipsi = 0; ipsi < exp_nr_psi; ipsi++)
            {
                for (int iover_rot = 0; iover_rot < exp_nr_over_rot; iover_rot++)
                {
                    for (int itrans = 0; itrans < exp_nr_trans; itrans++)
                    {
                        for (int iover_trans = 0; iover_trans < exp_nr_over_trans; iover_trans++)
                        {
                            if (!exp_Rot_significant[iclass].rptrAll()[ipsi])
                                continue;
                            
                            int iclass_rot = iclass * exp_nr_psi + ipsi;
                            int iclass_rot_over = iclass_rot * exp_nr_over_rot + iover_rot;
                            int iclass_rot_trans_over = iclass_rot_over*exp_nr_trans*exp_nr_over_trans + itrans*exp_nr_over_trans+iover_trans;
                            
                            int tid = omp_get_thread_num();
                            auto wsum_pdf_direction_tid 		=	thread_wsum_pdf_direction.wptr(tid, nr_classes);
                            auto wsum_pdf_class_tid 			=	thread_wsum_pdf_class[tid].wptr(nr_classes);
                            auto wsum_prior_offsetx_class_tid 	=	thread_wsum_prior_offsetx_class[tid].wptr(nr_classes);
                            auto wsum_prior_offsety_class_tid 	=	thread_wsum_prior_offsety_class[tid].wptr(nr_classes);
                            auto sumw_group_tid 				= 	thread_sumw_group[tid].wptr(exp_nr_images);
                            auto wsum_sigma2_offset_tid			=	thread_wsum_sigma2_offset[tid].wptr(1);
                            // Only deal with this sampling point if its weight was significant
                            double weight = exp_Mweight[iimage*exp_Mweight_xsize+iclass_rot_trans_over];
                            // Only sum weights for non-zero weights
                            if (weight >= exp_significant_weight.rptrAll()[iimage])
                            {
                                
                                // Normalise the weight (do this after the comparison with exp_significant_weight!)
                                weight /= exp_sum_weight.rptrAll()[iimage];
                                
                                // Store sum of weights for this group
                                sumw_group_tid[iimage] += weight;
                                
                                // Store weights for this class and orientation
                                wsum_pdf_class_tid[iclass] += weight;
                                
                                wsum_pdf_direction_tid[iclass] += weight;
                                
                                wsum_prior_offsetx_class_tid[iclass] += weight * (exp_old_offsetx.rptrAll()[iimage] + samplingGrid.exp_over_trans_x[itrans][iover_trans]);
                                wsum_prior_offsety_class_tid[iclass] += weight * (exp_old_offsety.rptrAll()[iimage] + samplingGrid.exp_over_trans_y[itrans][iover_trans]);
                                
                                auto offsetx = mlModel.prior_offsetx_class.rptrAll()[iclass] - exp_old_offsetx.rptrAll()[iimage] - samplingGrid.exp_over_trans_x[itrans][iover_trans];
                                auto offsety = mlModel.prior_offsety_class.rptrAll()[iclass] - exp_old_offsety.rptrAll()[iimage] - samplingGrid.exp_over_trans_y[itrans][iover_trans];
                                
                                //this may cause some false share!
                                wsum_sigma2_offset_tid[0] += weight * (offsetx*offsetx+offsety*offsety);
                                
                            } // end if weight >= exp_significant_weight
                        } // end loop iover_trans
                    } // end loop itrans
                }// end if iover_rot
            } // end loop iorient
        }// end of iclass
    } // end loop iimage
    //
    for (int thread = 0; thread < maxthreads; thread++)
    {
        auto wsum_pdf_class_tid = thread_wsum_pdf_class[thread].wptr(nr_classes);
        auto wsum_pdf_direction_tid = thread_wsum_pdf_direction.wptr(thread, nr_classes);
        auto wsum_prior_offsetx_class_tid = thread_wsum_prior_offsetx_class[thread].wptr(nr_classes);
        auto wsum_prior_offsety_class_tid = thread_wsum_prior_offsety_class[thread].wptr(nr_classes);
        for (int iclass = 0; iclass < nr_classes; iclass++)
        {
            mlModel.wsum_pdf_class.wptrAll()[iclass] += wsum_pdf_class_tid[iclass];
            mlModel.wsum_pdf_direction[iclass].wptrAll()[0] += wsum_pdf_direction_tid[iclass];
            mlModel.wsum_prior_offsetx_class.wptrAll()[iclass] += wsum_prior_offsetx_class_tid[iclass];
            mlModel.wsum_prior_offsety_class.wptrAll()[iclass] += wsum_prior_offsety_class_tid[iclass];
        }
        auto wsum_sigma2_offset_tid = thread_wsum_sigma2_offset[thread].wptr(1);
        mlModel.wsum_sigma2_offset += wsum_sigma2_offset_tid[0];
        
        auto thread_sumw_group_tid = thread_sumw_group[thread].wptr(exp_nr_images);
        for (int iimage = 0; iimage < exp_nr_images; iimage++) {
            int igroup = exp_metadata[iimage].GROUP_NO-1;
            mlModel.wsum_sumw_group.wptrAll()[igroup] += thread_sumw_group_tid[iimage];
        }
    }
    
#ifdef DATA_STREAM
    // global_data_stream.foutDouble(mlModel.wsum_sumw_group[0], "updateOtherParams()_wsum_sumw_group1", __FILE__, __LINE__);
    // global_data_stream.foutDouble(mlModel.wsum_sumw_group[mlModel.nr_groups-1], "updateOtherParams()_wsum_sumw_groupN", __FILE__, __LINE__);
    // global_data_stream.foutDouble(mlModel.wsum_sigma2_offset, "updateOtherParams()_wsum_sigma2_offset", __FILE__, __LINE__);
    // global_data_stream.foutDouble(mlModel.wsum_pdf_class[0], "updateOtherParams()_wsum_pdf_class0", __FILE__, __LINE__);
    // global_data_stream.foutDouble(mlModel.wsum_pdf_class[nr_classes-1], "updateOtherParams()_wsum_pdf_classN", __FILE__, __LINE__);
    // global_data_stream.foutDouble(mlModel.wsum_pdf_direction, exp_nr_dir, "updateOtherParams()_wsum_pdf_direction0", __FILE__, __LINE__);
    // global_data_stream.foutDouble(mlModel.wsum_pdf_direction+(nr_classes-1)*mlModel.nr_directions, exp_nr_dir, "updateOtherParams()_wsum_pdf_directionN", __FILE__, __LINE__);
    // global_data_stream.check();global_data_stream.flush();
#endif
    /* do in convertSquaredDifferencesToWeights() */
#ifdef DATA_STREAM
    global_data_stream.foutDouble(exp_metadata[0].ROT, "exp_metadata[0].ROT", __FILE__, __LINE__);
    global_data_stream.foutDouble(exp_metadata[0].TILT, "exp_metadata[0].TILT", __FILE__, __LINE__);
    global_data_stream.foutDouble(exp_metadata[0].PSI, "exp_metadata[0].PSI", __FILE__, __LINE__);
    global_data_stream.foutDouble(exp_metadata[0].XOFF, "exp_metadata[0].XOFF", __FILE__, __LINE__);
    global_data_stream.foutDouble(exp_metadata[0].YOFF, "exp_metadata[0].YOFF", __FILE__, __LINE__);
    global_data_stream.foutDouble(exp_metadata[0].CLASS, "exp_metadata[0].CLASS", __FILE__, __LINE__);
    global_data_stream.foutDouble(exp_metadata[0].PMAX, "exp_metadata[0].PMAX", __FILE__, __LINE__);
    global_data_stream.check();global_data_stream.flush();
#endif
}


void backProjection()
{
    int exp_current_Fsize2 = exp_current_size*(exp_current_size/2+1);
    
    // In doThreadPrecalculateShiftedImagesCtfsAndInvSigma2s() the origin of the exp_local_Minvsigma2s was omitted.
    // Set those back here
    for (int iimage = 0; iimage < exp_nr_images; iimage++)
    {
        int igroup = exp_metadata[iimage].GROUP_NO-1;
        auto exp_local_Minvsigma2s_iimage = exp_local_Minvsigma2s[iimage].wptr(exp_current_Fsize2);
        exp_local_Minvsigma2s_iimage[0] = 1. / (sigma2_fudge * mlModel.sigma2_noise[igroup][0]);
    }
    
#ifdef DATA_STREAM
    global_data_stream.foutDouble(exp_local_Minvsigma2s[0].wptr(exp_current_Fsize2),exp_nr_images*exp_current_Fsize2, "backProjection()_exp_local_Minvsigma2s", __FILE__, __LINE__);
    global_data_stream.check();global_data_stream.flush();
#endif
    
#pragma omp parallel for
    for (int iimage = 0; iimage < exp_nr_images; iimage++) {
        if (do_map == false) {
            auto exp_local_Minvsigma2s_iimage = exp_local_Minvsigma2s[iimage].wptr(exp_current_Fsize2);
            for (int n = 0; n < exp_current_Fsize2; n++) {
                exp_local_Minvsigma2s_iimage[n] = 1.0;
            }
        }
        //Apply CTF to reference
        if (do_ctf_correction == false) {
            auto exp_local_Fctfs_iimage = exp_local_Fctfs[iimage].wptr(exp_current_Fsize2);
            for (int n = 0; n < exp_current_Fsize2; n++) {
                exp_local_Fctfs_iimage[n] = 1.0;
            }
        }
    }
    // replace exp_local_Minvsigma2s with exp_local_Fctfs*exp_local_Minvsigma2s
    auto& Mctf_invsigma2 = exp_local_Minvsigma2s;
#pragma omp parallel for
    for (int iimage = 0; iimage < exp_nr_images; iimage++) {
        auto Mctf_invsigma2_iimage = Mctf_invsigma2[iimage].wptr(exp_current_Fsize2);
        auto exp_local_Fctfs_iimage = exp_local_Fctfs[iimage].wptr(exp_current_Fsize2);
        auto exp_local_Minvsigma2s_iimage = exp_local_Minvsigma2s[iimage].wptr(exp_current_Fsize2);
        for (int n = 0; n < exp_current_Fsize2; n++) {
            Mctf_invsigma2_iimage[n] = exp_local_Fctfs_iimage[n]*exp_local_Minvsigma2s_iimage[n];
        }
    }

    // NOTE : multiply nomask image by Minvsigma2s and Fctfs(exp_local_Minvsigma2s_X_Fctfs)
    // this will reduce some memory access in inner for loop
#pragma omp parallel for
    for (int iimage = 0; iimage < exp_nr_images; iimage++)
    {
        auto Fimages_nomask_real = particleModel.Fimages_nomask_real[iimage].wptr(exp_current_Fsize2);
        auto Fimages_nomask_imag = particleModel.Fimages_nomask_imag[iimage].wptr(exp_current_Fsize2);
        auto exp_local_Minvsigma2s_X_Fctfs_iimage = Mctf_invsigma2[iimage].wptr(exp_current_Fsize2);
        for (int n = 0; n < exp_current_Fsize2; n++) {
            Fimages_nomask_real[n] *= exp_local_Minvsigma2s_X_Fctfs_iimage[n];
            Fimages_nomask_imag[n] *= exp_local_Minvsigma2s_X_Fctfs_iimage[n];
        }
    }

    // replace exp_local_Minvsigma2s with exp_local_Fctfs^2*exp_local_Minvsigma2s
    auto& Mctf2_invsigma2 = exp_local_Minvsigma2s;
#pragma omp parallel for
    for (int iimage = 0; iimage < exp_nr_images; iimage++) {
        auto Mctf2_invsigma2_iimage = Mctf2_invsigma2[iimage].wptr(exp_current_Fsize2);
        auto exp_local_Fctfs_iimage = exp_local_Fctfs[iimage].wptr(exp_current_Fsize2);
        for (int n = 0; n < exp_current_Fsize2; n++) {
            Mctf2_invsigma2_iimage[n] = exp_local_Fctfs_iimage[n]*Mctf2_invsigma2_iimage[n];
        }
    }
    
    //
    particleModel.pregetShiftedImagesNomask(exp_Fimgs_shifted_real, exp_Fimgs_shifted_imag, exp_current_size, samplingGrid);
    
#pragma omp parallel for collapse(3)
    for (int iclass = exp_iclass_min; iclass <= exp_iclass_max; iclass++)//50 or 100
    {
        for (int ipsi = 0; ipsi < exp_nr_psi; ipsi++)//36
        {
            for (int iover_rot = 0; iover_rot < exp_nr_over_rot; iover_rot++)//2
            {
                int iclass_rot = iclass * exp_nr_psi + ipsi;
                int iclass_rot_over = iclass_rot * exp_nr_over_rot + iover_rot;
                
                auto Frefs_Rot_real = exp_Frefs_Rot_real.wptr(iclass, ipsi*exp_nr_over_rot+iover_rot, exp_current_Fsize2);
                auto Frefs_Rot_imag = exp_Frefs_Rot_imag.wptr(iclass, ipsi*exp_nr_over_rot+iover_rot, exp_current_Fsize2);
                auto Fweight_Rot = exp_Fweight_Rot.wptr(iclass, ipsi*exp_nr_over_rot+iover_rot, exp_current_Fsize2);
                
                for (int n = 0; n < exp_current_Fsize2; n++)
                    Frefs_Rot_real[n] = Frefs_Rot_imag[n] = Fweight_Rot[n] = 0.;

            }
        }
    }
    
#ifdef DATA_STREAM
    // turn off data stream inside openmp for loop
    global_data_stream.turnOff();
#endif
    
#pragma omp parallel for collapse(3) schedule(dynamic)
    for (int iclass = exp_iclass_min; iclass <= exp_iclass_max; iclass++)
    {
        for (int ipsi = 0; ipsi < exp_nr_psi; ipsi++)//36
        {
            for (int iover_rot = 0; iover_rot < exp_nr_over_rot; iover_rot++)
            {
                if (!exp_Rot_significant[iclass].rptrAll()[ipsi])
                    continue;
                
                int iclass_rot = iclass * exp_nr_psi + ipsi;
                int iclass_rot_over = iclass_rot * exp_nr_over_rot + iover_rot;
                
                auto Frefctf_real 	= exp_Frefs_Rot_real.wptr(iclass, ipsi*exp_nr_over_rot+iover_rot, exp_current_Fsize2);
                auto Frefctf_imag 	= exp_Frefs_Rot_imag.wptr(iclass, ipsi*exp_nr_over_rot+iover_rot, exp_current_Fsize2);
                auto Fweight_Rot	= exp_Fweight_Rot.wptr(iclass, ipsi*exp_nr_over_rot+iover_rot, exp_current_Fsize2);
                // Inside the loop over all translations and all part_id sum all shift Fimg's and their weights
                // Then outside this loop do the actual backprojection
                // Now that reference projection has been made loop over someParticles!
                for (int iimage = 0; iimage < exp_nr_images; iimage++)
                {
                    for (int itrans = 0; itrans < exp_nr_trans; itrans++)
                    {
                        for (int iover_trans = 0; iover_trans < exp_nr_over_trans; iover_trans++)
                        {
                            auto Mctf2_invsigma2 = exp_local_Minvsigma2s[iimage].wptr(exp_current_Fsize2);
                            
                            int iclass_rot_trans_over = iclass_rot_over*exp_nr_trans*exp_nr_over_trans + itrans*exp_nr_over_trans+iover_trans;
                            double weight = exp_Mweight[iimage*exp_Mweight_xsize+iclass_rot_trans_over];
                            
                            // Only deal with this sampling point if its weight was significant
                            if (weight >= exp_significant_weight.rptrAll()[iimage])
                            {
                                weight /= exp_sum_weight.rptrAll()[iimage];
                                // Get the shifted image
                                int itrans_over_trans = itrans*exp_nr_over_trans+iover_trans;
                                auto Fimg_shift_nomask_real = exp_Fimgs_shifted_real.wptr(iimage, itrans_over_trans, exp_current_Fsize2);
                                auto Fimg_shift_nomask_imag = exp_Fimgs_shifted_imag.wptr(iimage, itrans_over_trans, exp_current_Fsize2);
                                // Store sum of weight*SSNR*Fimg in data and sum of weight*SSNR in weight
                                // Use the FT of the unmasked image to back-project in order to prevent reconstruction artefacts! SS 25oct11
                                #pragma ivdep
                                for (int n = 0; n < exp_current_Fsize2; n++)
                                {

                                    Frefctf_real[n] += Fimg_shift_nomask_real[n] * weight;
                                    Frefctf_imag[n] += Fimg_shift_nomask_imag[n] * weight;
                                    // now Fweight stores sum of all w
                                    // Note that CTF needs to be squared in Fweight, weightxinvsigma2 already contained one copy
                                    Fweight_Rot[n] += Mctf2_invsigma2[n] * weight;
                                }
                            }// end if weight >= exp_significant_weight
#ifdef DATA_STREAM
                            // global_data_stream.foutInt(iclass, "backProjection()_iclass", __FILE__, __LINE__);
                            // global_data_stream.foutInt(idir, "backProjection()_idir", __FILE__, __LINE__);
                            // global_data_stream.foutInt(ipsi, "backProjection()_ipsi", __FILE__, __LINE__);
                            // global_data_stream.foutInt(iover_rot, "backProjection()_iover_rot", __FILE__, __LINE__);
                            global_data_stream.foutInt(iimage, "backProjection()_iimage", __FILE__, __LINE__);
                            global_data_stream.foutInt(itrans, "backProjection()_itrans", __FILE__, __LINE__);
                            global_data_stream.foutInt(iover_trans, "backProjection()_iover_trans", __FILE__, __LINE__);
                            global_data_stream.foutDouble(weight, "backProjection()_weight", __FILE__, __LINE__);
                            global_data_stream.foutDouble(exp_significant_weight.wptrAll()[iimage], "backProjection()_exp_significant_weight[iimage]", __FILE__, __LINE__);
                            global_data_stream.foutDouble(Frefctf_real, exp_current_Fsize2, "backProjection()_backProjection()_Fimg_real", __FILE__, __LINE__);
                            global_data_stream.foutDouble(Frefctf_imag, exp_current_Fsize2, "backProjection()_backProjection()_Fimg_imag", __FILE__, __LINE__);
                            global_data_stream.foutDouble(Fweight_Rot, exp_current_Fsize2, "backProjection()_backProjection()_Fweight", __FILE__, __LINE__);
                            global_data_stream.check();global_data_stream.flush();
#endif
                        }//end loop iover_trans
                    }// end loop itrans
                }// end loop iimage
            }// end loop iover_rot
        } // end loop iorient
    } //end loop iclass
    
#pragma omp parallel for
    for (int iclass = exp_iclass_min; iclass <= exp_iclass_max; iclass++)
    {
        for (int ipsi = 0; ipsi < exp_nr_psi; ipsi++)
        {
            for (int iover_rot = 0; iover_rot < exp_nr_over_rot; iover_rot++)
            {
                if (!exp_Rot_significant[iclass].rptrAll()[ipsi])
                    continue;
                
                int tid = omp_get_thread_num();
                
                auto Frefctf_real = exp_Frefs_Rot_real.wptr(iclass, ipsi*exp_nr_over_rot+iover_rot, exp_current_Fsize2);
                auto Frefctf_imag = exp_Frefs_Rot_imag.wptr(iclass, ipsi*exp_nr_over_rot+iover_rot, exp_current_Fsize2);
                auto Fweight = exp_Fweight_Rot.wptr(iclass, ipsi*exp_nr_over_rot+iover_rot, exp_current_Fsize2);
                
                FDOUBLE A[3][3];
                Euler_angles2matrix(0, 0, samplingGrid.exp_over_psi[ipsi][iover_rot], A);
                
                mapModel.set2DFourierTransform(tid, iclass, Frefctf_real, Frefctf_imag, exp_current_size, A, false, Fweight);
                
            }// end loop iover_rot
        } // end loop iorient
    } //end loop iclass
#ifdef DATA_STREAM
    global_data_stream.turnOn();
    int vol_size2 = mapModel.projector[0].pad_size;
    vol_size2 = vol_size2*(vol_size2/2+1);
    // global_data_stream.foutDoubleComplex((std::complex<double>*)mapModel.backprojector[0].data.wptr(), vol_size2, "backProjection()_backRef_0", __FILE__, __LINE__);
    // global_data_stream.foutDoubleComplex((std::complex<double>*)mapModel.backprojector[nr_classes-1].data.wptr(), vol_size2, "backProjection()_backRef_N-1", __FILE__, __LINE__);
    global_data_stream.foutDouble(mapModel.backprojector[0].weight.wptr(), vol_size2, "backProjection()_backRef_weight_0", __FILE__, __LINE__);
    global_data_stream.foutDouble(mapModel.backprojector[nr_classes-1].weight.wptr(), vol_size2, "backProjection()_backRef_weight_N-1", __FILE__, __LINE__);
    global_data_stream.check();global_data_stream.flush();
#endif
}

    
void maximization()
{

	NODE0ONLY std::cout << " Maximization ..." << std::endl;

	// First reconstruct the images for each class
#pragma omp parallel for
	for (int iclass = 0; iclass < nr_classes; iclass++)
	{
#ifdef DATA_STREAM
        global_data_stream.foutInt(iclass, "maximization()_iclass", __FILE__, __LINE__);
        global_data_stream.foutDouble(mapModel.Irefs[iclass].wptr(), mapModel.Irefs[iclass].dimzyx, "maximization()_mymodel.Iref[iclass]", __FILE__, __LINE__);
        global_data_stream.foutInt(gridding_nr_iter, "maximization()_gridding_nr_iter", __FILE__, __LINE__);
        global_data_stream.foutInt(do_map, "maximization()_do_map", __FILE__, __LINE__);
        global_data_stream.foutDouble(tau2_fudge_factor, "maximization()_mymodel.tau2_fudge_factor", __FILE__, __LINE__);
        global_data_stream.foutDouble(mlModel.tau2_class[iclass].wptr(ori_size/2+1), ori_size/2+1, "maximization()_mymodel.tau2_class[iclass]", __FILE__, __LINE__);
        global_data_stream.foutDouble(mlModel.sigma2_class[iclass].wptr(ori_size/2+1), ori_size/2+1, "maximization()_mymodel.sigma2_class[iclass]", __FILE__, __LINE__);
        global_data_stream.foutDouble(mlModel.data_vs_prior_class[iclass].wptr(ori_size/2+1), ori_size/2+1, "maximization()_mymodel.data_vs_prior_class[iclass]", __FILE__, __LINE__);
        global_data_stream.foutDouble(mlModel.fsc_halves_class[iclass].wptr(ori_size/2+1), ori_size/2+1, "maximization()_mymodel.fsc_halves_class[iclass]", __FILE__, __LINE__);
        global_data_stream.foutDouble(mlModel.wsum_pdf_class.wptrAll()[iclass], "maximization()_wsum_model.pdf_class[iclass]", __FILE__, __LINE__);
        global_data_stream.foutDouble(mapModel.minres_map, "maximization()_minres_map", __FILE__, __LINE__);
        global_data_stream.check();global_data_stream.flush();
#endif
        // TOOD : remember to remove this bug
        // TODO TODO TODO TODO TODO
        // if(mlModel.wsum_pdf_clsss[iclass] > 0.)
        if (mlModel.pdf_class.rptrAll()[iclass] > 0.)
		{
            // update
            auto sigma2_iclass = mlModel.sigma2_class[iclass].wptr(ori_size/2+1);
            mlModel.sigma2_class[iclass].zero();
            // update or not depend on update_tau2_with_fsc
            auto tau2_iclass = mlModel.tau2_class[iclass].wptr(ori_size/2+1);
            // update or not depend on update_tau2_with_fsc
            auto data_vs_prior_class_iclass = mlModel.data_vs_prior_class[iclass].wptr(ori_size/2+1);
            mlModel.data_vs_prior_class[iclass].zero();
            //
            auto fsc_halves_class_iclass = mlModel.fsc_halves_class[iclass].wptr(ori_size/2+1);
            
            // void reconstruction(int iclass,int max_iter_preweight,bool do_map,double tau2_fudge,double* tau2,double* sigma2,double* data_vs_prior,
            //                     const double* fsc, /* only input*/, double normalise = 1., bool update_tau2_with_fsc = false,
            //                     bool is_whole_instead_of_half = false,int minres_map = -1)
            
            mapModel.reconstruction(iclass, gridding_nr_iter, do_map, tau2_fudge_factor, tau2_iclass, sigma2_iclass,
                                    data_vs_prior_class_iclass, fsc_halves_class_iclass,
                                    mlModel.wsum_pdf_class.rptrAll()[iclass], false, false);
		}
		else
		{
            mapModel.Irefs[iclass].fill(0);
		}
#ifdef DATA_STREAM
        global_data_stream.foutInt(iclass, "maximization()_iclass", __FILE__, __LINE__);
        global_data_stream.foutDouble(mapModel.Irefs[iclass].wptr(), mapModel.Irefs[iclass].dimzyx, "maximization()_mymodel.Iref[iclass]", __FILE__, __LINE__);
        global_data_stream.foutInt(gridding_nr_iter, "maximization()_gridding_nr_iter", __FILE__, __LINE__);
        global_data_stream.foutInt(do_map, "maximization()_do_map", __FILE__, __LINE__);
        global_data_stream.foutDouble(tau2_fudge_factor, "maximization()_mymodel.tau2_fudge_factor", __FILE__, __LINE__);
        global_data_stream.foutDouble(mlModel.tau2_class[iclass].wptr(ori_size/2+1), ori_size/2+1, "maximization()_mymodel.tau2_class[iclass]", __FILE__, __LINE__);
        global_data_stream.foutDouble(mlModel.sigma2_class[iclass].wptr(ori_size/2+1), ori_size/2+1, "maximization()_mymodel.sigma2_class[iclass]", __FILE__, __LINE__);
        global_data_stream.foutDouble(mlModel.data_vs_prior_class[iclass].wptr(ori_size/2+1), ori_size/2+1, "maximization()_mymodel.data_vs_prior_class[iclass]", __FILE__, __LINE__);
        global_data_stream.foutDouble(mlModel.fsc_halves_class[iclass].wptr(ori_size/2+1), ori_size/2+1, "maximization()_mymodel.fsc_halves_class[iclass]", __FILE__, __LINE__);
        global_data_stream.foutDouble(mlModel.wsum_pdf_class.wptrAll()[iclass], "maximization()_wsum_model.pdf_class[iclass]", __FILE__, __LINE__);
        global_data_stream.foutDouble(mapModel.minres_map, "maximization()_minres_map", __FILE__, __LINE__);
        global_data_stream.check();global_data_stream.flush();
#endif
	}
    
	// Then perform the update of all other model parameters
	maximizationOtherParameters();
    
#ifdef DATA_STREAM
    global_data_stream.foutDouble(mlModel.avg_norm_correction, "maximizationOtherParameters()_model_avg_norm_correction", __FILE__, __LINE__);
    global_data_stream.foutDouble(mlModel.pdf_class.wptrAll()[0], "maximizationOtherParameters()_model_pdf_class[0]", __FILE__, __LINE__);
    global_data_stream.foutDouble(mlModel.pdf_class.wptrAll()[nr_classes-1], "maximizationOtherParameters()_model_pdf_class[nr_classes-1]", __FILE__, __LINE__);
    global_data_stream.foutDouble(mlModel.scale_correction.wptr(mlModel.nr_groups), mlModel.nr_groups, "maximizationOtherParameters()_model_scale_correction", __FILE__, __LINE__);
    global_data_stream.foutDouble(mlModel.pdf_direction[0].wptr(mlModel.nr_directions), mlModel.nr_directions, "maximizationOtherParameters()_model_pdf_direction[0]", __FILE__, __LINE__);
    global_data_stream.foutDouble(mlModel.pdf_direction[nr_classes-1].wptr(mlModel.nr_directions), mlModel.nr_directions, "maximizationOtherParameters()_model_pdf_direction[nr_classes-1]", __FILE__, __LINE__);
    global_data_stream.foutDouble(mlModel.sigma2_noise[0].wptr(mlModel.ori_Fsize), mlModel.ori_Fsize, "maximizationOtherParameters()_model_sigma2_noise1", __FILE__, __LINE__);
    global_data_stream.foutDouble(mlModel.sigma2_noise[mlModel.nr_groups-1].wptr(mlModel.ori_Fsize), mlModel.ori_Fsize, "maximizationOtherParameters()_model_sigma2_noiseN", __FILE__, __LINE__);
    global_data_stream.foutDouble(mlModel.sigma2_offset, "maximizationOtherParameters()_sigma2_offset", __FILE__, __LINE__);
    global_data_stream.foutDouble(mlModel.LL, "maximizationOtherParameters()_model_LL", __FILE__, __LINE__);
    global_data_stream.foutDouble(mlModel.ave_Pmax, "maximizationOtherParameters()_model_ave_Pmax", __FILE__, __LINE__);
    global_data_stream.foutDouble(mlModel.tau2_class[0].wptr(ori_size/2+1), ori_size/2+1, "maximizationOtherParameters()_model_tau2_class[0]", __FILE__, __LINE__);
    global_data_stream.foutDouble(mlModel.tau2_class[nr_classes-1].wptr(ori_size/2+1), ori_size/2+1, "maximizationOtherParameters()_model_tau2_class[nr_classes-1]", __FILE__, __LINE__);
    global_data_stream.foutDouble(mlModel.data_vs_prior_class[0].wptr(ori_size/2+1), ori_size/2+1, "maximizationOtherParameters()_data_vs_prior_class[0]", __FILE__, __LINE__);
    global_data_stream.foutDouble(mlModel.data_vs_prior_class[nr_classes-1].wptr(ori_size/2+1), ori_size/2+1, "maximizationOtherParameters()_data_vs_prior_class[nr_classes-1]", __FILE__, __LINE__);
    global_data_stream.foutDouble(mapModel.Irefs[0].wptr(), mapModel.Irefs[0].dimzyx, "maximizationOtherParameters()_classProjector.Irefs[0]", __FILE__, __LINE__);
    global_data_stream.foutDouble(mapModel.Irefs[nr_classes-1].wptr(), mapModel.Irefs[nr_classes-1].dimzyx, "maximizationOtherParameters()_classProjector.Irefs[nr_classes-1]", __FILE__, __LINE__);
    global_data_stream.check();global_data_stream.flush();
#endif
    
	// Keep track of changes in hidden variables
	hiddenVarMonitor.updateOverallChangesInHiddenVariables();
    
#ifdef DATA_STREAM
    global_data_stream.foutDouble(hiddenVarMonitor.current_changes_optimal_classes, "updateOverallChangesInHiddenVariables()_current_changes_optimal_classes", __FILE__, __LINE__);
    global_data_stream.foutDouble(hiddenVarMonitor.current_changes_optimal_orientations, "updateOverallChangesInHiddenVariables()_current_changes_optimal_orientations", __FILE__, __LINE__);
    global_data_stream.foutDouble(hiddenVarMonitor.current_changes_optimal_offsets, "updateOverallChangesInHiddenVariables()_current_changes_optimal_offsets", __FILE__, __LINE__);
    global_data_stream.foutDouble(hiddenVarMonitor.nr_iter_wo_large_hidden_variable_changes, "updateOverallChangesInHiddenVariables()_nr_iter_wo_large_hidden_variable_changes", __FILE__, __LINE__);
    global_data_stream.foutDouble(hiddenVarMonitor.smallest_changes_optimal_classes, "updateOverallChangesInHiddenVariables()_smallest_changes_optimal_classes", __FILE__, __LINE__);
    global_data_stream.foutDouble(hiddenVarMonitor.smallest_changes_optimal_offsets, "updateOverallChangesInHiddenVariables()_smallest_changes_optimal_offsets", __FILE__, __LINE__);
    global_data_stream.foutDouble(hiddenVarMonitor.smallest_changes_optimal_orientations, "updateOverallChangesInHiddenVariables()_smallest_changes_optimal_orientations", __FILE__, __LINE__);
    global_data_stream.check();global_data_stream.flush();
#endif
    
}


void maximizationOtherParameters()
{
	// Calculate total sum of weights, and average CTF for each class (for SSNR estimation)
	double sum_weight = 0.;
	for (int iclass = 0; iclass < nr_classes; iclass++)
        sum_weight += mlModel.wsum_pdf_class.rptrAll()[iclass];//sum_weight += wsum_model.pdf_class[iclass];

	// Update average norm_correction
	if (do_norm_correction)
	{
        mlModel.avg_norm_correction = mlModel.wsum_avg_norm_correction / sum_weight;
	}


	// Update model.pdf_class vector (for each k)
	for (int iclass = 0; iclass < nr_classes; iclass++)
	{
        mlModel.pdf_class.wptrAll()[iclass] = mlModel.wsum_pdf_class.rptrAll()[iclass] / sum_weight;

		// for 2D also update priors of translations for each class!
        if (mlModel.wsum_pdf_class.rptrAll()[iclass] > 0.){
            mlModel.prior_offsetx_class.wptrAll()[iclass] = mlModel.wsum_prior_offsetx_class.rptrAll()[iclass] / mlModel.wsum_pdf_class.rptrAll()[iclass];
            mlModel.prior_offsety_class.wptrAll()[iclass] = mlModel.wsum_prior_offsety_class.rptrAll()[iclass] / mlModel.wsum_pdf_class.rptrAll()[iclass];
        }
        else{
            mlModel.prior_offsetx_class.wptrAll()[iclass] = 0.;
            mlModel.prior_offsety_class.wptrAll()[iclass] = 0.;
        }

        // only one direction
        assert(mlModel.nr_directions==1);
        mlModel.pdf_direction[iclass][0] = mlModel.wsum_pdf_direction[iclass][0] / sum_weight;

	}


    mlModel.sigma2_offset = (mlModel.wsum_sigma2_offset) / (2. * sum_weight);

	// TODO: update estimates for sigma2_rot, sigma2_tilt and sigma2_psi!
    for (int igroup = 0; igroup < mlModel.nr_groups; igroup++)
    {
        auto sigma2_noise_igroup = mlModel.sigma2_noise[igroup].wptr(mlModel.ori_Fsize);
        auto wsum_sigma2_noise_igroup = mlModel.wsum_sigma2_noise[igroup].wptr(mlModel.ori_Fsize);
        // Factor 2 because of the 2-dimensionality of the complex-plane
        for (int n = 0; n < mlModel.ori_Fsize; n++)
        {
            sigma2_noise_igroup[n] = wsum_sigma2_noise_igroup[n] / (2. * mlModel.wsum_sumw_group.rptrAll()[igroup] * Npix_per_shell.rptrAll()[n]);
        }
    }


	// After the first iteration the references are always CTF-corrected
    if (do_ctf_correction)
    	refs_are_ctf_corrected = true;

	// Some statistics to output

    mlModel.LL = mlModel.wsum_LL;
    
    mlModel.ave_Pmax = mlModel.wsum_ave_Pmax / sum_weight;

}

void writeClassesAndMetadata(std::string fn_class,std::string fn_metadata){
    
    NODE0ONLY
    {
        // A. Write Classes
        mapModel.writeResult(fn_class);
        // B.Write Metadata
        metadata.writeToStar(fn_metadata);
    }
}

//TODO: could this inner for loop be vectorized
inline bool isSignificantAnyParticleAnyTranslation(int iorient)
{
    int ihidden_start = iorient * exp_nr_trans;
    int ihidden_end = iorient * exp_nr_trans + exp_nr_trans;
    for (int iimage = 0; iimage < exp_nr_images; iimage++)
    {
        for (int ihidden = ihidden_start; ihidden < ihidden_end; ihidden++)
        {
            if (exp_Mcoarse_significant[iimage].rptrAll()[ihidden])
                return true;
        }
    }
    return false;
}

void calspace(int ori_size,int current_size,int set_nr_pool)
{
    // TOOD TODO use the Heap
    // to statistics all the datastructure size
    // TODO TODO TODO...
    
//    double perGb = 1/1024./1024./1024.;
//    
//    double images_data_size               = nr_local_images*sizeof(double)*ori_size*(ori_size/2+1)*perGb;
//    double metadata_size                  = nr_global_images*sizeof(MetaDataElem)*perGb;
//    double images_id_size                 = nr_global_images*sizeof(int)*perGb;
//    double model_Irefs_size               = nr_classes*ori_size*ori_size*sizeof(double)*perGb;
//    double model_tau2_class_size          = nr_classes*(ori_size/2+1)*sizeof(double)*perGb;
//    double model_sigma2_noise_size        = (ori_size/2+1)*sizeof(double)*perGb;
//    double model_data_vs_prior_class_size = nr_classes*(ori_size/2+1)*sizeof(double)*perGb;
//    double model_sigma2_class_size        = nr_classes*(ori_size/2+1)*sizeof(double)*perGb;
//    double model_fsc_halves_class_size    = nr_classes*(ori_size/2+1)*sizeof(double)*perGb;
//    double model_pdf_class_size           = nr_classes*sizeof(double)*perGb;
//    double model_pdf_direction_size       = nr_classes*sizeof(double)*perGb;
//    double model_prior_offsetx_class_size = nr_classes*sizeof(double)*perGb;
//    double model_prior_offsety_class_size = nr_classes*sizeof(double)*perGb;
//    double wsum_sigma2_noise_size         = (ori_size/2+1)*sizeof(double)*perGb;
//    double wsum_pdf_class_size            = nr_classes*sizeof(double)*perGb;
//    double wsum_prior_offsetx_class_size  = nr_classes*sizeof(double)*perGb;
//    double wsum_prior_offsety_class_size  = nr_classes*sizeof(double)*perGb;
//    double wsum_pdf_direction_size        = nr_classes*sizeof(double)*perGb;
//    double fix_size = images_data_size + metadata_size + images_id_size + model_Irefs_size + model_tau2_class_size + model_sigma2_noise_size \
//                    +  model_data_vs_prior_class_size + model_sigma2_class_size + model_fsc_halves_class_size + model_pdf_class_size + model_pdf_direction_size \
//                    +  model_prior_offsetx_class_size + model_prior_offsety_class_size + wsum_sigma2_noise_size + wsum_pdf_class_size + wsum_prior_offsetx_class_size \
//                    +  wsum_prior_offsety_class_size + wsum_pdf_direction_size;
//    
//    std::cout<<"fixed size : "<<std::endl;
//    std::cout<<"images_data size               : "<<images_data_size<<" GB."<<std::endl;
//    std::cout<<"metadata size                  : "<<metadata_size<<" GB."<<std::endl;
//    std::cout<<"image_id size                  : "<<images_id_size<<" GB."<<std::endl;
//    std::cout<<"model_Irefs size               : "<<model_Irefs_size<<" GB."<<std::endl;
//    std::cout<<"model_tau2_class size          : "<<model_tau2_class_size<<" GB."<<std::endl;
//    std::cout<<"model_sigma2_noise size        : "<<model_sigma2_noise_size<<" GB."<<std::endl;
//    std::cout<<"model_data_vs_prior_class size : "<<model_data_vs_prior_class_size<<" GB."<<std::endl;
//    std::cout<<"model_sigma2_class size        : "<<model_sigma2_class_size<<" GB."<<std::endl;
//    std::cout<<"model_fsc_halves_class size    : "<<model_fsc_halves_class_size<<" GB."<<std::endl;
//    std::cout<<"model_pdf_class size           : "<<model_pdf_class_size<<" GB."<<std::endl;
//    std::cout<<"model_pdf_direction size       : "<<model_pdf_class_size<<" GB."<<std::endl;
//    std::cout<<"model_prior_offsetx_class size : "<<model_prior_offsetx_class_size<<" GB."<<std::endl;
//    std::cout<<"model_prior_offsety_class size : "<<model_prior_offsety_class_size<<" GB."<<std::endl;
//    std::cout<<"wsum_sigma2_noise size         : "<<wsum_sigma2_noise_size<<" GB."<<std::endl;
//    std::cout<<"wsum_pdf_class size            : "<<wsum_pdf_class_size<<" GB."<<std::endl;
//    std::cout<<"wsum_prior_offsetx_class size  : "<<wsum_prior_offsetx_class_size<<" GB."<<std::endl;
//    std::cout<<"wsum_prior_offsety_class size  : "<<wsum_prior_offsety_class_size<<" GB."<<std::endl;
//    std::cout<<"wsum_pdf_direction size        : "<<wsum_pdf_direction_size<<" GB."<<std::endl;
//    std::cout<<"total fixed size               : "<<fix_size<<" GB"<<std::endl;
//    std::cout<<"------------------------------------------------------"<<std::endl;
//    
//
//    int r_max = std::min(current_size / 2, ori_size / 2);
//    int model_Frefs_pad_dim = 2 * (padding_factor * r_max + 1) + 1;
//    
//    double model_Frefs_pad_size          = 2*nr_classes*model_Frefs_pad_dim*(model_Frefs_pad_dim/2+1)*sizeof(double)*perGb;
//    double wsum_data_size                = 2*nr_classes*model_Frefs_pad_dim*(model_Frefs_pad_dim/2+1)*sizeof(double)*perGb;
//    double wsum_weight_size              = nr_classes*model_Frefs_pad_dim*(model_Frefs_pad_dim/2+1)*sizeof(double)*perGb;
//    
//    double exp_over_rot_psi_size         = sampler2d.NrPsi(adaptive_oversampling)*sizeof(double)*perGb;
//    double exp_over_trans_x_size         = sampler2d.NrTrans(adaptive_oversampling)*sizeof(double)*perGb;
//    double exp_over_trans_y_size         = sampler2d.NrTrans(adaptive_oversampling)*sizeof(double)*perGb;
//    double exp_Frefs_Rot_size            = 2*nr_classes*sampler2d.NrPsi(adaptive_oversampling)*current_size*(current_size/2+1)*sizeof(double)*perGb;
//    double exp_Fweight_size              = nr_classes*sampler2d.NrPsi(adaptive_oversampling)*current_size*(current_size/2+1)*sizeof(double)*perGb;
//    double exp_Rot_significant_size      = nr_classes*sampler2d.NrPsi(adaptive_oversampling)*sizeof(bool)*perGb;
//    double unfixed_size = model_Frefs_pad_size + wsum_data_size + wsum_weight_size + exp_over_rot_psi_size + exp_over_trans_x_size \
//                        + exp_over_trans_y_size +  exp_Frefs_Rot_size + exp_Fweight_size + exp_Rot_significant_size;
//    
//    
//    set_nr_pool = 1;
//    double exp_metadata_size             = set_nr_pool*sizeof(MetaDataElem)*perGb;
//    double exp_imgs_size                 = set_nr_pool*ori_size*ori_size*sizeof(double)*perGb;
//    double exp_power_imgs_size           = set_nr_pool*(ori_size/2+1)*sizeof(double)*perGb;
//    double exp_highres_Xi2_imgs_size     = set_nr_pool*sizeof(double)*perGb;
//    double exp_min_diff2_size            = set_nr_pool*sizeof(double)*perGb;
//    double exp_old_offsetx_size          = set_nr_pool*sizeof(double)*perGb;
//    double exp_old_offsety_size          = set_nr_pool*sizeof(double)*perGb;
//    double exp_wsum_norm_correction_size = set_nr_pool*sizeof(double)*perGb;
//    double exp_Fimgs_size                = 2*set_nr_pool*current_size*(current_size/2+1)*sizeof(double)*perGb;
//    double exp_Fimgs_nomask_size         = 2*set_nr_pool*current_size*(current_size/2+1)*sizeof(double)*perGb;
//    double exp_Fctfs_size                = set_nr_pool*current_size*(current_size/2+1)*sizeof(double)*perGb;
//    double exp_significant_weight_size   = set_nr_pool*sizeof(double)*perGb;
//    double exp_sum_weight_size           = set_nr_pool*sizeof(double)*perGb;
//    double exp_Fimgs_shifted_size        = 2*set_nr_pool*sampler2d.NrTrans(adaptive_oversampling)*current_size*(current_size/2+1)*sizeof(double)*perGb;
//    double exp_Fimgs_shifted_nomask_size = 0;//2*set_nr_pool*sampling2d.NrTrans(adaptive_oversampling)*current_size*(current_size/2+1)*sizeof(double)*perGb;
//    double exp_local_Fctfs_size          = set_nr_pool*current_size*(current_size/2+1)*sizeof(double)*perGb;
//    double exp_local_Minvsigma2s_size    = set_nr_pool*current_size*(current_size/2+1)*sizeof(double)*perGb;
//    double exp_Mweight_size              = set_nr_pool*nr_classes*sampler2d.NrPoints(adaptive_oversampling)*sizeof(double)*perGb;
//    double exp_Mcoarse_significant_size  = set_nr_pool*nr_classes*sampler2d.NrPoints(0)*sizeof(bool)*perGb;
//    double nr_pool_size = exp_metadata_size + exp_imgs_size + exp_power_imgs_size + exp_highres_Xi2_imgs_size + exp_min_diff2_size \
//                        + exp_old_offsetx_size + exp_old_offsety_size + exp_wsum_norm_correction_size + exp_Fimgs_size + exp_Fimgs_nomask_size \
//                        + exp_Fctfs_size + exp_significant_weight_size + exp_sum_weight_size + exp_Fimgs_shifted_size + exp_Fimgs_shifted_nomask_size \
//                        + exp_local_Fctfs_size + exp_local_Minvsigma2s_size + exp_Mweight_size + exp_Mcoarse_significant_size;
//    
//    set_nr_pool = (6.5 - fix_size - unfixed_size)/nr_pool_size;
//    std::cout<<"------------------------------------------------------"<<std::endl;
//    std::cout<<"suggestion nr_pool             : "<<set_nr_pool<<std::endl;
//    std::cout<<"real nr_pool                   : "<<nr_pool<<std::endl;
//    std::cout<<"------------------------------------------------------"<<std::endl;
//    //nr_pool = set_nr_pool;
//    unfixed_size += nr_pool_size*nr_pool;
//    
//    std::cout<<"unfixed size(Expectation step) : "<<std::endl;
//    std::cout<<"wsum_data size                 : "<<wsum_data_size<<" GB."<<std::endl;
//    std::cout<<"exp_over_rot_psi size          : "<<exp_over_rot_psi_size<<" GB."<<std::endl;
//    std::cout<<"exp_over_trans_x size          : "<<exp_over_trans_x_size<<" GB."<<std::endl;
//    std::cout<<"exp_over_trans_y size          : "<<exp_over_trans_y_size<<" GB."<<std::endl;
//    std::cout<<"exp_Frefs_Rot size             : "<<exp_Frefs_Rot_size<<" GB."<<std::endl;
//    std::cout<<"exp_Fweight_size size          : "<<exp_Fweight_size<<" GB."<<std::endl;
//    std::cout<<"exp_Rot_significant size       : "<<exp_Rot_significant_size<<" GB."<<std::endl;
//    std::cout<<"------------------------------------------------------"<<std::endl;
//    std::cout<<"exp_metadata size              : "<<nr_pool*exp_metadata_size<<" GB."<<std::endl;
//    std::cout<<"exp_imgs size                  : "<<nr_pool*exp_imgs_size<<" GB."<<std::endl;
//    std::cout<<"exp_power_imgs size            : "<<nr_pool*exp_power_imgs_size<<" GB."<<std::endl;
//    std::cout<<"exp_highres_Xi2_imgs size      : "<<nr_pool*exp_highres_Xi2_imgs_size<<" GB."<<std::endl;
//    std::cout<<"exp_min_diff2 size             : "<<nr_pool*exp_min_diff2_size<<" GB."<<std::endl;
//    std::cout<<"exp_old_offsetx_size size      : "<<nr_pool*exp_old_offsetx_size<<" GB."<<std::endl;
//    std::cout<<"exp_old_offsety size           : "<<nr_pool*exp_old_offsety_size<<" GB."<<std::endl;
//    std::cout<<"exp_wsum_norm_correction size  : "<<nr_pool*exp_wsum_norm_correction_size<<" GB."<<std::endl;
//    std::cout<<"exp_Fimgs size                 : "<<nr_pool*exp_Fimgs_size<<" GB."<<std::endl;
//    std::cout<<"exp_Fimgs_nomask size          : "<<nr_pool*exp_Fimgs_nomask_size<<" GB."<<std::endl;
//    std::cout<<"exp_Fctfs size                 : "<<nr_pool*exp_Fctfs_size<<" GB."<<std::endl;
//    std::cout<<"exp_significant_weight size    : "<<nr_pool*exp_significant_weight_size<<" GB."<<std::endl;
//    std::cout<<"exp_sum_weight size            : "<<nr_pool*exp_sum_weight_size<<" GB."<<std::endl;
//    std::cout<<"exp_Fimgs_shifted size         : "<<nr_pool*exp_Fimgs_shifted_size<<" GB."<<std::endl;
//    std::cout<<"exp_Fimgs_shifted_nomask size  : "<<nr_pool*exp_Fimgs_shifted_nomask_size<<" GB."<<std::endl;
//    std::cout<<"exp_local_Fctfs size           : "<<nr_pool*exp_local_Fctfs_size<<" GB."<<std::endl;
//    std::cout<<"exp_local_Minvsigma2s size     : "<<nr_pool*exp_local_Minvsigma2s_size<<" GB."<<std::endl;
//    std::cout<<"exp_Mweight size               : "<<nr_pool*exp_Mweight_size<<" GB."<<std::endl;
//    std::cout<<"exp_Mcoarse_significant size   : "<<nr_pool*exp_Mcoarse_significant_size<<" GB."<<std::endl;
//    std::cout<<"total unfixed size             : "<<unfixed_size<<" GB."<<std::endl;
//    std::cout<<"------------------------------------------------------"<<std::endl;
//    
//    
//    std::cout<<"thread size(Expectation step) thread variable size : "<<std::endl;
//    
//    double thread_wsum_sigma2_offset_size = 0;//maxthreads*sizeof(double)*perGb;
//    double thread_sumw_group_size = 0;//maxthreads*sizeof(double)*perGb;
//    double thread_max_weight_size = maxthreads*nr_pool*sizeof(double)*perGb;
//    double thread_wsum_sigma2_noise_size = maxthreads*current_size*(current_size/2+1)*sizeof(double)*perGb;
//    double thread_wsum_norm_correction_size = maxthreads*current_size*(current_size/2+1)*sizeof(double)*perGb;
//    double thread_wsum_pdf_direction_size = maxthreads*nr_classes*sizeof(double)*perGb;
//    double thread_wsum_pdf_class_size = maxthreads*nr_classes*sizeof(double)*perGb;
//    double thread_wsum_prior_offsetx_class_size = maxthreads*nr_classes*sizeof(double)*perGb;
//    double thread_wsum_prior_offsety_class_size = maxthreads*nr_classes*sizeof(double)*perGb;
//    double thread_wsum_data_size = 0;//2*maxthreads*nr_classes*model_Frefs_pad_dim*(model_Frefs_pad_dim/2+1)*sizeof(double)*perGb;
//    double thread_wsum_weight_size = 0;//maxthreads*nr_classes*model_Frefs_pad_dim*(model_Frefs_pad_dim/2+1)*sizeof(double)*perGb;
//    double thread_size = thread_wsum_sigma2_offset_size + thread_sumw_group_size + thread_max_weight_size + thread_wsum_sigma2_noise_size \
//                       + thread_wsum_norm_correction_size + thread_wsum_pdf_direction_size + thread_wsum_pdf_class_size \
//                       + thread_wsum_prior_offsetx_class_size + thread_wsum_prior_offsety_class_size + thread_wsum_data_size \
//                       + thread_wsum_weight_size;
//    std::cout<<"thread size(Expectation step)  : "<<thread_size<<" GB."<<std::endl;
//    std::cout<<"------------------------------------------------------"<<std::endl;
//    std::cout<<"total size                     : "<<(fix_size+unfixed_size+thread_size)<<" GB."<<std::endl;
//    std::cout<<"------------------------------------------------------"<<std::endl;
    
}

void debugStoreWeightedSums(){
    
    std::cout << " WARNING: norm_correction : "<< exp_metadata[0].NORM  << " for particle " << 0 <<std::endl;
    std::cout << " mymodel.current_size : " << current_size << " mymodel.ori_size= " << ori_size <<std::endl;
    std::cout << " coarse_size : " << coarse_size << std::endl;
    std::cout << " DIRECT_A2D_ELEM(exp_metadata2, my_image_no, exp_nr_imagas-1) : " <<exp_metadata[exp_nr_images-1].NORM << std::endl;
    std::cout << " mymodel.avg_norm_correction : " << mlModel.avg_norm_correction << std::endl;
    std::cout << " exp_wsum_norm_correction[ipart] : " << exp_wsum_norm_correction.rptrAll()[0] << std::endl;
    // std::cout << " old_norm_correction : " << old_norm_correction << std::endl;
    std::cout << " wsum_model.avg_norm_correction : " << mlModel.wsum_avg_norm_correction << std::endl;
    // std::cout << " group_id : " << group_id << " mymodel.scale_correction[group_id] : " << mymodel.scale_correction[group_id] << std::endl;
    std::cout << " mymodel.sigma2_noise[group_id = 0] : " << mlModel.sigma2_noise[0][0] << std::endl;
    std::cout << " wsum_model.sigma2_noise[group_id = 0] : " << mlModel.wsum_sigma2_noise[0][0] << std::endl;
    std::cout << " exp_power_imgs[my_image_no = 0] : " << exp_power_imgs[0][0] << std::endl;
    // std::cout << " exp_wsum_scale_correction_XA[ipart] : " << exp_wsum_scale_correction_XA[ipart] << " exp_wsum_scale_correction_AA[ipart] : " << exp_wsum_scale_correction_AA[ipart] << std::endl;
    // std::cout << " wsum_model.wsum_signal_product_spectra[group_id] : " << wsum_model.wsum_signal_product_spectra[group_id] << " wsum_model.wsum_reference_power_spectra[group_id] : " << wsum_model.wsum_reference_power_spectra[group_id] << std::endl;
    std::cout << " exp_min_diff2[ipart = 0] : " << exp_min_diff2.rptrAll()[0] << std::endl;
    // std::cout << " ml_model.scale_correction[group_id] : " << mymodel.scale_correction[group_id] << std::endl;
    std::cout << " exp_significant_weight[ipart = 0] : " << exp_significant_weight.rptrAll()[0] << std::endl;
    // std::cout << " exp_max_weight[ipart = 0] : " << exp_max_weight[0] << std::endl;
    
    std::cerr << " part_id : " << 0 << std::endl;
    std::cerr << " ipart : " << 0 << std::endl;
    std::cerr << " exp_min_diff2[ipart = 0] : " << exp_min_diff2.rptrAll()[0] << std::endl;
    // std::cerr << " logsigma2 : " << logsigma2 << std::endl;
    int group_id = 0;//mydata.getGroupId(part_id, 0);
    std::cerr << " group_id : " << group_id << std::endl;
    // std::cerr << " ml_model.scale_correction[group_id] : " << mymodel.scale_correction[group_id] << std::endl;
    std::cerr << " exp_significant_weight[ipart = 0] : " << exp_significant_weight.rptrAll()[0] << std::endl;
    // std::cerr << " exp_max_weight[ipart = 0]= " << exp_max_weight[0] << std::endl;
    std::cerr << " ml_model.sigma2_noise[group_id = 0] : " << mlModel.sigma2_noise[0][0] << std::endl;
}

} // end namespace MLoptimizer