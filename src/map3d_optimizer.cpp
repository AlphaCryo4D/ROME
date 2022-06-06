/***************************************************************************
 *
 * Authors: "Yongbei(Glow) Ma,Jiayi (Timmy) Wu, Youdong (Jack) Mao"
 * Dana-Farber Cancer Institute, Harvard Medical School and Peking University
 * Bevin Brett
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

#include "./map3d_optimizer.h"

namespace Map3dOptimizer_old {
    //
    DataStream global_data_stream;
    
    const int thread_per_L2_cache = 4;
    // configuration
    // workq_bdwb0_hyperthread // boardwell
    // 36 core and hyperthread,so 4 threads will share one L2 cache
    // const int thread_per_L2_cache = 4;
    // knldevq_ekl // KNL
    // const int thread_per_L2_cache = 2;
    
    // model
    HealpixSampler sampler3d;
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
    int sampler3d_healpix_order;
    std::string sampler3d_fn_sym;
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
    Exp_Mweight_old exp_Mweight_coarse;
    Exp_Mweight_old exp_Mweight_fine;
    
    // ------------ thread variable ----------------- //
    Aligned3dArray<FDOUBLE> thread_Frefctf_real;
    Aligned3dArray<FDOUBLE> thread_Frefctf_imag;
    Aligned3dArray<FDOUBLE> thread_Fimg_real;
    Aligned3dArray<FDOUBLE> thread_Fimg_imag;
    Aligned3dArray<FDOUBLE> thread_Fweight;
    Aligned3dArray<FDOUBLE> thread_wsum_pdf_direction;
    Aligned2dArray< char > threadfake_do_scale_norm_class;
    std::vector<std::vector<GridIndex>> thread_exp_max_weight_index;
    //
    //
    StatusTracer statusTracer;
    
    void setupStatusTracer()
    {
        statusTracer.clear();
        for (int i = 0; i < mapModel.Irefs.size(); i++) {
#if defined(FLOAT_PRECISION)
            statusTracer.appendFloatPtr(mapModel.Irefs[i].wptr(), mapModel.Irefs[i].dimzyx, "mapmodel_Irefs", true);
#else
            statusTracer.appendDoublePtr(mapModel.Irefs[i].wptr(), mapModel.Irefs[i].dimzyx, "mapmodel_Irefs"); // model read
#endif
        }
        //
        statusTracer.appendDoublePtr(&mapModel.current_resolution, 1, "mapModel_current_resolution"); // model read
        //
        for (int igroup = 0; igroup < mlModel.nr_groups; igroup++) {
#if defined(FLOAT_PRECISION)
            statusTracer.appendFloatPtr(mlModel.sigma2_noise[igroup].wptr(mlModel.ori_Fsize), mlModel.ori_Fsize, "mlmodel_sigma2_noise", true);
            statusTracer.appendFloatPtr(mlModel.wsum_signal_product_spectra[igroup].wptr(mlModel.ori_Fsize), mlModel.ori_Fsize, "mlmodel_wsum_signal_product_spectra", true);
            statusTracer.appendFloatPtr(mlModel.wsum_reference_power_spectra[igroup].wptr(mlModel.ori_Fsize), mlModel.ori_Fsize, "mlmodel_wsum_reference_power_spectra", true);
#else
            statusTracer.appendDoublePtr(mlModel.sigma2_noise[igroup].wptr(mlModel.ori_Fsize),
                                         mlModel.ori_Fsize, "mlmodel_sigma2_noise_"+num2str(igroup,6)); // model read
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
                                         mlModel.ori_Fsize, "mlmodel_tau2_class_"+num2str(iclass,6)); // model read
            statusTracer.appendDoublePtr(mlModel.sigma2_class[iclass].wptr(mlModel.ori_Fsize),
                                         mlModel.ori_Fsize, "mlmodel_sigma2_class_"+num2str(iclass,6)); // model read
            statusTracer.appendDoublePtr(mlModel.data_vs_prior_class[iclass].wptr(mlModel.ori_Fsize),
                                         mlModel.ori_Fsize, "mlmodel_data_vs_prior_class_"+num2str(iclass,6)); // model read
            statusTracer.appendDoublePtr(mlModel.fsc_halves_class[iclass].wptr(mlModel.ori_Fsize),
                                         mlModel.ori_Fsize, "mlmodel_fsc_halves_class_"+num2str(iclass,6)); // model read
            statusTracer.appendDoublePtr(mlModel.pdf_direction[iclass].wptr(mlModel.nr_directions),
                                         mlModel.nr_directions, "mlmodel_pdf_direction_"+num2str(iclass,6));// will set 1 if increase heal_pix order
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
                                     mlModel.nr_groups, "mlmodel_scale_correction"); // model read
        statusTracer.appendDoublePtr(mlModel.prior_offsetx_class.wptr(mlModel.nr_classes),
                                     mlModel.nr_classes, "mlmodel_prior_offsetx_class");// only read in 2D
        statusTracer.appendDoublePtr(mlModel.prior_offsety_class.wptr(mlModel.nr_classes),
                                     mlModel.nr_classes, "mlmodel_prior_offsety_class");// only read in 2D
        statusTracer.appendDoublePtr(mlModel.pdf_class.wptr(mlModel.nr_classes),
                                     mlModel.nr_classes, "mlmodel_pdf_class"); // model read
        //
        statusTracer.appendDoublePtr(&mlModel.avg_norm_correction, 1, "mlmodel_avg_norm_correction"); // model read
        statusTracer.appendDoublePtr(&mlModel.ave_Pmax, 1, "mlmodel_ave_Pmax"); // model read
        statusTracer.appendDoublePtr(&mlModel.sigma2_offset, 1, "mlmodel_sigma2_offset"); // model read
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
        statusTracer.appendFloatPtr(&sampler3d.random_perturbation,1,"sampling3d.random_perturbation",true); // sampling read
#else
        statusTracer.appendDoublePtr(&sampler3d.random_perturbation,1,"sampling3d.random_perturbation"); // sampling read
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
    //
void setupMLoptimizer()
{
    
#ifdef USEMPI
    // MPI::Init();//init mpi outside
    nodes = MPI::COMM_WORLD.Get_size();
    node = MPI::COMM_WORLD.Get_rank();
    MPI::Get_processor_name(nodeName,nodeNameLen);
    // std::cout<<nodeName<<" "<<node<<" "<<nodes<<std::flush<<std::endl;
#else
    nodes = 1;
    node = 0;
#endif
    
    if (continue_fn!="NULL") {
        if (continue_fn.find("_optimiser.star")!=std::string::npos)
            star_fn = continue_fn.substr(0,continue_fn.find("_optimiser"))+"_data.star";
        else if(continue_fn.find("_backup.back")!=std::string::npos)
            star_fn = continue_fn.substr(0,continue_fn.find("_backup"))+"_data.star";
        else
            ERROR_REPORT("Wrong continue file name "+continue_fn+",use *_optimiser.star or *_backup.back");
        iter = atoi( continue_fn.substr(continue_fn.find("_it")+3,continue_fn.find_last_of('_')-continue_fn.find("_it")-3).c_str() );
        write_fn = write_fn + "_ct" + std::to_string((long long)iter);
        NODE0ONLY std::cout<<"-Continue starting from iter "<<iter<<std::endl;
    }
    
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
    global_data_stream.foutInt(&sampler3d_healpix_order, 1, "sampling.healpix_order", __FILE__, __LINE__);
    global_data_stream.foutDouble(&sampler3d.psi_step, 1, "sampling.psi_step", __FILE__, __LINE__);
    global_data_stream.foutDouble(-91, "sampling.limit_tilt", __FILE__, __LINE__);
    global_data_stream.foutDouble(&offset_range, 1, "sampling.offset_range", __FILE__, __LINE__);
    global_data_stream.foutDouble(&offset_step, 1, "sampling.offset_step", __FILE__, __LINE__);
    global_data_stream.foutDouble(0.5, "sampling.perturbation_factor", __FILE__, __LINE__);
    global_data_stream.check();global_data_stream.flush();
#endif
    
    // ------------ initialize sampling ---------- //
    // set sampling 3d
    sampler3d.initialize(offset_step,offset_range,-1,sampler3d_healpix_order, sampler3d_fn_sym, sigma2_angle==0?NOPRIOR:PRIOR_ROTTILT_PSI);
    
    maxthreads = omp_get_max_threads();
    
    NODE0ONLY std::cout<<"maxthreads = "<<maxthreads<<std::endl;
    
    // omp_set_num_threads(8);
    if(0 == node) std::cout<<"number of threads = "<<maxthreads<<std::endl;
    
    metadata.readFromStar(star_fn);
    nr_global_images = metadata.numberOfParticles();
    
    Mrcs::MrcsHead mrcsHead;
    Mrcs::readMrcsHead(metadata[0].IMAGE.NAME, mrcsHead);
    NODE0ONLY std::cout<<"x = "<<mrcsHead.NC<<",y = "<<mrcsHead.NC<<",z = "<<mrcsHead.NS<<",mode = "<<mrcsHead.MODE<<std::endl;
    
    ori_size = mrcsHead.NC;
    
    // randomly shuffle metadata
    metadata.shuffle(random_seed);
    
    // NOTE : SAME_IMAGE_DIVISION uses to produce the same result for different nodes
#ifdef SAME_IMAGE_DIVISION
    //divided like not mpi version,just to easily debug
    nr_local_images = divide_equally_pool(nr_global_images,nr_pool,nodes,node, first_local_image, last_local_image);
#else
    nr_local_images = divide_equally(nr_global_images,nodes,node, first_local_image, last_local_image);
#endif
    
    // read the local images data
    images.resize(nr_local_images);
    for (auto& image : images) image.init(ori_size*ori_size);
    
    // All nodes read the file at once rather than broadcasting
    float *buffer = (float*)aMalloc(sizeof(float)*ori_size*ori_size,64);
    std::string preFileName = metadata[first_local_image].IMAGE.NAME;
    FILE* mrcsFile = fopen((metadata[first_local_image].IMAGE.NAME).c_str(),"rb");
    
    for (int iimage = 0; iimage < nr_local_images; iimage++) {
        
        if (metadata[iimage+first_local_image].IMAGE.NAME != preFileName)
        {
            fclose(mrcsFile);
            mrcsFile = fopen((metadata[iimage+first_local_image].IMAGE.NAME).c_str(),"rb");
            ERROR_CHECK(mrcsFile==NULL, "read file failed.");
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
        auto image_data = images[iimage].wptr(ori_size*ori_size);
        for (int i = 0; i < ori_size*ori_size; i++) {
            image_data[i] = buffer[i];
        }
    }
    
    aFree(buffer);
    
    fclose(mrcsFile);
    
    // ------------ initialize sampling ---------- //
    // Also randomize random-number-generator for perturbations on the angles
    dontShare_Random_generator.init(random_seed);
    
}

void prepare()
{
    // --------------- check pixel size ------------------ //
    if (metadata[0].MAGNIFICATION!=0 && metadata[0].DETECTOR_PIXEL_SIZE!=0) {
        double new_pixel_size = 10000. * metadata[0].DETECTOR_PIXEL_SIZE / metadata[0].MAGNIFICATION;
        if (fabs(new_pixel_size-pixel_size)>0.01) {
            NODE0ONLY std::cout<<"MODIFYING pixel size from "<<pixel_size<<" to "<<new_pixel_size
            			<<" based on magnification information in the input STAR file"<<std::endl;
            pixel_size = new_pixel_size;
        }
        for (int iimage = 1; iimage < nr_local_images; iimage++) {
            double my_pixel_size = 10000. * metadata[iimage].DETECTOR_PIXEL_SIZE / metadata[iimage].MAGNIFICATION;
            ERROR_CHECK(my_pixel_size != new_pixel_size, "ERROR inconsistent magnification and detector pixel sizes in images in input STAR file");
        }
    }
    // --------------- initialize MAP Model ---------------- //
    if (particle_diameter < 0.)
        particle_diameter = (ori_size - width_mask_edge) * pixel_size;
    
    // --------------- initialize Particle Model ----------------- //
    particleModel.initialize(ori_size, pixel_size, particle_diameter, width_mask_edge,
                             sigma2_fudge, random_seed, do_norm_correction, do_zero_mask,
                             do_shifts_onthefly,maxthreads,&global_data_stream);
    
    // ini_high ,user set
    mapModel.initialize(nr_classes, ori_size, particle_diameter, pixel_size, ini_high, maxthreads,
                        width_mask_edge, width_fmask_edge, 5, true, do_map, 3, sampler3d_fn_sym);
    
    if(iter==0) {
        mapModel.initializeRef(ref_fn);
    }
    
    // ------------ initialize model and wsum --------------- //
    mlModel.initialize(ori_size, nr_classes, metadata.numberOfGroups(), sampler3d.NrDir());
    
    if (mlModel.nr_groups == 1) do_scale_correction = false;
    NODE0ONLY std::cout<<"do_scale_correction = "<<do_scale_correction<<std::endl;
    //
    // continue
    setupStatusTracer();
    //
    if (continue_fn!="NULL") {
        readResult();
        // After the first iteration the references are always CTF-corrected
        if (do_ctf_correction)
            refs_are_ctf_corrected = true;
    }
    // ------------ initialize model ----------- //
    
    // Calculate initial sigma noise model from power_class spectra of the individual images
    Image Mavg;Mavg.init(ori_size*ori_size);Mavg.zero();
    // check the data whether normalized
    // checkNormalize(images_data, ori_size, nr_local_images, particle_diameter, pixel_size);
    
    // initialize Mavg,wsum_sigma2_noise,wsum_sumw_group,
    // set nr_particles_group
    mlModel.calculateSumOfPowerSpectraAndAverageImage(Mavg, images, do_zero_mask, metadata, first_local_image, mapModel);
    
#ifdef USEMPI
    // TODO : some diff after 3 iteration when use "--scale"
    // reduce data
    if (iter==0) {
        int local_temp_size = std::max(mlModel.nr_groups,ori_size*ori_size);
        FDOUBLE *local_temp = (FDOUBLE*)aMalloc(sizeof(FDOUBLE)*local_temp_size,64);
        memcpy(local_temp, Mavg.rptrAll(), sizeof(FDOUBLE)*ori_size*ori_size);
        MPI::COMM_WORLD.Reduce(local_temp,Mavg.wptrAll(),ori_size*ori_size,MPI_FDOUBLE,MPI::SUM,0);
        
        // NOTE : because we only get average sigma2_noise,do not need to reduce all nr_groups's sigma2 noise
        memcpy(local_temp, mlModel.wsum_sigma2_noise[0].rptrAll(), sizeof(FDOUBLE)*(ori_size/2+1));
        MPI::COMM_WORLD.Reduce(local_temp,mlModel.wsum_sigma2_noise[0].wptrAll(),ori_size/2+1,MPI_FDOUBLE,MPI::SUM,0);
        
        memcpy(local_temp, mlModel.wsum_sumw_group.rptrAll(), sizeof(FDOUBLE)*mlModel.nr_groups);
        MPI::COMM_WORLD.Reduce(local_temp,mlModel.wsum_sumw_group.wptrAll(),mlModel.nr_groups,MPI_FDOUBLE,MPI::SUM,0);
        
        aFree(local_temp);
        
        int* local_temp2 = (int*)aMalloc(sizeof(int)*mlModel.nr_groups,64);
        memcpy(local_temp2, mlModel.nr_particles_group.rptrAll(), sizeof(int)*mlModel.nr_groups);
        MPI::COMM_WORLD.Allreduce(local_temp2,mlModel.nr_particles_group.wptrAll(),mlModel.nr_groups,MPI::INT,MPI::SUM);
        
        aFree(local_temp2);
        if (false/*node == 1*/) {
            for (int i = 0; i < mlModel.nr_groups; i++) {
                std::cerr<<mlModel.nr_particles_group.rptrAll()[i]<<" ";
            }
            std::cerr<<std::endl;
        }
    }
#endif
    
#ifdef DATA_STREAM
    // TODO : check why different with relion.
    global_data_stream.foutInt(metadata.numberOfGroups(), "numberOfGroups", __FILE__, __LINE__);
    global_data_stream.foutInt(metadata.numberOfMicrographs(), "numberOfMicrographs", __FILE__, __LINE__);
    global_data_stream.foutInt(metadata.numberOfParticles(), "numberOfOriginalParticles", __FILE__, __LINE__);
    global_data_stream.foutDouble(mapModel.Irefs[0].wptr(), ori_size*ori_size*ori_size, "ref_1", __FILE__, __LINE__);
    global_data_stream.foutDouble(mapModel.Irefs[nr_classes-1].wptr(), ori_size*ori_size*ori_size, "ref_1", __FILE__, __LINE__);
    global_data_stream.foutInt(1, "mymodel.orientational_prior_mode", __FILE__, __LINE__);
    global_data_stream.foutInt(mapModel.ref_dim, "mymodel.ref_dim", __FILE__, __LINE__);
    global_data_stream.foutInt(sampler3d.NrDir(), "sampling.NrDirections()", __FILE__, __LINE__);
    global_data_stream.foutDouble(mlModel.pdf_direction[0].wptr(sampler3d.NrDir()), sampler3d.NrDir(), "mymodel.pdf_direction[0]", __FILE__, __LINE__);
    global_data_stream.foutDouble(mlModel.pdf_direction[nr_classes-1].wptr(sampler3d.NrDir()), sampler3d.NrDir(), "mymodel.pdf_direction_nr_class", __FILE__, __LINE__);
    if (iter==0) {
        global_data_stream.foutDouble(Mavg.wptrAll(), ori_size*ori_size, "Mavg", __FILE__, __LINE__);
        global_data_stream.foutDouble(mlModel.wsum_sigma2_noise[0].wptr(ori_size/2+1), ori_size/2+1, "wsum_model_sigma2_noise", __FILE__, __LINE__);
        global_data_stream.foutDouble(mlModel.wsum_sumw_group.wptrAll()[0], "wsum_model_sumw_group1", __FILE__, __LINE__);
        global_data_stream.foutDouble(mlModel.wsum_sumw_group.wptrAll()[mlModel.nr_groups-1], "wsum_model_sumw_groupN", __FILE__, __LINE__);
        global_data_stream.foutInt(mlModel.nr_particles_group.wptrAll()[0], "mymodel_sumw_nr_particles_group1", __FILE__, __LINE__);
        global_data_stream.foutInt(mlModel.nr_particles_group.wptrAll()[mlModel.nr_groups-1], "mymodel_sumw_nr_particles_groupN", __FILE__, __LINE__);
    }
    global_data_stream.check();global_data_stream.flush();
#endif
    //
    if (iter==0) {
        for (int i = 0; i < ori_size*ori_size; i++)
            Mavg[i] /= mlModel.wsum_sumw_group.rptrAll()[0];
    }

    // Set model_sigma2_noise and model_Iref from averaged poser spectra and Mavg
    NODE0ONLY {
        if (iter==0) {
            mlModel.setSigmaNoiseEstimates(Mavg);
        }
    }
    
#ifdef DATA_STREAM
    if (iter==0) {
        global_data_stream.foutDouble(mlModel.sigma2_noise[0].wptr(ori_size/2+1), ori_size/2+1, "sigma2_noise1", __FILE__, __LINE__);
        global_data_stream.foutDouble(mlModel.sigma2_noise[mlModel.nr_groups-1].wptr(ori_size/2+1), ori_size/2+1, "sigma2_noiseN", __FILE__, __LINE__);
        global_data_stream.foutInt(ori_size/2+1, "sigma2_size", __FILE__, __LINE__);
        global_data_stream.check();global_data_stream.flush();
    }
#endif
    
    // First low-pass filter the initial references
    if(iter==0){
        mapModel.applyLowPassFilter();
    }
    
#ifdef DATA_STREAM
    global_data_stream.foutDouble(mapModel.Irefs[0].wptr(), ori_size*ori_size*ori_size, "ref_1", __FILE__, __LINE__);
    global_data_stream.foutDouble(mapModel.Irefs[nr_classes-1].wptr(), ori_size*ori_size*ori_size, "ref_N", __FILE__, __LINE__);
    global_data_stream.check();global_data_stream.flush();
#endif
    
    // Initialise the model_data_versus_prior ratio to get the initial current_size right
    // model_tau2_class,model_data_vs_prior_class
    NODE0ONLY {
        if (iter==0) {
            mlModel.initialiseDataVersusPrior(mapModel,tau2_fudge_factor); // fix_tau was set in initialiseGeneral
        }
    }
    
#ifdef DATA_STREAM
    global_data_stream.foutDouble(mlModel.data_vs_prior_class[0].wptr(ori_size/2+1), ori_size/2+1, "data_vs_prior_class_1", __FILE__, __LINE__);
    global_data_stream.foutDouble(mlModel.data_vs_prior_class[nr_classes-1].wptr(ori_size/2+1), ori_size/2+1, "data_vs_prior_class_N", __FILE__, __LINE__);
    global_data_stream.check();global_data_stream.flush();
#endif
    
#ifdef USEMPI
    for (int igroup = 0; igroup < mlModel.nr_groups; igroup++)
        MPI::COMM_WORLD.Bcast(mlModel.sigma2_noise[igroup].wptrAll(),mlModel.ori_Fsize,MPI_FDOUBLE,0);
    
    for (int iclass = 0; iclass < nr_classes; iclass++){
        // in 3D case,only read reference from mrc file
        // MPI::COMM_WORLD.Bcast(mapModel.Irefs[iclass].wptr(),mapModel.Irefs[iclass].dimzyx,MPI::DOUBLE,0);
        MPI::COMM_WORLD.Bcast(mlModel.tau2_class[iclass].wptrAll(),mlModel.ori_Fsize,MPI_FDOUBLE,0);
        MPI::COMM_WORLD.Bcast(mlModel.data_vs_prior_class[iclass].wptrAll(),mlModel.ori_Fsize,MPI_FDOUBLE,0);
    }
#endif
    
    Mavg.fini();
}

void destroyMLoptimizer()
{
    images			.resize(0);
    mapModel		.finalize();
    mlModel			.finalize();
    particleModel	.finalize();
}

// ------------------------- EM-Iteration  ------------------------- //

void iterate()
{
    // Update the current resolution and image sizes, and precalculate resolution pointers
    // The rest of the time this will be done after maximization and before writing output files,
    // so that current resolution is in the output files of the current iteration
    bool set_by_ini_high = ini_high > 0. && (iter == 0 || (iter == 1 && do_firstiter_cc) );
    static int nr_iter_wo_resol_gain = 0;
    mapModel.updateCurrentResolution(mlModel.data_vs_prior_class,set_by_ini_high,nr_iter_wo_resol_gain);
    
#ifdef DATA_STREAM
    global_data_stream.foutDouble(mapModel.current_resolution, "updateCurrentResolution()_current_resolution", __FILE__, __LINE__);
#endif
    
    bool has_already_reached_convergence = false;
    for (iter = iter + 1; iter <= nr_iter; iter++)
    {
        NODE0ONLY std::cout<<"current_resolution = "<<mapModel.current_resolution<<std::endl;
        
#ifdef DATA_STREAM
        global_data_stream.foutInt(iter, "iterate()_iter", __FILE__, __LINE__);
#endif
        
        const double starttime = dtime();
        // update coarse_size,current_size,Npix_per_shell,Mresol_coarse,Mresol_fine
        double angularSampler = sampler3d.getAngularSampling();
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
        
        //#define DO_RECONSTRUCT_EACH_NODE
#ifdef DO_RECONSTRUCT_EACH_NODE
        mapModel.reduceData(MPI::COMM_WORLD,false);
#else
        mapModel.reduceData(MPI::COMM_WORLD,true);
#endif
        gatherMetaDataToMaster(metadata);
        
        MPI::COMM_WORLD.Barrier();
        
#endif
        
        maximization();
        
        // Apply masks to the reference images
        // At the last iteration, do not mask the map for validation purposes
        if(do_solvent && !has_converged)
            mapModel.applySolventFlatten(mask_fn);
        
#ifdef DATA_STREAM
        global_data_stream.foutDouble(mapModel.Irefs[0].wptr(), mapModel.Irefs[0].dimzyx, "iterate()_Iref[0]", __FILE__, __LINE__);
        global_data_stream.foutDouble(mapModel.Irefs[nr_classes-1].wptr(), mapModel.Irefs[nr_classes-1].dimzyx, "iterate()_Iref[nr_classes-1]", __FILE__, __LINE__);
#endif
        
        // Re-calculate the current resolution, do this before writing to get the correct values in the output files
        bool set_by_ini_high = ini_high > 0. && (iter == 0 || (iter == 1 && do_firstiter_cc) );
        mapModel.updateCurrentResolution(mlModel.data_vs_prior_class,set_by_ini_high,nr_iter_wo_resol_gain);
        mapModel.printResolution(mlModel.data_vs_prior_class,set_by_ini_high);
        
#ifdef DATA_STREAM
        global_data_stream.foutDouble(mapModel.current_resolution, "updateCurrentResolution()_current_resolution", __FILE__, __LINE__);
#endif
        
        // Write output files
        writeAndCheckResult();
        
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
    
    NODE0ONLY std::cout<<"setFourierTransformMaps,"<<__LINE__<<" , "<<__FILE__<<std::endl;
    NODE0ONLY std::cout<<"Create temporary threads data for 3D Volume,"<<std::endl;
    int pad_size = 2*(2*std::min(current_size/2, ori_size/2)+1)+1;
    int pad_Fsize = pad_size/2+1;
    NODE0ONLY std::cout<<"Memory needed for threads 3D Volume data ~= "<<maxthreads*nr_classes*(pad_size*pad_size*pad_Fsize/1024./1024./1024.)<<" GB."<<std::endl;
    // Waiting to complete computeFourierTransformMap in Reference class
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
    sampler3d.resetRandomlyPerturbedSampling();
    
    NODE0ONLY {
        // check error
        int n_trials_acc = 10;
        n_trials_acc = std::min(n_trials_acc, (int)nr_local_images);
        // calculateExpectedAngularErrors(0, n_trials_acc-1);
    }
    // the first time we use small nr_pool to initialize model_Iref,
    NODE0ONLY printMem(nr_pool);
    
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
    exp_wsum_scale_correction_XA.init(nr_images, ori_size/2+1);
    exp_wsum_scale_correction_AA.init(nr_images, ori_size/2+1);
    exp_imgs					.init(nr_images, ori_size*ori_size);
    exp_power_imgs				.init(nr_images, ori_size/2+1);
    //
	exp_metadata				.init(nr_images);
    
    // initialize sampling
    samplingGrid.initialize(sampler3d, adaptive_oversampling);
    exp_nr_trans = samplingGrid.exp_nr_trans();
    exp_nr_psi = samplingGrid.exp_nr_psi();
    exp_nr_dir = samplingGrid.exp_nr_dir();
    int max_nr_over_rot = sampler3d.oversamplingFactorOrientations(adaptive_oversampling);
    int max_nr_over_trans = sampler3d.oversamplingFactorTranslations(adaptive_oversampling);
    
    // intialzie some shift images and its ctf
    int exp_nr_all_trans    = sampler3d.NrTrans(adaptive_oversampling);
    exp_Fimgs_shifted_real	.init(nr_images, exp_nr_all_trans, current_Fsize2);
    exp_Fimgs_shifted_imag	.init(nr_images, exp_nr_all_trans, current_Fsize2);
    
    exp_Rot_significant		.init(nr_classes, exp_nr_dir*exp_nr_psi);
    exp_Mcoarse_significant	.init(nr_images, nr_classes*sampler3d.NrPoints(0));
    // base on the density of this matrix
    // for high-density case(in coarse easrch step or previous iterations ) use c-array
    // for low-density case(in fine search step or later iterations ) use stack data-structure
    // for low-density case,also need to improve the nr_pool to increase the image-processing each time
    exp_Mweight_coarse		.init(	nr_images, nr_classes, exp_nr_dir, exp_nr_psi, exp_nr_trans);
    exp_Mweight_fine		.init(	nr_images, nr_classes, exp_nr_dir, exp_nr_psi,
                              		max_nr_over_rot*exp_nr_trans*max_nr_over_trans);
    // -------------   thread data     --------------- //
    thread_exp_max_weight_index.resize(maxthreads);
    for (auto &exp_max_weight_index : thread_exp_max_weight_index)
        exp_max_weight_index.resize(nr_images);
    thread_Frefctf_real				.init(maxthreads, exp_nr_psi, current_Fsize2);
    thread_Frefctf_real				.fill_with_first_touch(0.);
    thread_Frefctf_imag				.init(maxthreads, exp_nr_psi, current_Fsize2);
    thread_Frefctf_imag				.fill_with_first_touch(0.);
    thread_Fimg_real				.init(maxthreads, max_nr_over_trans, ori_size*ori_size);
    thread_Fimg_real				.fill_with_first_touch(0.);
    thread_Fimg_imag				.init(maxthreads, max_nr_over_trans, ori_size*ori_size);
    thread_Fimg_imag				.fill_with_first_touch(0.);
    thread_Fweight					.init(maxthreads, exp_nr_psi, ori_size*ori_size);
    thread_Fweight					.fill_with_first_touch(0.);
    // ------- //
    thread_exp_sum_weight			.init(maxthreads, nr_images);
    thread_exp_min_diff2			.init(maxthreads, nr_images);
    thread_exp_max_weight			.init(maxthreads, nr_images);
    thread_wsum_norm_correction		.init(maxthreads, nr_images);
    thread_wsum_pdf_class			.init(maxthreads, nr_classes);
    thread_sumw_group				.init(maxthreads, nr_images);
    threadfake_do_scale_norm_class	.init(nr_classes, current_Fsize2);
    thread_wsum_sigma2_offset		.init(maxthreads, 1);
    thread_wsum_pdf_direction		.init(maxthreads, nr_classes, exp_nr_dir);
    thread_wsum_sigma2_noise		.init(maxthreads, nr_images, current_Fsize2);
    thread_wsum_scale_correction_XA	.init(maxthreads, nr_images, current_Fsize2);
    thread_wsum_scale_correction_AA	.init(maxthreads, nr_images, current_Fsize2);
    
    //
    assert(do_shifts_onthefly==false);
    particleModel.setup(nr_images, current_size, coarse_size, 0 , 0);
    
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
    exp_wsum_scale_correction_XA	.fini();
    exp_wsum_scale_correction_AA	.fini();
    exp_imgs						.fini();
    exp_power_imgs					.fini();
    //
    exp_metadata					.fini();
    //
    exp_Fimgs_shifted_real			.fini();
    exp_Fimgs_shifted_imag			.fini();
	//
    exp_Rot_significant				.fini();
    exp_Mcoarse_significant			.fini();
    exp_Mweight_coarse				.fini();
    exp_Mweight_fine				.fini();
    //
    samplingGrid					.finalize();
    // some thread variable
    thread_exp_max_weight_index		.resize(0);
    //
    thread_Frefctf_real				.fini();
    thread_Frefctf_imag				.fini();
    thread_Fimg_real				.fini();
    thread_Fimg_imag				.fini();
    thread_Fweight					.fini();
    //
    thread_exp_sum_weight			.fini();
    thread_exp_min_diff2			.fini();
    thread_exp_max_weight			.fini();
    thread_wsum_norm_correction		.fini();
    thread_wsum_pdf_class			.fini();
    thread_sumw_group				.fini();
    threadfake_do_scale_norm_class	.fini();
    thread_wsum_sigma2_offset		.fini();
    //
    thread_wsum_pdf_direction		.fini();
    thread_wsum_sigma2_noise		.fini();
    thread_wsum_scale_correction_XA	.fini();
    thread_wsum_scale_correction_AA	.fini();
    //
    particleModel					.destroy();
}

void expectation()
{
    
    NODE0ONLY std::cerr << " Expectation iteration " << iter<< " of " << nr_iter<<std::endl;
    
    int my_first_image,my_last_image,nr_images_done = 0;
    
    NODE0ONLY showProgressBar(0, nr_local_images);
    
    while (nr_images_done < nr_local_images)
    {
        my_first_image = first_local_image + nr_images_done;
        
        if ((!do_firstiter_cc && iter == 1) || (do_firstiter_cc && iter == 2)){
            // equally divided like nr_pool=1 case
            int iclass = divide_equally_which_group(nr_global_images, nr_classes, my_first_image);
            int first,last;
            divide_equally(nr_global_images, nr_classes, iclass, first, last);
            int suitable_pool = std::min(nr_pool, last-my_first_image+1);
            my_last_image = std::min(first_local_image+nr_local_images - 1, my_first_image + suitable_pool - 1);
        }
        else
            my_last_image = std::min(first_local_image+nr_local_images - 1, my_first_image + nr_pool - 1);
        
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
        // low-pass filter again and generate the seeds
        if (true/*do_generate_seeds*/)
        {
            if (do_firstiter_cc && iter == 1)
            {
                // In first (CC) iter, use a single reference (and CC)
                exp_iclass_min = exp_iclass_max = 0;
            }
            else if ( (do_firstiter_cc && iter == 2) || (!do_firstiter_cc && iter == 1))
            {
                // In second CC iter, or first iter without CC: generate the seeds
                // Now select a single random class
                // exp_part_id is already in randomized order (controlled by -seed)
                // WARNING: USING SAME iclass_min AND iclass_max FOR SomeParticles!!
                exp_iclass_min = exp_iclass_max = divide_equally_which_group(nr_global_images, nr_classes, exp_first_image);
            }
        }
        
#ifdef DATA_STREAM
        global_data_stream.foutInt(do_firstiter_cc, "expectationOneParticle()_do_firstiter_cc", __FILE__, __LINE__);
        global_data_stream.foutInt(exp_iclass_min, "expectationOneParticle()_exp_iclass_min", __FILE__, __LINE__);
        global_data_stream.foutInt(exp_iclass_max, "expectationOneParticle()_exp_iclass_max", __FILE__, __LINE__);
        global_data_stream.check();global_data_stream.flush();
#endif
        
        // prepare exp_image_data and exp_metadata
        // each node keep nr_local_images's images data and all(nr_global_images) metadata
        for (int iimage = 0; iimage < exp_nr_images; iimage++) {
            exp_metadata[iimage] = metadata[iimage+my_first_image];
            ::copy(exp_imgs[iimage].wptrAll(), ori_size*ori_size, images[my_first_image+iimage-first_local_image].rptrAll(), ori_size*ori_size);
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
        hiddenVarMonitor.monitorHiddenVariableChanges(sampler3d, metadata, my_first_image, &exp_metadata[0], exp_nr_images);
        
        // update the metadata
        for (int iimage = 0; iimage < exp_nr_images; iimage++) metadata[iimage+my_first_image] = exp_metadata[iimage];
        
        //
        nr_images_done += my_last_image - my_first_image + 1;
        
        NODE0ONLY showProgressBar(nr_images_done, nr_local_images);
    }
    NODE0ONLY std::cerr<<std::endl;
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
        bool do_coarse_search = (exp_ipass == 0);
        // Use smaller images in the first pass, larger ones in the second pass
        exp_current_size = (do_coarse_search && adaptive_oversampling > 0) ? coarse_size : current_size;
        
        // Use coarse sampling in the first pass, oversampled one the second pass
        // and initialize sampling data
        exp_current_oversampling = do_coarse_search ? 0 : adaptive_oversampling;
        samplingGrid.computeGrid3D(sampler3d, exp_current_oversampling);
        exp_nr_over_rot = samplingGrid.exp_nr_over_rot();
        exp_nr_over_trans = samplingGrid.exp_nr_over_trans();
		// particleModel.shiftAssistorsSetCurrDims(exp_nr_trans, exp_nr_over_trans);
        //
        bool do_cross_correlation = (iter==1&&do_firstiter_cc)|| do_always_cc;
        particleModel.preShiftedImagesCtfsAndInvSigma2s(exp_Fimgs_shifted_real, exp_Fimgs_shifted_imag,
                                                        exp_local_Minvsigma2s, exp_local_sqrtXi2, exp_local_Fctfs,
                                                        exp_current_size, do_cross_correlation, do_coarse_search,
                                                        samplingGrid, mlModel, Mresol_coarse, Mresol_fine);
        
        // get all rotated reference  and the significant rotation
        getReferenceAllOrientations();
        
        // get all reference and images 's squared differences
        getAllSquaredDifferences(do_coarse_search);
        
        // convert all squared differences to weight,and find significant(maximum) weight for each image
        if (do_coarse_search){
            convertSquaredDifferencesToWeights(exp_Mweight_coarse);
            findAllSignificantPoints(exp_Mweight_coarse);
        }
        else{
            convertSquaredDifferencesToWeights(exp_Mweight_fine);
            findAllSignificantPoints(exp_Mweight_fine);
        }
        
        NODE0ONLY
        {
            std::string iterStr = num2str(iter);
            if (do_coarse_search) {
                std::string fn_exp_mweight = write_path+write_fn+"_it"+iterStr+"_exp_Mweight_coarse";
                std::string note = "coarse_search_"+std::to_string((long long)exp_first_image);
                exp_Mweight_coarse.analysis(fn_exp_mweight,note);
            }
            else {
                std::string fn_exp_mweight = write_path+write_fn+"_it"+iterStr+"_exp_Mweight_fine";
                std::string note = "fine_search_"+std::to_string((long long)exp_first_image);
                exp_Mweight_fine.analysis(fn_exp_mweight,note);
            }
        }
    }// end loop over 2 exp_ipass iterations
    
    // update some parameter for maximization
    updateOtherParams();
    
    //
    backProjection();
    
    //
    storeWeightedSums();
}


void getReferenceAllOrientations()
{
    
    int exp_current_Fsize2 = exp_current_size*(exp_current_size/2+1);
    
#pragma omp parallel for collapse(3)
    for (int iclass = exp_iclass_min; iclass <= exp_iclass_max; iclass++)
    {
        for (int idir = 0; idir < exp_nr_dir; idir++)
        {
            for (int ipsi = 0; ipsi < exp_nr_psi; ipsi++)
            {
                if (mlModel.pdf_class.rptrAll()[iclass] <= 0.)
                    continue;
                int iorient = idir*exp_nr_psi+ipsi;
                int iorientclass = iclass * exp_nr_dir * exp_nr_psi + iorient;
                assert(iorientclass>=0);
                // In the first pass, always proceed
                // In the second pass, check whether one of the translations for this orientation of any of
                // the particles had a significant weight in the first pass
                // if so, proceed with projecting the reference in that direction
                if (exp_ipass == 0) {
                    exp_Rot_significant[iclass].wptrAll()[iorient] = true;
                }
                else
                    exp_Rot_significant[iclass].wptrAll()[iorient] = isSignificantAnyParticleAnyTranslation(iorientclass);
            }// end loop of ipsi
        }// end loop of idir
    }// end loop of iclass
}

// TODO,better one
inline void decodeOverRotTrans(int ihidden,int& iover_rot,int& itrans,int& iover_trans)
{
    // decode the real position,psi->over_rot->trans->over_trans
    int denominator = exp_nr_over_rot*exp_nr_over_trans;
    itrans = ihidden / denominator;
    
    ihidden = ihidden % denominator;
    denominator = exp_nr_over_trans;
    iover_rot = ihidden / denominator;
    iover_trans = ihidden % denominator;
}

static void check_scale(FDOUBLE& myscale,int& igroup){
    static bool have_warned_small_scale = false;
    ERROR_CHECK(std::isnan(myscale) || std::isinf(myscale), "Scale is Nan or Inf.");
    if (myscale > 10000.)
    {
        std::cerr << " rlnMicrographScaleCorrection= " << myscale << " group= " << igroup + 1 << std::endl;
        ERROR_REPORT("ERROR: rlnMicrographScaleCorrection is very high. Did you normalize your data?");
    }
    else if (myscale < 0.001)
    {
        if (!have_warned_small_scale)
        {
            std::cout << " WARNING: ignoring group " << igroup + 1 << " with very small or negative scale (" << myscale <<
            "); Use larger groups for more stable scale estimates." << std::endl;
            have_warned_small_scale = true;
        }
        myscale = 0.001;
    }
};

namespace CoarseAndFineSearch
{
    int exp_current_Fsize2;
    int nr_image_block;
    int nr_psi_block;
    
    FDOUBLE COMPUTE_FLAG = (std::numeric_limits<FDOUBLE>::max)();
    
    //
    inline void get_projected_class2(int iclass,FDOUBLE over_rot,FDOUBLE over_tilt,FDOUBLE over_psi,
                                    FDOUBLE* &Frefctf_real,FDOUBLE* &Frefctf_imag,int n_start,int n_end)
    {
        FDOUBLE A[3][3];
        Euler_angles2matrix(over_rot, over_tilt, over_psi, A);
        mapModel.get2DFourierTransformOneTile(iclass, Frefctf_real, Frefctf_imag,
                                              n_start, n_end, exp_current_size, A, false);
        
        // mapModel.get2DFourierTransform(iclass, Frefctf_real, Frefctf_imag,exp_current_size, A, false);
        
    }
    //
    inline void applyCtfToClass(int iimage,FDOUBLE* &Frefctf_real,FDOUBLE* &Frefctf_imag,int n_start,int n_end)
    {
        // Apply CTF to reference projection
        // after first iteration refs_are_ctf_corrected = true
        if (do_ctf_correction && refs_are_ctf_corrected)
        {
            auto exp_local_Fctfs_iimage = exp_local_Fctfs[iimage].rptr(exp_current_Fsize2);
#pragma vector aligned
#pragma ivdep
            for (int n = n_start; n < n_end; n++) {
                Frefctf_real[n] = Frefctf_real[n]*exp_local_Fctfs_iimage[n];
                Frefctf_imag[n] = Frefctf_imag[n]*exp_local_Fctfs_iimage[n];
            }
        }
        if (do_scale_correction)
        {
            int igroup = exp_metadata[iimage].GROUP_NO-1;
            auto myscale = mlModel.scale_correction.rptrAll()[igroup];
            check_scale(myscale,igroup);
#pragma vector aligned
#pragma ivdep
            for (int n = n_start; n < n_end; n++)
            {
                Frefctf_real[n] = Frefctf_real[n]*myscale;
                Frefctf_imag[n] = Frefctf_imag[n]*myscale;
            }
        }
    }
    //
    inline void get_shifted_image(int tid,int iimage,int itrans,int iover_trans,
                                  FDOUBLE* &Fimg_shift_real,FDOUBLE* &Fimg_shift_imag,int n_start,int n_end)
    {
        // Get the shifted image from memory
        int itrans_over_trans = itrans * exp_nr_over_trans + iover_trans;
        Fimg_shift_real = exp_Fimgs_shifted_real.wptr(iimage, itrans_over_trans, exp_current_Fsize2);
        Fimg_shift_imag = exp_Fimgs_shifted_imag.wptr(iimage, itrans_over_trans, exp_current_Fsize2);
    }
    //
    inline void get_diff(int iimage,double& diff2,
                         const FDOUBLE* Frefctf_real,const FDOUBLE* Frefctf_imag,
                         const FDOUBLE* Fimg_shift_real,const FDOUBLE* Fimg_shift_imag,
                         const FDOUBLE* Minvsigma2_iimage,
                         int n_start,int n_end)
    {
        // do cross-correlation
        if ((iter == 1 && do_firstiter_cc) || do_always_cc)
        {
            // Do not calculate squared-differences, but signal product
            // Negative values because smaller is worse in this case
            diff2 = 0.;
            double suma2 = 0.;
#pragma vector aligned
#pragma ivdep
            for(int n = n_start;n < n_end;n++)
            {
                diff2 -= Frefctf_real[n] * Fimg_shift_real[n];
                diff2 -= Frefctf_imag[n] * Fimg_shift_imag[n];
                suma2 += Frefctf_real[n] * Frefctf_real[n] + Frefctf_imag[n] * Frefctf_imag[n];
            }
            // Normalised cross-correlation coefficient: divide by power of reference (power of image is a constant)
            diff2 /= sqrt(suma2) * exp_local_sqrtXi2.rptrAll()[iimage];
            
        }
        else // do maximum-likelihood
        {
            if(n_start==0) diff2 = 0;
#pragma vector aligned
#pragma ivdep
            for(int n = n_start;n < n_end;n++)
            {
                auto diff_real = Frefctf_real[n] - Fimg_shift_real[n];
                auto diff_imag = Frefctf_imag[n] - Fimg_shift_imag[n];
                
                diff2 += (diff_real * diff_real + diff_imag * diff_imag) * Minvsigma2_iimage[n];
            }
            
            // Calculate the actual squared difference term of the Gaussian probability function
            // If current_size < model ori_size diff2 is initialised to the sum of
            // all |Xij|2 terms that lie between current_size and ori_size
            // Factor two because of factor 2 in division below, NOT because of 2-dimensionality of the complex plane!
            if(n_end==exp_current_Fsize2) diff2 = (diff2 + exp_highres_Xi2_imgs.rptrAll()[iimage]) / 2;
        }
    }
    //
    void global_data_stream_before()
    {
#ifdef DATA_STREAM
        int tid = 0;
        for (int iimage = 0; iimage < exp_nr_images; iimage++){
            auto Minvsigma2_iimage 	= exp_local_Minvsigma2s[iimage].wptr(exp_current_Fsize2);
            for (int iclass = exp_iclass_min; iclass <= exp_iclass_max; iclass++){
                for (int idir = 0; idir < exp_nr_dir; idir++){
                    if (mlModel.pdf_direction[iclass].rptrAll()[idir] <= 0. ||
                        mlModel.pdf_class.rptrAll()[iclass] <= 0.) continue;
                    for (int ipsi = 0; ipsi < exp_nr_psi; ipsi++){
                        // exp_iclass loop does not always go from 0 to nr_classes!
                        int iorient = idir*exp_nr_psi+ipsi;
                        if (!exp_Rot_significant[iclass].rptrAll()[iorient]) continue;
                        FDOUBLE *exp_over_rot,*exp_over_tilt,*exp_over_psi;
                        samplingGrid.setOrientation(sampler3d, exp_current_oversampling, idir, ipsi,
                                                    exp_over_rot, exp_over_tilt, exp_over_psi);
                        for (int iover_rot = 0; iover_rot < exp_nr_over_rot; iover_rot++)
                        {
                            auto Frefctf_real 		= thread_Frefctf_real.wptr(tid, ipsi, exp_current_Fsize2);
                            auto Frefctf_imag 		= thread_Frefctf_imag.wptr(tid, ipsi, exp_current_Fsize2);
                            // get projected image
                            get_projected_class2(iclass,exp_over_rot[iover_rot],exp_over_tilt[iover_rot],exp_over_psi[iover_rot],
                                                 Frefctf_real,Frefctf_imag,0,exp_current_Fsize2);
                            // write out data stream
                            global_data_stream.foutInt(iclass, "getAllSquaredDifferences()_iclass", __FILE__, __LINE__);
                            global_data_stream.foutInt(iimage, "getAllSquaredDifferences()_iimage", __FILE__, __LINE__);
                            global_data_stream.foutInt(idir, "getAllSquaredDifferences()_idir", __FILE__, __LINE__);
                            global_data_stream.foutInt(ipsi, "getAllSquaredDifferences()_ipsi", __FILE__, __LINE__);
                            global_data_stream.foutInt(iover_rot, "getAllSquaredDifferences()_iover_rot", __FILE__, __LINE__);
                            global_data_stream.foutDouble(mlModel.pdf_class.rptrAll()[iclass], "getAllSquaredDifferences()_pdf_class[iclass]", __FILE__, __LINE__);
                            global_data_stream.foutDouble(mlModel.pdf_direction[iclass][idir], "getAllSquaredDifferences()_pdf_direction[iclass][idir]", __FILE__, __LINE__);
                            global_data_stream.foutDouble(exp_over_rot[iover_rot], "getAllSquaredDifferences()_exp_over_rot[iorient][iover_rot]", __FILE__, __LINE__);
                            global_data_stream.foutDouble(exp_over_tilt[iover_rot], "getAllSquaredDifferences()_exp_over_tilt[iorient][iover_rot]", __FILE__, __LINE__);
                            global_data_stream.foutDouble(exp_over_psi[iover_rot], "getAllSquaredDifferences()_exp_over_psi[iorient][iover_rot]", __FILE__, __LINE__);
                            global_data_stream.foutInt(do_firstiter_cc, "getAllSquaredDifferences()_do_firstiter_cc", __FILE__, __LINE__);
                            global_data_stream.foutInt(do_always_cc, "getAllSquaredDifferences()_do_always_cc", __FILE__, __LINE__);
                            global_data_stream.foutInt(do_scale_correction, "getAllSquaredDifferences()_do_scale_correction", __FILE__, __LINE__);
                            global_data_stream.foutInt(do_ctf_correction, "getAllSquaredDifferences()_do_ctf_correction", __FILE__, __LINE__);
                            global_data_stream.foutInt(refs_are_ctf_corrected, "getAllSquaredDifferences()_refs_are_ctf_corrected", __FILE__, __LINE__);
                            global_data_stream.foutInt(exp_metadata[iimage].GROUP_NO-1, "getAllSquaredDifferences()_GROUP_NO", __FILE__, __LINE__);
                            global_data_stream.foutDouble(mlModel.scale_correction.wptrAll()[exp_metadata[iimage].GROUP_NO-1], "getAllSquaredDifferences()_scale_correction", __FILE__, __LINE__);
                            global_data_stream.foutDouble(Frefctf_real, exp_current_Fsize2, "getAllSquaredDifferences()_Fref_real", __FILE__, __LINE__);
                            global_data_stream.foutDouble(Frefctf_imag, exp_current_Fsize2, "getAllSquaredDifferences_Fref_imag", __FILE__, __LINE__);
                            global_data_stream.foutDouble(Minvsigma2_iimage, exp_current_Fsize2, "getAllSquaredDifferences_Fref_imag()_Minvsigma2", __FILE__, __LINE__);
                            global_data_stream.check();global_data_stream.flush();
                        }
                    }
                }
            }
        }
#endif
    }
    //
    void global_data_stream_after()
    {
#ifdef DATA_STREAM
        for (int iimage = 0; iimage < exp_nr_images; iimage++)
        {
            if ((iter == 1 && do_firstiter_cc) || do_always_cc)
                global_data_stream.foutDouble(exp_local_sqrtXi2.wptrAll()[iimage], "getAllSquaredDifferences()_exp_local_sqrtXi2", __FILE__, __LINE__);
            else
                global_data_stream.foutDouble(exp_highres_Xi2_imgs.wptrAll()[iimage], "getAllSquaredDifferences()_exp_highres_Xi2_imgs", __FILE__, __LINE__);
        }
        // NOTE,exp_Mweight is checked in convertSquaredDifferencesToWeights()....
        global_data_stream.foutDouble(exp_min_diff2.wptrAll()[0], "getAllSquaredDifferences()_exp_min_diff2[0]", __FILE__, __LINE__);
        global_data_stream.check();global_data_stream.flush();
#endif
    }
    //
    void getAllSquaredDifferencesCoarse()
    {
#pragma omp parallel for collapse(3) schedule(dynamic)
        for (int iimage = 0; iimage < exp_nr_images; iimage++) // ~= 8(nr_pool)
        {
            for (int iclass = exp_iclass_min; iclass <= exp_iclass_max; iclass++) // 4~10
            {
                for (int idir = 0; idir < exp_nr_dir; idir++) // 192(healpix_order=1)~768(healpix_order=2)
                {
                    //
                    if (mlModel.pdf_direction[iclass].rptrAll()[idir] <= 0. ||
                        mlModel.pdf_class.rptrAll()[iclass] <= 0.) continue;
                    int tid = omp_get_thread_num();
                    int iover_rot = 0;
                    int iover_trans = 0;
                    // thread data
                    FDOUBLE *Fimg_shift_real,*Fimg_shift_imag;
                    auto Minvsigma2_iimage 	= exp_local_Minvsigma2s[iimage].wptr(exp_current_Fsize2);
                    auto thread_exp_min_diff2_tid = thread_exp_min_diff2[tid].wptr(exp_nr_images);
                    
                    for (int ipsi = 0; ipsi < exp_nr_psi; ipsi++)
                    {
                        FDOUBLE *exp_over_rot,*exp_over_tilt,*exp_over_psi;
                        samplingGrid.setOrientation(sampler3d, exp_current_oversampling, idir, ipsi,
                                                    exp_over_rot, exp_over_tilt, exp_over_psi);
                        auto Frefctf_real 	= thread_Frefctf_real.wptr(tid, ipsi, exp_current_Fsize2);
                        auto Frefctf_imag 	= thread_Frefctf_imag.wptr(tid, ipsi, exp_current_Fsize2);
                        // get projected image
                        get_projected_class2(iclass,exp_over_rot[0],exp_over_tilt[0],exp_over_psi[0],Frefctf_real,Frefctf_imag,0,exp_current_Fsize2);
                        //
                        applyCtfToClass(iimage, Frefctf_real, Frefctf_imag,0,exp_current_Fsize2);

                    }
                    
                    for (int ipsi = 0; ipsi < exp_nr_psi; ipsi++) // 24(healpix_order=1)~48(healpix_order=2)
                    {
                        auto Frefctf_real 		= thread_Frefctf_real.wptr(tid, ipsi, exp_current_Fsize2);
                        auto Frefctf_imag 		= thread_Frefctf_imag.wptr(tid, ipsi, exp_current_Fsize2);
                        auto& exp_Mweight_coarse_sub = exp_Mweight_coarse.wptr_sparse(iimage, iclass, idir, ipsi);
                        
                        for (int itrans = 0; itrans < exp_nr_trans; itrans++)// 100(offset_range=10,step=2) or 400(offset_range=10,step=1)
                        {
                            //  get shifted image
                            get_shifted_image(tid,iimage,itrans,iover_trans,Fimg_shift_real,Fimg_shift_imag,0,exp_current_Fsize2);
                            
                            // compute difference
                            double diff2 = 0;
                            get_diff(iimage,diff2,Frefctf_real,Frefctf_imag,
                                     Fimg_shift_real,Fimg_shift_imag,Minvsigma2_iimage,
                                     0,exp_current_Fsize2);
                            
                            // Store all diff2 in exp_Mweight_coarse
                            int ihidden = itrans*exp_nr_over_rot + iover_rot;
                            ihidden = ihidden*exp_nr_over_trans + iover_trans;
                            exp_Mweight_coarse_sub.push_back(std::pair<int,double>(ihidden,diff2));
                            //
                            if (diff2 < thread_exp_min_diff2_tid[iimage]) {
                                thread_exp_min_diff2_tid[iimage] = diff2;
                            }
                        } // end loop itrans
                    } // end loop ipsi
                } // end loop idir
            } // end loop iclass
        } // end loop iimage
    }
    //
    void getAllSquaredDifferencesFine()
    {
        // set tile
        int ipsiStep = exp_nr_psi/thread_per_L2_cache;
        
#pragma omp parallel for collapse(4) schedule(dynamic)
        for (int iimage = 0; iimage < exp_nr_images; iimage++) // =(nr_pool)
        {
            for (int iclass = exp_iclass_min; iclass <= exp_iclass_max; iclass++) // 4~10
            {
                for (int idir = 0; idir < exp_nr_dir; idir++) // 192(healpix_order=2)~768(healpix_order=3)
                {
                    for (int ipsi_start = 0; ipsi_start < exp_nr_psi; ipsi_start+=ipsiStep)
                    {
                        //
                        if (mlModel.pdf_direction[iclass][idir] <= 0. ||
                            mlModel.pdf_class.rptrAll()[iclass] <= 0.) continue;
                        int tid = omp_get_thread_num();
                        // thread data
                        FDOUBLE *Fimg_shift_real,*Fimg_shift_imag;
                        auto Minvsigma2_iimage 		  = exp_local_Minvsigma2s[iimage].wptr(exp_current_Fsize2);
                        auto thread_exp_min_diff2_tid = thread_exp_min_diff2[tid].wptr(exp_nr_images);
                        int ipsi_end = std::min(ipsi_start+ipsiStep, exp_nr_psi);
                        for (int ipsi = ipsi_start; ipsi < ipsi_end; ipsi++) // 24(healpix_order=2)~48(healpix_order=3)
                        {
                            int iorient = idir*exp_nr_psi+ipsi;
                            // exp_iclass loop does not always go from 0 to nr_classes!
                            int iorientclass = iclass * exp_nr_dir * exp_nr_psi + iorient;
                            assert(iorientclass >= 0);
                            
                            FDOUBLE *exp_over_rot,*exp_over_tilt,*exp_over_psi;
                            samplingGrid.setOrientation(sampler3d, exp_current_oversampling, idir, ipsi,
                                                        exp_over_rot, exp_over_tilt, exp_over_psi);
                            auto& exp_Mweight_fine_sub = exp_Mweight_fine.wptr_sparse(iimage, iclass, idir, ipsi);
                            
                            bool compute_project_class = true;
                            for (int itrans = 0; itrans < exp_nr_trans; itrans++)// 100(offset_range=10,step=2) or 400(offset_range=10,step=1)
                            {
                                int iclass_rot_trans = iorientclass * exp_nr_trans + itrans;
                                if (!exp_Mcoarse_significant[iimage].rptrAll()[iclass_rot_trans])
                                    continue;
                                // the projected class, only need compute once
                                if (compute_project_class)
                                {
                                    for (int iover_rot = 0; iover_rot < exp_nr_over_rot; iover_rot++) // 1(oversampling=0) or 8(oversampling=1)
                                    {
                                        assert(exp_nr_over_rot < exp_nr_psi);
                                        auto Frefctf_real 			= thread_Frefctf_real.wptr(tid, iover_rot, exp_current_Fsize2);
                                        auto Frefctf_imag 			= thread_Frefctf_imag.wptr(tid, iover_rot, exp_current_Fsize2);
                                        // get projected image
                                        get_projected_class2(iclass,exp_over_rot[iover_rot],exp_over_tilt[iover_rot],exp_over_psi[iover_rot],
                                                             Frefctf_real,Frefctf_imag, 0, exp_current_Fsize2);
                                        //
                                        applyCtfToClass(iimage, Frefctf_real, Frefctf_imag, 0, exp_current_Fsize2);
                                    }
                                    compute_project_class = false;
                                }
                                
                                // the shifted images
                                if (do_shifts_onthefly)
                                {
                                    for (int iover_trans = 0; iover_trans < exp_nr_over_trans; iover_trans++)
                                    {
                                        auto Fimg_shift_real 		= thread_Fimg_real.wptr(tid, iover_trans, exp_current_Fsize2);
                                        auto Fimg_shift_imag 		= thread_Fimg_imag.wptr(tid, iover_trans, exp_current_Fsize2);
                                        particleModel.getShiftedMaskImageOneTileFine(iimage, Fimg_shift_real, Fimg_shift_imag,
                                                                                 	 0, exp_current_Fsize2,itrans,iover_trans,itrans*exp_nr_over_trans+iover_trans);
                                    }
                                }
                                //
                                for (int iover_rot = 0; iover_rot < exp_nr_over_rot; iover_rot++) //8(oversampling=1)
                                {
                                    for (int iover_trans = 0; iover_trans < exp_nr_over_trans; iover_trans++)//4(oversampling=1)
                                    {
                                        auto Frefctf_real 		= thread_Frefctf_real.wptr(tid, iover_rot, exp_current_Fsize2);
                                        auto Frefctf_imag 		= thread_Frefctf_imag.wptr(tid, iover_rot, exp_current_Fsize2);
                                        if (do_shifts_onthefly)
                                        {
                                            Fimg_shift_real 	= thread_Fimg_real.wptr(tid, iover_trans, exp_current_Fsize2);
                                            Fimg_shift_imag 	= thread_Fimg_imag.wptr(tid, iover_trans, exp_current_Fsize2);
                                        }
                                        else
                                        {
                                            //  get shifted image
                                            get_shifted_image(tid,iimage,itrans,iover_trans,Fimg_shift_real,Fimg_shift_imag,0,exp_current_Fsize2);
                                        }
                                        
                                        // Store all diff2 in exp_Mweight
                                        double diff2 = 0;
                                        int ihidden = itrans*exp_nr_over_rot + iover_rot;
                                        ihidden = ihidden*exp_nr_over_trans + iover_trans;
                                        
                                        // compute difference
                                        get_diff(iimage,diff2,Frefctf_real,Frefctf_imag,
                                                 Fimg_shift_real,Fimg_shift_imag,
                                                 Minvsigma2_iimage,0,exp_current_Fsize2);
                                        //
                                        exp_Mweight_fine_sub.push_back( std::pair<int, double>(ihidden,diff2) );
                                        //
                                        if (diff2 < thread_exp_min_diff2_tid[iimage]) {
                                            thread_exp_min_diff2_tid[iimage] = diff2;
                                        }
                                    } // end loop iover_trans
                                }// end loop iover_rot
                            } // end loop itrans
                        } // end loop ipsi
                    }// end loop ipsi_tile
                } // end loop idir
            } // end loop iclass
        } // end loop iimage
        //
    } //
    //
    void getAllSquaredDifferences(bool do_coarse_search)
    {
        exp_current_Fsize2 = exp_current_size*(exp_current_size/2+1);
        
        thread_exp_min_diff2.fill_with_first_touch((std::numeric_limits<double>::max)());
        
#ifdef DATA_STREAM
        // turn off data stream inside openmp for loop
        // global_data_stream.turnOff();
        global_data_stream_before();
#endif
        
        if (do_coarse_search)
        {
            exp_Mweight_coarse.reset();
            getAllSquaredDifferencesCoarse();
        }
        else{
            if (/*tune*/false) {
                for (int i = 0; i < 10; i++) {
                    double t_start = dtime();
                    exp_Mweight_fine.reset();
                    getAllSquaredDifferencesFine();
                    double t_end = dtime();
                    std::cout<<"Iteration : "<<i<<",Time costs : "<<(t_end-t_start)<<std::endl;
                }
				EXIT_ABNORMALLY;
            }
            exp_Mweight_fine.reset();
            getAllSquaredDifferencesFine();
        }
        
#pragma omp parallel for
        for (int iimage = 0; iimage < exp_nr_images; iimage++) {
            exp_min_diff2.wptrAll()[iimage] = (std::numeric_limits<double>::max)();
            for (int thread = 0; thread < maxthreads; thread++) {
                auto thread_exp_min_diff2_tid = thread_exp_min_diff2[thread].wptr(exp_nr_images);
                if (thread_exp_min_diff2_tid[iimage] < exp_min_diff2.rptrAll()[iimage]) {
                    exp_min_diff2.wptrAll()[iimage] = thread_exp_min_diff2_tid[iimage];
                }
            }
        }
        
#ifdef DATA_STREAM
        //
        global_data_stream.turnOn();
        //
        global_data_stream_after();
#endif
        
    }
    
    inline void get_wdiff2(FDOUBLE weight,FDOUBLE& sum_wdiff2,
                           const FDOUBLE* Frefctf_real,const FDOUBLE* Frefctf_imag,
                           const FDOUBLE* Fimg_shift_real,const FDOUBLE* Fimg_shift_imag,
                           const char* threadfake_do_scale_norm_class_iclass,
                           FDOUBLE* thread_wsum_sigma2_noise_tid_iimage,
                           FDOUBLE* thread_wsum_scale_correction_XA_tid_iimage,
                           FDOUBLE* thread_wsum_scale_correction_AA_tid_iimage)
    {
#define SET_SIGMA2_NOISE(n) 													\
/* Use FT of masked image for noise estimation!	*/								\
FDOUBLE diff_real 	= Frefctf_real[n] - Fimg_shift_real[n];							\
FDOUBLE diff_imag 	= Frefctf_imag[n] - Fimg_shift_imag[n];							\
FDOUBLE wdiff2 		= weight * (diff_real * diff_real + diff_imag * diff_imag);	\
/* group-wise sigma2_noise */													\
thread_wsum_sigma2_noise_tid_iimage[n] += wdiff2;								\
/* For norm_correction */														\
sum_wdiff2 += wdiff2;															// end set_sigma2_noise
        //																				//
#define SET_SCALE_CORR(n)														\
FDOUBLE sumXA  = Frefctf_real[n] * Fimg_shift_real[n];								\
sumXA +=	Frefctf_imag[n] * Fimg_shift_imag[n]; 								\
FDOUBLE sumA2  = Frefctf_real[n] * Frefctf_real[n];								\
sumA2 += Frefctf_imag[n] * Frefctf_imag[n];										\
thread_wsum_scale_correction_XA_tid_iimage[n] += weight * sumXA;				\
thread_wsum_scale_correction_AA_tid_iimage[n] += weight * sumA2;				// end set_scale_corr
        //
        // Store weighted sum of squared differences for sigma2_noise estimation
        if (do_scale_correction)
        {
#pragma ivdep
            for (int n = 0; n < exp_current_Fsize2; n++)
            {
                auto do_scale_or_norm = threadfake_do_scale_norm_class_iclass[n];
                // do sigma_noise and scale_correction
                if (do_scale_or_norm == 'S')
                {
                    SET_SIGMA2_NOISE(n)
                    SET_SCALE_CORR(n)
                } // do sigma_noise
                else if (do_scale_or_norm == 'N')
                {
                    SET_SIGMA2_NOISE(n)
                }
            }
        }
        else
        {
#pragma ivdep
            for (int n = 0; n < exp_current_Fsize2; n++)
            {
                int ires = Mresol_fine.rptrAll()[n];
                if (ires > -1){
                    SET_SIGMA2_NOISE(n);
                }
            }
        }
    }
    
    void updateModel()
    {
        int ipsi_tile = exp_nr_psi/thread_per_L2_cache;
        
        for (int iclass = exp_iclass_min; iclass <= exp_iclass_max; iclass++)
        {
            auto data_vs_prior_iclass = mlModel.data_vs_prior_class[iclass].wptr(ori_size/2+1);
            auto threadfake_do_scale_norm_class_iclass  = threadfake_do_scale_norm_class.wptr(iclass, exp_current_Fsize2);
#pragma ivdep // flag do scale or norm(sigma_noise)
            for (int n = 0; n < exp_current_Fsize2; n++){
                int ires = Mresol_fine.rptrAll()[n];
                if (do_scale_correction && ires > -1 && data_vs_prior_iclass[ires] > 3.) // do scale and norm(sigma_noise)
                    threadfake_do_scale_norm_class_iclass[n] = 'S';
                else if (ires > -1)
                    threadfake_do_scale_norm_class_iclass[n] = 'N'; // only do norm(sigma_noise)
                else
                    threadfake_do_scale_norm_class_iclass[n] = 'O';// do nothing
            }
        }
#pragma omp parallel for collapse(4) schedule(dynamic)
        for (int iimage = 0; iimage < exp_nr_images; iimage++)
        {
            for (int iclass = exp_iclass_min; iclass <= exp_iclass_max; iclass++)
            {
                for (int idir = 0; idir < exp_nr_dir; idir++)
                {
                    for (int ipsi_tile_start = 0; ipsi_tile_start < exp_nr_psi; ipsi_tile_start+=ipsi_tile)
                    {
                        int ipsi_tile_end = std::min(ipsi_tile_start+ipsi_tile,exp_nr_psi);
                        for (int ipsi = ipsi_tile_start; ipsi < ipsi_tile_end; ipsi++)
                        {
                            const auto& exp_Mweight_fine_sub = exp_Mweight_fine.wptr_sparse(iimage, iclass, idir, ipsi);
                            if (exp_Mweight_fine_sub.size() == 0) continue;
                            // some thread data
                            int tid = omp_get_thread_num();
                            FDOUBLE *Fimg_shift_real,*Fimg_shift_imag;
                            FDOUBLE *Frefctf_real,*Frefctf_imag;
                            auto thread_wsum_sigma2_noise_tid_iimage 		= thread_wsum_sigma2_noise.wptr(tid, iimage, exp_current_Fsize2);
                            auto thread_wsum_scale_correction_XA_tid_iimage = thread_wsum_scale_correction_XA.wptr(tid, iimage, exp_current_Fsize2);
                            auto thread_wsum_scale_correction_AA_tid_iimage = thread_wsum_scale_correction_AA.wptr(tid, iimage, exp_current_Fsize2);
                            auto threadfake_do_scale_norm_class_iclass  	= threadfake_do_scale_norm_class.wptr(iclass, exp_current_Fsize2);
                            //
                            auto thread_wsum_norm_correction_tid 			= thread_wsum_norm_correction[tid].wptr(exp_nr_images);
                            FDOUBLE *exp_over_rot,*exp_over_tilt,*exp_over_psi;
                            samplingGrid.setOrientation(sampler3d, exp_current_oversampling, idir, ipsi,
                                                        exp_over_rot, exp_over_tilt, exp_over_psi);
                            // flag whether need to re-compute projected class
                            for (int iover_rot = 0; iover_rot < exp_nr_over_rot; iover_rot++)
                            {
                                Frefctf_real = thread_Frefctf_real.wptr(tid, iover_rot, exp_current_Fsize2);
                                Frefctf_real[0] = COMPUTE_FLAG;
                            }
                            // flag whether need to re-compute shifted image
                            int pre_itrans = 11111111;
                            // NOTE : the order of this for loop is same as getAllSquaredDifferencesFine()
                            for (auto const & ihidden_weight : exp_Mweight_fine_sub)
                            {
                                if (ihidden_weight.second < exp_significant_weight.rptrAll()[iimage])
                                    continue;
                                FDOUBLE weight = ihidden_weight.second;
                                //
                                int itrans,iover_rot,iover_trans;
                                decodeOverRotTrans(ihidden_weight.first, iover_rot, itrans, iover_trans);
                                // get projected class and shifted image
                                Frefctf_real = thread_Frefctf_real.wptr(tid, iover_rot, exp_current_Fsize2);
                                Frefctf_imag = thread_Frefctf_imag.wptr(tid, iover_rot, exp_current_Fsize2);
                                if (Frefctf_real[0]==COMPUTE_FLAG)
                                {
                                    get_projected_class2(iclass, exp_over_rot[iover_rot], exp_over_tilt[iover_rot], exp_over_psi[iover_rot],
                                                         Frefctf_real, Frefctf_imag, 0, exp_current_Fsize2);
                                    applyCtfToClass(iimage, Frefctf_real, Frefctf_imag, 0, exp_current_Fsize2);
                                }
                                if (do_shifts_onthefly)
                                {
                                    //
                                    if (pre_itrans != itrans)
                                    {
                                        // flag whether need to re-compute shifted image
                                        for (int iover_trans = 0; iover_trans < exp_nr_over_trans; iover_trans++)
                                        {
                                            Fimg_shift_real = thread_Fimg_real.wptr(tid, iover_trans, exp_current_Fsize2);
                                            Fimg_shift_real[0] = COMPUTE_FLAG;
                                        }
                                        pre_itrans = itrans;
                                    }
                                    Fimg_shift_real = thread_Fimg_real.wptr(tid, iover_trans, exp_current_Fsize2);
                                    Fimg_shift_imag = thread_Fimg_imag.wptr(tid, iover_trans, exp_current_Fsize2);
                                    if (Fimg_shift_real[0]==COMPUTE_FLAG)
                                    {
                                        particleModel.getShiftedMaskImageOneTileFine(iimage, Fimg_shift_real, Fimg_shift_imag,
                                                                                 	 0, exp_current_Fsize2, itrans, iover_trans, itrans*exp_nr_over_trans+iover_trans);
                                    }
                                }
                                else
                                {
                                    // Get the shifted image
                                    get_shifted_image(tid, iimage, itrans, iover_trans, Fimg_shift_real, Fimg_shift_imag, 0, exp_current_Fsize2);
                                }
                                
                                //
                                // Normalise the weight (do this after the comparison with exp_significant_weight!)
                                weight /= exp_sum_weight.rptrAll()[iimage];
                                
                                FDOUBLE sum_wdiff2 = 0.;
                                //
                                get_wdiff2(weight, sum_wdiff2, Frefctf_real, Frefctf_imag,
                                           Fimg_shift_real, Fimg_shift_imag, threadfake_do_scale_norm_class_iclass,
                                           thread_wsum_sigma2_noise_tid_iimage,
                                           thread_wsum_scale_correction_XA_tid_iimage,thread_wsum_scale_correction_AA_tid_iimage);
                                //
                                thread_wsum_norm_correction_tid[iimage] += sum_wdiff2;
                                //
                            }// end loop ihidden(itrans,iover_rot,iover_trans)
                        }// end loop ipsi
                    }// end loop ipsi block
                }// end loop idir
            }// end of iclass
        }// end loop iimage
    }//
    //
}

void getAllSquaredDifferences(bool do_coarse_search){
    CoarseAndFineSearch::getAllSquaredDifferences(do_coarse_search);
}

void convertSquaredDifferencesToWeights(Exp_Mweight_old& exp_Mweight)
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
    // do cross-correlation
    bool do_cross_correlation = (iter == 1 && do_firstiter_cc) || do_always_cc;
    if (do_cross_correlation)
    {
#pragma omp parallel for
        for (int iimage = 0; iimage < exp_nr_images; iimage++)
        {
            int tid = omp_get_thread_num();
            auto thread_exp_max_weight_index_tid = thread_exp_max_weight_index[tid].data();
            auto thread_exp_max_weight_tid = thread_exp_max_weight[tid].wptr(exp_nr_images);
            double mymin_diff2 = (std::numeric_limits<double>::max)();
            int mymin_ihidden;
            // Binarize the squared differences array to skip marginalisation
            // Find the smallest element in this row of exp_Mweight
            auto set_cc = [&](int iimage,int iclass,int idir,int ipsi,int ihidden,double cc){
                // ignore non-determined cc
                if (cc == -999.123456789)
                    return;
                
                // just search for the maximum
                if (cc < mymin_diff2)
                {
                    mymin_diff2 = cc;
                    mymin_ihidden = ihidden;
                    thread_exp_max_weight_tid[iimage] = 1.;
                    int iover_rot,itrans,iover_trans;
                    decodeOverRotTrans(ihidden, iover_rot, itrans, iover_trans);
                    thread_exp_max_weight_index_tid[iimage] = {iimage,iclass,idir,ipsi,iover_rot,itrans,iover_trans};
                }
            };
            for (int iclass = 0; iclass < nr_classes; iclass++) {
                for (int idir = 0; idir < exp_nr_dir; idir++) {
                    for (int ipsi = 0; ipsi < exp_nr_psi; ipsi++) {
                        for (auto const&ihidden_weight : exp_Mweight.wptr_sparse(iimage, iclass, idir, ipsi)){
                            set_cc(iimage,iclass,idir,ipsi,ihidden_weight.first,ihidden_weight.second);
                        }
                    }
                }
            }
            //
            auto thread_exp_sum_weight_tid = thread_exp_sum_weight[tid].wptr(exp_nr_images);
            thread_exp_sum_weight_tid[iimage] += thread_exp_max_weight_tid[iimage];
            // Set all except for the best hidden variable to zero and the smallest element to 1
            int mymin_iclass = thread_exp_max_weight_index_tid[iimage].iclass;
            int mymin_idir = thread_exp_max_weight_index_tid[iimage].idir;
            int mymin_ipsi = thread_exp_max_weight_index_tid[iimage].ipsi;
            exp_Mweight.clear_image(iimage);
            auto& exp_Mweight_sub = exp_Mweight.wptr_sparse(iimage, mymin_iclass, mymin_idir, mymin_ipsi);
            exp_Mweight_sub.push_back( std::pair<int, double>(mymin_ihidden,1) );
#ifdef DATA_STREAM
            global_data_stream.foutDouble(mymin_diff2, "convertSquaredDifferencesToWeights()_mymindiff2", __FILE__, __LINE__);
#endif
        }
    } // end do cross-correlation
    else // do maximum-likelihood
    {
#pragma omp parallel for collapse(3) schedule(dynamic)
        for (int iimage = 0; iimage < exp_nr_images;iimage++)
        {
            for (int iclass = exp_iclass_min; iclass <= exp_iclass_max; iclass++)
            {
                for (int idir = 0; idir < exp_nr_dir; idir++)
                {
                    for (int ipsi = 0; ipsi < exp_nr_psi; ipsi++)
                    {
                        int tid = omp_get_thread_num();
                        auto thread_exp_sum_weight_tid = thread_exp_sum_weight[tid].wptr(exp_nr_images);
                        auto thread_exp_max_weight_tid = thread_exp_max_weight[tid].wptr(exp_nr_images);
                        auto thread_exp_max_weight_index_tid = thread_exp_max_weight_index[tid].data();
                        auto pdf_orientation = mlModel.pdf_direction[iclass][idir];
                        auto& exp_Mweight_sub = exp_Mweight.wptr_sparse(iimage, iclass, idir, ipsi);
                        FDOUBLE offset_x,offset_y;
                        double pdf_offset = (std::numeric_limits<double>::max)();
                        int pre_itrans = (std::numeric_limits<int>::max)();
                        // for itrans,iover_rot,iover_trans
                        for (auto& ihidden_weight : exp_Mweight_sub)
                        {
                            int itrans,iover_rot,iover_trans;
                            decodeOverRotTrans(ihidden_weight.first, iover_rot, itrans, iover_trans);
                            if (itrans != pre_itrans) // recompute pdf_offset
                            {
                                pre_itrans = itrans;
                                // NOTE : To speed things up, only calculate pdf_offset at the coarse sampling.
                                // That should not matter much, and that way one does not need to calculate all the OversampledTranslations
                                sampler3d.getTranslation(itrans, offset_x, offset_y);
                                offset_x += exp_old_offsetx.rptrAll()[iimage];
                                offset_y += exp_old_offsety.rptrAll()[iimage];
                                // Convert offsets back to Angstroms to calculate PDF!
                                pdf_offset = mlModel.calculatePdfOffset(offset_x,offset_y,
                                                                        mlModel.prior_offsetx_class.rptr(nr_classes)[iclass],
                                                                        mlModel.prior_offsety_class.rptr(nr_classes)[iclass]);
                            }
                            
                            double weight = pdf_orientation * pdf_offset;
                            
                            // Sacrifice some performace to get same result as relion
                            double diff2 = ihidden_weight.second - exp_min_diff2.rptrAll()[iimage];
                            
                            // next line because of numerical precision of exp-function
                            if (diff2 > 700.) weight = 0.;
                            // TODO: use tabulated exp function?
                            else weight *= exp(-diff2);
                            
                            // Store the weight
                            ihidden_weight.second = weight;
                            //
                            thread_exp_sum_weight_tid[iimage] += weight;
                            if (weight > thread_exp_max_weight_tid[iimage])
                            {
                                thread_exp_max_weight_tid[iimage] = weight;
                                thread_exp_max_weight_index_tid[iimage] = {iimage,iclass,idir,ipsi,iover_rot,itrans,iover_trans};
                            }
#ifdef DATA_STREAM
                            global_data_stream.foutInt(iimage, "convertSquaredDifferencesToWeights()_iimage", __FILE__, __LINE__);
                            global_data_stream.foutInt(iclass, "convertSquaredDifferencesToWeights()_iclass", __FILE__, __LINE__);
                            global_data_stream.foutInt(idir, "convertSquaredDifferencesToWeights()_idir", __FILE__, __LINE__);
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
                        }// end loop ihidden(itrans,iover_rot,iover_trans)
                    }// end loop ipsi
                }// end loop idir
            }// end loop iclass
        }// end loop iimage
    } // end do maximum-likelihood
#ifdef DATA_STREAM
    global_data_stream.turnOn();
#endif
    //

#pragma omp parallel for
    for (int iimage = 0; iimage < exp_nr_images; iimage++)
    {
        double exp_max_weight = 0;
        GridIndex exp_max_weight_index;
        exp_sum_weight.wptrAll()[iimage] = 0.;
        for (int thread = 0; thread < maxthreads; thread++)
        {
            auto thread_exp_sum_weight_tid = thread_exp_sum_weight[thread].wptr(exp_nr_images);
            exp_sum_weight.wptrAll()[iimage] += thread_exp_sum_weight_tid[iimage];
            auto thread_exp_max_weight_tid = thread_exp_max_weight[thread].wptr(exp_nr_images);
            auto thread_exp_max_weight_index_tid = thread_exp_max_weight_index[thread].data();
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
        FDOUBLE *exp_over_rot,*exp_over_tilt,*exp_over_psi;
        samplingGrid.setOrientation(sampler3d, exp_current_oversampling, idir, ipsi,
                                    exp_over_rot, exp_over_tilt, exp_over_psi);
        //
        exp_metadata[iimage].PMAX = exp_max_weight/exp_sum_weight.rptrAll()[iimage];
        exp_metadata[iimage].CLASS = iclass+1;
        exp_metadata[iimage].ROT = exp_over_rot[iover_rot];
        exp_metadata[iimage].TILT = exp_over_tilt[iover_rot];
        exp_metadata[iimage].PSI = exp_over_psi[iover_rot];
        exp_metadata[iimage].XOFF = exp_old_offsetx.rptrAll()[iimage] + samplingGrid.exp_over_trans_x[itrans][iover_trans];
        exp_metadata[iimage].YOFF = exp_old_offsety.rptrAll()[iimage] + samplingGrid.exp_over_trans_y[itrans][iover_trans];
    }
#ifdef DATA_STREAM
    global_data_stream.foutDouble(exp_sum_weight.wptrAll()[0], "convertSquaredDifferencesToWeights()_exp_sum_weight[0]", __FILE__, __LINE__);
    global_data_stream.check();global_data_stream.flush();
#endif
}

void findAllSignificantPoints(Exp_Mweight_old& exp_Mweight)
{
    
#ifdef DATA_STREAM
    // turn off data stream inside openmp for loop
    // global_data_stream.turnOff();
#endif
    // Now, for each image,find the exp_significant_weight that encompasses adaptive_fraction of exp_sum_weight
#pragma omp parallel for schedule(dynamic)
    for (int iimage = 0;iimage < exp_nr_images; iimage++)
    {
        std::vector<double> sorted_weight;
        int sorted_weight_count = 0;
        
        double frac_weight = 0.;
        double my_significant_weight;
        int my_nr_significant_coarse_samples = 0;
        
        // get weight number
        for (int iclass = 0; iclass < nr_classes; iclass++)
            for (int idir = 0; idir < exp_nr_dir; idir++)
                for (int ipsi = 0; ipsi < exp_nr_psi; ipsi++)
                    sorted_weight_count += exp_Mweight.wptr_sparse(iimage, iclass, idir, ipsi).size();
        
        sorted_weight.resize(sorted_weight_count);
        sorted_weight_count = 0;
        // set
        for (int iclass = 0; iclass < nr_classes; iclass++)
            for (int idir = 0; idir < exp_nr_dir; idir++)
                for (int ipsi = 0; ipsi < exp_nr_psi; ipsi++)
                    for(auto const &ihidden_weight : exp_Mweight.wptr_sparse(iimage, iclass, idir, ipsi))
                        sorted_weight[sorted_weight_count++] = ihidden_weight.second;
        // thread-building block
        // tbb::parallel_sort(sorted_weight.begin(),sorted_weight.end());
        
        // MKL sort rountine
        int info = LAPACKE_dlasrt('I',sorted_weight_count,&sorted_weight[0]);
        if (info<0) {
            std::cerr<<"illegal value : "<<sorted_weight[-info]<<std::endl;
            ERROR_REPORT("the "+std::to_string((long long)-info)+"-th argument had an illegal value.");
        }
        
        // sort from small to large
        // std::sort(&sorted_weight[0], &sorted_weight[sorted_weight_count]);
        
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
                std::cout<<"sum of abs(sorted_weight) = "<<sumvec(&sorted_weight[0], sorted_weight_count)<<std::endl;
                std::cerr << "written sorted_weight.spi" << std::endl;
            }
            ERROR_REPORT("my_nr_significant_coarse_samples == 0");
        }
        
        if (exp_ipass==0)
        {
            exp_metadata[iimage].NR_SIGN = (double)my_nr_significant_coarse_samples;
            
            // Keep track of which coarse samplings were significant were significant for this particle
            // in ipass = 0,exp_Mcoarse_xsize eqaul to exp_Mweight_xsize
            // TODO,TODO better data structure for exp_Mcoarse_significant....
            // transform exp_Mweight coordinate to exp_Mcoarse_significant
            auto set_exp_Mcoarse_significant =[&](int iclass,int idir,int ipsi,int rot_trans_over,bool value)
            {
                int iover_rot,itrans,iover_trans;
                decodeOverRotTrans(rot_trans_over, iover_rot, itrans, iover_trans);
                assert(iover_rot==0);
                assert(iover_trans==0);
                size_t ihidden2 = iclass*exp_nr_dir + idir;
                ihidden2 = ihidden2*exp_nr_psi + ipsi;
                ihidden2 = ihidden2*exp_nr_over_rot + iover_rot;
                ihidden2 = ihidden2*exp_nr_trans + itrans;
                ihidden2 = ihidden2*exp_nr_over_trans + iover_trans;
                exp_Mcoarse_significant[iimage][ihidden2] = value;
            };
            
            exp_Mcoarse_significant[iimage].fill(false);
            //
            for (int iclass = 0; iclass < nr_classes; iclass++) {
                for (int idir = 0; idir < exp_nr_dir; idir++) {
                    for (int ipsi = 0; ipsi < exp_nr_psi; ipsi++) {
                        for (auto const & ihidden_weight : exp_Mweight.wptr_sparse(iimage, iclass, idir, ipsi))
                        {
                            if (ihidden_weight.second >= my_significant_weight)
                                set_exp_Mcoarse_significant(iclass,idir,ipsi,ihidden_weight.first,true);
                        }
                    }
                }
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
            // memory conflict
            wsum_sigma2_noise_igroup[ires] += exp_power_imgs_iimage[ires];
            // Also extend the weighted sum of the norm_correction
            
            // to remove outside????
            exp_wsum_norm_correction.wptrAll()[iimage] += exp_power_imgs_iimage[ires];
        }
        
        // Store norm_correction
        // Multiply by old value because the old norm_correction term was already applied to the image
        if (do_norm_correction)
        {
            double old_norm_correction = exp_metadata[iimage].NORM / mlModel.avg_norm_correction;
            // Now set the new norm_correction in the relevant position of exp_metadata22
            // The factor two below is because exp_wsum_norm_correctiom is similar to sigma2_noise, which is the variance for the real/imag components
            // The variance of the total image (on which one normalizes) is twice this value!
            exp_metadata[iimage].NORM = old_norm_correction * sqrt(exp_wsum_norm_correction.wptrAll()[iimage] * 2.);
            if (!(iter == 1 && do_firstiter_cc) && exp_metadata[iimage].NORM > 10.)
            {
#pragma omp critical
                std::cerr<<"Warning in storeWeightedSums(),please debug this function."<<std::endl;
            }
#pragma omp atomic
            mlModel.wsum_avg_norm_correction += old_norm_correction * sqrt(exp_wsum_norm_correction.wptrAll()[iimage] * 2.);
        }
        // Store weighted sums for scale_correction
        if (do_scale_correction)
        {
            // Divide XA by the old scale_correction and AA by the square of that, because was incorporated into Fctf
            for (int n = 0; n < mlModel.ori_Fsize; n++) {
                exp_wsum_scale_correction_XA[iimage].wptrAll()[n] /= mlModel.scale_correction.rptrAll()[igroup];
                exp_wsum_scale_correction_AA[iimage].wptrAll()[n] /= mlModel.scale_correction.rptrAll()[igroup] * mlModel.scale_correction.rptrAll()[igroup];
            }
            auto wsum_signal_product_spectra_group = mlModel.wsum_signal_product_spectra[igroup].wptr(mlModel.ori_Fsize);
            auto wsum_reference_power_spectra_group = mlModel.wsum_reference_power_spectra[igroup].wptr(mlModel.ori_Fsize);
            for (int n = 0; n < mlModel.ori_Fsize; n++) {
                wsum_signal_product_spectra_group[n]  += exp_wsum_scale_correction_XA[iimage].wptrAll()[n];
                wsum_reference_power_spectra_group[n] += exp_wsum_scale_correction_AA[iimage].wptrAll()[n];
            }
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
            // use previous step's sigma2_noise???
            if (ires > 0)
                logsigma2 += log( 2. * PI * sigma2_noise_igroup[ires]);
        }
        
        
        if (exp_sum_weight.rptrAll()[iimage]==0)
        {
            ERROR_REPORT("ERROR: exp_sum_weight[ipart]==0");
        }
        
        double dLL = 0.;
        if ((iter==1 && do_firstiter_cc) || do_always_cc)
            dLL = -exp_min_diff2.rptrAll()[iimage];
        else
            dLL = log(exp_sum_weight.rptrAll()[iimage]) - exp_min_diff2.rptrAll()[iimage] - logsigma2;
        
        // Also store dLL of each image in the output array
        exp_metadata[iimage].DLL = dLL;
        
#pragma omp atomic
        mlModel.wsum_LL += dLL;
#pragma omp atomic
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
    int ipsi_tile = exp_nr_psi/2;
    
    thread_wsum_scale_correction_XA	.fill_with_first_touch(0.);
    thread_wsum_scale_correction_AA	.fill_with_first_touch(0.);
    thread_wsum_sigma2_noise		.fill_with_first_touch(0.);
    thread_wsum_norm_correction		.fill_with_first_touch(0.);
    
    if (/*tune*/false)
    {
#define VTUNE_ANA 1
        
#ifdef VTUNE_ANA
        //
        for (int i = 0; i < 5; i++)
        {
            double time_start = dtime();
            CoarseAndFineSearch::updateModel();
            
            double time_end = dtime();
            std::cout	<<"count : "<<i<<"time : "<<(time_end-time_start)<<std::endl;
        }
#endif
        EXIT_ABNORMALLY;
    }
    
    CoarseAndFineSearch::updateModel();
    
    // reduce all the threads data
    exp_wsum_norm_correction		.fill(0.);
    exp_wsum_scale_correction_XA	.fill(0.);
    exp_wsum_scale_correction_AA	.fill(0.);
    
    for (int thread = 0; thread < maxthreads; thread++)
    {
        //
        auto thread_wsum_norm_correction_aux = thread_wsum_norm_correction[thread].wptr(exp_nr_images);
        for (int iimage = 0; iimage < exp_nr_images; iimage++)
            exp_wsum_norm_correction.wptrAll()[iimage] += thread_wsum_norm_correction_aux[iimage];
        //
        for (int iimage = 0; iimage < exp_nr_images; iimage++)
        {
            int igroup = exp_metadata[iimage].GROUP_NO-1;
            auto wsum_sigma2_noise_igroup = mlModel.wsum_sigma2_noise[igroup].wptr(mlModel.ori_Fsize);
            auto thread_wsum_sigma2_noise_iimage = thread_wsum_sigma2_noise.wptr(thread, iimage, exp_current_Fsize2);
            auto thread_wsum_scale_correction_XA_iimage = thread_wsum_scale_correction_XA.wptr(thread, iimage, exp_current_Fsize2);
            auto thread_wsum_scale_correction_AA_iimage = thread_wsum_scale_correction_AA.wptr(thread, iimage, exp_current_Fsize2);
            for (int n = 0; n < exp_current_Fsize2; n++)
            {
                // mapping the data to specturm
                int ires = Mresol_fine.rptrAll()[n];
                if (ires > -1){
                    wsum_sigma2_noise_igroup[ires] += thread_wsum_sigma2_noise_iimage[n];
                    if (do_scale_correction)
                    {
                        exp_wsum_scale_correction_XA[iimage].wptrAll()[ires] += thread_wsum_scale_correction_XA_iimage[n];
                        exp_wsum_scale_correction_AA[iimage].wptrAll()[ires] += thread_wsum_scale_correction_AA_iimage[n];
                    }
                }
            }
        }
    }
#ifdef DATA_STREAM
    global_data_stream.foutDouble(exp_wsum_norm_correction.wptrAll()[0], "updateOtherParams()_exp_wsum_norm_correction", __FILE__, __LINE__);
    // global_data_stream.foutDouble(mlModel.wsum_sigma2_noise, mlModel.ori_Fsize, "updateOtherParams()_wsum_sigma2_noise1", __FILE__, __LINE__);
    // global_data_stream.foutDouble(mlModel.wsum_sigma2_noise+(mlModel.nr_groups-1)*mlModel.ori_Fsize, mlModel.ori_Fsize, "updateOtherParams()_wsum_sigma2_noiseN", __FILE__, __LINE__);
    // todo,print different group
    if (do_scale_correction) {
        global_data_stream.foutDouble(exp_wsum_scale_correction_XA[0].wptrAll(), mlModel.ori_Fsize, "updateOtherParams()_exp_wsum_scale_correction_XA", __FILE__, __LINE__);
        global_data_stream.foutDouble(exp_wsum_scale_correction_AA[0].wptrAll(), mlModel.ori_Fsize, "updateOtherParams()_exp_wsum_scale_correction_AA", __FILE__, __LINE__);
    }
    global_data_stream.check();global_data_stream.flush();
#endif
    //
    thread_wsum_pdf_direction	.fill_with_first_touch(0.);
    thread_wsum_pdf_class		.fill_with_first_touch(0.);
    thread_wsum_sigma2_offset	.fill_with_first_touch(0.);
    thread_sumw_group			.fill_with_first_touch(0.);
    
#pragma omp parallel for collapse(3) schedule(dynamic)
    for (int iimage = 0; iimage < exp_nr_images; iimage++)
    {
        for (int iclass = exp_iclass_min; iclass <= exp_iclass_max; iclass++)
        {
            for (int idir = 0; idir < exp_nr_dir; idir++)
            {
                for (int ipsi = 0; ipsi < exp_nr_psi; ipsi++)
                {
                    //
                    int tid = omp_get_thread_num();
                    double thread_wsum_pdf_iclass = 0.;
                    double thread_wsum_sigma2_offset_i = 0.;
                    auto thread_wsum_pdf_direction_tid = thread_wsum_pdf_direction.wptr(tid, iclass, exp_nr_dir);
                    auto thread_sumw_group_tid = thread_sumw_group[tid].wptr(exp_nr_images);
                    auto& exp_Mweight_fine_sub = exp_Mweight_fine.wptr_sparse(iimage, iclass, idir, ipsi);
                    // TODO,vectorize this???
                    for (auto const ihidden_weight : exp_Mweight_fine_sub)
                    {
                        int ihidden = ihidden_weight.first;
                        if (ihidden_weight.second >= exp_significant_weight.rptrAll()[iimage])
                        {
                            FDOUBLE weight = ihidden_weight.second;
                            // Normalise the weight (do this after the comparison with exp_significant_weight!)
                            weight /= exp_sum_weight.rptrAll()[iimage];
                            
                            // Store sum of weights for this group
                            thread_sumw_group_tid[iimage] += weight;
                            
                            // Store weights for this class and orientation
                            thread_wsum_pdf_iclass += weight;
                            
                            thread_wsum_pdf_direction_tid[idir] += weight;
                            
                            int iover_rot,itrans,iover_trans;
                            decodeOverRotTrans(ihidden, iover_rot, itrans, iover_trans);
                            double offsetx = mlModel.prior_offsetx_class.rptrAll()[iclass] - exp_old_offsetx.rptrAll()[iimage] - samplingGrid.exp_over_trans_x[itrans][iover_trans];
                            double offsety = mlModel.prior_offsety_class.rptrAll()[iclass] - exp_old_offsety.rptrAll()[iimage] - samplingGrid.exp_over_trans_y[itrans][iover_trans];
                            
                            //this may cause some false share!
                            thread_wsum_sigma2_offset_i += weight * (offsetx*offsetx+offsety*offsety);
                        }
                    }
                    // Fix the NUMA-Issua for padding size small than cache-line
                    auto thread_wsum_pdf_class_tid = thread_wsum_pdf_class[tid].wptr(nr_classes);
                    auto thread_wsum_sigma2_offset_tid = thread_wsum_sigma2_offset[tid].wptr(1);
                    thread_wsum_pdf_class_tid[iclass] += thread_wsum_pdf_iclass;
                    thread_wsum_sigma2_offset_tid[0] += thread_wsum_sigma2_offset_i;
                }// end loop ipsi
            } // end loop idir
        }// end of iclass
    } // end loop iimage
    
    for (int thread = 0; thread < maxthreads; thread++)
    {
        auto thread_wsum_sigma2_offset_aux = thread_wsum_sigma2_offset[thread].wptr(1);
        mlModel.wsum_sigma2_offset += thread_wsum_sigma2_offset_aux[0];
        auto thread_wsum_pdf_class_aux = thread_wsum_pdf_class[thread].wptr(nr_classes);
        for (int iclass = 0; iclass < nr_classes; iclass++)
        {
            mlModel.wsum_pdf_class.wptr(nr_classes)[iclass] += thread_wsum_pdf_class_aux[iclass];
            auto thread_wsum_pdf_direction_tid = thread_wsum_pdf_direction.wptr(thread, iclass, exp_nr_dir);
            auto mlModel_wsum_pdf_direction = mlModel.wsum_pdf_direction[iclass].wptr(exp_nr_dir);
            for (int idir = 0; idir < exp_nr_dir; idir++){
                mlModel_wsum_pdf_direction[idir] += thread_wsum_pdf_direction_tid[idir];
            }
        }
        auto thread_sumw_group_tid = thread_sumw_group[thread].wptr(exp_nr_images);
        for (int iimage = 0; iimage < exp_nr_images; iimage++) {
            int igroup = exp_metadata[iimage].GROUP_NO-1;
            mlModel.wsum_sumw_group.wptr(mlModel.nr_groups)[igroup] += thread_sumw_group_tid[iimage];
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
    /*** do on convertSquaredDifferencesToWeights()
     ****do on convertSquaredDifferencesToWeights() ***///
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

// TODO : Delay calculating the backProjection()
// TOOD : unless we have enough images and its' significant weight
// TODO : than it will possible to add #omp for loop to each image.
// TODO : to fix the MAPModel.reduce() and large 3D map issue
namespace backProjectionBuffer
{
    int exp_current_Fsize2;
    
    void prepareNomaskImageAndMinvsigma2s()
    {
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
        for (int iimage = 0; iimage < exp_nr_images; iimage++)
        {
            auto exp_local_Fctfs_iimage = exp_local_Fctfs[iimage].wptr(exp_current_Fsize2);
            auto exp_local_Minvsigma2s_iimage = exp_local_Minvsigma2s[iimage].wptr(exp_current_Fsize2);
            
            if (!do_map)
            {
                for (int n = 0; n < exp_current_Fsize2; n++)
                    exp_local_Minvsigma2s_iimage[n] = 1.0;
            }
            // Apply CTF to reference
            if (!do_ctf_correction)
            {
                for (int n = 0; n < exp_current_Fsize2; n++)
                    exp_local_Fctfs_iimage[n] = 1.0;
            }
            if (do_scale_correction)
            {
                // For CTF-terms in BP
                int igroup = exp_metadata[iimage].GROUP_NO-1;
                auto myscale = mlModel.scale_correction.rptrAll()[igroup];
                check_scale(myscale, igroup);
                for (int n = 0; n < exp_current_Fsize2; n++)
                    exp_local_Fctfs_iimage[n] = exp_local_Fctfs_iimage[n] * myscale;
            }
            // NOTE : later code do some trick for Fimages and exp_local_Minvigma2s
            // Operation : exp_local_Minvsigma2s =  exp_local_Fctfs*exp_local_Minvsigma2s
            auto exp_local_Minvsigma2s_X_Fctfs_iimage = exp_local_Minvsigma2s[iimage].wptr(exp_current_Fsize2);
            for (int n = 0; n < exp_current_Fsize2; n++)
                exp_local_Minvsigma2s_X_Fctfs_iimage[n] = exp_local_Minvsigma2s_iimage[n]*exp_local_Fctfs_iimage[n];
            // NOTE : multiply nomask image by Minvsigma2s and Fctfs(exp_local_Minvsigma2s_X_Fctfs)
            // this will reduce some memory access in inner for loop
            auto Fimages_nomask_real = particleModel.Fimages_nomask_real[iimage].wptr(exp_current_Fsize2);
            auto Fimages_nomask_imag = particleModel.Fimages_nomask_imag[iimage].wptr(exp_current_Fsize2);
            for (int n = 0; n < exp_current_Fsize2; n++) {
                Fimages_nomask_real[n] *= exp_local_Minvsigma2s_X_Fctfs_iimage[n];
                Fimages_nomask_imag[n] *= exp_local_Minvsigma2s_X_Fctfs_iimage[n];
            }
            // NOTE : multiply Fctf to exp_local_Minvsigma2s
            // then exp_local_Minvsigma2s will be equal to Minvsigma2s*Fctf^2
            auto exp_local_Minvsigma2s_X_Fctfs2_iimage = exp_local_Minvsigma2s[iimage].wptr(exp_current_Fsize2);
            for (int n = 0; n < exp_current_Fsize2; n++)
                exp_local_Minvsigma2s_X_Fctfs2_iimage[n] = exp_local_Minvsigma2s_X_Fctfs_iimage[n]*exp_local_Fctfs_iimage[n];
        }
        
        //
        particleModel.pregetShiftedImagesNomask(exp_Fimgs_shifted_real, exp_Fimgs_shifted_imag, exp_current_size, samplingGrid);
    }
    
    inline void get_shifted_image(int tid,int iimage,int itrans,int iover_trans,
                                  FDOUBLE* &Fimg_shift_real,FDOUBLE* &Fimg_shift_imag,int exp_current_Fsize2)
    {
        int itrans_over_trans = itrans*exp_nr_over_trans+iover_trans;
        Fimg_shift_real = exp_Fimgs_shifted_real.wptr(iimage, itrans_over_trans, exp_current_Fsize2);
        Fimg_shift_imag = exp_Fimgs_shifted_imag.wptr(iimage, itrans_over_trans, exp_current_Fsize2);
    }
    
    inline void global_data_stream_inner(int iimage,int itrans,int iover_trans,double weight,
                                         FDOUBLE* Fimg_real,FDOUBLE* Fimg_imag,FDOUBLE* Fweight,int exp_current_Fsize2)
    {
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
        global_data_stream.foutDouble(Fimg_real, exp_current_Fsize2, "backProjection()_backProjection()_Fimg_real", __FILE__, __LINE__);
        global_data_stream.foutDouble(Fimg_imag, exp_current_Fsize2, "backProjection()_backProjection()_Fimg_imag", __FILE__, __LINE__);
        global_data_stream.foutDouble(Fweight, exp_current_Fsize2, "backProjection()_backProjection()_Fweight", __FILE__, __LINE__);
        global_data_stream.check();global_data_stream.flush();
#endif
    }
    
//    inline void addShiftedImages(int iimage,double weight,double* Fimg_shift_real,double* Fimg_shift_imag,
//                                 double* Fimg_real,double* Fimg_imag,double* Fweight,int n_start,int n_end)
//    {
//        assert(weight >= exp_significant_weight.rptrAll()[iimage]);
//        // Inside the loop over all translations and all part_id sum all shift Fimg's and their weights
//        // Then outside this loop do the actual backprojection
//        // Now that reference projection has been made loop over someParticles!
//        auto exp_local_Minvsigma2s_X_Fctfs2_iimage = exp_local_Minvsigma2s[iimage].wptr(exp_current_Fsize2);
//        weight /= exp_sum_weight.rptrAll()[iimage];
//        // Store sum of weight*SSNR*Fimg in data and sum of weight*SSNR in weight
//        // Use the FT of the unmasked image to back-project in order to prevent reconstruction artefacts! SS 25oct11
//#pragma ivdep
//        for (int n = n_start; n < n_end; n++)
//        {
//            Fimg_real[n] += Fimg_shift_real[n] * weight;
//            Fimg_imag[n] += Fimg_shift_imag[n] * weight;
//            // now Fweight stores sum of all w
//            // Note that CTF needs to be squared in Fweight, weightxinvsigma2 already contained one copy
//            Fweight[n] += exp_local_Minvsigma2s_X_Fctfs2_iimage[n] * weight;
//        }
//    }
    
    void doBackProjectionNew()
    {
        exp_current_Fsize2 = exp_current_size*(exp_current_size/2+1);
        //
        prepareNomaskImageAndMinvsigma2s();
        //
        int ipsi_tile = exp_nr_psi/thread_per_L2_cache;
#ifdef DATA_STREAM
        // turn off data stream inside openmp for loop
        global_data_stream.turnOff();
#endif
        for (int iclass = exp_iclass_min; iclass <= exp_iclass_max; iclass++)
        {
        // compare the shifted image with rotated direction
#pragma omp parallel for collapse(3) schedule(dynamic)
            for (int iimage = 0; iimage < exp_nr_images; iimage++)
            {
                for (int idir = 0; idir < exp_nr_dir; idir++)
                {
                    for (int ipsi_tile_start = 0; ipsi_tile_start < exp_nr_psi; ipsi_tile_start+=ipsi_tile)
                    {
                        int tid = omp_get_thread_num();
                        FDOUBLE *Frefctf_real,*Frefctf_imag,*Fweight;
                        FDOUBLE *Fimg_shift_real,*Fimg_shift_imag;
                        auto exp_local_Minvsigma2s_X_Fctfs2_iimage = exp_local_Minvsigma2s[iimage].wptr(exp_current_Fsize2);
                        //
                        int ipsi_tile_end = std::min(ipsi_tile_start+ipsi_tile, exp_nr_psi);
                        for (int ipsi = ipsi_tile_start; ipsi < ipsi_tile_end; ipsi++)
                        {
                            const auto& exp_Mweight_fine_sub = exp_Mweight_fine.wptr_sparse(iimage, iclass, idir, ipsi);
                            int exp_Mweight_fine_sub_size = exp_Mweight_fine_sub.size();
                            if (exp_Mweight_fine_sub_size == 0) continue;
                            FDOUBLE *exp_over_rot,*exp_over_tilt,*exp_over_psi;
                            samplingGrid.setOrientation(sampler3d, exp_current_oversampling, idir, ipsi,
                                                        exp_over_rot, exp_over_tilt, exp_over_psi);
                            // flag whether need to re-compute projected class
                            for (int iover_rot = 0; iover_rot < exp_nr_over_rot; iover_rot++)
                            {
                                Frefctf_real = thread_Frefctf_real.wptr(tid, iover_rot, exp_current_Fsize2);
                                Frefctf_real[0] = (std::numeric_limits<FDOUBLE>::max)();
                            }
                            
                            // NOTE : the order of this for loop is same as getAllSquaredDifferencesFine()
                            for (auto const & ihidden_weight : exp_Mweight_fine_sub)
                            {
                                if (ihidden_weight.second < exp_significant_weight.rptrAll()[iimage])
                                    continue;
                                assert(ihidden_weight.second >= exp_significant_weight.rptrAll()[iimage]);
                                FDOUBLE weight = ihidden_weight.second;
                                //
                                int itrans,iover_rot,iover_trans;
                                decodeOverRotTrans(ihidden_weight.first, iover_rot, itrans, iover_trans);
                                Frefctf_real 	= thread_Frefctf_real.wptr(tid, iover_rot, exp_current_Fsize2);
                                Frefctf_imag 	= thread_Frefctf_imag.wptr(tid, iover_rot, exp_current_Fsize2);
                                Fweight   		= thread_Fweight		.wptr(tid, iover_rot, exp_current_Fsize2);
                                if (Frefctf_real[0]==(std::numeric_limits<FDOUBLE>::max)())
                                {
                                    for (int n = 0; n < exp_current_Fsize2; n++)
                                        Frefctf_real[n] = Frefctf_imag[n] = Fweight[n] = 0.;
                                }
                                // get the shifted image
                                get_shifted_image(tid, iimage, itrans, iover_trans,Fimg_shift_real, Fimg_shift_imag, exp_current_Fsize2);
                                
                                // Inside the loop over all translations and all part_id sum all shift Fimg's and their weights
                                // Then outside this loop do the actual backprojection
                                // Now that reference projection has been made loop over someParticles!
                                weight /= exp_sum_weight.rptrAll()[iimage];
                                // Store sum of weight*SSNR*Fimg in data and sum of weight*SSNR in weight
                                // Use the FT of the unmasked image to back-project in order to prevent reconstruction artefacts! SS 25oct11
                                #pragma ivdep
                                for (int n = 0; n < exp_current_Fsize2; n++)
                                {
                                    Frefctf_real[n] += Fimg_shift_real[n] * weight;
                                    Frefctf_imag[n] += Fimg_shift_imag[n] * weight;
                                    // now Fweight stores sum of all w
                                    // Note that CTF needs to be squared in Fweight, weightxinvsigma2 already contained one copy
                                    Fweight[n] += exp_local_Minvsigma2s_X_Fctfs2_iimage[n] * weight;
                                }
                                //
                                global_data_stream_inner(iimage, itrans, iover_trans, weight, Frefctf_real, Frefctf_imag, Fweight, exp_current_Fsize2);
                                
                            }// end loop ihidden(itrans,iover_rot,iover_trans)
                            
                            // put all image back
                            for (int iover_rot = 0; iover_rot < exp_nr_over_rot; iover_rot++)
                            {
                                Frefctf_real	= thread_Frefctf_real.wptr(tid, iover_rot, exp_current_Fsize2);
                                Frefctf_imag	= thread_Frefctf_imag.wptr(tid, iover_rot, exp_current_Fsize2);
                                Fweight			= thread_Fweight.wptr(tid, iover_rot, exp_current_Fsize2);
                                if (Frefctf_real[0] != (std::numeric_limits<FDOUBLE>::max)())
                                {
                                    FDOUBLE A[3][3];
                                    Euler_angles2matrix(exp_over_rot[iover_rot], exp_over_tilt[iover_rot], exp_over_psi[iover_rot], A);
                                    mapModel.set2DFourierTransform(tid, iclass, Frefctf_real, Frefctf_imag, exp_current_size, A, false, Fweight);
                                }
                            }
                        }// end loop ipsi
                    }// end loop ipsi_tile
                } // end loop idir
            } // end loop iimage
            mapModel.reduceThreadMap(iclass);
        } //end loop iclass
        
        //
        // Reduce thread classes outside the backprojection()
#ifdef DATA_STREAM
        global_data_stream.turnOn();
        int vol_size3 = mapModel.projector[0].pad_size;
        vol_size3 = vol_size3*vol_size3*(vol_size3/2+1);
        // global_data_stream.foutDouble(mapModel.backprojector[0].data_real.wptr(), vol_size3, "backProjection()_backRef_real_0", __FILE__, __LINE__);
        // global_data_stream.foutDouble(mapModel.backprojector[0].data_imag.wptr(), vol_size3, "backProjection()_backRef_imag_0", __FILE__, __LINE__);
        // global_data_stream.foutDouble(mapModel.backprojector[nr_classes-1].data_real.wptr(), vol_size3, "backProjection()_backRef_real_N-1", __FILE__, __LINE__);
        // global_data_stream.foutDouble(mapModel.backprojector[nr_classes-1].data_imag.wptr(), vol_size3, "backProjection()_backRef_imag_N-1", __FILE__, __LINE__);
        // global_data_stream.foutDouble(mapModel.backprojector[0].weight.wptr(), vol_size3, "backProjection()_backRef_weight_0", __FILE__, __LINE__);
        // global_data_stream.foutDouble(mapModel.backprojector[nr_classes-1].weight.wptr(), vol_size3, "backProjection()_backRef_weight_N-1", __FILE__, __LINE__);
        global_data_stream.foutDouble(mapModel.backprojector[0].weight.wptr(), vol_size3, "backProjection()_backRef_weight_0", __FILE__, __LINE__);
        global_data_stream.foutDouble(mapModel.backprojector[(nr_classes-1)].weight.wptr(), vol_size3, "backProjection()_backRef_weight_N-1", __FILE__, __LINE__);
        global_data_stream.check();global_data_stream.flush();
#endif
    }
    
}

void backProjection()
{
    if (/*tune*/false)
    {
#define VTUNE_ANA 1
        
#ifdef VTUNE_ANA
        //
        for (int i = 0; i < 5; i++)
        {
            double time_start = dtime();
            backProjectionBuffer::doBackProjectionNew();
            
            double time_end = dtime();
            std::cout	<<"count : "<<i<<"time : "<<(time_end-time_start)<<std::endl;
        }
#endif
        EXIT_ABNORMALLY;
    }
    
    backProjectionBuffer::doBackProjectionNew();
}


void maximization()
{
    NODE0ONLY std::cout << " Maximization ..." << std::endl;
    
    int first_local_class,last_local_class;
    first_local_class = 0;last_local_class=nr_classes-1;
#ifdef DO_RECONSTRUCT_EACH_NODE
    int nr_local_classes = divide_equally(nr_classes, nodes, node, first_local_class, last_local_class);
#endif
    // First reconstruct the images for each class
#pragma omp parallel for
    for (int iclass = first_local_class; iclass <= last_local_class; iclass++)
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
        // if(mlModel.wsum_pdf_clsss[iclass] > 0.)
        if (mlModel.pdf_class.rptrAll()[iclass] > 0.)
        {
            // update
            mlModel.sigma2_class[iclass].zero();
            auto sigma2_iclass = mlModel.sigma2_class[iclass].wptr(ori_size/2+1);
            // update or not depend on update_tau2_with_fsc
            auto tau2_iclass = mlModel.tau2_class[iclass].wptr(ori_size/2+1);
            // update or not depend on update_tau2_with_fsc
            mlModel.data_vs_prior_class[iclass].zero();
            auto data_vs_prior_class_iclass = mlModel.data_vs_prior_class[iclass].wptr(ori_size/2+1);
            auto fsc_halves_class_iclass = mlModel.fsc_halves_class[iclass].wptr(ori_size/2+1);
            
            // void reconstruction(int iclass,int max_iter_preweight,bool do_map,double tau2_fudge,double* tau2,double* sigma2,double* data_vs_prior,
            //                     const double* fsc, /* only input*/, double normalise = 1., bool update_tau2_with_fsc = false,
            //                     bool is_whole_instead_of_half = false,int minres_map = -1)
            
            mapModel.reconstruction(iclass, gridding_nr_iter, do_map, tau2_fudge_factor, tau2_iclass, sigma2_iclass,
                                    data_vs_prior_class_iclass, fsc_halves_class_iclass,
                                    mlModel.wsum_pdf_class.rptrAll()[iclass], false, false, 1);
            
            
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
    
#if defined(USEMPI) && defined(DO_RECONSTRUCT_EACH_NODE)
    mlModel.broadcastData();
    mapModel.broadcastData();
#endif
    
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
    FDOUBLE sum_weight = 0.;
    for (int iclass = 0; iclass < nr_classes; iclass++)
        sum_weight += mlModel.wsum_pdf_class.rptrAll()[iclass];
    // Update average norm_correction
    if (do_norm_correction)
    {
        mlModel.avg_norm_correction = mlModel.wsum_avg_norm_correction / sum_weight;
    }
    
    auto sum=[&](FDOUBLE* V,int N){
        FDOUBLE s = 0;
        for (int i = 0; i < N; i++) s += V[i];
        return s;
    };
    
    if ( do_scale_correction && !(iter==1 && do_firstiter_cc) )
    {
        std::vector<FDOUBLE> sorted(mlModel.nr_groups);
        //
        for (int igroup = 0; igroup < mlModel.nr_groups; igroup++)
        {
            auto sumXA = sum(mlModel.wsum_signal_product_spectra[igroup].wptr(mlModel.ori_Fsize),mlModel.ori_Fsize);
            auto sumAA = sum(mlModel.wsum_reference_power_spectra[igroup].wptr(mlModel.ori_Fsize),mlModel.ori_Fsize);
            if (sumAA > 0.)
                mlModel.scale_correction.wptrAll()[igroup] = sumXA / sumAA;
            else
                mlModel.scale_correction.wptrAll()[igroup] = 1.;
            sorted[igroup] = mlModel.scale_correction.rptrAll()[igroup];
        }
        
        // TODO! Avoid extremities in scale estimates, because they lead to catastrophic events and instabilities in refinement
        // Let's exclude anything bigger than 5x the median or smaller than 1/5 of the median...
        // Use the median instead of the mean, because it is much more robust to outliers.
        std::sort(sorted.begin(), sorted.end());
        auto median = sorted[mlModel.nr_groups / 2];
        
        FDOUBLE avg_scale_correction = 0., nr_part = 0.;
        for (int igroup = 0; igroup < mlModel.nr_groups; igroup++)
        {
            
            if (mlModel.scale_correction.rptrAll()[igroup] > 5. * median)
                mlModel.scale_correction.wptrAll()[igroup] = 5. * median;
            else if (mlModel.scale_correction.rptrAll()[igroup] < median / 5.)
                mlModel.scale_correction.wptrAll()[igroup] =  median / 5.;
            
            avg_scale_correction += (FDOUBLE)(mlModel.nr_particles_group.rptrAll()[igroup]) * mlModel.scale_correction.rptrAll()[igroup];
            nr_part += (FDOUBLE)(mlModel.nr_particles_group.rptrAll()[igroup]);
        }
        
        // Constrain average scale_correction to one.
        avg_scale_correction /= nr_part;
        for (int igroup = 0; igroup < mlModel.nr_groups; igroup++)
        {
            mlModel.scale_correction.wptrAll()[igroup] /= avg_scale_correction;
            //            #define DEBUG_UPDATE_SCALE
#ifdef DEBUG_UPDATE_SCALE
            std::cerr<< "Group "<<igroup+1<<": scale_correction= "<<mlModel.scale_correction[igroup]<<std::endl;
            for (int i = 0; i < mlModel.ori_Fsize; i++)
                if (mlModel.wsum_reference_power_spectra[igroup*mlModel.ori_Fsize+i]> 0.)
                    std::cerr << " i= " << i << " XA= " << mlModel.wsum_signal_product_spectra[igroup*mlModel.ori_Fsize+i]
                    << " A2= " << mlModel.wsum_reference_power_spectra[igroup*mlModel.ori_Fsize+i]
                    << " XA/A2= " << mlModel.wsum_signal_product_spectra[igroup*mlModel.ori_Fsize+i]/mlModel.wsum_reference_power_spectra[igroup*mlModel.ori_Fsize+i] << std::endl;
#endif
        }
        
    }
    
    // Update model.pdf_class vector (for each k)
    for (int iclass = 0; iclass < nr_classes; iclass++)
    {
        mlModel.pdf_class.wptrAll()[iclass] = mlModel.wsum_pdf_class.rptrAll()[iclass] / sum_weight;
        
        // for 2D also update priors of translations for each class!
        // None
        // TODO
        
        for (int idir = 0; idir < sampler3d.NrDir(); idir++)
        {
            mlModel.pdf_direction[iclass][idir] = mlModel.wsum_pdf_direction[iclass][idir] / sum_weight;
        }
    }
    
    mlModel.sigma2_offset = (mlModel.wsum_sigma2_offset) / (2. * sum_weight);
    
    // TODO: update estimates for sigma2_rot, sigma2_tilt and sigma2_psi!
    if(!(iter == 1 && do_firstiter_cc))
    {
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
    }
    
    
    // After the first iteration the references are always CTF-corrected
    if (do_ctf_correction)
        refs_are_ctf_corrected = true;
    
    // Some statistics to output
    
    mlModel.LL = mlModel.wsum_LL;
    if ((iter==1 && do_firstiter_cc) || do_always_cc)
        mlModel.LL /= sum_weight; // this now stores the average ccf
    
    mlModel.ave_Pmax = mlModel.wsum_ave_Pmax / sum_weight;
    
    // After the first, special iteration, apply low-pass filter of -ini_high again
    if (iter == 1 && do_firstiter_cc)
    {
        mapModel.applyLowPassFilter();
        if (ini_high > 0.)
        {
            // Adjust the tau2_class and data_vs_prior_class, because they were calculated on the unfiltered maps
            // This is merely a matter of having correct output in the model.star file (these values are not used in the calculations)
            FDOUBLE radius = ori_size * pixel_size / ini_high;
            radius -= mapModel.width_fmask_edge / 2.;
            FDOUBLE radius_p = radius + mapModel.width_fmask_edge;
            
            for (int iclass = 0; iclass < nr_classes; iclass++)
            {
                auto tau2_aux = mlModel.tau2_class[iclass].wptr(ori_size/2+1);
                auto data_vs_prior_class_aux = mlModel.data_vs_prior_class[iclass].wptr(ori_size/2+1);
                for (int rr = 0; rr < ori_size/2+1; rr++)
                {
                    FDOUBLE r = (FDOUBLE)rr;
                    if (r < radius)
                        continue;
                    else if (r > radius_p)
                    {
                        tau2_aux[rr] = 0.;
                        data_vs_prior_class_aux[rr] = 0.;
                    }
                    else
                    {
                        FDOUBLE raisedcos = 0.5 - 0.5 * cos(rome_pi * (radius_p - r) / mapModel.width_fmask_edge);
                        tau2_aux[rr] *= raisedcos * raisedcos;
                        data_vs_prior_class_aux[rr] *= raisedcos * raisedcos;
                    }
                }
            }
        }
        
        if (true/*do_generate_seeds*/ && nr_classes > 1)
        {
            // In the first CC-iteration only a single reference was used
            // Now copy this one reference to all K references, for seed generation in the second iteration
            for (int iclass = 1; iclass < nr_classes; iclass++)
            {
                auto tau2_class_iclass = mlModel.tau2_class[iclass].wptr(ori_size/2+1);
                auto data_vs_prior_class_iclass = mlModel.data_vs_prior_class[iclass].wptr(ori_size/2+1);
                for (int i = 0; i < ori_size/2+1; i++) {
                    tau2_class_iclass[i] = mlModel.tau2_class[0][i];
                    data_vs_prior_class_iclass[i] = mlModel.data_vs_prior_class[0][i];
                }
                auto pdf_direction_iclass = mlModel.pdf_direction[iclass].wptr(mlModel.nr_directions);
                for (int i = 0; i < mlModel.nr_directions; i++) {
                    pdf_direction_iclass[i] = mlModel.pdf_direction[0][i];
                }
                
                auto& Iref_iclass = mapModel.Irefs[iclass];
                for (int i = 0; i < Iref_iclass.dimzyx; i++) {
                    Iref_iclass(0,0,i) = mapModel.Irefs[0](0,0,i);
                }
                mlModel.pdf_class.wptrAll()[iclass] = mlModel.pdf_class.rptrAll()[0] / nr_classes;
            }
            mlModel.pdf_class.wptrAll()[0] /= nr_classes;
        }
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
    
void readResult()
{
    if (continue_fn.find("_optimiser.star")!=std::string::npos) {
        std::string model_fn,sampling_fn;
        //
        std::tie(model_fn, sampling_fn) = readFromOptimizer(continue_fn);
        //
        mlModel.readFromModel(model_fn,mapModel);
        //
        sampler3d.readFromSampling(sampling_fn);
    }
    else if(continue_fn.find("_backup.back")) {
        std::string statusFn = continue_fn.substr(0,continue_fn.find("_backup"));
        statusTracer.recoveryStatus(statusFn);
    }
    else{
        ERROR_REPORT("Wrong continue file name "+continue_fn+",use *_optimiser.star or *_backup.back");
    }
}
    
void writeAndCheckResult()
{
    NODE0ONLY
    {
        // check result
        if(guidance_fn!="NULL" && guidance_fn!="BACK")
            statusTracer.checkStatus(guidance_fn+"_it"+num2str(iter));
        //
        if(guidance_fn=="BACK")
        {
            std::string fn 		= write_path+write_fn;
            std::string iterStr = num2str(iter);
            std::string statusFn= fn+"_it"+iterStr;
            statusTracer.backupStatus(statusFn);
        }
		// write result
        if(write_path != "NULL")
        {
            std::string fn = write_path+write_fn;
            std::string iterStr 	= num2str(iter);
            //
            std::string fn_optimiser = fn+"_it"+iterStr+"_optimiser";
            writeOutOptimizer(fn_optimiser);
            // Write Classes
            std::string fn_class    = fn+"_it"+iterStr+"_class";
            mapModel.writeResult(fn_class);
            // Write Metadata
            std::string fn_metadata = fn+"_it"+iterStr+"_data";
            metadata.writeToStar(fn_metadata);
            //
            std::string fn_model = fn+"_it"+iterStr+"_model";
            mlModel.writeOutModel(fn_model, mapModel, metadata);
            //
            mlModel.writeOutBild(fn_class, mapModel, sampler3d);
            //
            std::string fn_sampling = fn+"_it"+iterStr+"_sampling";
            sampler3d.writeOutSampling(fn_sampling);
        }
    }
}

void printMem(int set_nr_pool)
{
    std::string fn = write_path+write_fn;
    std::string iterStr	= num2str(iter);
    std::string fn_ram 	= fn+"_it"+iterStr+"_ram";
    std::ofstream ramFile;
    ramFile.open((fn_ram+".txt").c_str(), std::ios::out);

#define FOUT ramFile<<std::setw(40)<<std::left
    
    double perGb = 1/1024./1024./1024.;
    int ori_Fsize = (ori_size/2+1);
    int current_Fsize = (current_size/2+1);
    int current_Fsize2 = current_Fsize*current_size;
    //
    double images_data_size               = nr_local_images*sizeof(double)*ori_size*(ori_size/2+1)*perGb;
    double metadata_size                  = nr_global_images*sizeof(MetaDataElem)*perGb;
    double fix_size = images_data_size + metadata_size;
    FOUT<<"fixed size : "<<std::endl;
    FOUT<<"images_data size : "<<images_data_size<<" GB."<<std::endl;
    FOUT<<"metadata_size : "<<metadata_size<<" GB."<<std::endl;
    FOUT<<"Total : "<<fix_size<<" GB."<<std::endl;
    FOUT<<"------------------------------------------------"<<std::endl;
    mlModel.printSpaceInfo(ramFile);
    mapModel.printSpaceInfo(ramFile);
    //
    double exp_metadata_size             		=	set_nr_pool*sizeof(MetaDataElem)*perGb;
    double exp_imgs_size                 		=	set_nr_pool*ori_size*ori_size*sizeof(double)*perGb;
    double exp_power_imgs_size           		=	set_nr_pool*ori_Fsize*sizeof(double)*perGb;
    double exp_highres_Xi2_imgs_size     		=	set_nr_pool*sizeof(double)*perGb;
    double exp_min_diff2_size            		=	set_nr_pool*sizeof(double)*perGb;
    double exp_old_offsetx_size          		=	set_nr_pool*sizeof(double)*perGb;
    double exp_old_offsety_size					=	set_nr_pool*sizeof(double)*perGb;
    double exp_wsum_scale_correction_XA_size	=	set_nr_pool*ori_Fsize*sizeof(double)*perGb;
    double exp_wsum_scale_correction_AA_size	=	set_nr_pool*ori_Fsize*sizeof(double)*perGb;
    double exp_wsum_norm_correction_size		=	set_nr_pool*sizeof(double)*perGb;
    double exp_Fimgs_size                		= 	2*set_nr_pool*current_size*current_Fsize*sizeof(double)*perGb;
    double exp_Fimgs_nomask_size         		= 	2*set_nr_pool*current_size*current_Fsize*sizeof(double)*perGb;
    double exp_Fctfs_size                		= 	set_nr_pool*current_size*current_Fsize*sizeof(double)*perGb;
    double exp_significant_weight_size   		= 	set_nr_pool*sizeof(double)*perGb;
    double exp_sum_weight_size           		= 	set_nr_pool*sizeof(double)*perGb;
    double exp_Fimgs_shifted_nomask_size		= 	0;
    double exp_local_sqrtXi2_size				=	set_nr_pool*sizeof(double)*perGb;
    double exp_over_psi_size					=	sampler3d.NrPsi()*sampler3d.NrDir()*sampler3d.oversamplingFactorOrientations(adaptive_oversampling)*sizeof(double)*perGb;
    double exp_over_rot_size					=	sampler3d.NrPsi()*sampler3d.NrDir()*sampler3d.oversamplingFactorOrientations(adaptive_oversampling)*sizeof(double)*perGb;
    double exp_over_tilt_size					=	sampler3d.NrPsi()*sampler3d.NrDir()*sampler3d.oversamplingFactorOrientations(adaptive_oversampling)*sizeof(double)*perGb;
    double exp_over_trans_x						=	sampler3d.NrTrans()*sampler3d.oversamplingFactorTranslations(adaptive_oversampling)*sizeof(double)*perGb;
    double exp_over_trans_y						=	sampler3d.NrTrans()*sampler3d.oversamplingFactorTranslations(adaptive_oversampling)*sizeof(double)*perGb;
    double exp_Fimgs_shifted_size				=	2*set_nr_pool*sampler3d.NrTrans(adaptive_oversampling)*current_Fsize*sizeof(double)*perGb;
    double exp_local_Fctfs_size					=	set_nr_pool*current_Fsize2*sizeof(double)*perGb;
    double exp_local_Minvsigma2s_size			=	set_nr_pool*current_Fsize2*sizeof(double)*perGb;
    double exp_Rot_significant_size 			=	nr_classes*sampler3d.NrDir()*sampler3d.NrPsi()*sizeof(bool)*perGb;
    double exp_Mcoarse_significant_size			=	set_nr_pool*nr_classes*sampler3d.NrPoints(0)*sizeof(bool)*perGb;
    double exp_Mweight_coarse_size				=	set_nr_pool*nr_classes*sampler3d.NrPoints(0)*sizeof(double)*perGb;
    double exp_Mweight_fine_size				=	set_nr_pool*nr_classes*sampler3d.NrPoints(adaptive_oversampling)*sizeof(double)*perGb;
    //
    double unfixed_size = exp_metadata_size + exp_imgs_size + exp_power_imgs_size + exp_highres_Xi2_imgs_size + exp_min_diff2_size \
    + exp_old_offsetx_size + exp_old_offsety_size + exp_wsum_scale_correction_XA_size \
    + exp_wsum_scale_correction_AA_size + exp_wsum_norm_correction_size + exp_Fimgs_size + exp_Fimgs_nomask_size \
    + exp_Fctfs_size + exp_significant_weight_size + exp_sum_weight_size + exp_Fimgs_shifted_nomask_size \
    + exp_local_sqrtXi2_size + exp_over_psi_size + exp_over_rot_size + exp_over_tilt_size + exp_over_trans_x \
    + exp_over_trans_y + exp_Fimgs_shifted_size + exp_local_Fctfs_size + exp_local_Minvsigma2s_size + exp_Rot_significant_size \
    + exp_Mcoarse_significant_size + exp_Mweight_coarse_size + exp_Mweight_fine_size;
    int sugg_set_nr_pool = (6.5 - fix_size)/(unfixed_size/set_nr_pool);
    FOUT<<"------------------------------------------------------"<<std::endl;
    FOUT<<"suggestion nr_pool : "<<sugg_set_nr_pool<<std::endl;
    FOUT<<"real nr_pool : "<<set_nr_pool<<std::endl;
    FOUT<<"------------------------------------------------------"<<std::endl;
    // nr_pool = set_nr_pool;
    FOUT<<"unfixed size(Expectation step) : "<<std::endl;
    FOUT<<"exp_metadata_size : "<<exp_metadata_size<<" GB."<<std::endl;
    FOUT<<"exp_imgs_size : "<<exp_imgs_size<<" GB."<<std::endl;
    FOUT<<"exp_power_imgs_size : "<<exp_power_imgs_size<<" GB."<<std::endl;
    FOUT<<"exp_highres_Xi2_imgs_size : "<<exp_highres_Xi2_imgs_size<<" GB."<<std::endl;
    FOUT<<"exp_min_diff2_size : "<<exp_min_diff2_size<<" GB."<<std::endl;
    FOUT<<"exp_old_offsetx_size : "<<exp_old_offsetx_size<<" GB."<<std::endl;
    FOUT<<"exp_old_offsety_size : "<<exp_old_offsety_size<<" GB."<<std::endl;
    FOUT<<"exp_wsum_scale_correction_XA_size : "<<exp_wsum_scale_correction_XA_size<<" GB."<<std::endl;
    FOUT<<"exp_wsum_scale_correction_AA_size : "<<exp_wsum_scale_correction_AA_size<<" GB."<<std::endl;
    FOUT<<"exp_wsum_norm_correction_size  	: "<<exp_wsum_norm_correction_size<<" GB."<<std::endl;
    FOUT<<"exp_Fimgs_size : "<<exp_Fimgs_size<<" GB."<<std::endl;
    FOUT<<"exp_Fimgs_nomask_size : "<<exp_Fimgs_nomask_size<<" GB."<<std::endl;
    FOUT<<"exp_Fctfs_size : "<<exp_Fctfs_size<<" GB."<<std::endl;
    FOUT<<"exp_significant_weight_size : "<<exp_significant_weight_size<<" GB."<<std::endl;
    FOUT<<"exp_sum_weight_size : "<<exp_sum_weight_size<<" GB."<<std::endl;
    FOUT<<"exp_Fimgs_shifted_nomask_size : "<<exp_Fimgs_shifted_nomask_size<<" GB."<<std::endl;
    FOUT<<"exp_local_sqrtXi2_size : "<<exp_local_sqrtXi2_size<<" GB."<<std::endl;
    FOUT<<"exp_over_psi_size : "<<exp_over_psi_size<<" GB."<<std::endl;
    FOUT<<"exp_over_rot_size : "<<exp_over_rot_size<<" GB."<<std::endl;
    FOUT<<"exp_over_tilt_size : "<<exp_over_tilt_size<<" GB."<<std::endl;
    FOUT<<"exp_over_trans_x : "<<exp_over_trans_x<<" GB."<<std::endl;
    FOUT<<"exp_over_trans_y : "<<exp_over_trans_y<<" GB."<<std::endl;
    FOUT<<"exp_Fimgs_shifted_size : "<<exp_Fimgs_shifted_size<<" GB."<<std::endl;
    FOUT<<"exp_local_Fctfs_size : "<<exp_local_Fctfs_size<<" GB."<<std::endl;
    FOUT<<"exp_local_Minvsigma2s_size : "<<exp_local_Minvsigma2s_size<<" GB."<<std::endl;
    FOUT<<"exp_Rot_significant_size : "<<exp_Rot_significant_size<<" GB."<<std::endl;
    FOUT<<"exp_Mcoarse_significant_size : "<<exp_Mcoarse_significant_size<<" GB."<<std::endl;
    FOUT<<"exp_Mweight_coarse_size : "<<exp_Mweight_coarse_size<<" GB."<<std::endl;
    FOUT<<"exp_Mweight_fine_size : "<<exp_Mweight_fine_size<<" GB."<<std::endl;
    FOUT<<"Total : "<<unfixed_size<<" GB."<<std::endl;
    FOUT<<"------------------------------------------------------"<<std::endl;
    //
    double thread_Frefctf_size 					=	2*maxthreads*current_Fsize2*sizeof(double)*perGb;
    double thread_Fimg_nomask_size				=	2*maxthreads*ori_size*ori_size*sizeof(double)*perGb;
    double thread_Fimg_size						=	2*maxthreads*ori_size*ori_size*sizeof(double)*perGb;
    double thread_Fweight_size					=	maxthreads*ori_size*ori_size*sizeof(double)*perGb;
    double thread_max_weight_size				=	maxthreads*nr_pool*sizeof(double)*perGb;
    double thread_wsum_norm_correction_size		=	maxthreads*nr_pool*sizeof(double)*perGb;
    double thread_wsum_sigma2_noise_size		=	maxthreads*nr_pool*current_Fsize2*sizeof(double)*perGb;
    double thread_wsum_pdf_direction_size		=	maxthreads*nr_classes*exp_nr_dir*sizeof(double)*perGb;
    double thread_wsum_pdf_class_size			=	maxthreads*nr_classes*sizeof(double)*perGb;
    double thread_wsum_scale_correction_XA_size	=	maxthreads*nr_pool*current_Fsize2*sizeof(double)*perGb;
    double thread_wsum_scale_correction_AA_size	=	maxthreads*nr_pool*current_Fsize2*sizeof(double)*perGb;
    double thread_sumw_group_size				=	maxthreads*nr_pool*sizeof(double)*perGb;
    double thread_size = thread_Frefctf_size + thread_Fimg_nomask_size + thread_Fimg_size + thread_Fweight_size \
    + thread_max_weight_size + thread_wsum_norm_correction_size + thread_wsum_sigma2_noise_size \
    + thread_wsum_pdf_direction_size + thread_wsum_pdf_class_size + thread_wsum_scale_correction_XA_size \
    + thread_wsum_scale_correction_AA_size + thread_sumw_group_size;
    FOUT<<"thread size(Expectation step) : "<<std::endl;
    FOUT<<"thread_Frefctf_size : "<<thread_Frefctf_size<<" GB."<<std::endl;
    FOUT<<"thread_Fimg_nomask_size : "<<thread_Fimg_nomask_size<<" GB."<<std::endl;
    FOUT<<"thread_Fimg_size : "<<thread_Fimg_size<<" GB."<<std::endl;
    FOUT<<"thread_Fweight_size : "<<thread_Fweight_size<<" GB."<<std::endl;
    FOUT<<"thread_max_weight_size : "<<thread_max_weight_size<<" GB."<<std::endl;
    FOUT<<"thread_wsum_norm_correction_size : "<<thread_wsum_norm_correction_size<<" GB."<<std::endl;
    FOUT<<"thread_wsum_sigma2_noise_size : "<<thread_wsum_sigma2_noise_size<<" GB."<<std::endl;
    FOUT<<"thread_wsum_pdf_direction_size : "<<thread_wsum_pdf_direction_size<<" GB."<<std::endl;
    FOUT<<"thread_wsum_pdf_class_size : "<<thread_wsum_pdf_class_size<<" GB."<<std::endl;
    FOUT<<"thread_wsum_scale_correction_XA_size : "<<thread_wsum_scale_correction_XA_size<<" GB."<<std::endl;
    FOUT<<"thread_wsum_scale_correction_AA_size : "<<thread_wsum_scale_correction_AA_size<<" GB."<<std::endl;
    FOUT<<"thread_sumw_group_size : "<<thread_sumw_group_size<<" GB."<<std::endl;
    FOUT<<"Total : "<<thread_size<<" GB."<<std::endl;
    FOUT<<"------------------------------------------------------"<<std::endl;
    
    ramFile.close();
}


void debugStoreWeightedSums(){
    
    std::cout << " WARNING: norm_correction : "<< exp_metadata[0].NORM  << " for particle " << 0 <<std::endl;
    std::cout << " mymodel.current_size : " << current_size << " mymodel.ori_size= " << ori_size <<std::endl;
    std::cout << " coarse_size : " << coarse_size << std::endl;
    std::cout << " DIRECT_A2D_ELEM(exp_metadata2, my_image_no, exp_nr_imagas-1) : " <<exp_metadata[exp_nr_images-1].NORM << std::endl;
    // std::cout << " mymodel.avg_norm_correction : " << model_avg_norm_correction << std::endl;
    std::cout << " exp_wsum_norm_correction[ipart] : " << exp_wsum_norm_correction.rptrAll()[0] << std::endl;
    // std::cout << " old_norm_correction : " << old_norm_correction << std::endl;
    // std::cout << " wsum_model.avg_norm_correction : " << wsum_avg_norm_correction << std::endl;
    // std::cout << " group_id : " << group_id << " mymodel.scale_correction[group_id] : " << mymodel.scale_correction[group_id] << std::endl;
    // std::cout << " mymodel.sigma2_noise[group_id = 0] : " << model_sigma2_noise[0] << std::endl;
    // std::cout << " wsum_model.sigma2_noise[group_id = 0] : " << wsum_sigma2_noise[0] << std::endl;
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
    // std::cerr << " ml_model.sigma2_noise[group_id = 0] : " << model_sigma2_noise[0] << std::endl;
}
} // end namespace MLoptimizer old
