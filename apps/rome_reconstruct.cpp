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

#include "../src/option.h"
#include "../src/metadata.h"
#include "../src/mrcs.h"
#include "../src/reconstruct.h"
#include "../src/ctf.h"
#include "../src/debug.h"
#include "../src/mpi.h"
#include "../src/progressbar.h"

namespace Reconstruct{
    //
    std::string fn_mrc,fn_star,fn_fsc,fn_sub,fn_sym,fn_tmp;
    // Kaiser–Bessel function
    double blob_radius,blob_alpha,blob_order;
    //
    double pixel_size,maxres,angular_error,shift_error;
    //
    int nr_iters,padding_factor,r_max,r_min_nn,ref_dim,interpolator;
    //
    bool do_NN,do_ctf,ctf_phase_flipped,only_flip_phases,ctf_intact_first_peak;
    //
    MetaDataTable metadata;
    // image data
    int ori_size;
    // MPI
    int nr_global_images,nr_local_images,first_local_image,last_local_image;
    int nodes,node;
#ifdef USEMPI
    int nodeNameLen;
    char nodeName[MPI_MAX_PROCESSOR_NAME];
#endif
    int maxthreads;
    //
    void setupReconstruct()
    {
#ifdef USEMPI
        // MPI::Init();//init mpi outside
        nodes = MPI::COMM_WORLD.Get_size();
        node = MPI::COMM_WORLD.Get_rank();
        MPI::Get_processor_name(nodeName,nodeNameLen);
        NODE0ONLY std::cout<<nodeName<<" "<<node<<" "<<nodes<<std::flush<<std::endl;
#else
        nodes = 1;
        node = 0;
#endif

        // Read MetaData file, which should have the image names and their angles!
        metadata.readFromStar(fn_star);
        
        // TODO 2D
        assert(ref_dim==3);
        assert(metadata.containLabel(AngleRot)==true);
        assert(metadata.containLabel(AngleTilt)==true);
        assert(metadata.containLabel(AnglePsi)==true);
        assert(metadata.containLabel(OriginX)==true);
        assert(metadata.containLabel(OriginY)==true);
        
        nr_global_images = metadata.numberOfParticles();
        
        nr_local_images = divide_equally(nr_global_images,nodes,node, first_local_image, last_local_image);
        
        maxthreads = omp_get_max_threads();
        NODE0ONLY std::cout<<"maxthreads = "<<maxthreads<<",node = "<<node<<std::endl;
        NODE0ONLY std::cout<<"first_local_image = "<<first_local_image<<",last_local_image = "<<last_local_image<<std::endl;
        NODE0ONLY std::cout<<"nr_local_images = "<<nr_local_images<<std::endl;
        // Get dimension of the images

        Mrcs::MrcsHead mrcsHead;
        Mrcs::readMrcsHead(metadata[0].IMAGE.NAME, mrcsHead);
        ori_size = mrcsHead.NC;
        
        // srand(static_cast <unsigned> (time(NULL)) );
        
        dontShare_Random_generator.init(33);
        
        assert(do_NN==false);
        if (do_NN)
            interpolator = NEAREST_NEIGHBOUR;
        else
            interpolator = TRILINEAR;
    }
    
    //
    void reconstruct()
    {
        double rot, tilt, psi, trans_x, trans_y;
        FDOUBLE A3D[3][3];
        // Projector proj;
        int ori_Fsize = ori_size/2+1;
        int ori_Fsize2 = ori_size*ori_Fsize;
        
        if (metadata.containLabel(BeamTiltX) || metadata.containLabel(BeamTiltY)) {
            NODE0ONLY std::cerr<<"warning,not support beam.....todo..."<<std::endl;
        }
        
        if (maxres < 0.)
            r_max = ori_size/2; // default is Nyquist
        else
            r_max = std::ceil(ori_size * pixel_size / maxres);
        
        double memory_requirement = (ori_size*padding_factor*ori_size*padding_factor*ori_size*padding_factor*8./1024./1024./1024.)*2 + // projector.data && backprojector.data
        							(ori_size*padding_factor*ori_size*padding_factor*ori_size*padding_factor/2.*8./1024./1024./1024.); // backprojector.weight
        NODE0ONLY std::cerr << " Starting Back-projecting all images ..." << std::endl;
        NODE0ONLY std::cerr << " memory requirement : "<<memory_requirement<<" GB."<<std::endl;
        NODE0ONLY std::cerr << " Make sure your system has enough memory...."<<std::endl;
        
        MyBackProjector backprojector;
        backprojector.initialize(ori_size, 2*r_max, ref_dim, "C1", padding_factor, interpolator, r_min_nn, blob_order,blob_radius, blob_alpha);
        
        MyProjector projector;
        projector.initialize(ori_size, 2*r_max, ref_dim, padding_factor, interpolator, r_min_nn);
        
        if (fn_sub != "")
        {
            // sub.read(fn_sub);
            // projector.computeFourierTransformMap(sub(), dummy, 2 * r_max);
        }
#include "../src/util_heap_undefs.h"
        // initialize all threads data
        FFTWTransformer* transformer = new FFTWTransformer(ori_size,ori_size);
#include "../src/util_heap_defs.h"
        //
        CTF ctf;
        Aligned1dArray<FDOUBLE> Fctf(ori_size*(ori_size/2+1));
        Fctf.fill(1);
        
        //
        Aligned1dArray<float> floatBuffer(ori_size*ori_size);
        Aligned1dArray<double> originalImage(ori_size*ori_size);
        Aligned1dArray<double> centeredImage(ori_size*ori_size);
        Aligned1dArray<double> originalF2dReal(ori_Fsize2);
        Aligned1dArray<double> originalF2dImag(ori_Fsize2);
        Aligned1dArray<FDOUBLE> shiftedF2dReal(ori_Fsize2);
        Aligned1dArray<FDOUBLE> shiftedF2dImag(ori_Fsize2);
        //
        NODE0ONLY std::cout<<"ori_size = "<<ori_size<<std::endl;
        
        TIMEPOINT_INIT
        TIMEPOINT
        
        std::cout.precision(30);
        double read_data_time = 0,operate_time = 0;
		NODE0ONLY showProgressBar(0, nr_local_images);
        for (int iimage = first_local_image; iimage <= last_local_image; iimage++)
        {
            double t1 = dtime();
            
            auto floatBuffer_ptr = floatBuffer.wptr(ori_size*ori_size);
            auto originalImage_ptr = originalImage.wptr(ori_size*ori_size);
            auto centeredImage_ptr = centeredImage.wptr(ori_size*ori_size);
            auto originalF2dReal_ptr = originalF2dReal.wptr(ori_size*ori_size);
            auto originalF2dImag_ptr = originalF2dImag.wptr(ori_size*ori_size);
            auto shiftedF2dReal_ptr = shiftedF2dReal.wptr(ori_size*ori_size);
            auto shiftedF2dImag_ptr = shiftedF2dImag.wptr(ori_size*ori_size);
            auto Fctf_ptr = Fctf.wptr(ori_size*(ori_size/2+1));
            
            FILE* filehandle = fopen((metadata[iimage].IMAGE.NAME).c_str(),"rb");
            size_t image_index = metadata[iimage].IMAGE.INDEX;
            size_t offset = (256+(image_index-1)*ori_size*ori_size)*sizeof(float);
                        
            fseek(filehandle,offset,SEEK_SET);
            
            ERROR_CHECK(fread((char*)floatBuffer_ptr,ori_size*ori_size*sizeof(float),1,filehandle) == NULL, "read file failed.");
            
            fclose(filehandle);
            
            double t2 = dtime();
            
            read_data_time += (t2-t1);
            
            t1 = dtime();
            
            for (int i = 0; i < ori_size*ori_size; i++) {
                originalImage_ptr[i] = floatBuffer_ptr[i];
            }
	
            // NODE0ONLY std::cout<<sumvec(originalImage_ptr, ori_size*ori_size)<<std::endl;
            
            // Rotations
            if (ref_dim==2)
            {
                // TODO
                // rot = tilt = 0.;
            }
            else
            {
                rot = metadata[iimage].ROT;
                tilt = metadata[iimage].TILT;
            }
            psi = metadata[iimage].PSI;
            
            if (angular_error > 0.)
            {
                rot += dontShare_Random_generator.rnd_gaus(0., angular_error);
                tilt += dontShare_Random_generator.rnd_gaus(0., angular_error);
                psi += dontShare_Random_generator.rnd_gaus(0., angular_error);
                // NODE0ONLY std::cout << rnd_gaus(0., angular_error) << std::endl;
            }
            
            Euler_angles2matrix(rot, tilt, psi, A3D);
            
            // Translations (either through phase-shifts or in real space
            trans_x = metadata[iimage].XOFF;
            trans_y = metadata[iimage].YOFF;
            if (shift_error > 0.)
            {
                trans_x += dontShare_Random_generator.rnd_gaus(0., shift_error);
                trans_y += dontShare_Random_generator.rnd_gaus(0., shift_error);
            }
            
            // Use either selfTranslate OR shiftImageInFourierTransform!!
            //selfTranslate(img(), trans, WRAP);
            centerFFT2D(originalImage_ptr, centeredImage_ptr, ori_size, true);
            // NODE0ONLY std::cout<<sumvec(centeredImage_ptr, ori_size*ori_size)<<std::endl;
            
            transformer->FourierTransform(centeredImage_ptr, originalF2dReal_ptr, originalF2dImag_ptr);
            
            // NODE0ONLY std::cout<<sumvec(originalF2dReal_ptr,ori_Fsize2)<<" "<<sumVec(originalF2dImag_ptr,ori_Fsize2)<<std::endl;
            // NODE0ONLY std::cout<<sumvec(shiftedF2dReal_ptr,ori_Fsize2)<<" "<<sumVec(shiftedF2dImag_ptr,ori_Fsize2)<<std::endl;
            
            if (fabs(trans_x) > 0. || fabs(trans_y) > 0.)
            {
                shiftImageInFourierTransform(originalF2dReal_ptr, originalF2dImag_ptr, shiftedF2dReal_ptr, shiftedF2dImag_ptr,
                                             ori_size, trans_x, trans_y, ori_size);
            }
            
            // NODE0ONLY std::cout<<sumvec(shiftedF2dReal_ptr,ori_Fsize2)<<" "<<sumVec(shiftedF2dImag_ptr,ori_Fsize2)<<std::endl;
            
            // Apply CTF if necessary
            if (do_ctf)
            {
                ctf.setValues(	metadata[iimage].CTF_DEFOCUS_U,
                                metadata[iimage].CTF_DEFOCUS_V,
                                metadata[iimage].CTF_DEFOCUS_ANGLE,
                                metadata[iimage].CTF_VOLTAGE,
                                metadata[iimage].CTF_CS,
                                metadata[iimage].CTF_Q0,
                                metadata[iimage].CTF_BFAC);
                ctf.getFftwImage(Fctf_ptr, ori_size, ori_size, ori_size, pixel_size, ctf_phase_flipped, only_flip_phases, ctf_intact_first_peak, true);
                
                // NODE0ONLY std::cout<<sumvec(Fctf_ptr, ori_Fsize2)<<std::endl;
            }
            
            // Subtract reference projection
            if (fn_sub != "")
            {
                // TODO
                // Fsub.resize(F2D);
                // projector.get2DFourierTransform(Fsub, A3D, IS_NOT_INV);
                
                // Apply CTF if necessary
                // if (do_ctf)
                // {
                //     FOR_ALL_DIRECT_ELEMENTS_IN_MULTIDIMARRAY(Fsub)
                //     {
                //         DIRECT_MULTIDIM_ELEM(Fsub, n) *= DIRECT_MULTIDIM_ELEM(Fctf, n);
                //     }
                // }
                
                // FOR_ALL_DIRECT_ELEMENTS_IN_MULTIDIMARRAY(Fsub)
                // {
                //     DIRECT_MULTIDIM_ELEM(F2D, n) -= DIRECT_MULTIDIM_ELEM(Fsub, n);
                // }
                // Back-project difference image
                // backprojector.set2DFourierTransform(F2D, A3D, IS_NOT_INV);
            }
            else
            {
                // "Normal" reconstruction, multiply X by CTF, and W by CTF^2
                if (do_ctf)
                {
                    for (int n = 0; n < ori_Fsize2; n++) {
                        shiftedF2dReal_ptr[n] *= Fctf_ptr[n];
                        shiftedF2dImag_ptr[n] *= Fctf_ptr[n];
                        Fctf_ptr[n] *= Fctf_ptr[n];
                    }
                }
                
                //#define DEBUG_RECONSTRUCT_ONLY
#ifdef DEBUG_RECONSTRUCT_ONLY
                if (fn_img == "/lmb/home/scheres/data/betaGal_rh_withnoise/betaGal_2010_all_p_2x2_unflipped/img00001.win100")
                    //if (part_id == my_first_particle_id)
                {
                    std::cerr << " fn_img= " << fn_img << std::endl;
                    //std::cerr << " myscale= " << myscale << std::endl;
                    //std::cerr << " mymodel.avg_norm_correction= " << mymodel.avg_norm_correction << " normcorr= " << normcorr << std::endl;
                    //std::cerr << " sigma2_fudge= " << sigma2_fudge << " mymodel.tau2_fudge_factor= " << mymodel.tau2_fudge_factor<< std::endl;
                    //std::cerr << " A3D= " << A3D << std::endl;
                    std::cerr << " A3D= " << A3D << std::endl;
                    //std::cerr << " exp_R_mic= " << exp_R_mic << std::endl;
                    std::cerr << " rot= " << rot << " tilt= " << tilt << " psi= " << psi << " xoff= "<< XX(trans)<< " yoff= "<<YY(trans)<<std::endl;
                    //std::cerr << "mic_id= "<<mic_id<<" mymodel.sigma2_noise[mic_id]= " << mymodel.sigma2_noise[mic_id] << std::endl;
                    Image<DOUBLE> It;
                    It()=Fctf;
                    It.write("reconstruct_Fctf.spi");
                    It().resize(mysize, mysize);
                    MultidimArray<Complex > Faux = F2D;
                    FourierTransformer transformer;
                    transformer.inverseFourierTransform(Faux, It());
                    CenterFFT(It(), false);
                    It.write("reconstruct_Mimg.spi");
                }
#endif
            	// NODE0ONLY std::cout<<sumvec(shiftedF2dReal_ptr,shiftedF2dReal_ptr,ori_Fsize2)<<std::endl;
            	// NODE0ONLY std::cout<<sumVec(Fctf_ptr, ori_Fsize2)<<std::endl;
                backprojector.backproject(shiftedF2dReal_ptr, shiftedF2dImag_ptr, ori_size, A3D, IS_NOT_INV,
                                          backprojector.data, backprojector.weight, Fctf_ptr);
                
            }
            t2 = dtime();
            // operate_time += (t2-t1);
            // NODE0ONLY std::cout<<"iimage : "<<iimage<<", read_data_time : "<<read_data_time<<", operate_time : "<<operate_time<<std::endl;
            NODE0ONLY showProgressBar(iimage, nr_local_images);
        }
        //
        NODE0ONLY showProgressBar(nr_local_images, nr_local_images);
        //
        projector		.finalize();
        Fctf			.fini();
        floatBuffer		.fini();
        originalImage	.fini();
        centeredImage	.fini();
        originalF2dReal	.fini();	originalF2dImag	.fini();
        shiftedF2dReal	.fini();	shiftedF2dImag	.fini();
        
#include "../src/util_heap_undefs.h"
        delete transformer;
#include "../src/util_heap_defs.h"
        
        TIMEPOINT
        TIMEPOINT_FINA
        
        bool do_map = false;
        bool do_use_fsc = false;
        std::vector<FDOUBLE> fsc(ori_Fsize);
        if (fn_fsc != "")
        {
            // TODO
            // do_map = true;
            // do_use_fsc =true;
            // MetaDataTable MDfsc;
            // MDfsc.read(fn_fsc);
            // FOR_ALL_OBJECTS_IN_METADATA_TABLE(MDfsc)
            // {
            //     int idx;
            //     DOUBLE val;
            //     MDfsc.getValue(EMDL_SPECTRAL_IDX, idx);
            //     MDfsc.getValue(EMDL_MLMODEL_FSC_HALVES_REF, val);
            //     fsc(idx) =  val;
            // }
        }
        // reduce data
        int local_temp_size = backprojector.data.dimzyx*2;
        FDOUBLE* local_temp = (FDOUBLE*)aMalloc(sizeof(FDOUBLE)*local_temp_size,64);
        {
            memcpy(local_temp, (FDOUBLE*)backprojector.data.wptr(), sizeof(FDOUBLE)*local_temp_size);
            MPI::COMM_WORLD.Allreduce(local_temp,(FDOUBLE*)backprojector.data.wptr(),local_temp_size,MPI_FDOUBLE,MPI::SUM);
            memcpy(local_temp, backprojector.weight.wptr(), sizeof(FDOUBLE)*local_temp_size/2);
            MPI::COMM_WORLD.Allreduce(local_temp,backprojector.weight.wptr(),local_temp_size/2,MPI_FDOUBLE,MPI::SUM);
        }
        aFree(local_temp);
        
		//
       	NODE0ONLY{
            memory_requirement = (ori_size*padding_factor*ori_size*padding_factor*ori_size*padding_factor*8./1024./1024./1024.) + // backprojector.data
                                 (ori_size*padding_factor*ori_size*padding_factor*ori_size*padding_factor/2.*8./1024./1024./1024.) + // backprojector.weight
                                 (ori_size*padding_factor*ori_size*padding_factor*ori_size*padding_factor*8./1024./1024./1024.) + // reconstruct.Fconv
                                 (ori_size*padding_factor*ori_size*padding_factor*ori_size*padding_factor*8./1024./1024./1024.) + // reconstruct.Fweight,Fnewweight
                                 (ori_size*padding_factor*ori_size*padding_factor*ori_size*padding_factor*8./1024./1024./1024.*3/2);  // tmp data
            std::cout << " Starting the reconstruction..." << std::endl;
            std::cout << " For "<<ori_size<<"x"<<ori_size<<"x"<<ori_size<<" and pad=2 volume..."<<std::endl;
            std::cout << " Memory requirement : "<<memory_requirement*1.2<<" GB."<<std::endl;
            std::cout << " Make sure your system has enough memory...."<<std::endl;
            
            Vol<FDOUBLE> vol_out;
            vol_out.init(ori_size, ori_size, ori_size);
            std::vector<FDOUBLE> dummy(ori_Fsize);
            backprojector.reconstruct(vol_out, nr_iters, do_map, 1., &dummy[0], &dummy[0], &dummy[0], &fsc[0], 1., do_use_fsc, true, -1, maxthreads, fn_tmp);
            
            backprojector.finalize();
            
            Mrcs::MrcVolume oneVol(vol_out.wptr(),ori_size,pixel_size);
            //
            Mrcs::MrcsHead refHead;
            Mrcs::setMrcHead(refHead, pixel_size, ori_size);
            //
            oneVol.write(fn_mrc,refHead);
            
            std::cerr<<" Done writing map in "<<fn_mrc<<std::endl;
        }
        // free backprojector to decrease memory usage...
        if (node != 0) {
            backprojector.finalize();
        }
    }
    
    void destroyReconstruct()
    {
        //;
    }
};

int main(int argc, char * argv[]){
    
    Option option;
    // general option
    option.addOption("-i",                  "Input STAR file with the projection images and their orientations"                                                     );
    option.addOption("-o",                  "Name for output reconstruction(*.mrc)"                                                                                 );
    option.addOption("-angpix",             "Pixel size (in Angstroms)!!!!! please use -angpix instead of -pixel !!!!!"												);
    // advanced option
    option.addOption("-iter",               "Number of gridding-correction iterations",                                                                        "10" );
    option.addOption("-refdim",             "Dimension of the reconstruction (2D or 3D)",                                                                      "3"  );
    option.addOption("-sym",                "Symmetry group(TODO)",                                                                                            "c1" );
    option.addOption("-ctf", 				"Apply CTF correction",																							   "0"	);
    option.addOption("-ctf_intact_first_peak","Leave CTFs intact until first peak",																		   	   "0"	);
    option.addOption("-ctf_phase_flipped", 	"Images have been phase flipped",																				   "0"	);
    option.addOption("-only_flip_phases", 	"Do not correct CTF-amplitudes, only flip phases",																   "0"	);
    option.addOption("-maxres",             "Maximum resolution (in Angstrom) to consider in Fourier space (default Nyquist)",                                 "-1" );
    option.addOption("-subtract",           "Subtract projections of this map from the images used for reconstruction",                                        ""   );
    option.addOption("-pad",                "Padding factor",                                                                                                  "2"  );
    option.addOption("-r_min_nn",           "Minimum number of Fourier shells to perform linear Fourier-space interpolation",                                  "1"  );
    option.addOption("-NN",                 "Use nearest-neighbour instead of linear interpolation before gridding correction",                                "0"  );
    option.addOption("-blob_r",             "Radius of blob for gridding interpolation",                                                                       "1.9");
    option.addOption("-blob_m",             "Order of blob for gridding interpolation",                                                                        "0"  );
    option.addOption("-blob_a",             "Alpha-value of blob for gridding interpolation",                                                                  "15" );
    option.addOption("-angular_error",      "Apply random deviations with this standard deviation (in degrees) to each of the 3 Euler angles",                 "0"  );
    option.addOption("-shift_error",        "Apply random deviations with this standard deviation (in pixels) to each of the 2 translations",                  "0"  );
    option.addOption("-fsc",                "FSC-curve for regularized reconstruction",                                                                        "0"  );
    option.addOption("-tmp", 				"floder for temporary data,for small memory user.",															   	  "NULL");
    option.addOption("-refdim", 			"Dimension of the reconstruction (2D(TODO) or 3D)",																   "3"	);
    if (argc < 3) {
        option.printHelp();
        std::cerr << "example:" << std::endl;
        std::cerr << "../bin/rome_reconstruct -i class19.star -o class19_recon -angpix 1.72 "<<std::endl;
        std::cerr << std::endl;
        EXIT_ABNORMALLY;
    }
    
    option.readCommandLine(argc, argv);
    
#ifdef USEMPI
    MPI::Init();
    if(MPI::COMM_WORLD.Get_rank()  == 0)
        option.printValue();
#else
    option.printValue();
#endif
    
    // set parameters
    Reconstruct::fn_star        =   pathRemoveSuffix(option.getStrOption("-i"))+".star";
    Reconstruct::fn_mrc         =   pathRemoveSuffix(option.getStrOption("-o"))+".mrc";
    Reconstruct::fn_fsc         =   option.getStrOption("-fsc");
    Reconstruct::fn_sub         =   option.getStrOption("-subtract");
    Reconstruct::fn_sym         =   option.getStrOption("-sym");
    Reconstruct::fn_tmp			=	option.getStrOption("-tmp");
    
    Reconstruct::do_ctf			=	option.getBoolOption("-ctf");
    Reconstruct::ctf_intact_first_peak = option.getBoolOption("-ctf_intact_first_peak");
    Reconstruct::ctf_phase_flipped	= option.getBoolOption("-ctf_phase_flipped");
    Reconstruct::only_flip_phases = option.getBoolOption("-only_flip_phases");
    
    Reconstruct::pixel_size     =   option.getFloatOption("-angpix");
    Reconstruct::maxres         =   option.getFloatOption("-maxres");
    Reconstruct::blob_radius    =   option.getFloatOption("-blob_r");
    Reconstruct::blob_order     =   option.getFloatOption("-blob_m");
    Reconstruct::blob_alpha     =   option.getFloatOption("-blob_a");
    Reconstruct::angular_error  =   option.getFloatOption("-angular_error");
    Reconstruct::shift_error    =   option.getFloatOption("-shift_error");
    
    Reconstruct::nr_iters       =   option.getIntOption("-iter");
    Reconstruct::padding_factor =   option.getIntOption("-pad");
    Reconstruct::r_min_nn       =   option.getIntOption("-r_min_nn");
    Reconstruct::ref_dim        =   option.getIntOption("-refdim");
    
    Reconstruct::do_NN          =   option.getBoolOption("-NN");
    
    double t1 = dtime();
    //
    Reconstruct::setupReconstruct();
    Reconstruct::reconstruct();
    Reconstruct::destroyReconstruct();
    
    double t2 = dtime();
    if(Reconstruct::node  == 0)
    	std::cout<<"(reconstruction)Time cost : "<<(t2-t1)<<" seconds."<<std::endl;
    
#ifdef USEMPI
    MPI::Finalize();
#endif
    
    return 0;
}
