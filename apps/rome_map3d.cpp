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
#include <cstdio>
#include <iostream>

#include "../src/map_model.h"
#include "../src/time.h"
#include "../src/option.h"
#include "../src/map3d_optimizer.h"

int main(int argc, char * argv[]) {
        
    Option option;
    // general option
    option.addOption("-i",                  "Input metadata file with images needed to align"                                                                       );
    option.addOption("-o",                  "Output metadata"                                                                                                       );
    option.addOption("-ref",                "3d reference file name(*.mrc)"                                                                                         );
    option.addOption("-particle_diameter",  "Particle_diameter"                                                                                                     );
    option.addOption("-K",                  "Number of classes needed to classify!!!!! please use Uppercase —K instead of lowercase -k !!!!!" 						);
    option.addOption("-iter",               "Maximum number of iterations to perform"                                                                               );
    option.addOption("-angpix",             "Pixel size (in Angstroms)!!!!! please use -angpix instead of -pixel !!!!!"												);
    option.addOption("-ini_high",           "Resolution (in Angstroms) to which to limit refinement in the first iteration"                                         );
    option.addOption("-oversampling",       "Adaptive oversampling order to speed-up calculations (0=no oversampling, 1=2x, 2=4x, etc)"                             );
    option.addOption("-healpix_order",      "Healpix order for the angular sampling (before oversampling) on the (3D) sphere: hp2=15deg, hp3=7.5deg, etc"           );
    option.addOption("-tau2_fudge",         "Regularisation parameter (values higher than 1 give more weight to the data,4 for 3D)"                                 );
    option.addOption("-sym",                "Symmetry group"						                                                                                );
    // advanced option
    option.addOption("-pool",               "Particles processing together each EM expectation.",                                                              "50" );
    option.addOption("-random_seed",        "Seed of randomly shuffle the image data and noise generation",                                                    "33" );
    option.addOption("-offset_step",        "The offset step of image shift searching",                                                                        "2"  );
    option.addOption("-offset_range",       "The offset range of image shift searching,(-10~10)",                                                              "10" );
    //
    option.addOption("-ctf",                "Perform CTF correction?",                                                  		                               "0" 	);
    option.addOption("-only_flip_phases",	"Only perform CTF phase-flipping? (default is full amplitude-correction)",										   "0"	);
    option.addOption("-ctf_phase_flipped",	"Have the data been CTF phase-flipped?",																		   "0"	);
    option.addOption("-ctf_corrected_ref", 	"Have the input references been CTF-amplitude corrected?",														   "0"	);
    option.addOption("-intact_ctf_first_peak", "Ignore CTFs until their first peak?",																		   "0"	);
    option.addOption("-firstiter_cc",       "Perform CC-calculation in the first iteration (use this if references are not on the absolute intensity scale)",  "0"  );
    option.addOption("-always_cc",          "Perform CC-calculation in all iterations (useful for faster denovo model generation?)",                           "0"  );
    option.addOption("-scale",              "Perform intensity-scale corrections on image groups?",                                                            "0"  );
    option.addOption("-norm",               "Perform normalisation-error correction?",                                                                         "0"  );
    option.addOption("-zero_mask",          "Mask surrounding background in particles to zero (by default the solvent area is filled with random noise)",      "0"  );
    option.addOption("-flatten_solvent",    "Perform masking on the references as well?",                                                                      "0"  );
    option.addOption("-solvent_mask", 		"User-provided mask for the references (default is to use spherical mask with particle_diameter)",				  "NULL");
    option.addOption("-no_map",             "Do not use map estimitor",                                                                                        "0"  );
    //
    option.addOption("-datastream_in", 		"data stream for comparing",																					  "NULL");
    option.addOption("-datastream_out", 	"write data stream to file",																					  "NULL");
    option.addOption("-continue",        	"continue work from specified file.( from *_optimiser.star or *_backup.back file)",    						   	  "NULL");
    option.addOption("-testguidance",       "keep each iteration same as guidefile(use 'BACK' to backup the data)",     									  "NULL");
    //
    if (argc < 3) {
        option.printHelp();
        std::cerr << "example:" << std::endl;
        std::cerr << "./bin/rome_map3d -i dataset/class19.star -o dataset/class19_ml2d -K 10 -iter 30 -angpix 1.74 > dataset/ml2d_output.txt" << std::endl;
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
    
    double t1 = dtime();
    
    std::string set_star_fn         =   pathRemoveSuffix(option.getStrOption("-i"))+".star";
    std::string output_fn           =   pathRemoveSuffix(option.getStrOption("-o"));
    
    // set parameters
    Map3dOptimizer_old::star_fn                	=   set_star_fn;
    Map3dOptimizer_old::write_path             	=   pathGetDir(output_fn);
    Map3dOptimizer_old::write_fn               	=   pathGetFilename(output_fn);
    Map3dOptimizer_old::ref_fn                 	=   option.getStrOption("-ref");
    Map3dOptimizer_old::nr_classes             	=   option.getIntOption("-K");
    Map3dOptimizer_old::nr_pool                	=   option.getIntOption("-pool");
    Map3dOptimizer_old::pixel_size             	=   option.getFloatOption("-angpix");
    Map3dOptimizer_old::nr_iter                	=   option.getIntOption("-iter");
    Map3dOptimizer_old::random_seed            	=   option.getIntOption("-random_seed");
    Map3dOptimizer_old::offset_step            	=   option.getIntOption("-offset_step");
    Map3dOptimizer_old::offset_range           	=   option.getIntOption("-offset_range");
    //
    Map3dOptimizer_old::do_ctf_correction 	   	=   option.getBoolOption("-ctf");
    Map3dOptimizer_old::only_flip_phases		=   option.getBoolOption("-only_flip_phases");
    Map3dOptimizer_old::ctf_phase_flipped		=   option.getBoolOption("-ctf_phase_flipped");
    Map3dOptimizer_old::refs_are_ctf_corrected 	=   option.getBoolOption("-ctf_corrected_ref");
    Map3dOptimizer_old::intact_ctf_first_peak	=   option.getBoolOption("-intact_ctf_first_peak");
    //
    Map3dOptimizer_old::ini_high               	=   option.getFloatOption("-ini_high");
    Map3dOptimizer_old::adaptive_oversampling  	=   option.getIntOption("-oversampling");
    Map3dOptimizer_old::sampler3d_healpix_order =   option.getIntOption("-healpix_order");
    Map3dOptimizer_old::tau2_fudge_factor      	=   option.getFloatOption("-tau2_fudge");
    Map3dOptimizer_old::sampler3d_fn_sym        =   option.getStrOption("-sym");
    Map3dOptimizer_old::particle_diameter      	=   option.getFloatOption("-particle_diameter");
    Map3dOptimizer_old::do_firstiter_cc         =   option.getBoolOption("-firstiter_cc");
    Map3dOptimizer_old::do_always_cc           	=   option.getBoolOption("-always_cc");
    Map3dOptimizer_old::do_scale_correction    	=   option.getBoolOption("-scale");
    Map3dOptimizer_old::do_norm_correction     	=   option.getBoolOption("-norm");
    Map3dOptimizer_old::do_zero_mask           	=   option.getBoolOption("-zero_mask");
    Map3dOptimizer_old::do_solvent             	=   option.getBoolOption("-flatten_solvent");
    Map3dOptimizer_old::do_map                 	=   !option.getBoolOption("-no_map");
    Map3dOptimizer_old::mask_fn					=	option.getStrOption("-solvent_mask");
    
    Map3dOptimizer_old::data_stream_in_fn		=	option.getStrOption("-datastream_in");
    Map3dOptimizer_old::data_stream_out_fn		=	option.getStrOption("-datastream_out");
    Map3dOptimizer_old::continue_fn				=	option.getStrOption("-continue");
    Map3dOptimizer_old::guidance_fn       		=   option.getStrOption("-testguidance");
    //
    Map3dOptimizer_old::setupMLoptimizer();
    
    Map3dOptimizer_old::prepare();
    
    Map3dOptimizer_old::iterate();
    
    Map3dOptimizer_old::destroyMLoptimizer();
    
    
    double t2 = dtime();
    
    if(Map3dOptimizer_old::node  == 0)
        std::cout<<"ML2D costs : "<<(t2-t1)<<std::endl;
    
#ifdef USEMPI
    MPI::Finalize();
#endif
    
    return 0;
}
