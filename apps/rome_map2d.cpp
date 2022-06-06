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

#include "../src/map2d_optimizer.h"
#include "../src/time.h"
#include "../src/option.h"

int main(int argc, char * argv[]) {
    
    Option option;
    // general option
    option.addOption("-i",                  "Input metadata file with images needed to align"                               										);
    option.addOption("-o",                  "Output metadata"                                                               										);
    option.addOption("-K",                  "Number of classes needed to classify,!!!!! please use Uppercase —K instead of lowercase -k !!!!!" 						);
    option.addOption("-iter",               "Maximum number of iterations to perform"                                       										);
    option.addOption("-angpix",             "Pixel size (in Angstroms),!!!!! please use -angpix instead of -pixel !!!!!"   											);
    option.addOption("-tau2_fudge",         "Regularisation parameter (values higher than 1 give more weight to the data,2 for 2D)"                                 );
    option.addOption("-oversampling",       "Adaptive oversampling order to speed-up calculations (0=no oversampling, 1=2x, 2=4x, etc)"                             );
    // advanced option
    option.addOption("-particle_diameter",  "Particle_diameter",																							"-1"	);
    option.addOption("-pool",               "Number of images to be processed together for each EM step",              										"50" 	);
    option.addOption("-random_seed",        "Number of the random seed generator",                                     										"33" 	);
    option.addOption("-offset_step",        "Sampling rate (before oversampling) for origin offsets (in pixels)",      										"2"  	);
    option.addOption("-offset_range",       "Search range for origin offsets (in pixels)",                             										"10" 	);
    option.addOption("-psi_step",           "Sampling rate (before oversampling) for the in-plane angle",              										"10" 	);
    //
    option.addOption("-ctf",                "Perform CTF correction?",                                                  		                            "0" 	);
    option.addOption("-only_flip_phases",	"Only perform CTF phase-flipping? (default is full amplitude-correction)",										"0"		);
    option.addOption("-ctf_phase_flipped",	"Have the data been CTF phase-flipped?",																		"0"		);
    option.addOption("-ctf_corrected_ref", 	"Have the input references been CTF-amplitude corrected?",														"0"		);
    option.addOption("-intact_ctf_first_peak", "Ignore CTFs until their first peak?",																		"0"		);
    option.addOption("-norm",               "Do norm correction for image data",                                        									"0" 	);
    //
    option.addOption("-zero_mask", 			"Mask surrounding background in particles to zero (by default the solvent area is filled with random noise)",	"0"		);
    option.addOption("-flatten_solvent", 	"Perform masking on the references as well",																	"0"		);
    //
    option.addOption("-datastream_in", 		"data stream for comparing",																					"NULL"  );
    option.addOption("-datastream_out", 	"write data stream to file",																				    "NULL"  );
    option.addOption("-continue",           "continue work from specified iter.(remember to change -i input star)",    										"0"  	);
    option.addOption("-testguidance",       "keep each iteration same as guidefile",                                   										"NULL"	);
    
    if (argc < 3) {
        option.printHelp();
        std::cerr << "example:" << std::endl;
        std::cerr << "./bin/rome_map2d -i dataset/class19.star -o dataset/class19_ml2d -K 10 -iter 30 -angpix 1.74 > dataset/ml2d_output.txt" << std::endl;
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
    
    // set parameters
    // input file
    std::string star_fn 			      	=	pathRemoveSuffix(option.getStrOption("-i"))+".star";
    Map2dOptimizer_old::star_fn           	=   star_fn;
    Map2dOptimizer_old::mrcs_dir          	=   pathGetDir(star_fn);
    // output file
    std::string output_fn                  	=   pathRemoveSuffix(option.getStrOption("-o"));
    Map2dOptimizer_old::write_fn          	=   pathGetFilename(output_fn);
    Map2dOptimizer_old::write_path        	=   pathGetDir(output_fn);
    // basic
    Map2dOptimizer_old::nr_classes        	=   option.getIntOption("-K");
    Map2dOptimizer_old::nr_pool           	=   option.getIntOption("-pool");
    Map2dOptimizer_old::pixel_size        	=   option.getFloatOption("-angpix");
    Map2dOptimizer_old::nr_iter           	=   option.getIntOption("-iter");
    Map2dOptimizer_old::random_seed       	=   option.getIntOption("-random_seed");
    Map2dOptimizer_old::offset_step       	=   option.getIntOption("-offset_step");
    Map2dOptimizer_old::offset_range		=   option.getIntOption("-offset_range");
    Map2dOptimizer_old::psi_step          	=   option.getIntOption("-psi_step");
    Map2dOptimizer_old::tau2_fudge_factor   =   option.getFloatOption("-tau2_fudge");
    Map2dOptimizer_old::particle_diameter   =   option.getFloatOption("-particle_diameter");
    Map2dOptimizer_old::adaptive_oversampling=  option.getIntOption("-oversampling");
    // true or false
    Map2dOptimizer_old::do_ctf_correction 	=   option.getBoolOption("-ctf");
    Map2dOptimizer_old::only_flip_phases	=	option.getBoolOption("-only_flip_phases");
    Map2dOptimizer_old::ctf_phase_flipped	=	option.getBoolOption("-ctf_phase_flipped");
    Map2dOptimizer_old::refs_are_ctf_corrected= option.getBoolOption("-ctf_corrected_ref");
    Map2dOptimizer_old::intact_ctf_first_peak=	option.getBoolOption("-intact_ctf_first_peak");
    Map2dOptimizer_old::do_norm_correction	=   option.getBoolOption("-norm");
    Map2dOptimizer_old::do_zero_mask	  	=	option.getBoolOption("-zero_mask");
    Map2dOptimizer_old::do_solvent			=	option.getBoolOption("-flatten_solvent");
    // others
    Map2dOptimizer_old::data_stream_in_fn	=	option.getStrOption("-datastream_in");
    Map2dOptimizer_old::data_stream_out_fn	=	option.getStrOption("-datastream_out");
    Map2dOptimizer_old::iter              	=   option.getIntOption("-continue");
    Map2dOptimizer_old::guidance_fn       	=   option.getStrOption("-testguidance");
    Map2dOptimizer_old::ref_fn            	=	"NULL";
    
    //
    Map2dOptimizer_old::setupMap2dOptimizer();
    
    Map2dOptimizer_old::readImages();
    
    Map2dOptimizer_old::prepare();
    
    Map2dOptimizer_old::iterate();
    
    Map2dOptimizer_old::destroyMap2dOptimizer();

    
    double t2 = dtime();
    
    if(Map2dOptimizer_old::node  == 0)
        std::cout<<"ML2D costs : "<<(t2-t1)<<std::endl;
    
#ifdef USEMPI
    MPI::Finalize();
#endif
    
    return 0;
}
