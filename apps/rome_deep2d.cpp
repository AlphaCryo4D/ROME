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

#include "../src/gtm_optimizer.h"
#include "../src/option.h"
#include "../src/map2d_optimizer.h"

/*
void MLProgram() {
    
    Map2dOptimizer::setupMap2dOptimizer();
    Map2dOptimizer::readImages();
    Map2dOptimizer::prepare();
    Map2dOptimizer::iterate();
    Map2dOptimizer::writeClassesAndMetadata();
    Map2dOptimizer::destroyMap2dOptimizer();
    
}*/



int main(int argc, char * argv[]) {
    
    //---------- read input parameters  ----------//
    
    Option option;
    // general option
    option.addOption("-i"					,"Input metadata file with images needed to align"                               	);
    option.addOption("-o"					,"Output metadata"                                                               	);
    option.addOption("-map2d_K"				,"Number of classes needed to classify based on maximum likelihood method!!!!! please use Uppercase —map2d_K instead of lowercase -map2d_k !!!!!");
    option.addOption("-sml_K"				,"Number of classes needed to classify based on GTM!!!!! please use Uppercase —sml_K instead of lowercase -gtm_k !!!!!");
    option.addOption("-map2d_iter"			,"Maximum number of iterations to perform based on maximum likelihood method"    	);
    option.addOption("-sml_iter"			,"Maximum number of iterations to perform based on GTM"                          	);
    option.addOption("-angpix"				,"Pixel size (in Angstroms)!!!!! please use -angpix instead of -pixel !!!!!" 		);
    // add some options
    option.addOption("-tau2_fudge",         "Regularisation parameter (values higher than 1 give more weight to the data,2 for 2D)"                                 );
    option.addOption("-particle_diameter"   ,"Particle diameter"                                                        ,"-1"   );
    // advanced option
    option.addOption("-dimension"			,"GTM sampling grid dimension(Default 1)"									,"1"  	);
    option.addOption("-pca"					,"Using PCA to initialize,(Defalut is using Gauss-Circle template)"			,"0"  	);
    // option.addOption("-continue"			,"continue work from specified iter.(remember to change -i input star)"		,"0"	);
    // option.addOption("-testguidance"		,"keep each iteration same as guidefile"									,"NULL"	);
    option.addOption("-pool"				,"Number of images to be processed together for each EM step"				,"50" 	);
    option.addOption("-random_seed"			,"Number of the random seed generator"										,"33" 	);
    option.addOption("-offset_step"			,"Sampling rate (before oversampling) for origin offsets (in pixels)"		,"2"  	);
    option.addOption("-offset_range"		,"Search range for origin offsets (in pixels)"								,"10" 	);
    option.addOption("-psi_step"			,"Sampling rate (before oversampling) for the in-plane angle"				,"10" 	);
    option.addOption("-norm"				,"Do norm correction for image data"										,"0"  	);
    option.addOption("-ctf"					,"Perform CTF correction?"													,"0"  	);
    option.addOption("-only_flip_phases"	,"Only perform CTF phase-flipping? (default is full amplitude-correction)"	,"0"	);
    option.addOption("-ctf_phase_flipped"	,"Have the data been CTF phase-flipped?"									,"0"	);
    option.addOption("-ctf_corrected_ref",  "Have the input references been CTF-amplitude corrected?"                   ,"0"    );
    option.addOption("-intact_ctf_first_peak","Ignore CTFs until their first peak?"										,"0"	);
    option.addOption("-zero_mask",          "Mask surrounding background in particles to zero (by default the solvent area is filled with random noise)", "0" );
    option.addOption("-flatten_solvent",    "Perform masking on the references as well",                                                                  "0" );
    
    if (argc < 3) {
        option.printHelp();
        std::cerr << "example : " << std::endl;
        std::cerr << "../bin/rome_deep2d -i class19.star -o class19_deep2d -map2d_K 10 -sml_K 10 -map2d_iter 20 -sml_iter 30 -angpix 1.74 -ctf -norm > deep2d_output.txt" << std::endl;
        std::cerr << std::endl;
        EXIT_ABNORMALLY;
    }
    
    option.readCommandLine(argc, argv);
    
    //----------- ML2D main program  -------------------//
    
#ifdef USEMPI
    MPI::Init();
    if(MPI::COMM_WORLD.Get_rank() == 0)
        option.printValue();
#else
    option.printValue();
#endif

    // input file
    std::string star_fn 			      	=	pathRemoveSuffix(option.getStrOption("-i"))+".star";
    Map2dOptimizer_old::star_fn           	=   star_fn;
    Map2dOptimizer_old::mrcs_dir          	=   pathGetDir(star_fn);
    // output file
    std::string fn                  	  	=   pathRemoveSuffix(option.getStrOption("-o"));
    Map2dOptimizer_old::write_fn          	=   pathGetFilename(fn);
    Map2dOptimizer_old::write_path        	=   pathGetDir(fn);
    // basic
    Map2dOptimizer_old::nr_classes        	=   option.getIntOption("-map2d_K");
    Map2dOptimizer_old::nr_pool           	=   option.getIntOption("-pool");
    Map2dOptimizer_old::pixel_size        	=   option.getFloatOption("-angpix");
    Map2dOptimizer_old::nr_iter           	=   option.getIntOption("-map2d_iter");
    Map2dOptimizer_old::random_seed       	=   option.getIntOption("-random_seed");
    Map2dOptimizer_old::offset_step       	=   option.getIntOption("-offset_step");
    Map2dOptimizer_old::offset_range		=   option.getIntOption("-offset_range");
    Map2dOptimizer_old::psi_step          	=   option.getIntOption("-psi_step");
    //add some option
    Map2dOptimizer_old::particle_diameter   =   option.getIntOption("-particle_diameter");
    Map2dOptimizer_old::tau2_fudge_factor   =   option.getFloatOption("-tau2_fudge");
    //
    Map2dOptimizer_old::do_ctf_correction 	=   option.getBoolOption("-ctf");
    Map2dOptimizer_old::only_flip_phases	=	option.getBoolOption("-only_flip_phases");
    Map2dOptimizer_old::ctf_phase_flipped	=	option.getBoolOption("-ctf_phase_flipped");
    Map2dOptimizer_old::refs_are_ctf_corrected= option.getBoolOption("-ctf_corrected_ref");
    Map2dOptimizer_old::intact_ctf_first_peak=	option.getBoolOption("-intact_ctf_first_peak");
    Map2dOptimizer_old::do_norm_correction	=   option.getBoolOption("-norm");
    Map2dOptimizer_old::do_zero_mask	  	=	option.getBoolOption("-zero_mask");
    Map2dOptimizer_old::do_solvent			=	option.getBoolOption("-flatten_solvent");
    // others
    Map2dOptimizer_old::iter              	=   0;//option.getIntOption("-continue");
    Map2dOptimizer_old::guidance_fn       	=   "NULL";//option.getStrOption("-testguidance");
    Map2dOptimizer_old::ref_fn            	=	"NULL";
    
    double t1 = dtime();
    
    Map2dOptimizer_old::setupMap2dOptimizer();
    
    Map2dOptimizer_old::readImages();

    Map2dOptimizer_old::prepare();
    
    Map2dOptimizer_old::iterate();
    
    Map2dOptimizer_old::destroyMap2dOptimizer();
    
    double t2 = dtime();
    
    if(Map2dOptimizer_old::node == 0)
        std::cout<<"ML2D costs : "<<(t2-t1)<<std::endl;
    

#ifdef USEMPI
    MPI::COMM_WORLD.Barrier();
#endif
    
    //--------- GTM main program  ------------//
    
    
    int set_sampling_dim            = option.getIntOption("-dimension");
    if (set_sampling_dim > 3 || set_sampling_dim < 1) {std::cerr<<"Wrong sampling dimension."<<std::endl;EXIT_ABNORMALLY;}
    
    double set_precision            = 10e-12;//option.getFloatOption("-precision");
    double set_alpha                = 0.01;//option.getIntOption("-updateAlpha");
    bool set_update_beta            = 1;//option.getIntOption("-updateBeta");
    
    double probThreshold            = 1.0;//option.getFloatOption("-weightedsum");
    
    // ----- initialize Basis for multi-dimension
    auto split=[&](const std::string &s, char delim) {
        std::vector<int> elems;
        std::stringstream ss(s);
        std::string item;
        while (getline(ss, item, delim)) elems.push_back(atoi(item.c_str()));
        return elems;
    };
    
    double X_infos[9] = {0};
    double Mu_infos[9] = {0};
    int M = 0;
    std::vector<int> K_info;
    switch (set_sampling_dim) {
        case 1:
            K_info = split(option.getStrOption("-sml_K"),',');
            if (K_info.size() != 1) {std::cerr<<"Wrong for input -sml_K"<<std::endl;EXIT_ABNORMALLY;}
            X_infos[0]=1;X_infos[1]=K_info[0];X_infos[2]=K_info[0];
            M=X_infos[2]*0.8;Mu_infos[0]=1.*(double)M/(M-1);Mu_infos[1]=X_infos[2]*(double)M/(M-1);Mu_infos[2]=M;
            break;
        case 2:
            K_info = split(option.getStrOption("-sml_K"),',');
            if (K_info.size() != 2) {std::cerr<<"Wrong for input -sml_K"<<std::endl;EXIT_ABNORMALLY;}
            X_infos[0]=1;X_infos[1]=K_info[0];X_infos[2]=K_info[0];
            X_infos[3]=1;X_infos[4]=K_info[1];X_infos[5]=K_info[1];
            M=X_infos[2]*0.8;Mu_infos[0]=1.*(double)M/(M-1);Mu_infos[1]=X_infos[2]*(double)M/(M-1);Mu_infos[2]=M;
            M=X_infos[5]*0.8;Mu_infos[3]=1.*(double)M/(M-1);Mu_infos[4]=X_infos[5]*(double)M/(M-1);Mu_infos[5]=M;
            break;
        case 3:
            K_info = split(option.getStrOption("-sml_K"),',');
            if (K_info.size() != 3) {std::cerr<<"Wrong for input -sml_K"<<std::endl;EXIT_ABNORMALLY;}
            X_infos[0]=1;X_infos[1]=K_info[0];X_infos[2]=K_info[0];
            X_infos[3]=1;X_infos[4]=K_info[1];X_infos[5]=K_info[1];
            X_infos[6]=1;X_infos[7]=K_info[2];X_infos[8]=K_info[2];
            M = X_infos[2]*0.8;Mu_infos[0]=1.*(double)M/(M-1);Mu_infos[1]=X_infos[2]*(double)M/(M-1);Mu_infos[2]=M;
            M = X_infos[5]*0.8;Mu_infos[3]=1.*(double)M/(M-1);Mu_infos[4]=X_infos[5]*(double)M/(M-1);Mu_infos[5]=M;
            M = X_infos[8]*0.8;Mu_infos[6]=1.*(double)M/(M-1);Mu_infos[7]=X_infos[8]*(double)M/(M-1);Mu_infos[8]=M;
            break;
        default:
            break;
    }
    std::string set_star_fn             = 	pathRemoveSuffix(option.getStrOption("-o"))+"_iter"+num2str(option.getIntOption("-map2d_iter"))+".star";
    std::string set_result_fn       	= 	pathRemoveSuffix(option.getStrOption("-o"))+"_sml";
    GTMoptimizer::star_fn				=	set_star_fn;
    GTMoptimizer::result_fn				=	set_result_fn;
    GTMoptimizer::pixel_size			=	option.getFloatOption("-angpix");
    GTMoptimizer::nr_iter				=	option.getIntOption("-sml_iter");
    //
    GTMoptimizer::do_ctf_correction		=	option.getBoolOption("-ctf");
    GTMoptimizer::only_flip_phases		=	option.getBoolOption("-only_flip_phases");
    GTMoptimizer::ctf_phase_flipped		=	option.getBoolOption("-ctf_phase_flipped");
    GTMoptimizer::intact_ctf_first_peak	=	option.getBoolOption("-intact_ctf_first_peak");
    //
    GTMoptimizer::init_by_pca			=	option.getBoolOption("-pca");
    GTMoptimizer::numMIC				=	0;//option.getIntOption("-nummic");
    GTMoptimizer::loadMIC				= 	1;//option.getFloatOption("-loadmic");
    
    t1 = dtime();
    
    GTMoptimizer::setupGTMoptimizer(X_infos,Mu_infos,set_sampling_dim);
    
    GTMoptimizer::run(set_precision,probThreshold,set_alpha,set_update_beta);
    
    GTMoptimizer::destroyGTMoptimizer();
    
    
    t2 = dtime();
    
    if(GTMoptimizer::node == 0) std::cout<<"GTM costs : "<<(t2-t1)<<std::endl;
    
#ifdef USEMPI
    MPI::Finalize();
#endif
    
    return 0;
}
