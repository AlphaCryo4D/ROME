/***************************************************************************
 *
 * Authors: "Jiayi (Timmy) Wu, Yongbei(Glow) Ma, Youdong (Jack) Mao"
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

#include <ctime>
#include <cstdlib>

#include "../src/option.h"
#include "../src/gtm_optimizer.h"
#include "../src/string.h"

void fineClassify(Option& option)
{
    // do fine search
    double X_infos[3];
    double Mu_infos[3];
    
    std::string set_star_fn         = pathRemoveSuffix(option.getStrOption("-i"))+".star";
    std::string set_result_fn       = pathRemoveSuffix(option.getStrOption("-o"));
    
    double set_precision            = option.getFloatOption("-precision");
    double set_alpha                = option.getFloatOption("-alpha");
    bool set_update_beta            = option.getIntOption("-updatebeta");
    int set_sampling_dim            = option.getIntOption("-dimension");
    if (set_sampling_dim != 1) {std::cerr<<"Wrong sampling dimension for GTM fine searching."<<std::endl;EXIT_ABNORMALLY;}
    // turn off the mic offload mode for fine searching
    double probThreshold            = 1.0;//option.getFloatOption("-weightedsum");
    
    std::string config_fn           = option.getStrOption("-search");
    
    double t1 = dtime();

    GTMoptimizer::fineSearch.setup(set_star_fn,config_fn);
    int nr_selectedClass = GTMoptimizer::fineSearch.getSelectedClassNumber();
    
    for (int iclass = 0; iclass < nr_selectedClass; iclass++)
    {
        if(GTMoptimizer::fineSearch.searchNext(set_star_fn))
        {
        
            X_infos[0]                      	= 	1;
            X_infos[1]                      	= 	GTMoptimizer::fineSearch.getCurrentClassFineClassNumber();
            X_infos[2]                      	= 	GTMoptimizer::fineSearch.getCurrentClassFineClassNumber();
            
            int M                           	= 	X_infos[2]*0.8;
            Mu_infos[0]                     	= 	1.*(double)M/(M-1);
            Mu_infos[1]                    	 	= 	X_infos[2]*(double)M/(M-1);
            Mu_infos[2]                     	= 	M;
            
            GTMoptimizer::star_fn				=	set_star_fn;
            GTMoptimizer::result_fn				=	set_result_fn;
            GTMoptimizer::nr_iter				=	option.getIntOption("-iter");
            GTMoptimizer::pixel_size			=	option.getFloatOption("-angpix");
            //
            GTMoptimizer::do_ctf_correction		=	option.getBoolOption("-ctf");
            GTMoptimizer::only_flip_phases		=	option.getBoolOption("-only_flip_phases");
            GTMoptimizer::ctf_phase_flipped		=	option.getBoolOption("-ctf_phase_flipped");
            GTMoptimizer::intact_ctf_first_peak	=	option.getBoolOption("-intact_ctf_first_peak");
            //
            GTMoptimizer::numMIC				=	0;
            GTMoptimizer::loadMIC				=	0;
            GTMoptimizer::init_by_pca			=	false;
            
            GTMoptimizer::setupGTMoptimizer(X_infos,Mu_infos,set_sampling_dim);
            
            GTMoptimizer::run(set_precision,probThreshold,set_alpha,set_update_beta);
            
            GTMoptimizer::destroyGTMoptimizer();
            
        }
        if( iclass%10==0 || iclass==(nr_selectedClass-1) )
            GTMoptimizer::fineSearch.writeResult(set_result_fn+"_fine");
#ifdef USEMPI
        MPI::COMM_WORLD.Barrier();
#endif
    }
    
    double t2 = dtime();
    
    if(GTMoptimizer::node == 0) {
        std::cout<<"GTM costs : "<<(t2-t1)<<std::endl;
    }

}

void ordinaryClassify(Option& option)
{
    if (option.getBoolOption("-debug"))
    {
        /*** my baseline dataset to check the result ***/
        
        // Option option;
        //dim = 100,read data time : 7.81634,bcast data time : 2.22168
        //debug_refence_input = "../dataset/flower/flower_template_25";
        
        double X_infos[3];
        double Mu_infos[3];
        X_infos[0]  = 0;
        X_infos[1]  = 2*3.1415926535897;
        X_infos[2]  = 25;
        Mu_infos[0] = 0.*20./19.;
        Mu_infos[1] = 2*3.1415926535897*20./19.;
        Mu_infos[2] = 20;
        int set_sampling_dim = 1;
        
        double set_precision = 10e-12;
        double set_update_beta = true;
        double set_alpha = 0.01;
        double probThreshold = 1.0;
        
        double t1 = dtime();
        
        GTMoptimizer::star_fn				=	option.getStrOption("-i");
        GTMoptimizer::result_fn				=	pathRemoveSuffix(option.getStrOption("-o"));
        GTMoptimizer::pixel_size			=	option.getFloatOption("-angpix");
        GTMoptimizer::nr_iter				=	option.getIntOption("-iter");
        GTMoptimizer::ref_fn				=	option.getStrOption("-ref");
        //
        GTMoptimizer::do_ctf_correction		=	option.getBoolOption("-ctf");
        GTMoptimizer::only_flip_phases		=	option.getBoolOption("-only_flip_phases");
        GTMoptimizer::ctf_phase_flipped		=	option.getBoolOption("-ctf_phase_flipped");
        GTMoptimizer::intact_ctf_first_peak	=	option.getBoolOption("-intact_ctf_first_peak");
        //
        GTMoptimizer::init_by_pca			=	false;
        GTMoptimizer::numMIC				=	0;
        GTMoptimizer::loadMIC				=	0;
        
        GTMoptimizer::setupGTMoptimizer(X_infos,Mu_infos,set_sampling_dim);
        
        GTMoptimizer::run(set_precision,probThreshold,set_alpha,set_update_beta);
        
        GTMoptimizer::destroyGTMoptimizer();
        
        double t2 = dtime();
        
        if(GTMoptimizer::node == 0) std::cout<<"GTM costs : "<<(t2-t1)<<std::endl;
        
    }
    else
    {
        double set_precision            = option.getFloatOption("-precision");
        double set_alpha                = option.getFloatOption("-alpha");
        bool set_update_beta            = option.getIntOption("-updatebeta");
        int set_sampling_dim            = option.getIntOption("-dimension");
        if (set_sampling_dim > 3 || set_sampling_dim < 1) {std::cerr<<"Wrong sampling dimension."<<std::endl;EXIT_ABNORMALLY;}
        
        double probThreshold            = option.getFloatOption("-weightedsum");
        
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
                K_info = split(option.getStrOption("-K"),',');
                if (K_info.size() != 1) {std::cerr<<"Wrong for input -K"<<std::endl;EXIT_ABNORMALLY;}
                X_infos[0]=1;X_infos[1]=K_info[0];X_infos[2]=K_info[0];
                M=X_infos[2]*0.8;Mu_infos[0]=1.*(double)M/(M-1);Mu_infos[1]=X_infos[2]*(double)M/(M-1);Mu_infos[2]=M;
                break;
            case 2:
                K_info = split(option.getStrOption("-K"),',');
                if (K_info.size() != 2) {std::cerr<<"Wrong for input -K"<<std::endl;EXIT_ABNORMALLY;}
                X_infos[0]=1;X_infos[1]=K_info[0];X_infos[2]=K_info[0];
                X_infos[3]=1;X_infos[4]=K_info[1];X_infos[5]=K_info[1];
                M=X_infos[2]*0.8;Mu_infos[0]=1.*(double)M/(M-1);Mu_infos[1]=X_infos[2]*(double)M/(M-1);Mu_infos[2]=M;
                M=X_infos[5]*0.8;Mu_infos[3]=1.*(double)M/(M-1);Mu_infos[4]=X_infos[5]*(double)M/(M-1);Mu_infos[5]=M;
                break;
            case 3:
                K_info = split(option.getStrOption("-K"),',');
                if (K_info.size() != 3) {std::cerr<<"Wrong for input -K"<<std::endl;EXIT_ABNORMALLY;}
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
        
        double t1 = dtime();
        
        GTMoptimizer::star_fn				=	pathRemoveSuffix(option.getStrOption("-i"))+".star";
        GTMoptimizer::result_fn				=	pathRemoveSuffix(option.getStrOption("-o"));
        GTMoptimizer::pixel_size			=	option.getFloatOption("-angpix");
        GTMoptimizer::nr_iter				=	option.getIntOption("-iter");
        GTMoptimizer::ref_fn				=	option.getStrOption("-ref");
        //
        GTMoptimizer::do_ctf_correction		=	option.getBoolOption("-ctf");
        GTMoptimizer::only_flip_phases		=	option.getBoolOption("-only_flip_phases");
        GTMoptimizer::ctf_phase_flipped		=	option.getBoolOption("-ctf_phase_flipped");
        GTMoptimizer::intact_ctf_first_peak	=	option.getBoolOption("-intact_ctf_first_peak");
        //
        GTMoptimizer::init_by_pca			=	option.getBoolOption("-pca");
        GTMoptimizer::numMIC				=	option.getIntOption("-nummic");
        GTMoptimizer::loadMIC				= 	option.getFloatOption("-loadmic");
        
        GTMoptimizer::setupGTMoptimizer(X_infos,Mu_infos,set_sampling_dim);
        
        GTMoptimizer::run(set_precision,probThreshold,set_alpha,set_update_beta);
        
        GTMoptimizer::destroyGTMoptimizer();
        
        double t2 = dtime();
        
        if(GTMoptimizer::node == 0) std::cout<<"GTM costs : "<<(t2-t1)<<std::endl;
        
    }

}

int main(int argc,char *argv[])
{
    Option option;
    // general option
    option.addOption("-i"                   ,"Input metadata file with images needed to align"                                              );
    option.addOption("-o"                   ,"Output metadata"                                                                              );
    option.addOption("-K"                   ,"Number of classes needed to classify!!!!! please use Uppercase —K instead of lowercase -k !!!!!");
    option.addOption("-iter"                ,"Maximum number of iterations to perform"                                                      );
    option.addOption("-angpix"              ,"Pixel size (in Angstroms)!!!!! please use -angpix instead of -pixel !!!!!"		            );
    // advanced option
    option.addOption("-dimension"           ,"Sampling grid dimension(Default 1)"                                               ,"1"        );
    option.addOption("-pca"                 ,"Using PCA to initialize,(Defalut is using Gauss-Circle template)"                 ,"0"        );
    option.addOption("-ref"					,"2d reference file name(*.mrcs)"                                                   ,"NULL"		);
    option.addOption("-alpha"               ,"Value of the variance of prior model Gaussian distribution"                       ,"0.01"     );
    option.addOption("-updatebeta"          ,"Update the variance of noise whether or not, 0 or 1"                              ,"1"        );
    option.addOption("-precision"           ,"The condition of GTM stop(less value faster stop)"                                ,"10e-12"   );
    option.addOption("-nummic"              ,"The number of mic cards used to compute"                                          ,"0"       	);
    option.addOption("-loadmic"             ,"The percentage of job put to compute in mic card"                                 ,"0"     	);
    option.addOption("-weightedsum"         ,"Probability threshold (0~1) for weighted class averaging(Default maximum prob)"  	,"1"     	);
    option.addOption("-search"              ,"Follow the configure file to classify each class(-search default or -search class.config)","NULL");
    //
    option.addOption("-ctf"                 ,"Perform CTF correction?"                                                          ,"0"        );
    option.addOption("-only_flip_phases"	,"Only perform CTF phase-flipping? (default is full amplitude-correction)"			,"0"		);
    option.addOption("-ctf_phase_flipped"	,"Have the data been CTF phase-flipped?"											,"0"		);
    option.addOption("-intact_ctf_first_peak","Ignore CTFs until their first peak?"												,"0"		);
    //
    option.addOption("-debug"				,"debug"																			,"0"		);
    
    if (argc < 3) {
        option.printHelp();
        std::cout<<"example : ./bin/rome_sml -i dataset/class19.star -o dataset/class19_gtm -K 10 -iter 30 -angpix 1.74 > dataset/gtm_output.txt"<<std::endl;
        EXIT_ABNORMALLY;
    }

    option.readCommandLine(argc, argv);
    
#ifdef USEMPI
    MPI::Init();
    if(MPI::COMM_WORLD.Get_rank() == 0)
        option.printValue();
#else
    option.printValue();
#endif
    
    if (option.getStrOption("-search") == "NULL")
        ordinaryClassify(option);
    else
        fineClassify(option);
    
#ifdef USEMPI
    MPI::Finalize();
#endif
}
