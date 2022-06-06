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

#include "resmap_util.h"		// used for building precompiled headers on Windows

#include "resmap_map_optimizer_base.h"

namespace MapOptimizerBase
{
#ifdef USEMPI
    int nodeNameLen;
    char nodeName[MPI_MAX_PROCESSOR_NAME];
#endif
    //
#define SEP
#define ELT(T,N,V) T N = V;
    MAPOPTIMIZER_BASE_VARS
#undef SEP
#undef ELT
    //
    //   ------------ track the change of hidden variable ---------------   //
    void HiddenVariableMonitor::monitorHiddenVariableChanges(HealpixSampler& sampling,const MetaDataTable& metadata_old,int first_image,
                                                             const MetaDataElem* metadata_new,int nr_images)
    {
        for (int iimage = 0; iimage < nr_images; iimage++)
        {
            // Old optimal parameters
            double old_rot 	= metadata_old[first_image+iimage].ROT;
            double old_tilt = metadata_old[first_image+iimage].TILT;
            double old_psi 	= metadata_old[first_image+iimage].PSI;
            double old_xoff = metadata_old[first_image+iimage].XOFF;
            double old_yoff = metadata_old[first_image+iimage].YOFF;
            int old_iclass 	= metadata_old[first_image+iimage].CLASS;
            // New optimal parameters
            double new_rot 	= metadata_new[iimage].ROT;
            double new_tilt = metadata_new[iimage].TILT;
            double new_psi 	= metadata_new[iimage].PSI;
            double new_xoff = metadata_new[iimage].XOFF;
            double new_yoff = metadata_new[iimage].YOFF;
            int new_iclass 	= metadata_new[iimage].CLASS;
            // Some orientational distance....
            sum_changes_optimal_orientations += sampling.calculateAngularDistance(new_rot, new_tilt, new_psi, old_rot, old_tilt, old_psi);
            sum_changes_optimal_offsets += (new_xoff-old_xoff)*(new_xoff-old_xoff) + (new_yoff-old_yoff)*(new_yoff-old_yoff);
            if (new_iclass != old_iclass) sum_changes_optimal_classes += 1.;
            sum_changes_count += 1.;
        }
    }
    //
    void HiddenVariableMonitor::reduceData()
    {
#ifdef USEMPI
        double temp = sum_changes_optimal_orientations;
        MPI::COMM_WORLD.Allreduce(&temp, &sum_changes_optimal_orientations, 1, MPI::DOUBLE, MPI::SUM);
        temp = sum_changes_optimal_offsets;
        MPI::COMM_WORLD.Allreduce(&temp, &sum_changes_optimal_offsets, 1, MPI::DOUBLE, MPI::SUM);
        temp = sum_changes_optimal_classes;
        MPI::COMM_WORLD.Allreduce(&temp, &sum_changes_optimal_classes, 1, MPI::DOUBLE, MPI::SUM);
        temp = sum_changes_count;
        MPI::COMM_WORLD.Allreduce(&temp, &sum_changes_count, 1, MPI::DOUBLE, MPI::SUM);
#endif
    }
    //
    void HiddenVariableMonitor::bcastData(int bcastNode)
    {
#ifdef USEMPI
        MPI::COMM_WORLD.Bcast(&current_changes_optimal_classes, 1, MPI::DOUBLE, bcastNode);
        MPI::COMM_WORLD.Bcast(&current_changes_optimal_orientations, 1, MPI::DOUBLE, bcastNode);
        MPI::COMM_WORLD.Bcast(&current_changes_optimal_offsets, 1, MPI::DOUBLE, bcastNode);
        MPI::COMM_WORLD.Bcast(&nr_iter_wo_large_hidden_variable_changes, 1, MPI::INT, bcastNode);
        MPI::COMM_WORLD.Bcast(&smallest_changes_optimal_classes, 1, MPI::INT, bcastNode);
        MPI::COMM_WORLD.Bcast(&smallest_changes_optimal_offsets, 1, MPI::DOUBLE, bcastNode);
        MPI::COMM_WORLD.Bcast(&smallest_changes_optimal_orientations, 1, MPI::DOUBLE, bcastNode);
#endif
    }
    void HiddenVariableMonitor::updateOverallChangesInHiddenVariables()
    {
        // Calculate hidden variable changes
        current_changes_optimal_classes = sum_changes_optimal_classes / sum_changes_count;
        current_changes_optimal_orientations = sum_changes_optimal_orientations / sum_changes_count;
        current_changes_optimal_offsets = sqrt(sum_changes_optimal_offsets / (2. * sum_changes_count));
        
        // Reset the sums
        sum_changes_optimal_classes = 0.;
        sum_changes_optimal_orientations = 0.;
        sum_changes_optimal_offsets = 0.;
        sum_changes_count = 0.;
        
        // Update nr_iter_wo_large_hidden_variable_changes if all three assignment types are within 3% of the smallest thus far
        if (1.03 * current_changes_optimal_classes >= smallest_changes_optimal_classes &&
            1.03 * current_changes_optimal_offsets >= smallest_changes_optimal_offsets &&
            1.03 * current_changes_optimal_orientations >= smallest_changes_optimal_orientations)
            nr_iter_wo_large_hidden_variable_changes++;
        else
            nr_iter_wo_large_hidden_variable_changes = 0;
        
        // Update smallest changes in hidden variables thus far
        if (current_changes_optimal_classes < smallest_changes_optimal_classes)
            smallest_changes_optimal_classes = round(current_changes_optimal_classes);
        if (current_changes_optimal_offsets < smallest_changes_optimal_offsets)
            smallest_changes_optimal_offsets = current_changes_optimal_offsets;
        if (current_changes_optimal_orientations < smallest_changes_optimal_orientations)
            smallest_changes_optimal_orientations = current_changes_optimal_orientations;
        
    }
    //  ---------------   setup   TODO   ---------------------- //
    //
    // void setupMLoptimizer();
    
    //
    // void destroyMLoptimizer();
    
    // Interpret command line for the initial start of a run
    // void prepare();
    
    // Perform expectation-maximization iterations
    // void iterate();
    //
    //
    void writeOutOptimizer(std::string fn_optimiser)
    {
        std::ofstream optimiserFile;
        optimiserFile.open((fn_optimiser+".star").c_str(), std::ios::out);
        optimiserFile << "# ROME optimiser" << std::endl << "#" <<std::endl;
        //
        {
            optimiserFile << std::endl;
            optimiserFile << "data_optimiser_general" <<std::endl;
            optimiserFile << std::endl;
#define COUTMETADATA(v1,v2) optimiserFile << "_rln" << std::left<<std::setw(30) << v1 << std::right<<std::setw(18) << v2 <<std::endl;
            //
            COUTMETADATA("OutputRootName"				, fn_optimiser.substr(0,fn_optimiser.find("_it")) 							)
            if (do_split_random_halves) {
                COUTMETADATA("ModelStarFile"			, fn_optimiser.substr(0,fn_optimiser.find("_optimiser"))+"_half1_model.star")
                COUTMETADATA("ModelStarFile2"			, fn_optimiser.substr(0,fn_optimiser.find("_optimiser"))+"_half2_model.star")
            }
            else {
                COUTMETADATA("ModelStarFile"			, fn_optimiser.substr(0,fn_optimiser.find("_optimiser"))+"_model.star" 		)
            }
            COUTMETADATA("ExperimentalDataStarFile"		, fn_optimiser.substr(0,fn_optimiser.find("_optimiser"))+"_data.star" 		)
            COUTMETADATA("OrientSamplingStarFile"		, fn_optimiser.substr(0,fn_optimiser.find("_optimiser"))+"_sampling.star"	)
            COUTMETADATA("CurrentIteration"				, iter																		)
            COUTMETADATA("NumberOfIterations"			, nr_iter																	)
            COUTMETADATA("DoSplitRandomHalves"			, 0																			)
            COUTMETADATA("JoinHalvesUntilThisResolution", -1.0																		)
            COUTMETADATA("AdaptiveOversampleOrder"		, adaptive_oversampling														)
            COUTMETADATA("AdaptiveOversampleFraction"	, adaptive_fraction															)
            COUTMETADATA("RandomSeed"					, random_seed																)
            COUTMETADATA("ParticleDiameter"				, 999																		) // TODO
            COUTMETADATA("WidthMaskEdge"				, width_mask_edge															)
            COUTMETADATA("DoZeroMask"					, do_zero_mask																)
            COUTMETADATA("DoSolventFlattening"			, do_solvent																)
            // TODO.......
            //
#undef COUTMETADATA
            optimiserFile << std::endl << std::endl;
        }
        optimiserFile.close();
    }
    //
    std::pair<std::string, std::string> readFromOptimizer(std::string fn_optimiser,bool debug_flag /*= false*/)
    {
        std::string model_fn,sampling_fn;
        ifstreamCheckingExistence optimiserFile(fn_optimiser.c_str());
        {
            bool startingRead = false;
            std::string line;
            while (true) {
                if (startingRead) {
                    double doubleTmp;std::string stringTmp;
#define CINMETADATADOUBLE(V) optimiserFile >> line;optimiserFile >> doubleTmp;NODE0ONLY if(debug_flag) {NODE0ONLY std::cout<<std::setw(30)<<line<<" "<<doubleTmp<<std::endl;}
#define CINMETADATASTR(V) optimiserFile >> line;optimiserFile >> stringTmp;NODE0ONLY if(debug_flag) {NODE0ONLY std::cout<<std::setw(30)<<line<<" "<<stringTmp<<std::endl;}
                    //
                    CINMETADATASTR(		"OutputRootName"				)
                    if(do_split_random_halves) {
                        CINMETADATASTR(		"ModelStarFile"					);model_fn = stringTmp;
                        CINMETADATASTR(		"ModelStarFile2"				);
                    }
                    else {
                        CINMETADATASTR(		"ModelStarFile"					);model_fn = stringTmp;
                    }
                    CINMETADATASTR(		"ExperimentalDataStarFile"		);//assert(star_fn==stringTmp);
                    CINMETADATASTR(		"OrientSamplingStarFile"		);sampling_fn = stringTmp;
                    CINMETADATADOUBLE(	"CurrentIteration"				);iter = doubleTmp;
                    CINMETADATADOUBLE(	"NumberOfIterations"			)
                    CINMETADATADOUBLE(	"DoSplitRandomHalves"			)
                    CINMETADATADOUBLE(	"JoinHalvesUntilThisResolution"	)
                    CINMETADATADOUBLE(	"AdaptiveOversampleOrder"		)
                    CINMETADATADOUBLE(	"AdaptiveOversampleFraction"	)
                    CINMETADATADOUBLE(	"RandomSeed"					)
                    CINMETADATADOUBLE(	"ParticleDiameter"				)
                    CINMETADATADOUBLE(	"WidthMaskEdge"					)
                    CINMETADATADOUBLE(	"DoZeroMask"					)
                    CINMETADATADOUBLE(	"DoSolventFlattening"			)
                    //
#undef CINMETADATADOUBLE
#undef CINMETADATASTR
                    break;
                }
                else{
                    getline(optimiserFile,line);
                    if(debug_flag) NODE0ONLY std::cout<<line<<std::endl;
                    if ((line.find("data_optimiser_general") !=std::string::npos) ){
                        startingRead = true;
                        getline(optimiserFile,line);assert(line=="");// escape a empty line
                    }
                    ERROR_CHECK(optimiserFile.eof(), "end of optimiser file,can not find data_optimiser_general.");
                }
            }
        }
        optimiserFile.close();
        return std::make_pair(model_fn, sampling_fn);
    }
    //
    // ---------------------------- below is all useful data structure    ----------------------------- //
    //
    // ----------------    IndexCoder   ----------------------- //
    void IndexCoder::init(std::initializer_list<int> _index_size)
    {
        fini();
        for (auto const & __index_size : _index_size)
            index_size[index_num++] = __index_size;
    }
    void IndexCoder::fini(){
        index_num=0;
        for (int i = 0; i < max_index_num; i++) index_size[i] = 0;
    }
    int IndexCoder::code(std::initializer_list<int> _index)
    {
        int index = 0;int i = 0;
        for (auto const & __index : _index){
            assert(__index < index_size[i]);
            index = index*index_size[i++]+__index;
        }
        return index;
    }
    void IndexCoder::unitTest()
    {
        IndexCoder indexCoder;
        int dimx = 6,dimy = 5,dimz = 4,dimk = 7;
        indexCoder.init({dimx,dimy,dimz,dimk});
        for (int x = 0; x < dimx; x++)
            for (int y = 0; y < dimy; y++)
                for (int z = 0; z < dimz; z++)
                    for (int k = 0; k < dimk; k++)
                    {
                        int index = indexCoder.code({x,y,z,k});
                        std::cout<<index<<" ";
                        int x2,y2,z2,k2;
                        indexCoder.decode(index,k2, z2, y2, x2);
                        std::cout<<x2<<" "<<y2<<" "<<z2<<" "<<k2<<std::endl;
                        assert(x2==x);assert(y2==y);assert(z2==z);assert(k2==k);
                    }
    }
    
    
    
    
#ifdef USEMPI
    /*******  MPI function  *********/
    
    // reduce the metadata (ROT,TILT,PSI,XOFF,YOFF,ZOFF,CLASS,DLL,PMAX,NR_SIGN,NORM) to node0
    // or maybe donot reduce this,directly write metadata to disk
    void gatherMetaDataToMaster(MetaDataTable& metadata,std::vector<int> masterNodes,bool do_split_random_halves)
    {
        int nr_global_images = metadata.numberOfParticles();
        //
#define METADATA_ELEM \
        ELT(0,ROT) ELT(1,TILT) ELT(2,PSI) ELT(3,XOFF) ELT(4,YOFF) ELT(5,ZOFF) \
        ELT(6,CLASS) ELT(7,DLL) ELT(8,PMAX) ELT(9,NR_SIGN) ELT(10,NORM)
#define METADATA_LEN 11
        FDOUBLE* metadataSendBuf = (FDOUBLE*)aMalloc(sizeof(FDOUBLE)*nr_local_images*METADATA_LEN,64);
        FDOUBLE* metadataRecvBuf;
        int displs[nodes],rcounts[nodes];
        MPI::Datatype oneMetadataElem = MPI_FDOUBLE.Create_contiguous( METADATA_LEN );
        oneMetadataElem.Commit();
        metadataRecvBuf = (FDOUBLE*)aMalloc(sizeof(FDOUBLE)*nr_global_images*METADATA_LEN,64);
        for	(int iimage = 0; iimage < nr_local_images; iimage++)
        {
#define ELT(I,T) metadataSendBuf[iimage*METADATA_LEN+I] = metadataRecvBuf[iimage*METADATA_LEN+I] = metadata[iimage+first_local_image].T;
            METADATA_ELEM
#undef ELT
        }
        int dummy;
        for	(int rank = 0; rank < nodes; rank++){
            int nr_local_images,first_local_image,last_local_image;
            setImagesDistribution(dummy, nr_local_images, nodes, rank, first_local_image, last_local_image, metadata, do_split_random_halves);
            if (rank == 0){
                rcounts[rank] = nr_local_images;
                displs[rank] = 0;
            }else{
                rcounts[rank] = nr_local_images;
                displs[rank] = displs[rank-1]+rcounts[rank-1];
            }
        }
        
        // if(masterNode==node){
        //     std::cout<<"rcounts : "<<std::endl;
        //     for	(int i = 0; i < nodes; i++)
        //         std::cout<<std::setw(10)<<rcounts[i]<<" ";
        //     std::cout<<std::endl;
        //     std::cout<<"displs : "<<std::endl;
        //     for	(int i = 0; i < nodes; i++)
        //         std::cout<<std::setw(10)<<displs[i]<<" ";
        //     std::cout<<std::endl;
        // }
        
        for (auto masterNode : masterNodes) {
            MPI::COMM_WORLD.Gatherv(metadataSendBuf, nr_local_images, oneMetadataElem, metadataRecvBuf, rcounts, displs, oneMetadataElem, masterNode);
        }
        
        if(std::find(masterNodes.begin(), masterNodes.end(), node)!=masterNodes.end()){
            for (int iimage = 0;iimage < nr_global_images;iimage++){
                // set the metadata
#define ELT(I,T) metadata.accessAll(iimage).T = metadataRecvBuf[iimage*METADATA_LEN+I];
                METADATA_ELEM
#undef ELT
            }
        }
        oneMetadataElem.Free();
        aFree(metadataSendBuf);
        aFree(metadataRecvBuf);
#undef METADATA_ELEM
#undef METADATA_LEN
    }
#endif
    
}