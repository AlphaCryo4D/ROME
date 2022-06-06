/***************************************************************************
 *
 * Dana-Farber Cancer Institute, Harvard Medical School and Peking University
 * Intel® Parallel Computing Center for Structural Biology
 *
 * Authors: "Yong Bei Ma(galowma@gmail.com)"
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

#include "tileconfigure.h"


void Map3DTileCfg::setConstTile(int _iimage_tile,int _ipsi_tile,int _iover_rot_tile)
{
    const_iimage_tile = _iimage_tile/2;
    const_ipsi_tile = _ipsi_tile;
    const_iover_rot_tile = _iover_rot_tile;
}

bool Map3DTileCfg::setBestTileLearing(bool _tile_size_select,int images_batch_number)
{
    // check wheather need do tile size learning
    // global_search_number of batch of images : use random tile size
    // local_search_number of batch of images : use best tile add perturbation
    // batch_of_images-global_search_number-local_search_number : use best tile
    // assume the worst time and best time range is 1~10(gauss distribution)
    // when images_batch_number*5 > (1*local_search_number)+(5*global_search_number)
    // 								+1*(images_batch_number-global_search_number-local_search_number)+delta,it worths trying tile learning
    tile_size_select = _tile_size_select;
    if (tile_size_select) {
        double delta = 10;
        double best_time = (2*local_search_number)+(5*global_search_number)
        					+2*(images_batch_number-global_search_number-local_search_number)+delta;
        if(images_batch_number*5 > best_time){
            return true;
        }
        else {
            tile_size_select = false;
            return false;
        }
    }
    else {
        return false;
    }
}

Map3DTileCfg::SingleTileConfig Map3DTileCfg::askTileConfigure(const Bounds& bounds, BestTileConfigure& bestTileConfigure)
{
    auto setRandomTileSize = [&](SingleTileConfig& tilecfg) {
        for (auto& bound : bounds) {
            //int choice = (bound.end-bound.begin)/bound.interval;
            int choice = (bound.end/2-bound.begin)/bound.interval;
            int tile_size = choice? ( (rand() % choice+1) * bound.interval + bound.begin ) : bound.begin;
            tile_size = std::min(tile_size, bound.end);
            tilecfg.push_back(tile_size);
        }
    };
    auto addPerturbationToTileSize = [&](SingleTileConfig& tilecfg) {
        assert(tilecfg.size()==bounds.size());
        int tileNumber = bounds.size();
        for (int tile_i = 0; tile_i < tileNumber; tile_i++) {
            int _min = -(bounds[tile_i].interval-1);
            int _max =  (bounds[tile_i].interval-1);
            int perturbation = (rand()%(_max-_min+1))+_min;
            tilecfg[tile_i] = (tilecfg[tile_i]+perturbation<=0)?tilecfg[tile_i]:tilecfg[tile_i]+perturbation;
        }
    };
    
    SingleTileConfig tilecfg;
    if (tile_size_select)
    {
        switch (bestTileConfigure.getStatus())
        {
            case random_tile:
                setRandomTileSize(tilecfg);
                break;
            case perturbation_tile:
                tilecfg = bestTileConfigure.get();
                addPerturbationToTileSize(tilecfg);
                break;
            case best_tile:
                tilecfg = bestTileConfigure.get();
                break;
            default:
                break;
        }
    }
    else
    {
        setRandomTileSize(tilecfg);
    }

    askedTileTime = timeInMicroseconds();
    return tilecfg;
}

void Map3DTileCfg::recoredTileConfigure(const SingleTileConfig& cfg,MultiTileConfigure& tileConfigure,BestTileConfigure& bestTileConfigure)
{
    double costTime = timeInMicroseconds() - askedTileTime;
    tileConfigure.push_back(std::pair<SingleTileConfig,Microseconds>(cfg,costTime));
    if(tile_size_select) bestTileConfigure.insert(costTime, cfg);
}

void Map3DTileCfg::writeOut(std::string tile_size_fn)
{
    
    std::ofstream out;
    
    std::multimap<Microseconds,SingleTileConfig> tileConfigureByTime;
    
    auto writeAll = [&](MultiTileConfigure& tileConfigure,int tile_num,std::vector<std::string> tile_names,std::string file_name)
    {
        // open file and write out header
        //for (const auto& tile_name : tile_names) file_name += "_"+tile_name;
        file_name += "_all.txt";
        out.open(file_name.c_str());
        out<<std::setw(15)<<"I"<<'\t';
        for (const auto& tile_name : tile_names) out<<std::setw(15)<<tile_name<<'\t';
        out<<std::setw(15)<<"time"<<std::endl;
        //
        int i = 1;
        for (const auto& tile_time: tileConfigure) {
            assert(tile_time.first.size()==tile_num);
            out<<std::setw(15)<<(i++)<<'\t';
            for (int tile_i = 0; tile_i < tile_num; tile_i++) {
                out<<std::setw(15)<<tile_time.first[tile_i]<<'\t';
            }
            out<<std::setw(15)<<std::setprecision(10)<<tile_time.second<<std::endl;
        }
        out.close();
    };

    auto writeByAverage = [&](MultiTileConfigure& tileConfigure,int tile_num,std::vector<std::string> tile_names,std::string file_name)
    {
        //
    };
    
    //
    auto writeByTime = [&](MultiTileConfigure& tileConfigure,int tile_num,std::vector<std::string> tile_names,std::string file_name)
    {
        //sorted by time
        for (const auto& tile_time: tileConfigure)
            tileConfigureByTime.insert(std::pair<Microseconds,SingleTileConfig>(tile_time.second,tile_time.first));
        // open file and write out header
        file_name += "_sorted.txt";
        out.open(file_name.c_str());
        out<<std::setw(15)<<"I"<<'\t';
        for (const auto& tile_name : tile_names) out<<std::setw(15)<<tile_name<<'\t';
        out<<std::setw(15)<<"time"<<std::endl;
        //
        int i = 1;
        for (const auto& tile_time: tileConfigureByTime) {
            assert(tile_time.second.size()==tile_num);
            out<<std::setw(15)<<(i++)<<'\t';
            for (int tile_i = 0; tile_i < tile_num; tile_i++) {
                out<<std::setw(15)<<tile_time.second[tile_i]<<'\t';
            }
            out<<std::setw(15)<<std::setprecision(10)<<tile_time.first<<std::endl;
        }
        out.close();
        tileConfigureByTime.clear();
    };
    
#ifdef FIX_N_OVER_ROT_TILE
    // write out all tile configure and its time
    writeAll(tileConfigureForCoarseSearch,3,{"image_tile","imagesub_tile","itranssub_tile"},tile_size_fn+"_coarse");
    writeAll(tileConfigureForFineSearch,3,{"image_tile","iimagesub_tile","itranssub_tile"},tile_size_fn+"_fine");
    writeAll(tileConfigureForUpdateModel,3,{"image_tile","iimagesub_tile","itranssub_tile"},tile_size_fn+"_updatemodel");
    writeAll(tileConfigureForBackproject,3,{"image_tile","iimagesub_tile","itranssub_tile"},tile_size_fn+"_backproject");
    
    // write out sorted time and its configure
    writeByTime(tileConfigureForCoarseSearch,3,{"image_tile","imagesub_tile","itranssub_tile"},tile_size_fn+"_coarse");
    writeByTime(tileConfigureForFineSearch,3,{"image_tile","iimagesub_tile","itranssub_tile"},tile_size_fn+"_fine");
    writeByTime(tileConfigureForUpdateModel,3,{"image_tile","iimagesub_tile","itranssub_tile"},tile_size_fn+"_updatemodel");
    writeByTime(tileConfigureForBackproject,3,{"image_tile","iimagesub_tile","itranssub_tile"},tile_size_fn+"_backproject");
#else
    // write out all tile configure and its time
    writeAll(tileConfigureForCoarseSearch,4,{"N_tile","image_tile","imagesub_tile","itranssub_tile"},tile_size_fn+"_coarse");
    writeAll(tileConfigureForFineSearch,5,{"N_tile","image_tile","ioverrot_tile","iimagesub_tile","itranssub_tile"},tile_size_fn+"_fine");
    writeAll(tileConfigureForUpdateModel,5,{"N_tile","image_tile","ioverrot_tile","iimagesub_tile","itranssub_tile"},tile_size_fn+"_updatemodel");
    writeAll(tileConfigureForBackproject,5,{"N_tile","image_tile","ioverrot_tile","iimagesub_tile","itranssub_tile"},tile_size_fn+"_backproject");
    
    // write out sorted time and its configure
    writeByTime(tileConfigureForCoarseSearch,4,{"N_tile","image_tile","imagesub_tile","itranssub_tile"},tile_size_fn+"_coarse");
    writeByTime(tileConfigureForFineSearch,5,{"N_tile","image_tile","ioverrot_tile","iimagesub_tile","itranssub_tile"},tile_size_fn+"_fine");
    writeByTime(tileConfigureForUpdateModel,5,{"N_tile","image_tile","ioverrot_tile","iimagesub_tile","itranssub_tile"},tile_size_fn+"_updatemodel");
    writeByTime(tileConfigureForBackproject,5,{"N_tile","image_tile","ioverrot_tile","iimagesub_tile","itranssub_tile"},tile_size_fn+"_backproject");
#endif
    //
    
    tileConfigureByTime.clear();
    
    tileConfigureForCoarseSearch.clear();
    tileConfigureForFineSearch.clear();
    tileConfigureForUpdateModel.clear();
    tileConfigureForBackproject.clear();
}

void Map3DTileCfg::resetAll()
{
#define ELT(N, I) \
    tileConfigureFor##N.clear();	\
    bestTileConfigureFor##N.clear();
#define SEP
    ALL_MODULE
#undef ELT
#undef SEP
}

void Map3DTileCfg::printCompareMatrix(std::string fn,std::string title,const char* Mat,int M,int N)
{
    static std::ofstream out(fn);
    out<<title<<std::endl;
    for (int m = 0; m < M; m++) {
        for (int n = 0; n < N; n++){
            out<<bool(Mat[m*N+n])<<" ";
        }
        out<<std::endl;
    }
}