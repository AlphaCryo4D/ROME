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
#ifndef TILECONFIGURE_H_
#define TILECONFIGURE_H_

#include "../src/resmap/resmap_time.h"

#define FIX_N_OVER_ROT_TILE

class Map3DTileCfg : public NoCopy {
private:
    int const_iimage_tile;//=2*maxthreads
    int const_ipsi_tile;//=8
    int const_iover_rot_tile;//=8
    bool tile_size_select;// = false;
    static const int global_search_number = 12;
    static const int local_search_number = 6;
    static const int total_search_number = 18;
public:
    Map3DTileCfg():const_iimage_tile(-1),const_ipsi_tile(-1),const_iover_rot_tile(-1),tile_size_select(false){/*srand(time(NULL));*/}
    ~Map3DTileCfg(){}
    void setConstTile(int _iimage_tile,int _ipsi_tile,int _iover_rot_tile);
    int get_max_iimage_tile(){return const_iimage_tile;}
    int get_max_ipsi_tile(){return const_ipsi_tile;}
    int get_max_iover_rot_tile(){return const_iover_rot_tile;}
    
    void getTileForCoarseSearch(int N, bool can_N_tile, int& _N_tile,
        int& _ipsi_tile,int& _iimage_tile,int& _iimage_sub_tile,
        int& _itrans_sub_tile,bool do_local_searching)
    {
        _N_tile				= can_N_tile ? N : N;
        _ipsi_tile			= const_ipsi_tile;
        _iimage_tile		= can_N_tile ? const_iimage_tile : const_iimage_tile;
        _iimage_sub_tile	= can_N_tile ?  4 : 4;
        _itrans_sub_tile	= can_N_tile ?  8 : 8;

        // the ispi idir iimage matrix will be sparse in local search,donot set tile
		if (do_local_searching)
        {
            _ipsi_tile = 1;
            _iimage_tile = 1;
            _iimage_sub_tile = 1;
        }
    }
    
    void getTileForFineSearch(int N,bool can_N_tile,int& _N_tile,
        int& _ipsi_tile,int& _iover_rot_tile,
        int& _iimage_tile,int& _iimage_sub_tile,
        int& _itrans_sub_tile,bool do_local_searching)
    {
        _N_tile				= can_N_tile ? N : N;
        _ipsi_tile			= const_ipsi_tile;
        _iover_rot_tile  	= const_iover_rot_tile;	// other values lead to changes in output		std::min(4,exp_nr_over_rot);
        _iimage_tile		= can_N_tile ? const_iimage_tile : const_iimage_tile;
        _iimage_sub_tile	= can_N_tile ?  4 : 4;
        _itrans_sub_tile	= can_N_tile ?  8 : 8;
        
        // the ispi idir iimage matrix will be sparse in local search,donot set tile
        if (do_local_searching) {
            _ipsi_tile = 1;
            _iimage_tile = 1;
            _iimage_sub_tile = 1;
        }
    }
    
    void getTileForUpdateModel(int N,int& _N_tile,
        int& _ipsi_tile,int& _iover_rot_tile,
        int& _iimage_tile,int& _iimage_sub_tile,
        int& _itrans_sub_tile,bool do_local_searching)
    {
        _N_tile = N;
        _ipsi_tile = const_ipsi_tile;
        _iover_rot_tile = const_iover_rot_tile;
        _iimage_tile = const_iimage_tile;
        _iimage_sub_tile = 4;
        _itrans_sub_tile = 8;
        // the ispi idir iimage matrix will be sparse in local search,donot set tile
        if (do_local_searching) {
            _ipsi_tile = 1;
            _iimage_tile = 1;
            _iimage_sub_tile = 1;
        }
    }
    
    void getTileForBackproject(int N,int& _N_tile,int& _ipsi_tile,int& _iover_rot_tile,
        int& _iimage_tile,int& _iimage_sub_tile,int& _itrans_sub_tile,
        bool do_local_searching)
    {
        _N_tile = N;
        _ipsi_tile = const_ipsi_tile;
        _iover_rot_tile = const_iover_rot_tile;
        _iimage_tile = const_iimage_tile;
        _iimage_sub_tile = 4;
        _itrans_sub_tile = 8;
        // the ispi idir iimage matrix will be sparse in local search,donot set tile
        if (do_local_searching) {
            _ipsi_tile = 1;
            _iimage_tile = 1;
            _iimage_sub_tile = 1;
        }
    }
    
    typedef std::vector<int> SingleTileConfig;
    struct Bound { int begin; int end; int interval;};
    typedef std::vector<Bound> Bounds;
    Microseconds askedTileTime;
    
    // tile configure
    typedef std::vector<std::pair<SingleTileConfig, Microseconds> > MultiTileConfigure;
    // ordered map keep best tile configure
    enum BestTileConfigureStatus {random_tile,perturbation_tile,best_tile};
    struct TileConfigRunningTime {
        SingleTileConfig tileCfg;
        std::vector<double> running_times;
    };
    class BestTileConfigure : public NoCopy {
    private:
        std::multimap<Microseconds,TileConfigRunningTime> avgTimeAndTileRunningTimes;
        BestTileConfigureStatus status = random_tile;
    public:
        BestTileConfigure():status(random_tile){}
        ~BestTileConfigure(){avgTimeAndTileRunningTimes.clear();}
        void insert(Microseconds time,const SingleTileConfig cfg)
        {
            bool theTileAlreadyHas = false;
            for (auto& v: avgTimeAndTileRunningTimes)
            {
                if(v.second.tileCfg == cfg) {
                    theTileAlreadyHas = true;
                    std::pair<Microseconds,TileConfigRunningTime> new_v;
                    new_v.first = (time + v.first*v.second.running_times.size())/(v.second.running_times.size()+1);
                    new_v.second.tileCfg = cfg;
                    new_v.second.running_times = v.second.running_times;
                    new_v.second.running_times.push_back(time);
                    avgTimeAndTileRunningTimes.erase(v.first);
                    avgTimeAndTileRunningTimes.insert(new_v);
                    break;
                }
            }
            if (!theTileAlreadyHas) {
                std::pair<Microseconds,TileConfigRunningTime> new_v;
                new_v.first = time;
                new_v.second.tileCfg = cfg;
                new_v.second.running_times.push_back(time);
                avgTimeAndTileRunningTimes.insert(new_v);
            }
//            std::cout<<"------------------------------"<<std::endl;
//            int i = 1;
//            for (auto& v: avgTimeAndTileRunningTimes)
//            {
//                std::cout<<std::setw(10)<<i++;
//                std::cout<<std::setw(10)<<v.first;
//                for (const auto& v: v.second.tileCfg)
//                    std::cout<<std::setw(10)<<v;
//                for (const auto& v: v.second.running_times)
//                    std::cout<<std::setw(10)<<v;
//                std::cout<<std::endl;
//            }
            if (avgTimeAndTileRunningTimes.size()>global_search_number+local_search_number) {
                avgTimeAndTileRunningTimes.erase(avgTimeAndTileRunningTimes.end()->first);
                return;
            }
            if (avgTimeAndTileRunningTimes.size()>=global_search_number) status = perturbation_tile;
            if (avgTimeAndTileRunningTimes.size()>=global_search_number+local_search_number) status = best_tile;
        }
        SingleTileConfig get() {
            assert(status == perturbation_tile || status == best_tile);
            auto it = avgTimeAndTileRunningTimes.begin();
            if(status == perturbation_tile) {
                int count = rand() % 3 + 1; // get 1,2,3 elememt to add perturbation
                for (int i = 0; i < count; i++) it++;
            }
            return it->second.tileCfg;
        }
        void clear() {
            avgTimeAndTileRunningTimes.clear();
            status = random_tile;
        }
        BestTileConfigureStatus getStatus(){return status;}
    };
    
    bool setBestTileLearing(bool _learn_tile_size,int images_batch_number);
    SingleTileConfig askTileConfigure(const Bounds& bounds, BestTileConfigure& bestTileConfigure);
    void recoredTileConfigure(const SingleTileConfig& cfg,MultiTileConfigure& tileConfigure,BestTileConfigure& bestTileConfigure);
    
#define ALL_MODULE	\
    ELT(CoarseSearch	, 0		)	SEP \
	ELT(FineSearch		, 1		)	SEP \
	ELT(UpdateModel		, 2		)	SEP \
	ELT(Backproject		, 3		)
//
#define ELT(N, I) \
	MultiTileConfigure tileConfigureFor##N;	\
	BestTileConfigure bestTileConfigureFor##N;	\
	SingleTileConfig askTileConfigureFor##N(Bounds bounds) {return askTileConfigure(bounds,bestTileConfigureFor##N);}	\
    void recoredTileConfigureFor##N(SingleTileConfig cfg) {recoredTileConfigure(cfg,tileConfigureFor##N,bestTileConfigureFor##N);}
#define SEP
    ALL_MODULE
#undef ELT
#undef SEP
    
    void writeOut(std::string tile_size_fn);
    void resetAll();
    
    // Mat <- R(MxN)
    void printCompareMatrix(std::string fn,std::string title,const char* Mat,int M,int N);
};

#endif

