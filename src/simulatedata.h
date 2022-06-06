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

#ifndef SIMULATE_DATA_
#define SIMULATE_DATA_

#include "reconstruct.h"
#include "metadata.h"

class SimulatedData{
    
public:
    bool from_star_file = true;
    std::vector<MetaDataElem> metaDataElems;
    int nr_images;
    std::string output_fn;
    int nr_classes;
public:
    SimulatedData(int _nr_images){
        nr_images = _nr_images;
        std::vector<MetaDataElem> metaDataElems(nr_images);
    }
    ~SimulatedData(){}
    void initialize(std::string _input_mrc_fn,std::string _output_fn){
        
    }
    void randomShiftAndRotate(int offset_range){
        for (int iimage = 0; iimage < nr_images; iimage++) {
            // metaDataElems[iimage].ROT =
        }
    };
    void project(Vol<double>& ,double rot,double tilt,double psi){
        
    }
    void shift(double xoff,double yoff){
        
    }
    void applyCTF(){
        
    }
    void applyNoise(){
        
    }
    
    void unitTest(){
        
    }
};

#endif /* defined(____simulatedata__) */
