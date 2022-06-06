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

#include "util.h"		// used for building precompiled headers on Windows

#include "mrcs.h"

namespace Mrcs {
    
    // read Mrcshead from *.mrcs file
    void readMrcsHead(std::string fn_mrcs,MrcsHead& mrcsHead){
        
        std::ifstream inFile;
        inFile.open(fn_mrcs.c_str(),std::ios::in | std::ios::binary);
        
        // pointer to the MrcsHead structure
        int *mrcsHead_aux = &mrcsHead.NC;
        
        inFile.read((char*)mrcsHead_aux,256*sizeof(int));
        
        inFile.close();
    }
    
    // read MrcHead from *.mrc
    bool readMrcHead(std::string fn_mrc,MrcsHead& refHead,int size){
        FILE* refFile = fopen((fn_mrc).c_str(),"rb");
        ERROR_CHECK(refFile==NULL, "cannot open reference file.");
        bool is_ieee_le = true;
        int *refHead_aux = &refHead.NC;
        fread((char*)refHead_aux,256*sizeof(float),1,refFile);
        
        if (refHead.NC != size || refHead.NR != size || refHead.NS != size) {
            // the difference maybe caused by ieee-be or ieee-le format
            for (int i = 0; i < 256; i++)// ieee-be to ieee-le
                refHead_aux[i] = SWAP32(refHead_aux[i]);
            is_ieee_le = false;
            if (refHead.NC != size || refHead.NR != size || refHead.NS != size){
                std::cerr<<"reference file wrong: size:" << size << " diffs from "
                << "ref_head.NC:" << refHead.NC << " "
                << "ref_head.NR:" << refHead.NR << " "
                << "ref_head.NS:" << refHead.NS << std::endl;
                ERROR_REPORT("");
            }
        }
        
        fclose(refFile);
        return is_ieee_le;
    }
    
    // set mrcHead
    void setMrcHead(MrcsHead& refHead,double anpix,int size){
        refHead.NC = refHead.NR = refHead.NS = size;
        refHead.MODE = 2;
        refHead.NX = refHead.NY = refHead.NZ = size;
        refHead.Alpha = refHead.Beta = refHead.Gamma = 90;
        refHead.MAPC = 1;refHead.MAPR = 2;refHead.MAPS = 3;
        refHead.X_length = refHead.Y_length = refHead.Z_length = anpix*size;
    }
    
    // read Mrcs image Data
    void readMrcsData(std::string fn_mrcs,const MrcsHead& mrcsHead,FloatImages &mrcsImages){
        
        FILE* filehandle = fopen(fn_mrcs.c_str(),"r");
        
        MrcsHead head;
        int* head_aux = &head.NC;
        fread((char*)head_aux,sizeof(int)*256,1,filehandle);
        
        int nx = head.NC;
        int ny = head.NR;
        int N = head.NS;
        int mode = head.MODE;
        int numberOfImages = mrcsImages.nr_images();
        int imageSize = mrcsImages.imageSide();
        
        // check data
        ERROR_CHECK(nx != mrcsHead.NC || ny != mrcsHead.NR || N != mrcsHead.NS || mode != mrcsHead.MODE,
                   "mrcs file head not match.");
        ERROR_CHECK(mode != 2, "mode is not 2,ROME only support float mrcs data type.");
        ERROR_CHECK(N != numberOfImages, "not enough images in mrcs file.");
        ERROR_CHECK(nx != imageSize || ny != imageSize, "listofmrcsImage size wrong.");
        
        long dataSize = ny*nx*sizeof(float);
        std::vector<float> buffer(ny*nx);
        // NOTICE : sometimes the mrcs data from RELION may has nan number datatype
        // this is a bug in RELION
        int nan_count = 0;
        // read mrcs image one by none
        for (int iimage = 0; iimage < numberOfImages; iimage++) {
            
            long offset   = (256+iimage*nx*ny)*sizeof(int);
            
            ERROR_CHECK(fseek(filehandle,offset,SEEK_SET) == -1, "fseek mrcs data error.");
            
            fread((char*)&buffer[0],dataSize,1,filehandle);
            
            auto* image_data = mrcsImages.image_ptr(iimage);
            for (int i = 0; i < ny*nx; i++) {
                if (std::isnan(buffer[i])) nan_count++;
                image_data[i] = buffer[i];
            }
        }
        if (nan_count > 0) {
            std::cerr<<"########NOTICE : some nan number in your"<<fn_mrcs<<" file.##########"<<std::endl;
        }
        fclose(filehandle);
    }
    
    // write imagedata to mrcs file
    void writeMrcsData(std::string fn_mrcs,FloatImages &mrcsImages){
        
        std::string filename_mrcs = pathRemoveSuffix(fn_mrcs)+".mrcs";
        
        std::ofstream outFile;
        outFile.open(filename_mrcs.c_str(),std::ios::out | std::ios::binary);
        
        MrcsHead mrcsHead;
        int numberOfImages = mrcsImages.nr_images();
        int imageSize = mrcsImages.imageSide();
        mrcsHead.NC = imageSize;
        mrcsHead.NR = imageSize;
        mrcsHead.NS = numberOfImages;
        mrcsHead.MODE = 2;
        
        // set the mrcsHeader
        int* mrcsHead_aux = &mrcsHead.NC;
        outFile.write((char*)mrcsHead_aux,256*sizeof(int));
        
        // write image one by one
        int imageSize2 = imageSize*imageSize;
        std::vector<float> buffer(imageSize*imageSize);
        for (int iimage = 0; iimage < numberOfImages; iimage++) {
            auto* image_data = mrcsImages.image_ptr(iimage);
            for(int i = 0;i < imageSize2;i++)
                buffer[i] = image_data[i];
            if(!outFile.write((char*)&buffer[0],imageSize2*sizeof(float)))
                std::cout<<"write data failed."<<std::endl;
        }

        outFile.close();
        
    }
    
    // read mrc file
    void readMrcData(std::string fn_mrc,FloatVolume &mrcVolume){
        //
        int volumeSize = mrcVolume.volumeSide();
        float* volumeData = mrcVolume.ptr();
        MrcsHead ref_head;
        bool is_ieee_le = readMrcHead(fn_mrc, ref_head, volumeSize);
        //
        FILE* refFile = fopen((fn_mrc).c_str(),"rb");
        ERROR_CHECK(refFile==NULL, "cannot open reference file.");
        fseek(refFile, 256*4+ref_head.NSYMBT, SEEK_SET);
        fread((char*)&volumeData[0],volumeSize*volumeSize*volumeSize*sizeof(float),1,refFile);
        
        int volumeSize3 = volumeSize*volumeSize*volumeSize;

        if (is_ieee_le) {
            //
        }
        else{
            float temp;
            for (int i = 0; i < volumeSize3; i++) {
                ((int*)&temp)[0] = SWAP32(((int*)&volumeData[0])[i]);
                volumeData[i] = temp;
            }
        }
        
        fclose(refFile);
    }
    
    // write 3D volume to mrc file
    void writeMrcData(std::string fn_mrcs,FloatVolume &mrcVolume,MrcsHead mrcHead,double anpix){
        //
        int volumeSize = mrcVolume.volumeSide();
        mrcHead.NSYMBT = 0;
        assert(mrcHead.MODE == 2);
        assert(mrcHead.NX == mrcHead.NY);
        assert(mrcHead.NY == mrcHead.NZ);
        assert(mrcHead.NZ == volumeSize);
        
        int* mrcHead_aux = &mrcHead.NC;
        
        // write image one by one
        size_t volumeSize3 = volumeSize*volumeSize*volumeSize;
        std::vector<float> buffer(volumeSize*volumeSize*volumeSize);

        std::string filename_mrcs = pathRemoveSuffix(fn_mrcs)+".mrc";
        
        std::ofstream outFile;
        outFile.open(filename_mrcs.c_str(),std::ios::out | std::ios::binary);
        outFile.write((char*)mrcHead_aux,256*sizeof(int));
        
        auto* volume_data = mrcVolume.ptr();
        for(int i = 0;i < volumeSize3;i++)
            buffer[i] = volume_data[i];
        if(!outFile.write((char*)&buffer[0],volumeSize3*sizeof(float)))
            std::cout<<"write data failed."<<std::endl;
        
        outFile.close();
        
    }
    
    
    void UnitTest(std::string fn_mrcs){
        
        std::cout<<"testing Mrcs."<<std::endl;
        MrcsHead mrcsHead;
        std::cout<<"print initialize info for mrcs head.--------"<<std::endl;
        std::cout<<"sizeof(mrcsHead) = "<<sizeof(mrcsHead)<<",should be 1025=256*4."<<std::endl;
        mrcsHead.print(std::cout);
        
        std::cout<<"read mrcs head.--------"<<std::endl;
        readMrcsHead(fn_mrcs, mrcsHead);
        mrcsHead.print(std::cout);
        
        std::cout<<"testing mrcs head copy function--------"<<std::endl;
        MrcsHead mrcsHead2;
        mrcsHead2 = mrcsHead;
        mrcsHead2.print(std::cout);
        
        std::cout<<"read mrcs data.--------"<<fn_mrcs<<std::endl;
        int N = mrcsHead.NS;
        int size = mrcsHead.NC;
        MrcsImages listOfImages(size,N);
        
        readMrcsData(fn_mrcs, mrcsHead, listOfImages);
        
        std::cout<<"write mrcs data : "<<fn_mrcs<<".mrcs"<<"--------"<<std::endl;
        std::cout<<"using rome_vier.py to see the different."<<std::endl;
        writeMrcsData(fn_mrcs, listOfImages);
    }
} // namespace Mrcs
