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


#include "../src/spider.h"
#include "../src/option.h"
#include "../src/metadata.h"
#include "../src/fft_fftw3.h"
#include "../src/ctf.h"
#include "../src/string.h"

template<typename T>
T* alignedMalloc(size_t size,int align){return (T*)aMalloc(sizeof(T)*size,align);}
template<typename T>
void alignedFree(T* &ptr){aFree(ptr);}

/*** define some functionality module ***/

// write class average
void writeClassAverage(std::string fn_star,std::string fn_mrcs,int K,double pixel_size,double AverageBeta,double AverageAlpha)
{
    std::string mrcs_dir = fn_star.substr(0,fn_star.find_last_of("/")+1);

    MetaDataTable metadata;
    metadata.readFromStar(fn_star);
    
    int N = metadata.numberOfParticles();
    Mrcs::MrcsHead mrcshead;
    Mrcs::readMrcsHead(mrcs_dir+metadata[0].IMAGE.NAME, mrcshead);
    int ori_size = mrcshead.NC;
    int D = ori_size*(ori_size/2+1);
    // prepare data
    float* buffer   = alignedMalloc<float>(ori_size*ori_size, 64);
    float* buffer2  = alignedMalloc<float>(ori_size*ori_size, 64);
    float* TReal    = alignedMalloc<float>(N*D, 64);
    float* TImag    = alignedMalloc<float>(N*D, 64);
    
    FFTWFTransformer fftford_worker(ori_size,ori_size);
    SOAComplexArray<float> fft_out(D);
    CTF ctf;
    float* CTFvalue = alignedMalloc<float>((size_t)N*ori_size*(ori_size/2+1),64);
    double* CTFbuffer = alignedMalloc<double>(ori_size*(ori_size/2+1),64);
    // read image one by one
    for (size_t n = 0; n < N; n++) {
        
        FILE* filehandle = fopen((mrcs_dir+metadata[n].IMAGE.NAME).c_str(),"rb");
        
        // change type to int may cause bug,avoid offset out of range
        long image_id = metadata[n].IMAGE.INDEX;
        
        long offset = (256+(image_id-1)*ori_size*ori_size)*sizeof(float);
        
        fseek(filehandle,offset,SEEK_SET);
        
        if(fread((char*)(buffer),ori_size*ori_size*sizeof(float),1,filehandle) == NULL){
            std::cerr<<"read file failed."<<std::endl;
            EXIT_ABNORMALLY;
        }
        
        Spider::RT_SF(buffer, buffer2, ori_size, ori_size, 0, metadata[n].XOFF, metadata[n].YOFF);
        Spider::RT_SF(buffer2, buffer, ori_size, ori_size, -metadata[n].PSI, 0,0);
        
        fclose(filehandle);
        
        fftford_worker.FourierTransform(buffer, fft_out,false);
        
        // copy the data
        for (int i = 0; i < ori_size*(ori_size/2+1); i++) {
            TReal[n*D+i] = fft_out.real[i];
            TImag[n*D+i] = fft_out.imag[i];
        }
        
        // get CTF
        ctf.setValues(metadata[n].CTF_DEFOCUS_U,metadata[n].CTF_DEFOCUS_V,metadata[n].CTF_DEFOCUS_ANGLE,metadata[n].CTF_VOLTAGE,
                      metadata[n].CTF_CS,metadata[n].CTF_Q0,metadata[n].CTF_BFAC);
        
        // Have the data been CTF phase-flipped?
        bool ctf_phase_flipped = false;
        // Only perform CTF phase-flipping? (default is full amplitude-correction)
        bool only_flip_phases = false;
        // Ignore CTFs until their first peak
        bool intact_ctf_first_peak = false;
        // double pixel_size = 1.74;
        ctf.getFftwImage(CTFbuffer,ori_size,ori_size, ori_size,pixel_size,ctf_phase_flipped, only_flip_phases, intact_ctf_first_peak, true);
        
        for (int i = 0; i < D; i++) {
            CTFvalue[n*D+i] = CTFbuffer[i];
        }
    }
    
    // calculate the classaverage
    // int K = 50;
    // double AverageBeta = 1.8935207e-05;//1;
    // double AverageAlpha = 0.01;
    float* avTReal = alignedMalloc<float>(K*D,64);
    float* avTImag = alignedMalloc<float>(K*D,64);
    
    memset(avTReal, 0, sizeof(float)*K*D);
    memset(avTImag, 0, sizeof(float)*K*D);

    for(int d = 0;d < D;d++){
        for(int k = 0;k < K;k++){
            
            for(int n = 0;n < N;n++){
                if (metadata[n].CLASS == k+1) {
                    avTReal[k*D+d] += metadata[n].PMAX*AverageBeta*CTFvalue[n*D+d]*TReal[n*D+d];
                    avTImag[k*D+d] += metadata[n].PMAX*AverageBeta*CTFvalue[n*D+d]*TImag[n*D+d];
                }
            }
            
            double temp = 0;
            for(int n = 0;n < N;n++)
                if (metadata[n].CLASS == k+1) {
                    temp += metadata[n].PMAX*AverageBeta*CTFvalue[n*D+d]*CTFvalue[n*D+d];
                }
            
            temp += AverageAlpha;
            avTReal[k*D+d] = avTReal[k*D+d]/temp;
            avTImag[k*D+d] = avTImag[k*D+d]/temp;
        }
    }
    
    // write out the result
    float* image_data = alignedMalloc<float>(K*ori_size*ori_size, 64);
    for (int k = 0; k < K; k++) {
        for (int i = 0; i < D; i++) {
            fft_out.real[i] = avTReal[k*D+i];
            fft_out.imag[i] = avTImag[k*D+i];
        }
        fftford_worker.inverseFourierTransform(fft_out, image_data+k*ori_size*ori_size);
    }
    Mrcs::MrcsImages writeImages(image_data,ori_size,K);
    Mrcs::writeMrcsData(fn_mrcs, writeImages);
    std::cout<<"----------- complete write classAverage -----------"<<std::endl;
    
    alignedFree(buffer);
    alignedFree(buffer2);
    alignedFree(TReal);
    alignedFree(TImag);
    alignedFree(CTFvalue);
    alignedFree(CTFbuffer);
    alignedFree(avTReal);
    alignedFree(avTImag);
    alignedFree(image_data);
}

// convert *.mrcs file to *.dat
void convertMrcsToDat(std::string fn_mrcs,std::string fn_dat)
{
    Mrcs::MrcsHead mrcsHead;
    Mrcs::readMrcsHead(fn_mrcs, mrcsHead);
    int N = mrcsHead.NS;
    int image_size = mrcsHead.NC;
    Mrcs::MrcsImages listOfImages(image_size,N);
    std::cout<<"--------- starting reading mrcs image data ---------"<<std::endl;
    Mrcs::readMrcsData(fn_mrcs, mrcsHead, listOfImages);
    std::cout<<"--------- completing reading mrcs image data ---------"<<std::endl;
    std::cout<<"--------- starting writing dat image data ---------"<<std::endl;
    fn_dat = fn_dat.substr(0,fn_dat.find_first_of("."));
    Spider::writeDatData(fn_dat, listOfImages);
    std::cout<<"--------- completing writing dat image data ---------"<<std::endl;
}

// convert *.dat file to *.mrcs
void convertDatToMrcs(std::string fn_dat,std::string fn_mrcs)
{
    Spider::DatHead datHead;
    Spider::readDatHead(fn_dat, datHead);
    int N = (datHead.istack > 0)?datHead.maxim:datHead.nz;
    int image_size = datHead.ny;
    Spider::DatImages listOfImages(image_size,N);
    std::cout<<"--------- starting reading dat image data ---------"<<std::endl;
    Spider::readDatData(fn_dat, datHead, listOfImages);
    std::cout<<"--------- completing reading dat image data ---------"<<std::endl;
    std::cout<<"--------- starting writing mrcs image data ---------"<<std::endl;
    fn_mrcs = fn_mrcs.substr(0,fn_mrcs.find_first_of("."));
    Mrcs::writeMrcsData(fn_mrcs, listOfImages);
    std::cout<<"--------- completing writing mrcs image data ---------"<<std::endl;
}

// select the mrcs file from star file
void selectImages(std::string fn_star,std::string fn_out)
{

    std::string mrcs_dir = fn_star.substr(0,fn_star.find_last_of("/")+1);
    
    // read star file
    MetaDataTable metadata;
    metadata.readFromStar(fn_star);
    
    // read mrcs head
    Mrcs::MrcsHead mrcsHead;
    readMrcsHead(mrcs_dir+metadata[0].IMAGE.NAME, mrcsHead);
    size_t N = metadata.numberOfParticles();
    int ori_size = mrcsHead.NC;

    float *selectImageData = alignedMalloc<float>((size_t)N*ori_size*ori_size, 64);
    // read image one by one
    for (size_t n = 0; n < N; n++) {
        
        FILE* filehandle = fopen((mrcs_dir+metadata[n].IMAGE.NAME).c_str(),"rb");
        
        // change type to int may cause bug,avoid offset out of range
        long image_id = metadata[n].IMAGE.INDEX;
        // change the image id and fn
        metadata[n].IMAGE.INDEX = n+1;
        metadata[n].IMAGE.NAME = fn_out+".mrcs";
        
        long offset = (256+(image_id-1)*ori_size*ori_size)*sizeof(float);
        
        fseek(filehandle,offset,SEEK_SET);
        
        if(fread((char*)(selectImageData+n*ori_size*ori_size),ori_size*ori_size*sizeof(float),1,filehandle) == NULL){
            std::cerr<<"read file failed."<<std::endl;
            EXIT_ABNORMALLY;
        }
        
        fclose(filehandle);
    }
    
    Mrcs::MrcsImages listOfImages(selectImageData,ori_size,N);
    
    Mrcs::writeMrcsData(fn_out, listOfImages);
    metadata.writeToStar(fn_out);
    
    alignedFree(selectImageData);
}

// adjust the mrcs file from star file
void adjustImages(std::string fn_star,std::string fn_out)
{
    
    std::string mrcs_dir = fn_star.substr(0,fn_star.find_last_of("/")+1);
    
    // read star file
    MetaDataTable metadata;
    metadata.readFromStar(fn_star);
    // read mrcs head
    Mrcs::MrcsHead mrcsHead;
    readMrcsHead(mrcs_dir+metadata[0].IMAGE.NAME, mrcsHead);
    int ori_size = mrcsHead.NC;
    int N = metadata.numberOfParticles();

    float* selectImageData = alignedMalloc<float>((size_t)N*ori_size*ori_size,64);
    float* buffer = alignedMalloc<float>(ori_size*ori_size,64);
    // read image one by one
    for (size_t n = 0; n < N; n++) {
        
        FILE* filehandle = fopen((mrcs_dir+metadata[n].IMAGE.NAME).c_str(),"rb");
        
        // change type to int may cause bug,avoid offset out of range
        long image_id = metadata[n].IMAGE.INDEX;
        // change the image id and fn
        metadata[n].IMAGE.INDEX = n+1;
        metadata[n].IMAGE.NAME = fn_out+".mrcs";
        
        long offset = (256+(image_id-1)*ori_size*ori_size)*sizeof(float);
        
        fseek(filehandle,offset,SEEK_SET);
        
        if(fread((char*)(selectImageData+n*ori_size*ori_size),ori_size*ori_size*sizeof(float),1,filehandle) == NULL){
            std::cerr<<"read file failed."<<std::endl;
            EXIT_ABNORMALLY;
        }
        
        fclose(filehandle);
        
        Spider::RT_SF(selectImageData+n*ori_size*ori_size, buffer, ori_size, ori_size, 0, metadata[n].XOFF, metadata[n].YOFF);
        Spider::RT_SF(buffer, selectImageData+n*ori_size*ori_size, ori_size, ori_size, -metadata[n].PSI, 0,0);
        metadata[n].XOFF = 0;
        metadata[n].YOFF = 0;
        metadata[n].PSI = 0;
        if (n % 2000 == 0) std::cout<<"complete adjust "<<n<<std::endl;
    }
    
    Mrcs::MrcsImages listOfImages(selectImageData,ori_size,N);
    
    Mrcs::writeMrcsData(fn_out, listOfImages);
    metadata.writeToStar(fn_out);
    
    alignedFree(selectImageData);
    alignedFree(buffer);
}

// apply low pass filter to images
void applyFilter2Images(std::string fn_mrcs,std::string fn_out,double filter)
{
    Mrcs::MrcsHead mrcsHead;
    Mrcs::readMrcsHead(fn_mrcs, mrcsHead);
    int N = mrcsHead.NS;
    int image_size = mrcsHead.NC;
    Mrcs::MrcsImages listOfImages(image_size,N);
    Mrcs::readMrcsData(fn_mrcs, mrcsHead, listOfImages);
    
    Spider::applyFilter(listOfImages, filter);
    
    Mrcs::writeMrcsData(fn_out, listOfImages);
}

// clean the empty classaverage
void cleanEmptyClass(std::string fn_in,std::string fn_out)
{
    // get the old metadata table
    MetaDataTable metaDataTable;
    std::vector<int> imageNumberPerClass;
    metaDataTable.readFromStar(fn_in+".star");
    metaDataTable.statClassImageNumber(imageNumberPerClass,true);
    int nonEmptyClassNumber = 0;
    for (int iclass = 0; iclass < imageNumberPerClass.size(); iclass++) {
        if (imageNumberPerClass[iclass] != 0) nonEmptyClassNumber++;
    }
    // get the image size
    std::string fn_mrcs = fn_in+".mrcs";
    Mrcs::MrcsHead mrcsHead;
    Mrcs::readMrcsHead(fn_mrcs, mrcsHead);
    int image_size = mrcsHead.NC;
    // set the new data
    std::vector<MetaDataElem> metaDataElems(metaDataTable.numberOfParticles());
    int particleIndex = 0;
    float* nonEmptyClassAverage = (float*)aMalloc(sizeof(float)*nonEmptyClassNumber*image_size*image_size,64);
    int nonEmptyClassIndex = 0;
	//
    float buffer[image_size*image_size];
    FILE* mrcsFile = fopen(fn_mrcs.c_str(),"rb");
    ERROR_CHECK(NULL == mrcsFile,"previous classaverage file should be put in same folder and has same name as star file.");
    
    for (int iclass = 0; iclass < imageNumberPerClass.size(); iclass++) {
        if (imageNumberPerClass[iclass] == 0) continue;
        long offset = (256+iclass*image_size*image_size)*sizeof(float);
        
        fseek(mrcsFile,offset,SEEK_SET);
        if(fread((char*)buffer,image_size*image_size*sizeof(float),1,mrcsFile) == NULL)
            ERROR_REPORT("read mrcs data failed.");
        
        // copy the metadata
        for (int i = 0; i < metaDataTable.numberOfParticles(); i++) {
            if (metaDataTable[i].CLASS == iclass+1) {
                metaDataElems[particleIndex] = metaDataTable[i];
                metaDataElems[particleIndex].CLASS = nonEmptyClassIndex+1;
                particleIndex++;
            }
        }
        // copy the classaverage
        memcpy(nonEmptyClassAverage+nonEmptyClassIndex*image_size*image_size, buffer, image_size*image_size*sizeof(float));
        nonEmptyClassIndex++;
    }
    
    fclose(mrcsFile);
	
    // write the new star file
    MetaDataTable newMetaDataTable;
    newMetaDataTable.readFromMetaDataElements(metaDataElems);
    newMetaDataTable.writeToStar(fn_out);
    
    // write the new mrcs file
    Mrcs::MrcsImages listOfImages(nonEmptyClassAverage,image_size,nonEmptyClassIndex);
    listOfImages.write(fn_out);
    
    
}

void otherModules(){
    
#if 0 //subtract data
    int nr_global_images = 31990;
    std::string mrcs_dir="../../fc_sub/";
    int ori_size = 256;
    
    std::cout<<"read data1."<<std::endl;
    MyMetaData* metadata1 = (MyMetaData*)aMalloc(sizeof(MyMetaData)*nr_global_images,64);
    std::vector<std::string> metadata_fn1;
    std::string star_fn1="../../fc_sub/sub2.star";
    
    readStarData(star_fn1,metadata1,metadata_fn1,nr_global_images);
    
    float *buffer = (float*)aMalloc(sizeof(float)*ori_size*ori_size,64);
    
    float* data1 = (float*)aMalloc(sizeof(float)*nr_global_images*ori_size*ori_size,64);
    for (int n = 0; n < nr_global_images; n++) {
        
        FILE* filehandle = fopen((mrcs_dir+metadata_fn1[n]).c_str(),"rb");
        
        int image_id = metadata1[n].IMAGE_NAME;
        
        long offset = (256+(image_id-1)*ori_size*ori_size)*sizeof(float);
        
        fseek(filehandle,offset,SEEK_SET);
        
        if(fread((char*)buffer,ori_size*ori_size*sizeof(float),1,filehandle) == NULL){
            std::cerr<<"read file failed."<<std::endl;
            EXIT_ABNORMALLY;
        }
        
        memcpy(data1+n*ori_size*ori_size,buffer,sizeof(float)*ori_size*ori_size);
        
        fclose(filehandle);
    }
    
    std::cout<<"read data2."<<std::endl;
    MyMetaData* metadata2 = (MyMetaData*)aMalloc(sizeof(MyMetaData)*nr_global_images,64);
    std::vector<std::string> metadata_fn2;
    std::string star_fn2="../../fc_sub/unfil2.star";
    
    readStarData(star_fn2,metadata2,metadata_fn2,nr_global_images);
    
    float* data2 = (float*)aMalloc(sizeof(float)*nr_global_images*ori_size*ori_size,64);
    
    for (int n = 0; n < nr_global_images; n++) {
        
        FILE* filehandle = fopen((mrcs_dir+metadata_fn2[n]).c_str(),"rb");
        
        int image_id = metadata2[n].IMAGE_NAME;
        
        long offset = (256+(image_id-1)*ori_size*ori_size)*sizeof(float);
        
        fseek(filehandle,offset,SEEK_SET);
        
        if(fread((char*)buffer,ori_size*ori_size*sizeof(float),1,filehandle) == NULL){
            std::cerr<<"read file failed."<<std::endl;
            EXIT_ABNORMALLY;
        }
        
        memcpy(data2+n*ori_size*ori_size,buffer,sizeof(float)*ori_size*ori_size);
        
        fclose(filehandle);
    }
    
    for (int i = 0; i < nr_global_images*ori_size*ori_size; i++)
    {
        data2[i] -= data1[i];
    }
    
    int mrcsHead[256];
    mrcsHead[0]=mrcsHead[1]=ori_size;mrcsHead[2]=nr_global_images;mrcsHead[3]=2;
    writeMrcsData(mrcs_dir+"fc_sub2",mrcsHead,data2);
    
#endif
    
}


int main(int argc, char * argv[])
{
    
    if (argc < 2) {
        std::cout<<"you can choice this functionality : "<<std::endl;
        std::cout<<"-classaverage : Compute class averaging from a given file."<<std::endl;
        std::cout<<"-convert      : Convert image data file from given format (SPIDER form .dat or RELION form .mrcs) to formats (SPIDER form .dat or RELION form .mrcs) you want."<<std::endl;
        std::cout<<"-select       : Gather images from plenty of images into one stack."<<std::endl;
        std::cout<<"-adjust       : Shift and rotate images based on translations and rotation angle in star file."<<std::endl;
        std::cout<<"-applyfilter  : Perform low-pass filtering."<<std::endl;
        std::cout<<"-clean        : clean the empty class."<<std::endl;
        std::cout<<"for each functionality,using -help to see how to use,"<<std::endl;
        std::cout<<"example : rome_tool -classaverage -help."<<std::endl;
        EXIT_ABNORMALLY;
    }
    // write classaverage
    if (std::string(argv[1]) == "-classaverage")
    {
        Option option;
        option.addOption("-i", "Input metadata file with images(*.star)");
        option.addOption("-o", "Output metadata");
        option.addOption("-K", "Number of classes!!!!! please use Uppercase —K instead of lowercase -k !!!!!");
        option.addOption("-angpix","Pixel size (in Angstroms)!!!!! please use -angpix instead of -pixel !!!!!");
        option.addOption("-averageBeta", "The variance of noise when doing weighted class averaging","1");
        option.addOption("-averageAlpha", "The variance of prior model Gaussian distribution when doing weighted class averaging","0.01");
        if (argc < 4) {
            std::cout<<"rome_tool -classaverage usage : "<<std::endl;
            option.printHelp();
            EXIT_ABNORMALLY;
        }

        // remove the first argument and read the remaining argv
        option.readCommandLine(argc-1, argv+1);
        
        std::string fn_star = pathRemoveSuffix(option.getStrOption("-i"))+".star";
        std::string fn_mrcs = pathRemoveSuffix(option.getStrOption("-o"));
        int K = option.getIntOption("-K");
        double pixel_size = option.getFloatOption("-angpix");
        double averageBeta = option.getFloatOption("-averageBeta");
        double averageAlpha = option.getFloatOption("-averageAlpha");
        option.printValue();
        
        writeClassAverage(fn_star,fn_mrcs,K,pixel_size,averageBeta,averageAlpha);
    }
    
    // convert *.mrcs to *.dat or *.dat to *.mrcs
    if (std::string(argv[1]) == "-convert")
    {
        Option option;
        option.addOption("-i", "Input file name(*.mrcs or *.dat)");
        option.addOption("-o", "Output file name(*.mrcs or *.dat)");
        if (argc < 4) {
            std::cout<<"rome_tool -convert usage : "<<std::endl;
            option.printHelp();
            EXIT_ABNORMALLY;
        }
        
        // remove the first argument and read the remaining argv
        option.readCommandLine(argc-1, argv+1);
        
        std::string fn_input = option.getStrOption("-i");
        std::string fn_output = option.getStrOption("-o");
        option.printValue();
        std::string fn_input_suffix = pathGetSuffix(fn_input);
        std::string fn_output_suffix = pathGetSuffix(fn_output);

        if ((fn_input_suffix == ".mrcs") && (fn_output_suffix == ".dat")) {
            //
            convertMrcsToDat(fn_input, fn_output);
        }
        else if((fn_input_suffix == ".dat") && (fn_output_suffix == ".mrcs")){
            //
            convertDatToMrcs(fn_input, fn_output);
        }
        else
            std::cout<<"donnot support "<<fn_input<<" to "<<fn_output<<std::endl;
    }
    
    // select mrcs file from star file
    if (std::string(argv[1]) == "-select")
    {
        Option option;
        option.addOption("-i", "Input metadata file with images(*.star)");
        option.addOption("-o", "Output metadata");
        if (argc < 4){
            std::cout<<"rome_tool -select usage : "<<std::endl;
            option.printHelp();
            EXIT_ABNORMALLY;
        }
        // remove the first argument and read the remaining argv
        option.readCommandLine(argc-1, argv+1);
        std::string fn_star = pathRemoveSuffix(option.getStrOption("-i"))+".star";
        std::string fn_out  = pathRemoveSuffix(option.getStrOption("-o"));
        option.printValue();
        selectImages(fn_star, fn_out);
    }
    
    // adjust the image file
    if (std::string(argv[1]) == "-adjust") {
        Option option;
        option.addOption("-i", "Input metadata file with images(*.star)");
        option.addOption("-o", "Output metadata");
        if (argc < 4){
            std::cout<<"rome_tool -adjust usage : "<<std::endl;
            option.printHelp();
            EXIT_ABNORMALLY;
        }
        // remove the first argument and read the remaining argv
        option.readCommandLine(argc-1, argv+1);
        std::string fn_star = pathRemoveSuffix(option.getStrOption("-i"))+".star";
        std::string fn_out  = pathRemoveSuffix(option.getStrOption("-o"));
        option.printValue();
        adjustImages(fn_star, fn_out);
    }
    
    // add low pass filter
    if (std::string(argv[1]) == "-applyfilter") {
        Option option;
        option.addOption("-i", "Input file name(*.mrcs)");
        option.addOption("-o", "Output file name");
        option.addOption("-filter", "Filter radius in frequency");
        if (argc < 4){
            std::cout<<"rome_tool -filter usage : "<<std::endl;
            option.printHelp();
            EXIT_ABNORMALLY;
        }
        // remove the first argument and read the remaining argv
        option.readCommandLine(argc-1, argv+1);
        std::string fn_mrcs = pathRemoveSuffix(option.getStrOption("-i"))+".mrcs";
        std::string fn_out  = pathRemoveSuffix(option.getStrOption("-o"));
        double filter       = option.getFloatOption("-filter");
        option.printValue();
        
        applyFilter2Images(fn_mrcs,fn_out,filter);
    }
    
    // clean the empty class
    if (std::string(argv[1]) == "-clean") {
        Option option;
        option.addOption("-i", "Input file name(*.star)");
        option.addOption("-o", "Output file name(*.star)");
        if (argc < 4){
            std::cout<<"rome_tool -clean usage : "<<std::endl;
            option.printHelp();
            EXIT_ABNORMALLY;
        }
        // remove the first argument and read the remaining argv
        option.readCommandLine(argc-1, argv+1);
        option.printValue();
        
        std::string fn_in = pathRemoveSuffix(option.getStrOption("-i"));
        std::string fn_out  = pathRemoveSuffix(option.getStrOption("-o"));
        
        cleanEmptyClass(fn_in,fn_out);
    }
    
    return 0;
}

