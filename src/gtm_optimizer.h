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

#ifndef GTM_OPTIMIZER_H_
#define GTM_OPTIMIZER_H_

#include "./util.h"

#define DONT_USE_OFFLOAD
#ifdef DONT_USE_OFFLOAD
#define MIC_OFFLOAD_ATTRIBUTES 
#else
#include "offload.h"
#define MIC_OFFLOAD_ATTRIBUTES __attribute__ (( target (mic)))
#endif

#include "./pca_optimizer.h"
#include "./spider.h"
#include "./initialize.h"
#include "./time.h"
#include "./mpi.h"
#include "./mrcs.h"
#include "./metadata.h"
#include "./progressbar.h"
#include "./sampling.h"
#include "./ctf.h"
#include "./basis.h"
#include "./string.h"

#define DEBUGGTM

// We have array indices > 2^31.  This will adversely effect vectorization
#define INDEXSZ(x) ((size_t)(x))
// We have array indices < 2^31
// #define INDEXSZ(x) ((int)(x))

namespace GTMoptimizer
{
    // MPI Variable -------------------
#ifdef USEMPI
    extern char nodeName[MPI_MAX_PROCESSOR_NAME];
#else
    extern char nodeName[4096];
#endif
    
    // node,total number of nodes in cluster
	extern int node,nodes,nodeNameLen,messageTag;
    
    // GTM divides the problem up in the cluster in subsets of the dimension D,
    // with each rank getting D/number of ranks in this job of the work,
    // as specified by the following three parameters
    extern MIC_OFFLOAD_ATTRIBUTES int d_start, d_end, d_interval;
    
    // Maximum number of threads allowed on the host nodes
    extern MIC_OFFLOAD_ATTRIBUTES int maxthreads;
    
    // Number of MIC cards connected to this host node (not all
    // ranks need to have Intel(R) Xeon Phi(TM) coprocessors - this is the
    // number we will use if we find them)
    extern int numMIC;
    
    // The percentage of job put to compute in mic card
    extern double loadMIC;
    
    // Memory space we are allowed to use on the MIC card
    extern size_t offloadlimit;
    
    // How much space we have for everything except R on the coprocessor
    extern size_t micavail;
    
    // The maximum fraction of d_interval the coprocessor can traverse 
    // given the values of N, M, and K
    extern int estnodesize;
    
    // Maximum number of MIC cards the user allows to use on any rank
    // that has attached Intel(R) Xeon Phi(TM) coprocessors
    extern int maxXeonPhi;
    
    // Maximum Number of threads to use on the MIC card
    extern MIC_OFFLOAD_ATTRIBUTES int numXeonPhithreads;
    
    // TODO
    extern MIC_OFFLOAD_ATTRIBUTES
        __attribute__((aligned(64))) double* PHI;
    
    // Maximum number of conformations
    extern MIC_OFFLOAD_ATTRIBUTES int K;
    
    // Maximum number of basis functions
    extern MIC_OFFLOAD_ATTRIBUTES int M;
    
    // Total number of images
    extern MIC_OFFLOAD_ATTRIBUTES int N;
    
    // Total number of "pixels" in each image represented in
    // Fourier space (elements in Real or Imaginary components)
    extern MIC_OFFLOAD_ATTRIBUTES int D;
    
    // The real part of the image data in Fourier space
    extern MIC_OFFLOAD_ATTRIBUTES
        __attribute__((aligned(64))) float* TReal;
    
    // The imaginary part of the image data in Fourier space
    extern MIC_OFFLOAD_ATTRIBUTES
        __attribute__((aligned(64))) float* TImag;

    // CTF function for each pixel of each image
    extern MIC_OFFLOAD_ATTRIBUTES
        __attribute__((aligned(64))) float* CTFvalue;
    
    // The real part of the weights of mapping
    extern MIC_OFFLOAD_ATTRIBUTES
        __attribute__((aligned(64))) float* WReal;
    
    // The imaginary part of the weights of mapping
    extern MIC_OFFLOAD_ATTRIBUTES
        __attribute__((aligned(64))) float* WImag;
    
    // The real part of the points mapped from latent space to
    // a point in the data space
    extern MIC_OFFLOAD_ATTRIBUTES
        __attribute__((aligned(64))) double* YReal;
    
    // The imaginary part of the points mapped from the latent space 
    // to a point in the data space
    extern MIC_OFFLOAD_ATTRIBUTES
        __attribute__((aligned(64))) double* YImag;
    
    // The posterior probability matrix
    extern MIC_OFFLOAD_ATTRIBUTES
        __attribute__((aligned(64))) double* R;
    
    // The variance of the isotropic Gaussian prior distribution over W
    extern MIC_OFFLOAD_ATTRIBUTES
        __attribute__((aligned(64))) float* Alpha;
    
    // The mean variance of noise
    extern MIC_OFFLOAD_ATTRIBUTES double AverageBeta;
    
    // The mean variance of the isotropic Gaussian prior distribution over W
    extern MIC_OFFLOAD_ATTRIBUTES double AverageAlpha;
    
    // The real part of classaverage
    extern MIC_OFFLOAD_ATTRIBUTES
        __attribute__((aligned(64))) float *avTReal;
    
    // The imaginary part of classaverage
    extern MIC_OFFLOAD_ATTRIBUTES
        __attribute__((aligned(64))) float *avTImag;
    
    // Temporary Matrix used in the mkl_solve,size M x sizeof(int)
    extern MIC_OFFLOAD_ATTRIBUTES
        __attribute__((aligned(64))) int *ipiv;
    
    // The number of basis exp functions
    extern int Mnl;
    
    // The distance between data point and mapping point
    extern double deviation;
    
    // The pixel size for each image (angstroms)
    extern double pixel_size;
    
    //
    extern bool do_ctf_correction;
    // Have the data been CTF phase-flipped?
    extern bool only_flip_phases;
    // Only perform CTF phase-flipping? (default is full amplitude-correction)
    extern bool ctf_phase_flipped;
    // Ignore CTFs until their first peak
    extern bool intact_ctf_first_peak;
    
    // metadata for each collection of images
    extern MetaDataTable metadata;
    
    // The header for the image data file
	extern int mrcsHead[256];
    
    // Various strings used for reporting
    extern char timerstring[16384];
    extern char buffer[256];
    
    // The size of the images in each dimension (D*2 is the total number 
    // of pixels in each image)
    extern int picdim;
    
    //
    extern std::string star_fn;
    
    // output file name
    extern std::string result_fn;
    
    // reference file name 
    extern std::string ref_fn;
    
    //
    extern bool init_by_pca;
    
    // Total number iterations and current iteration
    extern int iter, nr_iter;
    
    // do fine search for each class
    class FineSearch{
//#define DEBUG_FINE_SEARCH
    public:
        bool doFineSearch;
    private:
        // image and metadata data
        float* classAverage;
        int image_size;
        std::vector<MetaDataElem> metaDataElems;
        int classIndexAccumulator;
        int imageIndexAccumulator;
        // the number of image to each class
        std::vector<int> imageNumberPerClass;
        int imageNumber;
        int classNumber;
        // the number of class to classify for each class
        std::vector<int> fineClassNumberPerClass;
        int fineClassNumber;
        // selected the class which has more than 20 images
        std::vector<int> selectedClass;
        int selectedClassNumber;
        int selectedClassIndex;
        int currentClass;
    public:
        FineSearch()
        {
            doFineSearch = false;
            fineClassNumber = 0;
            classIndexAccumulator = 0;
            imageIndexAccumulator = 0;
            classNumber = 0;
            imageNumber = 0;
            selectedClassNumber = 0;
            selectedClassIndex = -1;
            currentClass = 0;
        }
        ~FineSearch(){destroyFineSearch();}
        
        void setup(std::string star_fn,std::string config_fn)
        {
            doFineSearch = true;
            MetaDataTable metaDataTable;
            metaDataTable.readFromStar(star_fn);
            // set the image number for each class
            metaDataTable.statClassImageNumber(imageNumberPerClass);
            // set the class number
            classNumber = imageNumberPerClass.size();
            // set the fine class number for each class
            fineClassNumberPerClass.resize(classNumber,0);
            if (config_fn == "default")
            {
                for (int iclass = 0; iclass < classNumber; iclass++)
                {
                    // NOTE : donot do classification if image number < 20
                    if (imageNumberPerClass[iclass] < 20){
                        fineClassNumberPerClass[iclass] = 1;
                        imageNumber += imageNumberPerClass[iclass];
                        fineClassNumber += fineClassNumberPerClass[iclass];
                        selectedClass.push_back(iclass+1);
                    }
                    else{
                        fineClassNumberPerClass[iclass] = imageNumberPerClass[iclass]/10 + 2;
                        imageNumber += imageNumberPerClass[iclass];
                        fineClassNumber += fineClassNumberPerClass[iclass];
                        selectedClass.push_back(iclass+1);
                    }
#ifdef DEBUG_FINE_SEARCH
                    MASTERNODE
                    {
                        std::cout<<" set up iclass : "<<iclass<<std::endl;
                        std::cout<<" fineClassNumberPerClass[iclass] : "<<fineClassNumberPerClass[iclass]<<std::endl;
                        std::cout<<" imageNumberPerClass[iclass] : "<<imageNumberPerClass[iclass]<<std::endl;
                        std::cout<<" fineClassNumber : "<<fineClassNumber<<std::endl;
                    }
#endif
                }
            }
            else
            {
                // read config file
                std::ifstream configFile(config_fn.c_str());
                std::string line;
                while (!configFile.eof()) {
                    getline(configFile,line);
                    if(line.find(',') == std::string::npos)
                        continue;
                    if(line.find('#') == std::string::npos){
                        std::string class_str = line.substr(0,line.find_last_of(","));
                        std::string fineClassNumber_str = line.substr(line.find_last_of(",")+1);
                        int iclass = atoi(class_str.c_str())-1;
                        int fineClassNumberIClass = atoi(fineClassNumber_str.c_str());
                        fineClassNumberPerClass[iclass] = fineClassNumberIClass;
                        selectedClass.push_back(iclass+1);
                        fineClassNumber += fineClassNumberPerClass[iclass];
                        imageNumber += imageNumberPerClass[iclass];
                    }
                }
            }
            // alloc the data
            Mrcs::MrcsHead mrcsHead;
            Mrcs::readMrcsHead(pathGetDir(star_fn)+metaDataTable[0].IMAGE.NAME, mrcsHead);
            image_size = mrcsHead.NC;
            classAverage = (float*)aMalloc(sizeof(float)*fineClassNumber*image_size*image_size,64);
            memset(classAverage, 0, sizeof(float)*fineClassNumber*image_size*image_size);
            //
            metaDataElems.resize(imageNumber);
            //
            MASTERNODE std::cout<<" done.... set up fine searching. "<<std::endl;
        }
        void destroyFineSearch(){
            if(doFineSearch) {
                aFree(classAverage);
                metaDataElems.resize(0);
                imageNumberPerClass.resize(0);
                fineClassNumberPerClass.resize(0);
                selectedClass.resize(0);
                doFineSearch = false;
            }
        }
        bool searchNext(std::string star_fn)
        {
            if (doFineSearch)
            {
                // get the next searching class
                selectedClassIndex ++;
                currentClass = selectedClass[selectedClassIndex];
                ERROR_CHECK(currentClass != selectedClassIndex+1, "throw some class...");
                // accumulate the previous selected class fine class number
                if(selectedClassIndex > 0){
                    classIndexAccumulator += fineClassNumberPerClass[selectedClass[selectedClassIndex-1]-1];
                    imageIndexAccumulator += imageNumberPerClass[selectedClass[selectedClassIndex-1]-1];
                }
                assert(classIndexAccumulator < fineClassNumber);
                assert(imageIndexAccumulator < imageNumber);
                // escaping classify if fineClass equal 1
                if (fineClassNumberPerClass[currentClass-1] == 1)
                {
                    MASTERNODE std::cout<<"escaping doing fineSearch,selected class : "<<currentClass<<" , "
                    					<<"image number : "<<imageNumberPerClass[currentClass-1]<<" , "
                    					<<"fine class number : "<<fineClassNumberPerClass[currentClass-1]<<" , "
                    					<<"classIndexAccumulator : "<<classIndexAccumulator<<" , "
                    					<<"imageIndexAccumulator : "<<imageIndexAccumulator<<std::endl<<std::flush;
                    // append metadata
                    MetaDataTable metaDataTable;
                    metaDataTable.readFromStar(star_fn);
                    metaDataTable.fliterByClass(currentClass);
                    for (int iimage = 0; iimage < metaDataTable.numberOfParticles(); iimage++)
                        metaDataTable[iimage].CLASS = 1;
                    append(metaDataTable);
                    // append classAverage
                    float buffer[image_size*image_size];
                    std::string mrcs_fn = pathRemoveSuffix(star_fn)+".mrcs";
                    FILE* mrcsFile = fopen(mrcs_fn.c_str(),"rb");
                    ERROR_CHECK(NULL == mrcsFile,"previous classaverage file should be put in same folder and has same name as star file.");
                    long offset = (256+(currentClass-1)*image_size*image_size)*sizeof(float);
                    fseek(mrcsFile,offset,SEEK_SET);
                    if(fread((char*)buffer,image_size*image_size*sizeof(float),1,mrcsFile) == NULL)
                        ERROR_REPORT("read mrcs data failed.");
                    fclose(mrcsFile);
                    append(buffer);
                    return false;
                }
                else{
                    MASTERNODE std::cout<<"starting doing fineSearch,selected class : "<<currentClass<<" , "
                            		<<"image number : "<<imageNumberPerClass[currentClass-1]<<" , "
                            		<<"fine class number : "<<fineClassNumberPerClass[currentClass-1]<<" , "
                            		<<"classIndexAccumulator : "<<classIndexAccumulator<<" , "
                            		<<"imageIndexAccumulator : "<<imageIndexAccumulator<<std::endl<<std::flush;
                    return true;
                }
            }
            
            ERROR_REPORT("Cannot call fineSearch without init.");
            return false;
        }
        
        int getCurrentClass(){return currentClass;}
        int getCurrentClassFineClassNumber(){return fineClassNumberPerClass[currentClass-1];}
        int getSelectedClassNumber(){return selectedClass.size();}
        
        void append(float* class_data)
        {
            if (doFineSearch)
            {	// directly append the fine classaverage data
                for (int iclass = 0; iclass < fineClassNumberPerClass[currentClass-1]; iclass++)
                {
                    for (int i = 0; i < image_size*image_size; i++) {
                        classAverage[(classIndexAccumulator+iclass)*image_size*image_size+i] = class_data[iclass*image_size*image_size+i];
                    }
                }
            }
            else
            	ERROR_REPORT("Cannot call fineSearch without init.");
        }
        
        void append(MetaDataTable& fineClassMetaData)
        {
            if (doFineSearch)
            {
                for (int iimage = 0; iimage < imageNumberPerClass[currentClass-1]; iimage++)
                {
                    fineClassMetaData[iimage].CLASS += classIndexAccumulator;
                    metaDataElems[imageIndexAccumulator+iimage] = fineClassMetaData[iimage];
                    fineClassMetaData[iimage].CLASS -= classIndexAccumulator;
                }
            }
            else
            	ERROR_REPORT("Cannot call fineSearch without init.");
        }
        
        void writeResult(std::string fn)
        {
            if (doFineSearch)
            {
                MASTERNODE
                {
                    int nextImageIndexAccumulator = imageIndexAccumulator + imageNumberPerClass[currentClass-1];
                    int nextClassIndexAccumulator = classIndexAccumulator + fineClassNumberPerClass[currentClass-1];
                    std::string fn_out = fn+"_"+num2str(nextClassIndexAccumulator,5);
                    
                    std::vector<MetaDataElem> metaDataElemsTmp;
                    metaDataElemsTmp.resize(nextImageIndexAccumulator);
                    for (int i = 0; i < nextImageIndexAccumulator; i++){
                        metaDataElemsTmp[i] = metaDataElems[i];
                    }
                    
                    MetaDataTable metaDataTable;
                    metaDataTable.readFromMetaDataElements(metaDataElemsTmp);
                    
                    metaDataTable.writeToStar(fn_out);
                    Mrcs::MrcsImages listOfImages(classAverage,image_size,nextClassIndexAccumulator);
                    Mrcs::writeMrcsData(fn_out,listOfImages);
                    std::cout<<"fineSearch finish,write result,selected class : "<<currentClass<<" , "
                            <<"image number : "<<imageNumberPerClass[currentClass-1]<<" , "
                    		<<"fine class number : "<<fineClassNumberPerClass[currentClass-1]<<" , "
                    		<<"next classIndexAccumulator : "<<nextClassIndexAccumulator<<" , "
                    		<<"next imageIndexAccumulator : "<<nextImageIndexAccumulator<<std::endl<<std::flush;
                }
            }
            else
                ERROR_REPORT("Cannot call fineSearch without init.");
        }
            
    };
    
    extern FineSearch fineSearch;
    
    // setup,prepare and destroy the data for GTM Optimizer -------------------------------
    
    void setupGTMoptimizer(double X_infos[],double Mu_infos[],int set_sampling_dim);
    
    double prepare();
    
    void destroyGTMoptimizer();
    
    // Functions for main GTM algorithm  -----------------------------------
    
	void updateR();

	void iterateUpdateW();
    
    MIC_OFFLOAD_ATTRIBUTES
	void getWd(const int &d,double *A,double *b,int maxthreads,double *AA,double *bRealArray,double *bImagArray);
    
	void updateY();
    
	double updateBeta();
    
	double updateAlpha();

	double update_P_Delta();

	// run GTM
    void run(double preci = 10e-10,double probThreshold = 1.0,double alpha = 0.01,bool update_beta = true);
    
	//   write result -----------------------------------
    void writeClassAverageCTF(std::string filename_mrcs,std::string fn_metadata,double probThreshold);
    
    void writeClassAverageCTFsub(double probThreshold = 1.0);
    
    void writeClassAverageCTFsub_maxPro();
    
    // Intel(R) MKL function to solve A*x=b linear equation
    MIC_OFFLOAD_ATTRIBUTES
	inline void mkl_solve(double *A,int N,double *b,int M);

    // Shifts Matrix by column mean
    void shiftMat(double *Mat,int M,int N);
    
	void analysisR(double threshold = 0.98);

	// convergence detection functions -----------------------------------
    
	// If iterated 3 times and each change smaller than precision,then convergence
	bool checkConvergeBeta(double precision = 10e-20);

	bool checkConvergeAlpha(double precision = 10e-20);

	bool checkConvergeDeviation(double precision = 10e-20);

	// some debug functions -----------------------------------
    
	void dumpT();
    
	void dumpR();
    
	void dumpW();
    
	void dumpY();
    
	void dumpBelta();
    
    void TestOffloadMemory(int mic_index);
    
    double calSpace();
    
    double printDataSize();
};


#endif