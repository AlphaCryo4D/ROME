/***************************************************************************
 *
 * Intel® Parallel Computing Center for Structural Biology
 * Principal Investigator : Youdong (Jack) Mao (Youdong_Mao@dfci.harvard.edu)
 * Dana-Farber Cancer Institute, Harvard Medical School and Peking University
 *
 * Authors: "Yong Bei Ma(galowma@gmail.com) Jian Wang(wj_hust08@hust.edu.cn)"
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

#ifndef MRCS_H_
#define MRCS_H_

#include "resmap_util.h"

#include "resmap_error.h"
#include "resmap_string.h"
#include "resmap_macros.h"
#include "resmap_platform.h"
#include "resmap_mpi.h"

class FloatImages {
public:
    // TBD : why this virtual deconstructor cannot work on my Mac
    // virtual ~ListOfImages() {}
    virtual int nr_images() = 0;
    virtual int imageSide() = 0;		// images are assumed to be square
    virtual float* image_ptr(size_t i) = 0;
};

class FloatVolume{
public:
    virtual int volumeSide() = 0;
    virtual float* ptr() = 0;
};

namespace Mrcs {
    // CCP4 formation like RELION,from : http://www.ccp4.ac.uk/html/maplib.html
    // only NC,NR,NS,MODE is used for *.mrcs
    // NSYMBT will be used in *.mrc
    struct MrcsHead{
        // TODO : make sure this mrcsHead file can be readed by RELION
#define MRCSHEAD_DATA \
        ELT(int,    NC,       0,  1, "of Columns    (fastest changing in map)"                             ) SEP \
        ELT(int,    NR,       0,  2, "of Rows"                                                             ) SEP \
        ELT(int,    NS,       0,  3, "of Sections   (slowest changing in map)"                             ) SEP \
        ELT(int,    MODE,     2,  4, "Data type,0,1,2,3,4,5 we use 2(float)"                               ) SEP \
        ELT(int,    NCSTART,  0,  5, "Number of first COLUMN  in map"                                      ) SEP \
        ELT(int,    NRSTART,  0,  6, "Number of first ROW     in map"                                      ) SEP \
        ELT(int,    NSSTART,  0,  7, "Number of first SECTION in map"                                      ) SEP \
        ELT(int,    NX,       0,  8, "Number of intervals along X"                                         ) SEP \
        ELT(int,    NY,       0,  9, "Number of intervals along Y"                                         ) SEP \
        ELT(int,    NZ,       0,  10, "Number of intervals along Z"                                        ) SEP \
        ELT(float,  X_length, 0,  11, "Cell Dimensions (Angstroms)"                                        ) SEP \
        ELT(float,  Y_length, 0,  12, "\""                                                                 ) SEP \
        ELT(float,  Z_length, 0,  13, "\""                                                                 ) SEP \
        ELT(float,  Alpha,    0,  14, "Cell Angles     (Degrees)"                                          ) SEP \
        ELT(float,  Beta,     0,  15, "\""                                                                 ) SEP \
        ELT(float,  Gamma,    0,  16, "\""                                                                 ) SEP \
        ELT(int,    MAPC,     0,  17, "Which axis corresponds to Cols.  (1,2,3 for X,Y,Z)"                 ) SEP \
        ELT(int,    MAPR,     0,  18, "Which axis corresponds to Rows   (1,2,3 for X,Y,Z)"                 ) SEP \
        ELT(int,    MAPS,     0,  19, "Which axis corresponds to Sects. (1,2,3 for X,Y,Z)"                 ) SEP \
        ELT(float,  AMIN,     0,  20, "Minimum density value"                                              ) SEP \
        ELT(float,  AMAX,     0,  21, "Maximum density value"                                              ) SEP \
        ELT(float,  AMEAN,    0,  22, "Mean    density value    (Average)"                                 ) SEP \
        ELT(int,    ISPG,     0,  23, "Space group number"                                                 ) SEP \
        ELT(int,    NSYMBT,   0,  24, "Number of bytes used for storing symmetry operators"                ) SEP \
        ELT(int,    LSKFLG,   0,  25, "Flag for skew transformation, =0 none, =1 if foll"                  ) SEP \
        ELT(int,    SKWMAT11, 0,  26, "Skew matrix S (in order S11, S12, S13, S21 etc) if LSKFLG .ne. 0."  ) SEP \
        ELT(int,    SKWMAT12, 0,  27, "Skew matrix S (in order S11, S12, S13, S21 etc) if LSKFLG .ne. 0."  ) SEP \
        ELT(int,    SKWMAT13, 0,  28, "Skew matrix S (in order S11, S12, S13, S21 etc) if LSKFLG .ne. 0."  ) SEP \
        ELT(int,    SKWMAT21, 0,  29, "Skew matrix S (in order S11, S12, S13, S21 etc) if LSKFLG .ne. 0."  ) SEP \
        ELT(int,    SKWMAT22, 0,  30, "Skew matrix S (in order S11, S12, S13, S21 etc) if LSKFLG .ne. 0."  ) SEP \
        ELT(int,    SKWMAT23, 0,  31, "Skew matrix S (in order S11, S12, S13, S21 etc) if LSKFLG .ne. 0."  ) SEP \
        ELT(int,    SKWMAT31, 0,  32, "Skew matrix S (in order S11, S12, S13, S21 etc) if LSKFLG .ne. 0."  ) SEP \
        ELT(int,    SKWMAT32, 0,  33, "Skew matrix S (in order S11, S12, S13, S21 etc) if LSKFLG .ne. 0."  ) SEP \
        ELT(int,    SKWMAT33, 0,  34, "Skew matrix S (in order S11, S12, S13, S21 etc) if LSKFLG .ne. 0."  ) SEP \
        ELT(int,    SKWTRN_X, 0,  35, "Skew translation t if LSKFLG .ne. 0.Skew transformation is from"    ) SEP \
        ELT(int,    SKWTRN_Y, 0,  36, "standard orthogonal coordinate frame (as used for atoms) to"        ) SEP \
        ELT(int,    SKWTRN_Z, 0,  37, "orthogonal map frame, as Xo(map) = S * (Xo(atoms) - t)"             ) SEP \
        ELT(int,    NONE1,    0,  38, "(some of these are used by the MSUBSX routines"                     ) SEP \
        ELT(int,    NONE2,    0,  39, "in MAPBRICK, MAPCONT and FRODO)"                                    ) SEP \
        ELT(int,    NONE3,    0,  40, " (all set to zero by default)"                                      ) SEP \
        ELT(int,    NONE4,    0,  41, " (all set to zero by default)"                                      ) SEP \
        ELT(int,    NONE5,    0,  42, " (all set to zero by default)"                                      ) SEP \
        ELT(int,    NONE6,    0,  43, " (all set to zero by default)"                                      ) SEP \
        ELT(int,    NONE7,    0,  44, " (all set to zero by default)"                                      ) SEP \
        ELT(int,    NONE8,    0,  45, " (all set to zero by default)"                                      ) SEP \
        ELT(int,    NONE9,    0,  46, " (all set to zero by default)"                                      ) SEP \
        ELT(int,    NONE10,   0,  47, " (all set to zero by default)"                                      ) SEP \
        ELT(int,    NONE11,   0,  48, " (all set to zero by default)"                                      ) SEP \
        ELT(int,    NONE12,   0,  49, " (all set to zero by default)"                                      ) SEP \
        ELT(int,    NONE13,   0,  50, " (all set to zero by default)"                                      ) SEP \
        ELT(int,    NONE14,   0,  51, " (all set to zero by default)"                                      ) SEP \
        ELT(int,    NONE15,   0,  52, " (all set to zero by default)"                                      ) SEP \
        ELT(int,    MAP,      0,  53, "Character string 'MAP ' to identify file type"                      ) SEP \
        ELT(int,    MACHST,   0,  54, "Machine stamp indicating the machine type which wrote file"         ) SEP \
        ELT(int,    ARMS,     0,  55, "Rms deviation of map from mean density"                             ) SEP \
        ELT(int,    NLABL,    0,  56, "Number of labels being used"                                        ) // end of macro
        
#define ELT(T,N,V,I,C) T N = V;
#define SEP
        MRCSHEAD_DATA
#undef ELT
#undef SEP
        // 10  80 character text labels (ie. A4 format)
        int LABEL[200] = {0};
        
        void print(std::ostream& os = std::cout) const{
            os<<"MrcsHead{"
#define ELT(T,N,V,I,C) <<" "<< #N <<" : "<<N
#define SEP
            MRCSHEAD_DATA
#undef ELT
#undef SEP
            <<" }"<<std::endl;
        }
        void operator=(const MrcsHead& rhs){
#define ELT(T,N,V,I,C) N = rhs.N;
#define SEP
            MRCSHEAD_DATA
#undef ELT
#undef SEP
            for(int i = 0;i < 200;i++) LABEL[i] = rhs.LABEL[i];
        }
        
    };
    
    // typedef FloatImages MrcsImages;
    
    // read Mrcshead from *.mrcs file
    void readMrcsHead(std::string fn_mrcs,MrcsHead& mrcsHead);
    
    // read MrcHead from *.mrc
    bool readMrcHead(std::string fn_mrc,MrcsHead& refHead,int size);
    
    // set mrc Head
    void setMrcHead(MrcsHead& refHead,double anpix,int size);
    
    // read Mrcs image Data
    void readMrcsData(std::string fn_mrcs,const MrcsHead& mrcsHead,FloatImages &mrcsImages);
    
    // write imagedata to mrcs file
    void writeMrcsData(std::string fn_mrcs,FloatImages &mrcsImages);
    
    // read 3D volume from mrc file
    void readMrcData(std::string fn_mrc,FloatVolume &mrcVolume);
    
    // write 3D volume to mrc file
    void writeMrcData(std::string fn_mrcs, float *volume_data, int volumeSize, MrcsHead &mrcHead);
    void writeMrcData(std::string fn_mrc, FloatVolume &mrcVolume, MrcsHead mrcHead, double anpix);
    
    // 2D images stack
    class MrcsImages : public FloatImages{
        float* image_data;
        int N,size;
        bool is_alloc;
    public:
        MrcsImages(int _size, int _N) : size(_size), N(_N) {
            is_alloc = true;
            image_data = (float*)aMalloc(sizeof(float)*N*size*size,64);
            for (int i = 0; i < N*size*size; i++) image_data[i] = 0.;
        }
        MrcsImages(float* data, int _size, int _N) : size(_size), N(_N) {
            is_alloc = false;
            image_data = data;
        }
        MrcsImages(const double* data, int _size, int _N) : size(_size), N(_N) {
            is_alloc = true;
            image_data = (float*)aMalloc(sizeof(float)*N*size*size,64);
            for (int i = 0; i < N*size*size; i++) image_data[i] = float(data[i]);
        }
        virtual ~MrcsImages() {if(is_alloc) aFree(image_data); image_data = nullptr;}
        virtual int nr_images() { return N; }
        virtual int imageSide() { return size; }
        virtual float* image_ptr(size_t i) { return image_data + size*size*i; }
        void write(std::string filename){
            writeMrcsData(filename,(*this));
        }
    };
    
    // 2D Fourier-image stack
    class MrcsFImages : public FloatImages{
        float* image_data;
        float* image_data_real;
        float* image_data_imag;
        int N,size,Fsize;
        bool is_alloc;
    public:
        MrcsFImages(int _size, int _N) : size(_size), N(_N), Fsize(_size/2+1) {
            is_alloc = true;
            image_data = (float*)aMalloc(sizeof(float)*N*size*size,64);
            for (int i = 0; i < N*size*size; i++) image_data[i] = 0.;
            image_data_real = (float*)aMalloc(sizeof(float)*N*size*Fsize,64);
            image_data_imag = (float*)aMalloc(sizeof(float)*N*size*Fsize,64);
            for (int i = 0; i < N*size*Fsize; i++) image_data_real[i] = image_data_imag[i] = 0.;
        }
        template<typename T>
        MrcsFImages(T* data_real, T* data_imag, int _size, int _N) : size(_size), N(_N), Fsize(_size/2+1) {
            is_alloc = false;
            int size2 = size*size;
            int Fsize2 = size*Fsize;
            image_data = (float*)aMalloc(sizeof(float)*N*size*size,64);
            for (int n = 0; n < N; n++) {
                //
                for (int i = 0; i < size; i++) {
                    for (int j = 0; j < Fsize; j++) {
                        image_data[n*size2+i*size+j] = sqrt(data_real[n*Fsize2+i*Fsize+j]*data_real[n*Fsize2+i*Fsize+j] + \
                                                            data_imag[n*Fsize2+i*Fsize+j]*data_imag[n*Fsize2+i*Fsize+j]);
                        image_data[n*size2+i*size+size-j] = image_data[n*Fsize2+i*size+j];
                    }
                }
            }
        }
        virtual ~MrcsFImages() {
            if(is_alloc) {
                aFree(image_data_real);
                aFree(image_data_imag);
            }
            aFree(image_data);
            image_data_real = nullptr;
            image_data_imag = nullptr;
            image_data = nullptr;
        }
        virtual int nr_images() { return N; }
        virtual int imageSide() { return size; }
        virtual float* image_ptr(size_t i) { return image_data + size*size*i; }
        void write(std::string filename){
            writeMrcsData(filename,(*this));
        }
    };
    
    // 3D Image Volume
    class MrcVolume : FloatVolume {
    public:
        using value_type = float;

        value_type* data = nullptr;
        int size;
        value_type anpix;
        bool is_alloc = false;
        std::shared_ptr<MrcsHead> head;

        MrcVolume() { init(); }

        MrcVolume(std::string fn_mrc) { init(); read(fn_mrc); }

        MrcVolume(int _size,value_type _anpix) : size(_size), anpix(_anpix) {
            is_alloc = true;
            data = (value_type*)aMalloc(sizeof(value_type)*size*size*size,64);
            for (int i = 0; i < size*size*size; i++) data[i] = 0.;
        }

        MrcVolume(float* _data, int _size, float _anpix) : size(_size), anpix(_anpix) {
            is_alloc = false;
            data = _data;
        }

        MrcVolume(const double* _data, int _size, double _anpix) : size(_size), anpix(_anpix) {
            is_alloc = true;
            data = (value_type*)aMalloc(sizeof(value_type)*size*size*size,64);
            for (int i = 0; i < size*size*size; i++) data[i] = value_type(_data[i]);
        }

        ~MrcVolume() {
            clear();
        }

        virtual int volumeSide() { return size; }

        virtual value_type* ptr() {return data;}

        void write(std::string filename) {
            writeMrcData(filename, (*this), *head, anpix);
        }

        void write(std::string filename, MrcsHead& mrcHead){
            writeMrcData(filename, (*this), mrcHead, anpix);
        }

        void read(std::string fn_mrc) {
            clear();
            bool ieee_le = true;
            head = std::make_shared<MrcsHead>();
            std::ifstream ifile(fn_mrc.c_str(), std::ios::in|std::ios::binary);
            ERROR_CHECK(!ifile.is_open(), strMerge("cannot open mrc file: ", fn_mrc));
            int *p = &(head->NC);
            ifile.read(reinterpret_cast<char *>(p), sizeof(char)*1024);
            if (head->MODE > 10) ieee_le = false;
            if (!ieee_le) {
                for (int i = 0; i < sizeof(char)*1024/sizeof(int); i++) {
                    p[i] = SWAP32(p[i]);
                }
            }
            size = head->NC;
            MPI_LOG << fn_mrc << std::endl;
            if(MPI_IS_ROOT) head->print();
            anpix = head->X_length / (value_type)size;
            is_alloc = true;
            data = (value_type*)aMalloc(sizeof(value_type)*size*size*size, 64);
            ifile.seekg(head->NSYMBT, std::ios::cur);
            ifile.read(reinterpret_cast<char *>(data), sizeof(value_type)*size*size*size);
            p = (int*)data;
            if (!ieee_le) {
                for (int i = 0; i < sizeof(value_type)*size*size*size/sizeof(int); i++) {
                    p[i] = SWAP32(p[i]);
                }
            }
            ifile.close();
        }

        void init() {
            size = 0;
            anpix = 0;
            data = nullptr;
            is_alloc = false;
        }

        void clear() {
            if (is_alloc) {
                aFree(data);
            }
            init();
        }

        //void read(std::string filename){
        //    readMrcData(filename,(*this));
        //}
    };
    
    void UnitTest(std::string fn_mrcs);
    
}
#endif
