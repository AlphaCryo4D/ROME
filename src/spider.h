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

#ifndef SPIDER_H_
#define SPIDER_H_

#include "./util.h"

#include "fftw/fftw3.h"

#include "mrcs.h"
#include "error.h"
#include "spider_dep.h"

namespace Spider{
    // spider *.dat header,from : http://spider.wadsworth.org/spider_doc/spider/docs/image_doc.html
    struct DatHead
    {
#define DATHEAD_DATA \
        ELTONE(float,  nz,       0, 1,"Number of slices (planes) in volume (=1 for an image) In some ancient 2D images this may be -1)"         ) SEP \
        ELTONE(float,  ny,       0, 2,"Number of rows per slice"                                                                                ) SEP \
        ELTONE(float,  irec,     0, 3,"Total number of records (including header records) in each image of a simple image or stacked image file") SEP \
        ELTONE(float,  unused,   0, 4,"Unassigned"                                                                                              ) SEP \
        ELTONE(float,  iform,    1, 5,"File type specifier,(1.R.@D image),(3.R.3D volume),(-11.FO.2D Fourier,odd),(-12.FE.2D Fourier,even),(-21.FO.3D Fourier,odd),(-22.FE.3D Fourier,even)"         ) SEP \
        ELTONE(float,  imami,    0, 6,"Maximum/minimum flag = 0 when the file is created, and = 1 when the maximum, minimum, average, and standard deviation have been computed and stored into this header record (see following locations)"         ) SEP \
        ELTONE(float,  fmax,     0, 7,"Maximum data value"                                                                                      ) SEP \
        ELTONE(float,  fmin,     0, 8,"Minimum data value"                                                                                      ) SEP \
        ELTONE(float,  av,       0, 9,"Average data value"                                                                                      ) SEP \
        ELTONE(float,  sig,      0, 10,"Standard deviation of data. A value of -1.0 or 0.0 indicates that SIG has not been computed"            ) SEP \
        ELTONE(float,  unused2,  0, 11,"No longer used"                                                                                         ) SEP \
        ELTONE(float,  nx,       0, 12,"Number of pixels (samples) per line"                                                                    ) SEP \
        ELTONE(float,  labrec,   0, 13,"Number of records in file header (label)"                                                               ) SEP \
        ELTONE(float,  iangle,   0, 14,"Flag that following three tilt angles are present"                                                      ) SEP \
        ELTONE(float,  phi,      0, 15,"Tilt angle: phi (See note #2 below)"                                                                    ) SEP \
        ELTONE(float,  theta,    0, 16,"Tilt angle: theta"                                                                                      ) SEP \
        ELTONE(float,  gamma,    0, 17,"Tilt angle: gamma (also called psi)"                                                                    ) SEP \
        ELTONE(float,  xoff,     0, 18,"X translation"                                                                                          ) SEP \
        ELTONE(float,  yoff,     0, 19,"Y translation"                                                                                          ) SEP \
        ELTONE(float,  zoff,     0, 20,"Z translation"                                                                                          ) SEP \
        ELTONE(float,  scale,    0, 21,"scale"                                                                                                  ) SEP \
        ELTONE(float,  labbyt,   0, 22,"Total number of bytes in header"                                                                        ) SEP \
        ELTONE(float,  lenbyt,   0, 23,"Record length in bytes"                                                                                 ) SEP \
        ELTONE(float,  istack,   0, 24,"Position has a value of 0 in simple 2D or 3D (non-stack) files. In an "image stack" there is one overall stack header followed by a stack of images, in which each image has its own image header. A value of >0 in this position in the overall stack header indicates a stack of images. A value of <0 in this position in the overall stack header indicates an indexed stack of images and gives the maximum image number (MAXINDX) allowed in the index"         ) SEP \
        ELTONE(float,  unused3,  0, 25,"Unused now. Prior to release 9.0, a '-1' at this location in an overall stack indicated a valid stack and in the stacked images,a value of 1 indicated that this image was in use (existed)"         ) SEP \
        ELTONE(float,  maxim,    0, 26,"Position is only used in the overall header for a stacked image file. There, this position contains the number of the highest image currently used in the stack. This number is updated, if necessary, when an image is added or deleted from the stack"         ) SEP \
        ELTONE(float,  imgnum,   0, 27,"Position is only used in a stacked image header. There, this position contains the number of the current image or zero if this image is unused"         ) SEP \
        ELTONE(float,  lastindx, 0, 28,"Position is only used in overall header of indexed stacks. There, this position is the highest index location currently in use"         ) SEP \
        ELTONE(float,  unused4,  0, 29,"Unassigned"                                                                                             ) SEP \
        ELTONE(float,  unused5,  0, 30,"Unassigned"                                                                                             ) SEP \
        ELTONE(float,  kangle,   0, 31,"Flag that additional rotation angles follow in header. 1 = one additional angle set is present, 2 = two additional angle sets"         ) SEP \
        ELTONE(float,  phi1,     0, 32,"Angle"                                                                                                  ) SEP \
        ELTONE(float,  theta1,   0, 33,"Angle"                                                                                                  ) SEP \
        ELTONE(float,  psi1,     0, 34,"Angle"                                                                                                  ) SEP \
        ELTONE(float,  phi2,     0, 35,"Angle"                                                                                                  ) SEP \
        ELTONE(float,  theta2,   0, 36,"Angle"                                                                                                  ) SEP \
        ELTONE(float,  psi2,     0, 37,"Angle"                                                                                                  ) SEP \
        ELTONE(float,  pixsiz,   0, 38,"Pixel size (Angstroms)"                                                                                 ) SEP \
        ELTONE(float,  ev,       0, 39,"Electron voltage used"                                                                                  ) SEP \
        ELTONE(float,  proj,     0, 40,"Project number"                                                                                         ) SEP \
        ELTONE(float,  mic,      0, 41,"Micrograph number"                                                                                      ) SEP \
        ELTONE(float,  num,      0, 42,"Micrograph number"                                                                                      ) SEP \
        ELTONE(float,  glonum,   0, 43,"Global image number"                                                                                    ) SEP \
        ELTVEC(float,  unused6,  0, 44,4,"Unassigned"                                                                                           ) SEP \
        ELTVEC(float,  forxmipp, 0, 48,29,"Reserved for XMIPP or other local transforms"                                                        ) SEP \
        ELTVEC(float,  unused7,  0, 77,24,"Unassigned"                                                                                          ) SEP \
        ELTONE(float,  psi3,     0, 101,"Projection angle: Psi (From 'PJ 3Q')"                                                                  ) SEP \
        ELTONE(float,  theta3,   0, 102,"Projection angle: Theta (From 'PJ 3Q')"                                                                ) SEP \
        ELTONE(float,  phi3,     0, 103,"Projection angle: Phi (From 'PJ 3Q')"                                                                  ) SEP \
        ELTONE(float,  langle,   0, 104,"If = 1 then projection angles: PSI3, THETA3 & PHI3 are present in header"                              ) SEP \
        ELTVEC(float,  unused8,  0, 105,107,"Unassigned"                                                                                        ) SEP \
        ELTVEC(char,   cdat,     0, 212,12,"Character *11 Creation date e.g. 27-MAY-1999"                                                       ) SEP \
        ELTVEC(char,   ctim,     0, 215,8,"Character *8 Creation time e.g. 09:43:19"                                                            ) SEP \
        ELTVEC(char,   ctit,     0, 217,160,"Character *160 Title"                                                                              ) // end of macro
        
#define ELTONE(T,N,V,I,C) T N = V;
#define ELTVEC(T,N,V,I,L,C) T N[L] = {V};
#define SEP
        DATHEAD_DATA
#undef ELTONE
#undef ELTVEC
#undef SEP
        
        void print(std::ostream& os = std::cout) const{
            os<<"DatHead{"
#define ELTONE(T,N,V,I,C) <<" "<< #N <<" : "<<N
#define ELTVEC(T,N,V,I,L,C) <<" "<< #N <<" : "<<N
#define SEP
            DATHEAD_DATA
#undef ELTONE
#undef ELTVEC
#undef SEP
            <<std::endl;
        }
        
        void operator=(const DatHead& rhs){
#define ELTONE(T,N,V,I,C) N = rhs.N;
#define ELTVEC(T,N,V,I,L,C) for(int i=0;i<L;i++) N[i] = rhs.N[i];
#define SEP
            DATHEAD_DATA
#undef ELTONE
#undef ELTVEC
#undef SEP
        }
        
    };
    
    typedef Mrcs::MrcsImages DatImages;
    
    // reading and writing *.dat,the ordering is big-endian(ieee-bg) or
    // little-endian(ieee-le),the spider using ieee-be by default
    void readDatHead(std::string fn_dat,DatHead& dathead,std::string ordering = "ieee-be");

    void readDatData(std::string fn_dat,DatHead& dathead,FloatImages& datImages,std::string ordering = "ieee-be");
    
    void writeDatData(std::string filename,float* data,int nz,int ny,int nx,bool is_stack = false,std::string ordering = "ieee-be");
    void writeDatData(std::string filename,FloatImages& datImages, std::string ordering = "ieee-be");
    
    // Image operations
    
    template<typename T>
    inline T modulo(T A,T P){
        return (double)((double)A - floor((double)A/P) * P);
    }
    
    // apply filter to the images
    void applyFilter(FloatImages& datImages,double filter);
    
    // anticlockwise rotation(positive),shift from left(negative) to righ(positive)t,top(negative) to bottom(positive)
    // void applyRotationAndShift(ListOfImages& listOfImages,float angle,float shiftx,float shifty);
    // for advanced user use,rotate and shift image one by one
    // this is only used in GTM setup step to rotate and shift the image to right position
    void RT_SF(const float* image_data,float *image_out_data,int nx,int ny,float angle,float shiftx,float shifty);
    
    // for *.dat read and write
    void UnitTest(std::string fn_mrcs);
    
    // for image rotate and shift
    void UnitTest2(std::string fn_mrcs);
}

#endif
