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

#ifndef METADATA_H_
#define METADATA_H_

#include "./util.h"
#include "./memory.h"
#include "./error.h"
#include "./string.h"
#include "./macros.h"
#include "./mpi.h"

// metadata name format number@filename ,
// like : 3@particles.mrcs
struct MetaName{
public:
    size_t INDEX;
    std::string NAME;
    MetaName():INDEX(0),NAME(""){}
    MetaName(const char* full_name){set(full_name);}
    MetaName(const std::string& full_name){set(full_name);}
    ~MetaName(){}
    bool operator!=(const MetaName& rhs) const {return !(rhs.INDEX==INDEX && rhs.NAME==NAME);}
    bool operator==(const MetaName& rhs) const {return  (rhs.INDEX==INDEX && rhs.NAME==NAME);}
    bool operator>(const MetaName& rhs) const {return  NAME > rhs.NAME ? true : INDEX > rhs.INDEX;}
    bool operator<(const MetaName& rhs) const {return  NAME < rhs.NAME ? true : INDEX < rhs.INDEX;}
private:
    void set(const std::string& full_name){
        INDEX = atoi(full_name.substr(0,full_name.find_last_of('@')).c_str());
        NAME = full_name.substr(full_name.find_last_of('@')+1);
    }
};
inline std::ostream& operator<<(std::ostream& os, MetaName const & rhs) {os<<num2str(rhs.INDEX,7)<<'@'<<rhs.NAME; return os; }

//  metadata elements
//          datatype   , var                , default   ,  enum (star file's header)
#define METADATA_ELTS \
    ELTDOU(FDOUBLE     , CTF_VOLTAGE		, 	0.0     ,   Voltage                     ) SEP \
    ELTDOU(FDOUBLE     , CTF_DEFOCUS_U		,   0.0     ,   DefocusU                    ) SEP \
    ELTDOU(FDOUBLE     , CTF_DEFOCUS_V		,   0.0     ,   DefocusV                    ) SEP \
    ELTDOU(FDOUBLE     , CTF_DEFOCUS_ANGLE	,   0.0     ,   DefocusAngle                ) SEP \
    ELTDOU(FDOUBLE     , CTF_CS             ,   0.0     ,   SphericalAberration         ) SEP \
    ELTDOU(FDOUBLE     , CTF_Q0             ,   0.0     ,   AmplitudeContrast           ) SEP \
    ELTDOU(FDOUBLE     , CTF_BFAC			,   0.0     ,   Bfactor                     ) SEP \
    \
    ELTSTR(MetaName    , IMAGE              , ""        ,   ImageName                   ) SEP \
    ELTSTR(std::string , MICROGRAPH         , ""        ,   MicrographName              ) SEP \
    /* for each group a different noise spectrum and signal scale factor is estimated */      \
    /* similar defocus values often have similar scale intensity factors and noise */         \
    /* power spectra,so if some groups image number is less,you can group some Micrograph */  \
    /* only the group_name for reading */                                                     \
    ELTINT(int         , GROUP_NO			, 1         ,   GroupNumber                 ) SEP \
    ELTSTR(std::string , GROUP_NAME			, ""        ,   GroupName                   ) SEP \
	/* Random subset this particle belongs to(only used in 'autorefine') */					  \
	ELTINT(int		   , SUBSET				, 0			,	RandomSubset				) SEP \
    ELTDOU(FDOUBLE     , ROT				, 0.0       ,   AngleRot                    ) SEP \
    ELTDOU(FDOUBLE     , TILT				, 0.0       ,   AngleTilt                   ) SEP \
    ELTDOU(FDOUBLE     , PSI				, 0.0       ,   AnglePsi                    ) SEP \
    ELTDOU(FDOUBLE     , XOFF				, 0.0       ,   OriginX                     ) SEP \
    ELTDOU(FDOUBLE     , YOFF				, 0.0       ,   OriginY                     ) SEP \
    ELTDOU(FDOUBLE     , ZOFF				, 0.0       ,   OriginZ                     ) SEP \
    ELTINT(int         , CLASS				, 0         ,   ClassNumber                 ) SEP \
    ELTDOU(FDOUBLE     , NORM				, 1.0       ,   NormCorrection              ) SEP \
    ELTDOU(FDOUBLE     , SCALE				, 1.0       ,   ScaleCorrection             ) SEP \
    ELTDOU(FDOUBLE     , DLL				, 0.0       ,   LogLikeliContribution       ) SEP \
    ELTDOU(FDOUBLE     , PMAX				, 0.0       ,   MaxValueProbDistribution    ) SEP \
    ELTINT(int         , NR_SIGN			, 0         ,   NrOfSignificantSamples      ) SEP \
    \
    ELTDOU(FDOUBLE     , ROT_PRIOR			, 999.0     ,   AngleRotPrior               ) SEP \
    ELTDOU(FDOUBLE     , TILT_PRIOR         , 999.0     ,   AngleTiltPrior              ) SEP \
    ELTDOU(FDOUBLE     , PSI_PRIOR			, 999.0     ,   AnglePsiPrior               ) SEP \
    ELTDOU(FDOUBLE     , XOFF_PRIOR         , 999.0     ,   OriginXPrior                ) SEP \
    ELTDOU(FDOUBLE     , YOFF_PRIOR         , 999.0     ,   OriginYPrior                ) SEP \
    ELTDOU(FDOUBLE     , ZOFF_PRIOR         , 999.0     ,   OriginZPrior                ) SEP \
    \
    ELTDOU(FDOUBLE     , MAGNIFICATION      ,   0.0     ,   Magnification               ) SEP \
    ELTDOU(FDOUBLE     , DETECTOR_PIXEL_SIZE,   0.0     ,   DetectorPixelSize           ) SEP \
    ELTDOU(FDOUBLE     , BEAMTILT_X         ,   0.0     ,   BeamTiltX                   ) SEP \
    ELTDOU(FDOUBLE     , BEAMTILT_Y         ,   0.0     ,   BeamTiltY                   ) SEP \
    \
    ELTDOU(FDOUBLE     , MAT_0_0			,   1.0     ,   Matrix_1_1                  ) SEP \
    ELTDOU(FDOUBLE     , MAT_0_1			,   0.0     ,   Matrix_1_2                  ) SEP \
    ELTDOU(FDOUBLE     , MAT_0_2			,   0.0     ,   Matrix_1_3                  ) SEP \
    ELTDOU(FDOUBLE     , MAT_1_0			,   0.0     ,   Matrix_2_1                  ) SEP \
    ELTDOU(FDOUBLE     , MAT_1_1			,   1.0     ,   Matrix_2_2                  ) SEP \
    ELTDOU(FDOUBLE     , MAT_1_2			,   0.0     ,   Matrix_2_3                  ) SEP \
    ELTDOU(FDOUBLE     , MAT_2_0			,   0.0     ,   Matrix_3_1                  ) SEP \
    ELTDOU(FDOUBLE     , MAT_2_1			,   0.0     ,   Matrix_3_2                  ) SEP \
    ELTDOU(FDOUBLE     , MAT_2_2			,   1.0     ,   Matrix_3_3                  )


enum MetaDataLabel{
#define SEP ,
#define ELTINT(T,N,I,C) C
#define ELTDOU(T,N,I,C) C
#define ELTSTR(T,N,I,C) C
    METADATA_ELTS
#undef ELTINT
#undef ELTDOU
#undef ELTSTR
#undef SEP
    ,UnDefinedLabel
};

namespace LabelParser{
    MetaDataLabel stringToLabel(const std::string& name);
    std::string labelToString(MetaDataLabel label);
};

// one metadata element
struct MetaDataElem{
    size_t INNERID = -1;
#define SEP
#define ELTINT(T,N,I,C) T N = I;
#define ELTDOU(T,N,I,C) T N = I;
#define ELTSTR(T,N,I,C) T N = I;
    METADATA_ELTS
#undef ELTINT
#undef ELTDOU
#undef ELTSTR
#undef SEP
    
    void set(MetaDataLabel label,std::string val){
        switch(label)
        {
#define SEP
#define ELTINT(T,N,I,C) case C : N = atoi(val.c_str()); break;
#define ELTDOU(T,N,I,C) case C : N = atof(val.c_str()); break;
#define ELTSTR(T,N,I,C) case C : N = val            ; break;
            METADATA_ELTS
#undef ELTINT
#undef ELTDOU
#undef ELTSTR
#undef SEP
        }
    }
    
    void get(MetaDataLabel label,std::ostream& os) const {
        // os.precision(6);
        auto printDouble = [&](double v){
            if((fabs(v) > 0. && fabs(v) < 0.001) || fabs(v) > 100000.)
                os<<std::setw(12)<<std::scientific<<v<<"\t";
            else
                os<<std::setw(12)<<std::fixed<<v<<"\t";
        };
        auto printInt = [&](int v){
            os<<std::setw(6)<<v<<"\t";
        };
        
        switch(label)
        {
#define SEP
#define ELTINT(T,N,I,C) case C : printInt(N)        ; break;
#define ELTDOU(T,N,I,C) case C : printDouble(N)     ; break;
#define ELTSTR(T,N,I,C) case C : os<<N<<"\t";  ; break;
            METADATA_ELTS
#undef ELTINT
#undef ELTDOU
#undef ELTSTR
#undef SEP
        }
    }
    
    void operator=(MetaDataElem const & rhs) {
        INNERID = rhs.INNERID;
#define SEP
#define ELTINT(T,N,I,C) N = rhs.N;
#define ELTDOU(T,N,I,C) N = rhs.N;
#define ELTSTR(T,N,I,C) N = rhs.N;
        METADATA_ELTS
#undef ELTINT
#undef ELTDOU
#undef ELTSTR
#undef SEP
    }

    bool operator==(MetaDataElem const & rhs) const {
        const char* diff = NULL;
		//
		// NR_SIGN is a challenge, because it is an integer but the exact number (number of significant samples) is imprecise
		//
		auto compareInts = [&](const char* name, int lhs, int rhs) {
			if (lhs == rhs) return;
			if (strcmp(name, "NR_SIGN") == 0 && std::abs(lhs-rhs) < std::max(lhs,rhs)/100) return;
			diff = name;
			std::cerr << "MetaDataElem:" << name << " lhs:" << lhs << " != " << rhs<< std::endl;
		};
#define SEP
#define ELTINT(T,N,I,C) compareInts(#N, N, rhs.N);
#define ELTDOU(T,N,I,C) if (!nearEnoughTemplate(N, rhs.N)) { diff = #N;	std::cerr << "MetaDataElem:" << #N << " lhs:" << N << " != " << rhs.N << std::endl; }
#define ELTSTR(T,N,I,C) if (!nearEnoughTemplate(N, rhs.N)) { diff = #N;	std::cerr << "MetaDataElem:" << #N << " lhs:" << N << " != " << rhs.N << std::endl; }
    METADATA_ELTS
        if (INNERID != rhs.INNERID) { diff = "innerID"; std::cerr << "MetaDataElem:INNERID lhs:" << INNERID << " != " << rhs.INNERID << std::endl; }
        if (!diff) {
			return true;
		}
        return false;
#undef ELTINT
#undef ELTDOU
#undef ELTSTR
#undef SEP
    }

    void print(std::ostream& os) const {
        os << "MetaDataElems{"
#define SEP
#define ELTINT(T,N,I,C) << " " << #N << ":" << N
#define ELTDOU(T,N,I,C) << " " << #N << ":" << N
#define ELTSTR(T,N,I,C) << " " << #N << ":" << N
    METADATA_ELTS
        << " }";
#undef ELTINT
#undef ELTDOU
#undef ELTSTR
#undef SEP
    }
    
};

inline std::ostream& operator<<(std::ostream& os, MetaDataElem const & rhs) { rhs.print(os); return os; }

class MetaDataTable : public NoCopy {
public:
    MetaDataTable():N(0),latestLabelStatus(nullptr) {}
    ~MetaDataTable() { freeMetaData(); }
    void clear()     { freeMetaData(); }
    
    // initialize the metadata table
    void readFromMetaDataElements(std::vector<MetaDataElem>& _metaDataElems);

	// Read/Write *.star file from/to metadata
    void readFromStar(std::string fn_star);
    void writeToStar (std::string fn_star) const;

	//
    void printTable(std::ostream& os = std::cout) const;
    
    // get the metadata elements size
    int numberOfParticles(int random_subset = -1) const {
        if (random_subset==1) return nr_ori_particles_subset1;
        else if (random_subset==2) return nr_ori_particles_subset2;
        else return N;
    }
    void selectHalfMetadata(int random_subset) {
        if (random_subset==1) {
            subset_start = 0;subset_end = nr_ori_particles_subset1;
        }
        else if (random_subset==2) {
            subset_start = nr_ori_particles_subset1;subset_end = nr_ori_particles_subset1+nr_ori_particles_subset2;
        }
        else {
            subset_start = 0;subset_end = N;
        }
    }
    int numberOfGroups()		const 	{return nr_groups;}
    int numberOfMicrographs()	const	{return nr_micrograph;}
    //
    void append(std::vector<MetaDataElem>& _metadataElems);
    
    inline MetaDataElem       & operator[](int i)       {assert(i>-1);assert(i<(subset_end-subset_start)); return metaDataElems[metadata_order[i+subset_start]];}
    inline MetaDataElem const & operator[](int i) const {assert(i>-1);assert(i<(subset_end-subset_start)); return metaDataElems[metadata_order[i+subset_start]];}
    inline MetaDataElem 	  & accessAll(int i)		{assert(i>-1);assert(i<N); return metaDataElems[metadata_order[i]];}
    inline MetaDataElem const & accessAll(int i) const  {assert(i>-1);assert(i<N); return metaDataElems[metadata_order[i]];}
    // randomly shuffle the metadata elements
    // and split metadata to two random subset if no randomsubset read from *.star(only used in autorefine)
    void shuffle(int random_seed = -1,bool do_split_random_halves = false);
    
    // some statistics function
    void statClassImageNumber(std::vector<int> &imageNumberPerClass,bool showHist = false);
    // some fliter function
    void fliterByClass(int selectedClass);
    //
    bool containLabel(MetaDataLabel label);
    
    void unitTest(std::string fn);
    
private:
    void freeMetaData();

	// all
	class LabelStatus;
	LabelStatus*			  latestLabelStatus;
    std::vector<MetaDataElem> metaDataElems;

    // other metadata elements
    std::map<int,std::string>              undefinedMetaDataLabels; // label:column_index
    std::map<int,std::vector<std::string>> undefinedMetaDataElems;	// inner_id:....

    // metadata element number
    int N;
    int nr_groups		= 0;
    int nr_micrograph	= 0;
    // bool do_split_random_halves = false; // split data to two random subset(in autorefine)
    int nr_ori_particles_subset1 = 0;
    int nr_ori_particles_subset2 = 0;
    int subset_start,subset_end;
	//
    std::vector<int> metadata_order;

	int getImageNumber(std::string starFileName);

    void printTable(std::ostream& os, LabelStatus const& labelStatus) const;

};

enum ElemType {ElemTypeChar,ElemTypeInt,ElemTypeDouble};

class SmallMetataDataTable : public NoCopy {
private:
    std::string tableName;
    std::vector<std::string> MetaDataElemsName;
    std::vector<ElemType> MetaDataElemsType;
    std::vector<void*> MetaDataElems;
    int MetaDataElemsNumber;
public:
    SmallMetataDataTable(std::string _tableName){tableName = _tableName;}
    ~SmallMetataDataTable(){
        MetaDataElemsName.resize(0);
        MetaDataElemsType.resize(0);
        MetaDataElems.resize(0);
        MetaDataElemsNumber = 0;
    }
    void appendName(std::initializer_list<std::string> _MetaDataElemsName){
        for (auto const &__name : _MetaDataElemsName)
            MetaDataElemsName.push_back(__name);
    }
    void appendType(std::initializer_list<ElemType> _MetaDataElemsType){
        for (auto const &__type : _MetaDataElemsType)
            MetaDataElemsType.push_back(__type);
    }
    void appendElem(std::initializer_list<void*> _MetaDataElems,int _MetaDataElemsNumber){
        MetaDataElemsNumber = _MetaDataElemsNumber;
        for (auto const &__elem : _MetaDataElems)
            MetaDataElems.push_back(__elem);
    }
    void print(std::ostream& os)
    {
        assert(MetaDataElemsName.size()==MetaDataElemsType.size());
        assert(MetaDataElemsType.size()==MetaDataElems.size());
        os << tableName <<std::endl<<std::endl<<"loop_"<<std::endl;
		int ScaleCorrectionIndex = 0; // remove rlnScaleCorrection
        for (int i = 0; i < MetaDataElemsName.size(); i++) {
			if (MetaDataElemsName[i].compare("ScaleCorrection") == 0 ) { // remove rlnScaleCorrection
				ScaleCorrectionIndex = i;
				continue;
			}
			else if (i > ScaleCorrectionIndex)
				os << "_rln" << MetaDataElemsName[i] << " #" << std::to_string((long long)i) << std::endl;
            else
				os << "_rln"<<MetaDataElemsName[i]<<" #"<<std::to_string((long long)i+1)<<std::endl;
        }
        for (int metadataIndex = 0; metadataIndex < MetaDataElemsNumber; metadataIndex++) {
            for (int i = 0; i < MetaDataElems.size(); i++) {
				if (i == ScaleCorrectionIndex)  continue;
                switch (MetaDataElemsType[i]) {
                    case ElemTypeChar: {
                        auto string_ptr = (std::string**)MetaDataElems[i];
                        os << std::setw(12) << *string_ptr[metadataIndex] << " ";
                    }   break;
                    case ElemTypeInt: {
                        auto int_ptr = (int*)MetaDataElems[i];
                        os << std::setw(12) << int_ptr[metadataIndex] << " ";
                    }   break;
                    case ElemTypeDouble: {
                        auto* double_ptr = (double*)MetaDataElems[i];
                        auto var = double_ptr[metadataIndex];
                        if (var < 1 && var != 0) os << std::setw(12) << std::scientific << std::setprecision(7) << var << " ";
                        else os << std::setw(12) << std::fixed << std::setprecision(7) << var << " ";
                    }   break;
                    default:
                        break;
                }
            }
            os << std::endl;
        }
        os << std::endl << std::endl;
    }
    bool read(std::istream& is)
    {
        assert(MetaDataElemsName.size()==MetaDataElemsType.size());
        assert(MetaDataElemsType.size()==MetaDataElems.size());
        // reset to beginning
        is.clear();is.seekg(0, std::ios::beg);
        bool startingRead = false;
        std::string line;
        int metadataIndex = 0;
        //
        while (true)
        {
            getline(is,line);
            if ( (line.find(tableName) !=std::string::npos) ){
                startingRead = true;
                getline(is,line);assert(line=="");// escape a empty line
                getline(is,line);assert(line.find("loop_")!=std::string::npos);// escape loop_ line
                // head
                for (int headIndex = 0; headIndex < MetaDataElems.size(); headIndex++) {
                    getline(is,line);
                    if( line.find(MetaDataElemsName[headIndex]) == std::string::npos )
                        std::cerr<<"cannot find "+MetaDataElemsName[headIndex]+" in _model.star."<<std::endl;
                }
                break;
            }
            ERROR_CHECK(is.eof(), "end of model file,can not find "+tableName+".");
        }
        //
        bool metadata_num_fit = true;
        for (int metadataIndex = 0; metadataIndex < MetaDataElemsNumber; metadataIndex++)
        {
            getline(is,line);
            if (line.empty()) {
                metadata_num_fit = false;
                break;
            }
            std::istringstream lineStream(line);
            for (int i = 0; i < MetaDataElems.size(); i++)
            {
                switch (MetaDataElemsType[i]) {
                    case ElemTypeChar:{
                        std::string lineTmp;lineStream >> lineTmp;
                        ((std::string**)MetaDataElems[i])[metadataIndex] = 
#include "./util_heap_undefs.h"
							new std::string(lineTmp);
#include "./util_heap_defs.h"
					}	break;
                    case ElemTypeInt:
                        lineStream >> ((int*)MetaDataElems[i])[metadataIndex];
                    	break;
                    case ElemTypeDouble:
                    	lineStream >> ((double*)MetaDataElems[i])[metadataIndex];
                        break;
                    default:
                        break;
                }
            }
        }
        //if (!metadata_num_fit) {MASTERNODE std::cout<<"_rln"<<MetaDataElemsName[0]<<" has less elements to read!!!!!"<<std::endl;}
        return metadata_num_fit;
    }
    bool contain(std::istream& is,std::string head)
    {
        assert(MetaDataElemsName.size()==MetaDataElemsType.size());
        assert(MetaDataElemsType.size()==MetaDataElems.size());
        // reset to beginning
        is.clear();is.seekg(0, std::ios::beg);
        bool startingRead = false;
        std::string line;
        bool contain = false;
        //
        while (true)
        {
            getline(is,line);
            if ( (line.find(tableName) != std::string::npos) ){
                startingRead = true;
                getline(is,line);assert(line=="");// escape a empty line
                getline(is,line);assert(line.find("loop_")!=std::string::npos);// escape loop_ line
                // head
                while (startingRead) {
                    getline(is,line);
                    if(line.find("_rln") == std::string::npos) break;
                    if(line.find(head) != std::string::npos) contain = true;
                }
                break;
            }
            ERROR_CHECK(is.eof(), "end of model file,can not find "+tableName+".");
        }
        return contain;
    }
};

#endif
