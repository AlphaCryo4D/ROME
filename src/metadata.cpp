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

#include "metadata.h"

#ifdef _WIN32
#include <Windows.h>
#undef GROUP_NAME
#endif

namespace LabelParser{
    static std::map<std::string,MetaDataLabel> str2label = {
        {
#define SEP },{
#define ELTINT(T,N,I,C) "rln"#C,C
#define ELTDOU(T,N,I,C) "rln"#C,C
#define ELTSTR(T,N,I,C) "rln"#C,C
            METADATA_ELTS
#undef ELTINT
#undef ELTDOU
#undef ELTSTR
#undef SEP
        }
    };
    static std::vector<std::string> label2str = {
#define SEP ,
#define ELTINT(T,N,I,C) "rln"#C
#define ELTDOU(T,N,I,C) "rln"#C
#define ELTSTR(T,N,I,C) "rln"#C
        METADATA_ELTS
#undef ELTINT
#undef ELTDOU
#undef ELTSTR
#undef SEP
    };
    MetaDataLabel stringToLabel(const std::string& name) {
        auto it = str2label.find(name);
        if (it == str2label.end()) return UnDefinedLabel;
        return str2label[name];
    };
    std::string labelToString(MetaDataLabel label){
        if (label >= label2str.size()) {
            return "UnDefinedLabel";
        }
        return label2str[label];
    };
}


class MetaDataTable::LabelStatus : public NoCopy {
public:
	std::vector<bool> metadataLabelStatus;

	static LabelStatus* make(LabelStatus const & src) {
#include "./util_heap_undefs.h"
		return sNewA(LabelStatus, (src));
#include "./util_heap_defs.h"
	}
	LabelStatus(LabelStatus const & src) : metadataLabelStatus(src.metadataLabelStatus) {}

	static LabelStatus* make(size_t size) {
#include "./util_heap_undefs.h"
		return sNewA(LabelStatus, (size));
#include "./util_heap_defs.h"
	}
	LabelStatus(size_t size) : metadataLabelStatus(size) {
		for (auto const &elt : metadataLabelStatus) {
			assert(!elt);
		}
	}
	~LabelStatus() {}
	// TODO : better way to set label status
	void flagMoreIfNecessary(std::vector<MetaDataElem> const & metaDataElems) {
	    // for each label,if it has different value with each other
	    // or different value with default one
	    // flag it as requiring writing

		MetaDataElem default_one;

	    for (int i = 0; i < metadataLabelStatus.size(); i++){

			// if this label already flagged, don't change
	        if (metadataLabelStatus[i] == true) continue;
	        
			// if this label has a different value, flag
	        std::stringstream first_one,second_one;
	        metaDataElems[0].get(MetaDataLabel(i),first_one);
	        default_one     .get(MetaDataLabel(i),second_one);

			if (first_one.str() != second_one.str()) {
	            metadataLabelStatus[i] = true;
	            continue;
	        }

	        // compare this with other elements
	        for (auto& metaDataElem : metaDataElems) {
				second_one.str("");
	            metaDataElem.get(MetaDataLabel(i),second_one);
	            if (first_one.str()!= second_one.str()) {
	                metadataLabelStatus[i] = true;
	                break;
	            }
	        }
	    }

	    if (metadataLabelStatus[DefocusU] == true)
	        metadataLabelStatus[DefocusAngle] = true;

	    // must be writed out
	    // metadataLabelStatus[NormCorrection] = true;
	}
    bool containLabel(MetaDataLabel label){
        return metadataLabelStatus[label];
    }
};


// read the metadata table
void MetaDataTable::readFromMetaDataElements(std::vector<MetaDataElem>& _metaDataElems)
{
    // TODO : it will remove some undefined label,try to keep it
	freeMetaData();

	assert(!latestLabelStatus);
	latestLabelStatus = LabelStatus::make(LabelParser::label2str.size());

	metadata_order.resize(N);
    for (int i = 0; i < N; i++) metadata_order[i] = i;
    //
    N = _metaDataElems.size();
    metaDataElems = _metaDataElems;

    latestLabelStatus->flagMoreIfNecessary(metaDataElems);
}


void MetaDataTable::readFromStar(std::string fn_star)
{
    static const bool debug = false;//true;
    
	freeMetaData();

	assert(!latestLabelStatus);
	latestLabelStatus = LabelStatus::make(LabelParser::label2str.size());

    N = getImageNumber(fn_star);
    metaDataElems .resize(N);
    metadata_order.resize(N);
    for (int i = 0; i < N; i++) metadata_order[i] = i;
    
	ifstreamCheckingExistence starFile(fn_star.c_str());
    std::string starLine;
    // column's metadata label
    std::vector<MetaDataLabel> columnLabels(124);
    for (auto &v : columnLabels) v = UnDefinedLabel;
    //
    int columnNumber = 0;
    int starHeadCount = 0;
    // read the metadata first
    while (true) {
        // get label
        getline(starFile,starLine);
        if (debug) std::cout << ".star line: " << starLine << " ,length: " << starLine.length() << std::endl;
        
        if (starLine.find("_rln") != std::string::npos) {
            // get the label
            std::string label_str;
            int i = starLine.find("_")+1;
            while (starLine[i] != ' ' && starLine[i] != '\0')
                label_str += starLine[i++];
            
            auto label = LabelParser::stringToLabel(label_str);
            if (label != UnDefinedLabel){
                latestLabelStatus->metadataLabelStatus[label] = true;  // record the label status
            }
            else {
                // std::cerr<<"warning !!!!!!! undefined label : "<<label_str<<" in star file."<<std::endl;
            }

            // get the column index
            int column_index_num;
            if(starLine.find("#") != std::string::npos){
                std::string column_index_str;
                i = starLine.find("#")+1;
                while (starLine[i] != '\0')
                    column_index_str += starLine[i++];
                column_index_num = atoi(column_index_str.c_str())-1;
            }
            else{
                column_index_num = starHeadCount;
            }

            
            if (label != UnDefinedLabel){
                columnLabels[column_index_num] = label;
            }
            else{
                undefinedMetaDataLabels[column_index_num] = label_str;
            }
            // NOTE : support some column doesnot include in head _rln
            columnNumber = std::max(columnNumber, column_index_num+1);
            
            if (debug) std::cout<<"--- find label: "<<LabelParser::labelToString(label)<<" , column index: "
                                <<" column number : "<<columnNumber<<std::endl;
            starHeadCount++;
        }
        else if(starHeadCount > 1)// skip some header line
            break;
        
        if (starFile.eof()) {
            break;
        }
    }
    ERROR_CHECK(starHeadCount == 0, "star file format error.");
    
    int bodyLine = 0;
    MetaDataElem metadata_value;
    std::string tempString;
    while (starLine.length() > 10 && bodyLine < N)
    {
        bodyLine++;
        if (debug) std::cout << starLine << ".star line: " << bodyLine << " " << starLine.length() << std::endl;
        
        std::istringstream iString(starLine);
        std::vector<std::string> undefinedMetaDataElem;
        
        for (int column_index = 0; column_index < columnNumber; column_index++) {
            iString>>tempString;
            auto& label = columnLabels[column_index];
            if (label != UnDefinedLabel) {
                metadata_value.set(label,tempString);
                if (debug) std::cout<<LabelParser::labelToString(label)<<":"<<tempString<<" ";
            }
            else if(undefinedMetaDataLabels.find(column_index) != undefinedMetaDataLabels.end()){
                undefinedMetaDataElem.push_back(tempString);
            }
        }
        if(debug) std::cout<<std::endl;
        
        metadata_value.INNERID = bodyLine-1;
        metaDataElems[bodyLine-1] = metadata_value;
        
        undefinedMetaDataElems[metadata_value.INNERID] = undefinedMetaDataElem;
        
        if (debug) std::cout<<metaDataElems[bodyLine-1]<<std::endl;
        
        if (starFile.eof()) {
            break;
        }
        //
        getline(starFile,starLine);
    }
    
    
    if (bodyLine < N)
    {
        std::stringstream message;
        message << "input star file have not enough pictures, instead of " << N << " it only has " << bodyLine;
        ERROR_REPORT(message.str());
    }
    
    latestLabelStatus->flagMoreIfNecessary(metaDataElems);
    
    // print result
    if (debug) printTable();
    
    // sort metadata by micrograph
    auto compare = [&](MetaDataElem lhs,MetaDataElem rhs){
        return lhs.MICROGRAPH < rhs.MICROGRAPH;
    };
    std::stable_sort(metaDataElems.begin(), metaDataElems.end(), compare);
    
    // set group No and micrograph
    nr_groups = 0;
    nr_micrograph = 0;
    std::map<std::string, int> groupName;
    std::map<std::string,bool> micrograph;
    for (auto& element : metaDataElems) {
        // set group number
        // TODO,if only have micrograph metadata
        auto group_no = groupName.find(element.GROUP_NAME);
        if (group_no == groupName.end()) {
            nr_groups++;
            groupName[element.GROUP_NAME]=nr_groups;
            element.GROUP_NO = nr_groups;
        }
        else{
            element.GROUP_NO = group_no->second;
        }
        // set micrograph
        auto mic = micrograph.find(element.MICROGRAPH);
        if (mic == micrograph.end()) {
            nr_micrograph++;
            micrograph[element.MICROGRAPH] = true;
        }
    }
    // if user does not set groupName
    // set group name as group
    // then same micrograph particles will belong to one group
    if (nr_groups == 1 && nr_micrograph != 1) {
        //
        groupName.clear();
        nr_groups = 0;
        for (auto& element : metaDataElems) {
            element.GROUP_NAME = element.MICROGRAPH;
            auto group_no = groupName.find(element.GROUP_NAME);
            if (group_no == groupName.end()) {
                nr_groups++;
                groupName[element.GROUP_NAME]=nr_groups;
                element.GROUP_NO = nr_groups;
            }
            else{
                element.GROUP_NO = group_no->second;
            }
        }
    }
    // set group index
    if(debug) printTable();
    //
    subset_start = 0;
    subset_end = N;
}

void MetaDataTable::writeToStar(std::string fn_star) const
{
    std::ofstream  fh;
    fh.open((fn_star+".star").c_str(), std::ios::out);
    
    printTable(fh);
    
    fh.close();
}

void MetaDataTable::printTable(std::ostream& os) const {
    latestLabelStatus->flagMoreIfNecessary(metaDataElems);
	printTable(os, *latestLabelStatus);
}

/* write _data.star file*/
void MetaDataTable::printTable(std::ostream& os, LabelStatus const & labelStatus) const {
	auto & metadataLabelStatus = labelStatus.metadataLabelStatus;
    int column_index = 0;
    // write out header
    os<<"\ndata_\n\nloop_\n";
	int ScaleCorrection = -1;
    for (int i = 0; i < metadataLabelStatus.size(); i++){
        if (metadataLabelStatus[i] == true){
			if (LabelParser::labelToString(MetaDataLabel(i)).compare("rlnScaleCorrection") == 0) {
				ScaleCorrection = i;
				continue;
			}
            os<<"_"<<LabelParser::labelToString(MetaDataLabel(i))<<" #"<<(column_index+1)<<std::endl;
            column_index++;
        }
    }
    // other undefined labels
    for (auto& undefinedLabel : undefinedMetaDataLabels)
    {
        os<<"_"<<undefinedLabel.second<<" #"<<(column_index+1)<<std::endl;
        column_index++;
    }
    
    for (auto& metaDataElem : metaDataElems) {
        for (int i = 0; i < metadataLabelStatus.size(); i++)
            if (metadataLabelStatus[i] == true && i != ScaleCorrection)
                metaDataElem.get(MetaDataLabel(i),os);
        // other undefined label
		auto& map = undefinedMetaDataElems.find(metaDataElem.INNERID)->second;
        for (auto& other : map) {
            os<<std::setw(12)<<other<<"\t";
        }
        os<<std::endl;
    }
    os<<std::endl;
}


// some statistics function
void MetaDataTable::statClassImageNumber(std::vector<int> &imageNumberPerClass,bool showHist){
    int nr_classes = 0;
    for (int img_id = 0; img_id < N; img_id++) {
        if (metaDataElems[img_id].CLASS > nr_classes)
            nr_classes = metaDataElems[img_id].CLASS;
    }
    imageNumberPerClass.resize(nr_classes,0);
    for (int img_id = 0; img_id < N; img_id++) {
        imageNumberPerClass[metaDataElems[img_id].CLASS-1] ++;
    }
    if (showHist) {
        // histogram of per-class image number distribution
        int classHist[11] = {0};
        for (int iclass = 0; iclass < nr_classes; iclass++) {
            int index = imageNumberPerClass[iclass]/10;
            if (index > 10) index = 10;
            classHist[index]++;
        }
        std::cout<<" Total class number : "<<nr_classes<<std::endl;
        for (int i = 0; i < 10; i++)
            std::cout<< " Class number (Image number between "<<i*10<<"~"<<(i+1)*10<<" ) : "<<classHist[i]<<std::endl;
        std::cout<<" Class number (Image number larger than 100 ) : "<<classHist[10]<<std::endl;
    }
}

// some fliter funtion
void MetaDataTable::fliterByClass(int selectedClass){
    std::vector<MetaDataElem> metaDataElems_fliter;
    for (int img_id = 0; img_id < N; img_id++) {
        if (metaDataElems[img_id].CLASS == selectedClass)
            metaDataElems_fliter.push_back(metaDataElems[img_id]);
    }
    N = metaDataElems_fliter.size();
    metaDataElems.resize(0);
    metaDataElems.resize(N);
    for (int img_id = 0; img_id < N; img_id++) {
        metaDataElems[img_id] = metaDataElems_fliter[img_id];
    }
    metadata_order.resize(N);
    for (int i = 0; i < N; i++) metadata_order[i] = i;
}

// randomly shuffle the Metadata
// and split metadata to two random subset if no randomsubset read from *.star(only used in autorefine)
void MetaDataTable::shuffle(int random_seed,bool do_split_random_halves)
{
    //
    ERROR_CHECK(metadata_order.size()!=N, "metadata_order size error.");
    //
    if (do_split_random_halves)
    {
        nr_ori_particles_subset1 = 0;
        nr_ori_particles_subset2 = 0;
        if (latestLabelStatus->containLabel(RandomSubset))
        {
            for (auto& metaDataElem : metaDataElems) {
                int random_subset = metaDataElem.SUBSET;
                if (random_subset == 1) nr_ori_particles_subset1++;
                else if (random_subset == 2) nr_ori_particles_subset2++;
                else ERROR_REPORT("ERROR MetaDataTable::shuffleAndSplit: invalid number for random subset (i.e. not 1 or 2): " + std::to_string((long long)random_subset));
            }
        }
        else
        {
            srand(random_seed);
            for (auto& metaDataElem : metaDataElems) {
                int random_subset = rand() % 2 + 1;
                metaDataElem.SUBSET = random_subset; // randomly 1 or 2
                if (random_subset == 1) nr_ori_particles_subset1++;
                else if (random_subset == 2) nr_ori_particles_subset2++;
                else ERROR_REPORT("ERROR MetaDataTable::shuffleAndSplit: invalid number for random subset (i.e. not 1 or 2): " + std::to_string((long long)random_subset));
            }
        }
    }
    //
    // Randomise
    if(random_seed!=-1) srand(random_seed);
    if (do_split_random_halves)
    {
        std::vector<int> metadata_order1, metadata_order2;
        assert(metaDataElems.size()==nr_ori_particles_subset1+nr_ori_particles_subset2);
        assert(nr_ori_particles_subset1>5);
        assert(nr_ori_particles_subset2>5);
        // Fill the two particle lists
        for (int i = 0; i < metaDataElems.size(); i++)
        {
            int random_subset = metaDataElems[i].SUBSET;
            if (random_subset == 1) metadata_order1.push_back(i);
            else if (random_subset == 2) metadata_order2.push_back(i);
            else ERROR_REPORT("ERROR MetaDataTable::shuffleAndSplit: invalid number for random subset (i.e. not 1 or 2): " + std::to_string((long long)random_subset));
        }
        
        // Just a silly check for the sizes of the ori_particle_lists (to be sure)
        assert(metadata_order1.size()==nr_ori_particles_subset1);
        assert(metadata_order2.size()==nr_ori_particles_subset2);
        
        // Randomise the two particle lists
        std::random_shuffle(metadata_order1.begin(), metadata_order1.end());
        std::random_shuffle(metadata_order2.begin(), metadata_order2.end());
        
        // First fill metadata_order with the first subset, then with the second
        metadata_order.clear();
        for (const auto& i : metadata_order1) metadata_order.push_back(i);
        for (const auto& i : metadata_order2) metadata_order.push_back(i);
        assert(metadata_order.size()==metadata_order1.size()+metadata_order2.size());
    }
    else
    {
        std::random_shuffle(metadata_order.begin(), metadata_order.end());
    }
    
}


void MetaDataTable::append(std::vector<MetaDataElem>& _metadataElems){
    for (auto& v : _metadataElems) {
        metaDataElems.push_back(v);
        N++;
    }
}

int MetaDataTable::getImageNumber(std::string starFileName)
{
    int N = 0;
    ifstreamCheckingExistence starFile(starFileName.c_str());
    std::string starLine;
    while (getline(starFile,starLine)) {
        if (starLine.find(".mrcs") != std::string::npos) N++;
    }
    starFile.close();
    return N;
}

void MetaDataTable::freeMetaData(){
    if (N == 0) return;
    N = 0;
	sDelete(latestLabelStatus);
    undefinedMetaDataLabels.clear();
    undefinedMetaDataElems.clear();
    metadata_order.resize(0);
    metaDataElems.resize(0);
}

bool MetaDataTable::containLabel(MetaDataLabel label){
    return latestLabelStatus->containLabel(label);
}

void MetaDataTable::unitTest(std::string fn){
    readFromStar(fn+".star");
    writeToStar(fn+"_result");
    std::cout<<"hello world."<<std::endl;
}

