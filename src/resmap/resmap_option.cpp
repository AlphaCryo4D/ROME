/***************************************************************************
 *
 * Intel® Parallel Computing Center for Structural Biology
 * Principal Investigator : Youdong (Jack) Mao (Youdong_Mao@dfci.harvard.edu)
 * Dana-Farber Cancer Institute, Harvard Medical School and Peking University
 *
 * Authors: "Yong Bei Ma(galowma@gmail.com) Bevin R Brett(bevin_brett@hotmail.com)"
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

#include "resmap_util.h"		// used for building precompiled headers on Windows

#include "resmap_option.h"

static inline std::string trimKey(std::string key) {
	// skip leading whitespace
	auto nws = key.find_first_not_of(" \t");
	if (nws != std::string::npos) key = key.substr(nws, key.size()-nws);
	// see if it is an option, if not return the empty string
    if (key[0] != '-') return std::string();
	auto nh = key.find_first_not_of("-");
    return (nh != std::string::npos) ? key.substr(nh-1, key.size()-nh+1) : key;
}

// add all options
void Option::addOption(std::string key,std::string comment,std::string default_value){
    Elem oneOption;
    key = trimKey(key);
    maxKeyLength = std::max(maxKeyLength,key.length());
    size_t commentLength = comment.find('\n');
	if (commentLength == std::string::npos) commentLength = comment.size();
    maxCommentLength  = std::max(maxCommentLength, commentLength);	// assumes the longest line in a comment with \n is the first line
    oneOption.key     = key;
	oneOption.value   = default_value;
	oneOption.comment = comment;
    AllOptions.push_back(oneOption);
    IOParser[key] = default_value;
}

int Option::readArg(const char* arg0, const char* arg1) {

	auto key = trimKey(arg0);
	if (IOParser.find(key) == IOParser.end()) { // if it does not have this option, reject it
		std::cerr << "command line option '" << key << "' is not known." << std::endl;
		return 1;
	}
	
	if (arg1) {
		std::string value = arg1;
		if (IOParser.find(trimKey(value)) == IOParser.end()) { // if the value is not an option it must be the argument - weird, but supports negative numbers, sigh
			IOParser[key] = value;
			return 2;
		}
	}

	IOParser[key] = "1";
	return 1;
}

// read command line
void Option::readCommandLine(int argc, char * argv[]) {
    // add options one by one
    for (int i = 1; i < argc;) {
		i += readArg(argv[i], (i+1<argc) ? argv[i+1] : nullptr);
    }
}

void Option::readIncludeFile(std::string fnm) {
	ifstreamCheckingExistence is(fnm.c_str());
	while (!is.eof()) {
		std::string key,value;
		std::getline(is,key);
		if (is.fail()||is.bad()) break;
		auto nws = key.find_first_not_of(" \t");
		if (nws != std::string::npos) key = key.substr(nws,key.length()-nws);
		if (key.size() > 0 && key[0] == '#') continue;
		auto ws = key.find_first_of(" \t");
		if (ws != std::string::npos) {
			value = key.substr(ws, key.size()-ws);
			key   = key.substr(0, ws);
			auto nws = value.find_first_not_of(" \t");
			value    = value.substr(nws, value.size()-nws);
		}
		readArg(key.c_str(), value.size() ? value.c_str() : nullptr);
	}
}

std::string Option::getOption(std::string key){
    // trim the key to format : "-v"
    key = trimKey(key);
    if (IOParser.find(key) == IOParser.end()) {
        std::cerr<<"make sure you have added "<<key<<" option."<<std::endl;
        ERROR_REPORT("missing some key in option.")
    }
    std::string value = IOParser[key];
    if (value == Option::unspecified()) {
        std::cerr<<"you need input "<<key<<" option , which means : "<<std::endl;
        std::for_each(AllOptions.begin(), AllOptions.end(),
                      [&](Elem& elem){if(elem.key == key) std::cout<<elem.comment<<std::endl;});
        ERROR_REPORT("missing some input parameter.")
    }
    return value;
}

// read option
int Option::getIntOption(std::string key){
    std::string value = getOption(key);
    return atoi(value.c_str());
}

double Option::getFloatOption(std::string key){
    std::string value = getOption(key);
    float retval;
    sscanf(value.c_str(), "%f", &retval);
    return retval;
}

std::string Option::getStrOption(std::string key){
    std::string value = getOption(key);
    return value;
}

bool Option::getBoolOption(std::string key){
    std::string value = getOption(key);
    return bool(atoi(value.c_str()));
}

void Option::printValue(){
    if (!MPI_IS_ROOT) return;
    std::cout<<"--------------------------------  user's options    ------------------------------------------------------------"<<std::endl;
    std::for_each(AllOptions.begin(), AllOptions.end(),
                  [&](Elem& elem){
                      std::cout<<std::setw(maxKeyLength)<<std::left<<elem.key<<" , ";
                      std::cout<<std::setw(maxCommentLength)<<elem.comment;
                      std::cout<<" , "<<IOParser[elem.key]<<std::endl;});
    std::cout<<"---------------------------------------------------------------------------------------------------------------"<<std::endl;
}

void Option::printHelp(){
    if (!MPI_IS_ROOT) return;
    std::cerr<<"----------------------  general option(need to set)      -------------------------------------------------------"<<std::endl;
    std::for_each(AllOptions.begin(), AllOptions.end(),
                  [&](Elem& elem){if(elem.value == Option::unspecified())
                      std::cerr<<std::setw(maxKeyLength)<<std::left<<elem.key<<" , "<<elem.comment<<std::endl;});
    std::cerr<<"--------------------  advanced option(with default value)  -----------------------------------------------------"<<std::endl;
    std::for_each(AllOptions.begin(), AllOptions.end(),
                  [&](Elem& elem){if(elem.value != Option::unspecified())
                      std::cerr<<std::setw(maxKeyLength)<<elem.key<<" , "<<std::setw(maxCommentLength)<<elem.comment<<" , default : "<<elem.value<<std::endl;});
    std::cerr<<"----------------------------------------------------------------------------------------------------------------"<<std::endl;
}

