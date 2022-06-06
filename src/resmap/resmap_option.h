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

#ifndef OPTION_H_
#define OPTION_H_

#include "resmap_util.h"
#include "resmap_mpi.h"
#include "resmap_error.h"
#include "resmap_string.h"

class Option {
private:
    // all options
    std::map<std::string, std::string> IOParser;
    // all optinal options
    typedef struct {
        std::string key;
        std::string value;
        std::string comment;
    } Elem;
    std::vector<Elem> AllOptions;
    //
    size_t maxKeyLength,maxCommentLength;
    
public:
	Option() : maxKeyLength(0), maxCommentLength(0) {}
    ~Option() { IOParser.clear(); AllOptions.clear(); }
    
    // add all options
	static std::string unspecified() { return "<unspecified>"; }
    void addOption(std::string key,std::string commit,std::string default_value = unspecified());
    // read command line
    void readCommandLine(int argc, char * argv[]);
	void readIncludeFile(std::string fnm);
    
    // read option
    std::string getOption(std::string key);
    int getIntOption(std::string key);
    double getFloatOption(std::string key);
    std::string getStrOption(std::string key);
    bool getBoolOption(std::string key);
    
    //
    void printValue();
    //
    void printHelp();
private:
	int readArg(const char* arg0, const char* arg1);
};

#endif
