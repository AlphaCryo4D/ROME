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

#include "resmap_util.h"		// used for building precompiled headers on Windows

#include "resmap_string.h"

std::string pathRemoveSuffix(std::string path)
{
    /** 
    auto pos = path.find_last_of(".");
    // only work if path is a filename!
    if ((path.substr(pos)).find("/")==std::string::npos)
        return path.substr(0,pos);
    else
        return path;
     **/
    int pos = -1;
    if (path.find(".star") != std::string::npos) pos = path.find(".star");
    else if (path.find(".mrcs") != std::string::npos) pos = path.find(".mrcs");
    else if (path.find(".dat") != std::string::npos) pos = path.find(".dat");
    else if (path.find(".mrc") != std::string::npos) pos = path.find(".mrc");
    if (pos != -1)
        return path.substr(0,pos);
    else
        return path;
}

std::string pathGetFilename(std::string path)
{
    return path.substr(path.find_last_of("/")+1);
}

std::string pathGetSuffix(std::string path){
    return path.substr(path.find_last_of("."));
}

std::string pathGetDir(std::string path)
{
	if (path == "NULL") return path;
    return path.substr(0,path.find_last_of("/")+1);
}

std::string num2str(int number,int len)
{
	static const int maxLen = 255;
	len = std::min(len, maxLen);
    char numStr[maxLen+1];
    // format : 00x
    for (int i = 0; i < len; i++) {
        numStr[len-i-1] = '0' + number % 10;
        number = number / 10;
    }
    numStr[len] = 0;
    return (numStr);
}

namespace os { namespace path {

std::pair<std::string, std::string> splitext(std::string name) {
    int pos = name.find_last_of('.');
    if (pos == std::string::npos) {
        return {name, ""};
    }
    else {
        return {name.substr(0, pos), name.substr(pos)};
    }
//    std::smatch match;
//    if (std::regex_match(name, match, std::regex("^(.*)(\\.[^.]+)$"))) {
//        return {match[1], match[2]};
//    }
//    else {
//        return {name, ""};
//    }
}

std::string basename(std::string name) {
    int pos = name.find_last_of('/');
    if (pos == std::string::npos) {
        return name;
    }
    else {
        if (name.size() == pos+1) return "";
        else return name.substr(pos+1);
    }
//    std::smatch match;
//    if (std::regex_match(name, match, std::regex("^.*/([^/]+)$"))) {
//        return match[1];
//    }
//    else {
//        return name;
//    }
}

}} // namespace os::path

