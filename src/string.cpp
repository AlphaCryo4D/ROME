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

#include "string.h"
#include <algorithm>

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