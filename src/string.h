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

#ifndef STRING_H_
#define STRING_H_

#include <string>

// remove the suffix of a path
// example : ../Document/_iter10.star    ==> ./root/Document/_iter10
// example : /root/Document/_iter10.mrcs ==> /root/Document/_iter10
std::string pathRemoveSuffix(std::string path);

// get the filename
// example : /root/Document/_iter10.mrcs ==> _iter10.mrcs
// example : /root/Document/_iter10      ==> _iter10
std::string pathGetFilename(std::string path);

// get path extension
// example : /root/Document/_iter10.mrcs ==> .mrcs
std::string pathGetSuffix(std::string path);

// get path directory
// example : /root/Document/_iter10.mrcs ==> /root/Document/
std::string pathGetDir(std::string path);

// convert number to string
// example : 23 ==> "023"
std::string num2str(int number,int len = 3);

#endif
