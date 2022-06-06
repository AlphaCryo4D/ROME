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

#ifndef STRING_H_
#define STRING_H_

#include "resmap_util.h"

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

/**
 * Push nothing to a stream.
 */
template<typename Stream_>
void stream_push(Stream_ && stream) {
}

/**
 * Push sth to a stream.
 */
template<typename Stream_, typename Value_, typename... Values_>
void stream_push(Stream_ && stream, Value_ && value, Values_ && ...values) {
    stream << value;
    stream_push(stream, values...);
}

/**
 * Merge strings.
 */
template<typename... Strs_>
std::string strMerge(Strs_ && ...strs) {
    std::stringstream stream;
    stream_push(stream, strs...);
    return stream.str();
}

namespace os { namespace path {

/**
 * Split a path name from the extension.
 */
std::pair<std::string, std::string> splitext(std::string name);

/**
 * Get the base name of a path.
 */
std::string basename(std::string name);

}} // namespace os::path

#endif
