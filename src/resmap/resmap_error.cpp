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

#include "resmap_error.h"

rome_error::rome_error( const char* _what_arg ,const char* _file, int _line) : what_arg(_what_arg),file(_file),line(_line) {
	std::cout<<"Errors encountered, "<<what_arg<<" @ "<<file<<":"<<line<<std::endl;
	std::cerr<<"Errors encountered, "<<what_arg<<" @ "<<file<<":"<<line<<std::endl;
}

rome_error::rome_error( const std::string& _what_arg ,const char* _file, int _line) : what_arg(_what_arg.c_str()),file(_file),line(_line) {
	std::cout<<"Errors encountered, "<<what_arg<<" @ "<<file<<":"<<line<<std::endl;
	std::cerr<<"Errors encountered, "<<what_arg<<" @ "<<file<<":"<<line<<std::endl;
}

void rome_error::assertFalse(const char*        _what_arg, const char* _file, int _line) {
    std::cerr<<_what_arg<<" "<<_file<<" "<<_line<<std::endl;
    assert(false);
	throw rome_error(_what_arg, _file, _line);
}

void rome_error::assertFalse(const std::string& _what_arg, const char* _file, int _line) {
    std::cerr<<_what_arg<<" "<<_file<<" "<<_line<<std::endl;
	assert(false);
	throw rome_error(_what_arg, _file, _line);
}
