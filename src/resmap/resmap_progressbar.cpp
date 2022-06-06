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

#include "resmap_util.h"		// used for building precompiled headers on Windows

#include "resmap_progressbar.h"
#include "./resmap_time.h"

void showProgressBar(int curr, int max)
{
    int maxlen = max;
    if (maxlen > 59) maxlen = 59;
    int len, rate;
    char buf[maxlen+1];

    if (curr >= max) {
        len = maxlen;
        rate = 100;
    } else {
        len = (int)(curr * maxlen / max);
        rate = (int)(curr * 100 / max);
    }

	static int          prevLen;
	if (prevLen == len) return;
	prevLen = len;

	static bool			laterTime = false;
	static Microseconds prevNow;
	Microseconds now = timeInMicroseconds();
	if (laterTime && (now - prevNow < 500000) && (len < maxlen-4)) return;
	prevNow = now;
	laterTime = true;

	memset(buf, ' ', maxlen);
    memset(buf, '=', len);
    if (len>0) {
        buf[len-1] = '>';
    }
    buf[maxlen] = 0;

    std::cerr<<"\r"<<std::setfill(' ')<<std::setw(3)<<rate<<"% ["<<buf<<"]";
}

