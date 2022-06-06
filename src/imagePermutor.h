/***************************************************************************
 *
 * Authors: "Brett, Bevin"
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

#ifndef IMAGE_PERMUTOR_H_
#define IMAGE_PERMUTOR_H_

// Puts the iimage into an order that maximizes the reuse of the classes and assistors that are needed to do the various kernels with minimal memory traffic and maximal reuse
//
#include "./util.h"

class ImagePermutor {
public:
	ImagePermutor();
	virtual ~ImagePermutor();

	int numberOfImages();
	int image(int index) const;
	bool keyAlreadySet(int key);

	void setImageRange(int begin, int end);
	void setKey(int index, int key);
	void permute();

private:
	class Impl;
	Impl* _impl;
};



#endif
