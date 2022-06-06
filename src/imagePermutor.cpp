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

#include "util.h"		// used for building precompiled headers on Windows

#include "imagePermutor.h"

class ImagePermutor::Impl : NoCopy {
public:
	Impl() : significantKeys(0) {}
	~Impl() {}
	int significantKeys;
	int image_begin;
	struct Pair { 
		int key; 
		int value; 
		Pair() {}
		Pair(int key, int value) : key(key),value(value) {}
	};
	std::vector<Pair> v;
};

ImagePermutor::ImagePermutor() : _impl(sNew(Impl)) {}
ImagePermutor::~ImagePermutor() { sDelete(_impl); }

int ImagePermutor::numberOfImages() {
	return _impl->v.size();
}

int ImagePermutor::image(int index) const {
	auto & v = _impl->v;
	return v[index - _impl->image_begin].key;
}

void ImagePermutor::setImageRange(int begin, int end) {
	auto & v = _impl->v;
	_impl->image_begin = begin;
	auto size = end - begin;
	v.resize(size);
	for (int i = 0; i < size; i++) { v[i] = Impl::Pair(i,std::numeric_limits<int>::max()); }	// every index must be in it
}

bool ImagePermutor::keyAlreadySet(int key) {
	auto & vi = _impl->v[key - _impl->image_begin];
	return vi.value != std::numeric_limits<int>::max();
}

void ImagePermutor::setKey(int key, int value) {
	if (keyAlreadySet(key)) return;
	_impl->significantKeys++;
	assert(key >= 0);
	auto & vi = _impl->v[key - _impl->image_begin];
	assert(vi.key == key);
	vi.value = value;
}

void ImagePermutor::permute() {
	auto & v = _impl->v;
	std::sort(v.begin(), v.end(),
		[&](Impl::Pair const & lhs, Impl::Pair const & rhs)->bool {
		return lhs.value < rhs.value; 
	});
	// TODO - truncate the unwritten values

	if (0) {
		static int count = 0;
		if (count++ < 10) {
			std::cerr << "ImagePermutor::permute, v.size:" << v.size() << " with significantKeys:" << _impl->significantKeys << std::endl;
			int i = 0;
			for (; i < std::min(10, int(v.size())); i++) {
				std::cerr << "v[" << i << "].key:" << v[i].key << "->" << v[i].value << std::endl;
			}
			i = std::max(i, int(v.size()) - 10);
			for (; i < v.size(); i++) {
				std::cerr << "v[" << i << "].key:" << v[i].key << "->" << v[i].value << std::endl;
			}
		}
	}
}
