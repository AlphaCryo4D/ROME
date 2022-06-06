/***************************************************************************
*
* Authors: "Yongbei(Glow) Ma,Jiayi (Timmy) Wu, Youdong (Jack) Mao"
* Dana-Farber Cancer Institute, Harvard Medical School and Peking University
* Bevin Brett
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

#include "./exp_mweight.h"

#ifndef DONT_INCLUDE_SAMPLING
#include "./sampling.h"
#endif

class Exp_Mweight_old::Sparse_Data : public NestLoopStack<std::pair<int, double>> {
};

Exp_Mweight_old::Exp_Mweight_old() :
	sparse_data(sNew(Sparse_Data)),
	nr_images  (_nr_images ),
	nr_classes (_nr_classes),
	nr_dir     (_nr_dir	   ),
	nr_psi     (_nr_psi	   ),
	_nr_images (0),
	_nr_classes(0),
	_nr_dir    (0),
	_nr_psi    (0),
	capacity   (_capacity  ),
	current_fn("NULL"),
	os_file	  (),
	max_nr_significant_ihidden(0),
	max_density(0),
	max_size   (0)
{
}

Exp_Mweight_old::~Exp_Mweight_old() {
	fini();
	sDeleteConst(sparse_data);
}

void Exp_Mweight_old::init(int _nr_images, int _nr_classes, int _nr_dir, int _nr_psi, int _capacity)
{
	this->_nr_images = _nr_images;
	this->_nr_classes = _nr_classes;
	this->_nr_dir = _nr_dir;
	this->_nr_psi = _nr_psi;
	this->_capacity = _capacity;
	sparse_data->init({ _nr_images,_nr_classes,_nr_dir,_nr_psi });
}

void Exp_Mweight_old::fini() {
	sparse_data->fini();
}

void Exp_Mweight_old::clear_image(int iimage) {
	for (int iclass = 0; iclass < nr_classes; iclass++)
		for (int idir = 0; idir < nr_dir; idir++)
			for (int ipsi = 0; ipsi < nr_psi; ipsi++) {
				auto& _data = sparse_data->wptr({ iimage,iclass,idir,ipsi });
				_data.clear();
			}
}


void Exp_Mweight_old::reset() {
	for (int iimage = 0; iimage < nr_images; iimage++)
		for (int iclass = 0; iclass < nr_classes; iclass++)
			for (int idir = 0; idir < nr_dir; idir++)
				for (int ipsi = 0; ipsi < nr_psi; ipsi++) {
					auto& _data = sparse_data->wptr({ iimage,iclass,idir,ipsi });
					_data.clear();
				}
}

std::vector< std::pair<int, double> > & Exp_Mweight_old::wptr_sparse(int iimage, int iclass, int idir, int ipsi) {
	return sparse_data->wptr({ iimage,iclass,idir,ipsi });
}

void Exp_Mweight_old::analysis(std::string output_fn, std::string note)
{
	double per_vector_size = 3.5;
	auto resetStatic = [&]() {
		max_nr_significant_ihidden = 0;
		max_density = 0;
		max_size = 0;
	};
	// write the head
	auto writeHead = [&]() {
		os_file << std::setw(20) << " " << std::setw(20) << "nr_signif_ihidden"
			<< std::setw(20) << "density" << std::setw(20) << "memory_size(GB)" << std::endl;
	};
	// write the tail
	auto writeTail = [&]() {
		os_file << std::setw(20) << "max above : " << std::setw(20) << max_nr_significant_ihidden
			<< std::setw(20) << max_density
			<< std::setw(20) << max_size << std::endl;
		os_file << std::setw(20) << " ---------- " << std::setw(20) << " --------------- "
			<< std::setw(20) << " ----------- "
			<< std::setw(20) << " ----------- " << std::endl;
	};
	//
	if (current_fn == "NULL") {
		writeTail();
		current_fn = output_fn;
		os_file.open((current_fn + ".txt").c_str(), std::ios::out);
		writeHead();
		resetStatic();
	}
	else if (current_fn != output_fn) {
		writeTail();
		os_file.close();
		current_fn = output_fn;
		os_file.open((current_fn + ".txt").c_str(), std::ios::out);
		writeHead();
		resetStatic();
	}
	// write out the info
	double nr_significant_ihidden = 0;
	for (int iimage = 0; iimage < nr_images; iimage++)
		for (int iclass = 0; iclass < nr_classes; iclass++)
			for (int idir = 0; idir < nr_dir; idir++)
				for (int ipsi = 0; ipsi < nr_psi; ipsi++)
					nr_significant_ihidden += sparse_data->wptr({ iimage,iclass,idir,ipsi }).size();

	double density = nr_significant_ihidden / double(capacity*nr_classes*nr_dir*nr_images);
	double size = per_vector_size*nr_significant_ihidden*8. / 1024. / 1024. / 1024.;
	os_file << std::setw(20) << note << std::setw(20) << nr_significant_ihidden
		<< std::setw(20) << density
		<< std::setw(20) << size << std::endl;
	// compare with the max
	if (nr_significant_ihidden > max_nr_significant_ihidden) {
		max_nr_significant_ihidden = nr_significant_ihidden;
		max_density = density;
		max_size = size;
	}
}

#ifndef DONT_INCLUDE_SAMPLING
void Exp_Mweight_old::printSamplingInfo() {
	std::cout << "//" << std::setw(15) << "healpix_order," << std::setw(15) << "oversampling," << std::setw(10) << "nr_images," << std::setw(10) << "nr_class," \
		<< std::setw(15) << "nr_direction," << std::setw(15) << "nr_rotation," << std::setw(15) << "nr_over_rot" << std::setw(15) << "nr_trans," \
		<< std::setw(15) << "nr_over_trans," << std::setw(15) << "angle_precision," << std::setw(20) << "size," << std::setw(18) << "size(GB)" << std::endl;
	for (int i = 1; i < 11; i++) {
		HealpixSampler sampling3d;
		sampling3d.initialize(2, 10, -1, i);
		int adaptive_oversampling = 1;
		size_t nr_images = 8;
		size_t nr_classes = 8;
		size_t nr_dir = sampling3d.NrDir();
		size_t nr_rot = sampling3d.NrPsi();
		size_t nr_trans = sampling3d.NrTrans();
		size_t nr_over_rot = sampling3d.oversamplingFactorOrientations(adaptive_oversampling);
		size_t nr_over_trans = sampling3d.oversamplingFactorTranslations(adaptive_oversampling);
		size_t size = nr_images*nr_classes*nr_dir*nr_rot*nr_trans*nr_over_rot*nr_over_trans;
		double size_gb = size*8. / 1024. / 1024. / 1024.;
		std::cout << "//" << std::setw(14) << i << "," << std::setw(14) << adaptive_oversampling << "," << std::setw(9) << nr_images << "," << std::setw(9) << nr_classes << "," \
			<< std::setw(14) << nr_dir << "," << std::setw(14) << nr_rot << "," << std::setw(14) << nr_over_rot << "," << std::setw(14) << nr_trans << "," \
			<< std::setw(14) << nr_over_trans << "," << std::setw(14) << 360. / sampling3d.NrPsi(adaptive_oversampling) << "," << std::setw(20) << size << "," \
			<< std::setw(14) << size_gb << "(GB)" << std::endl;
	}
}
#endif
