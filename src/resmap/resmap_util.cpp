/***************************************************************************
 *
 * Intel® Parallel Computing Center for Structural Biology
 * Principal Investigator : Youdong (Jack) Mao (Youdong_Mao@dfci.harvard.edu)
 * Dana-Farber Cancer Institute, Harvard Medical School and Peking University
 *
 * Authors: "Bevin R Brett(bevin_brett@hotmail.com)"
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

#include <algorithm>
#include <ostream>
#include <sstream>

static int initial_omp_get_max_threads_maxThreads;
int initial_omp_get_max_threads() {
	if (!initial_omp_get_max_threads_maxThreads)
#undef omp_get_max_threads
		initial_omp_get_max_threads_maxThreads = omp_get_max_threads();
#define omp_get_max_threads initial_omp_get_max_threads
	return initial_omp_get_max_threads_maxThreads;
}
static std::string goSerialLabel, goParallelLabel;
static bool goSerialLabelFound,   goParallelLabelFound;
void setGoSerial(std::string label) {
	std::cerr << "setGoSerial looking for " << label << std::endl;
	goSerialLabel = label;
	goSerialLabelFound = false;
}
void setGoParallel(std::string label) {
	std::cerr << "setGoParallel looking for " << label << std::endl;
	goParallelLabel = label;
	goParallelLabelFound = false;
}
static std::map<std::string,bool> goSerialLabels;
void maybeGoSerial(const char* label) {
	if (goSerialLabel.size() == 0) return;
	if (goSerialLabels.find(label) == goSerialLabels.end()) {
		goSerialLabels.insert(std::make_pair<std::string,bool>(label, true));
		std::cerr << "maybeGoSerial reached " << label << std::endl;
		std::cout << "maybeGoSerial reached " << label << std::endl;
	}
	if (label == goParallelLabel) {
		goSerialLabelFound = false;
		if (goParallelLabelFound) return;
		goParallelLabelFound = true;
		auto maxThreads = omp_get_max_threads();
		std::cerr << "maybeGoSerial switching to parallel " << maxThreads << " at " << label << std::endl;
		std::cout << "maybeGoSerial switching to parallel " << maxThreads << " at " << label << std::endl;
		omp_set_num_threads(maxThreads);
	}
	if (label == goSerialLabel) {
		goParallelLabelFound = false;
		if (goSerialLabelFound) return;
		goSerialLabelFound = true;
		std::cerr << "maybeGoSerial switching to serial at " << label << std::endl;
		std::cout << "maybeGoSerial switching to serial at " << label << std::endl;
		omp_set_num_threads(1);
	}
}


Lock& PerformanceCounter_lock() {
	static Lock* lock;
	if (!lock) lock = sNew(Lock);	// and initialize it
	return *lock;
}


static PerformanceCounter* PerformanceCounter_head;

PerformanceCounter::PerformanceCounter() {
	ScopedAcquire acquire(PerformanceCounter_lock(), __FILE__, __LINE__);
	if (!PerformanceCounter_head) {
		PerformanceCounter_head = _prev = _next = this;
	} else {
		_next = PerformanceCounter_head;
		_prev = _next->_prev;
		_prev->_next = this;
		_next->_prev = this;
	}
}

PerformanceCounter::~PerformanceCounter() {
	ScopedAcquire acquire(PerformanceCounter_lock(), __FILE__, __LINE__);
	assert(!!PerformanceCounter_head);
	_prev->_next = _next;
	_next->_prev = _prev;
	PerformanceCounter_head = _next;
	if (PerformanceCounter_head == this) PerformanceCounter_head = nullptr;
}

void PerformanceCounter::showAll(std::ostream & os, bool andReset) {
	ScopedAcquire acquire(PerformanceCounter_lock(), __FILE__, __LINE__);
	if (!PerformanceCounter_head) return;
	auto p = PerformanceCounter_head; 
	for (;;) {
		p->show(os);
		if (andReset) p->reset();
		p = p->_next;
		if (p == PerformanceCounter_head) break;
	}
}

void IntPerformanceCounter::show(std::ostream & os) {
	if (count.v != 0) os << name << ":" << count.v << std::endl;
}

void IntPerformanceCounter::reset() {
	count.v = 0;
}


static std::atomic<int> next_uid = 0;

NoCopy::NoCopy() : _uid(++next_uid) {
	if (interesting()) std::cerr << "interesting NoCopy made _uid:" << uid() << std::endl;
}

void NoCopy::changeUid() {
	_uid = (++next_uid);
}

bool NoCopy::interesting() const {
	return false
		;
}

void exitAbnormally(const char* file, int line) {
	std::cout << "exitAbnormally called at " << file << ":" << line << std::endl; 	
	std::cerr << "exitAbnormally called at " << file << ":" << line << std::endl; 	
#include "./resmap_util_heap_undefs.h"
	exit(1);
#include "./resmap_util_heap_defs.h"
}


DoublePair::DoublePair(std::string s) {
	std::istringstream is(s);
	is >> lo;
	char colon; is >> colon;
	is >> hi;
}

std::string DoublePair::string() const {
	return std::to_string((long double)lo) + ":" + std::to_string((long double)hi);
}

ifstreamCheckingExistence::ifstreamCheckingExistence(const char* fileName) : std::ifstream(fileName) {
	if (!fail()) return;
	std::cerr << "Failed to open " << fileName << std::endl;
	std::cerr << "Current directory: " << currentDirectory() << std::endl;
	EXIT_ABNORMALLY;
}

ofstreamCheckingCreated::ofstreamCheckingCreated(const char* fileName) : std::ofstream(fileName) {
	if (!fail()) return;
	std::cerr << "Failed to create " << fileName << std::endl;
	std::cerr << "Current directory: " << currentDirectory() << std::endl;
	EXIT_ABNORMALLY;
}

template <> bool nearEnoughTemplate(double const & tentative, double const & known, bool info) {
	auto max = std::max(std::abs(tentative), std::abs(known));
	if (info) {
		std::cerr
			<< "nearEnoughTemplate<float> " << tentative << " ?= " << known
			<< ": Rel err (" << std::abs(tentative - known) / max << ") too large ( > 1.0e-1)"
			<< std::endl;
		std::cerr << std::endl;
	}
	return (max <= 1.0e-2) || std::abs(tentative - known) / max < 1.0e-1;
}

template <> bool nearEnoughTemplate(float const & tentative, float const & known, bool info) {
	double t = tentative;
	double k = known;
	return nearEnoughTemplate(t, k, info);
}

void MeanAndStdDeviation::compute() {
	if (_resultComputed) return;
	_resultComputed = true;
	_mean = _sum / _n;
	_stdDev = 0;
	if (_n <= 1) return;
	_stdDev = _sumOfSquares / _n - square(_mean);
	_stdDev *= _n / (_n - 1);						// See https://en.wikipedia.org/wiki/Standard_deviation#Corrected_sample_standard_deviation
	_stdDev = sqrt(std::abs(_stdDev));				// abs needed when rounding causes _stddev to be slightly negative
}
