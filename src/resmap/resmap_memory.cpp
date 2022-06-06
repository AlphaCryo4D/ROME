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

#include "resmap_memory.h"

static size_t populationCount;


//#define DEBUG_LEAKS
#if defined(DEBUG_LEAKS)

static omp_lock_t* _debug_lock;
static void acqDebugLock() {
#include "./resmap_util_heap_undefs.h"
	if (!_debug_lock) omp_init_lock(_debug_lock = new omp_lock_t);
#include "./resmap_util_heap_defs.h"
	omp_set_lock(_debug_lock);
}
static void relDebugLock() {
	omp_unset_lock(_debug_lock);
}

static size_t populationShowTrigger = 10000;
struct BirthPlace {
    const char* file; size_t line;
    BirthPlace(const char* file, size_t line) : file(file), line(line) {}
    bool operator<(BirthPlace const & rhs) const {
        auto fd = strcmp(file, rhs.file);
        return (fd != 0) ? (fd < 0) : (line < rhs.line);
    }
};
struct Summary {
    size_t count;
    size_t size;
    Summary() : count(0), size(0) {}
};

struct CensusCountable {
	bool isVector; size_t size;
    const char* file; int line;  size_t tick;
    CensusCountable(bool isVector, size_t size, const char* file, int line, size_t tick) : isVector(isVector), size(size), file(file), line(line), tick(tick) {}
};
typedef std::map<void*, CensusCountable> Population;
Population* population = nullptr;
#endif

size_t census() {
#if defined(DEBUG_LEAKS)
	size_t result;
	acqDebugLock();
	result = population ? population->size() : 0;
	relDebugLock();
	return result;
#else
    return populationCount;
#endif
}

static CensusDetailsProvider* censusDetailsProviders = nullptr;

CensusDetailsProvider::CensusDetailsProvider() : _next(censusDetailsProviders) {
	censusDetailsProviders = this;
}

CensusDetailsProvider::~CensusDetailsProvider() {
	auto prevPtr = &censusDetailsProviders;
	while (auto next = *prevPtr) {
		if (next == this) {
			*prevPtr = next->_next;
			return;
		}
		prevPtr = &next->_next;
	}
	assert(false);
}

void showPopulation() {
#if defined(DEBUG_LEAKS)
	std::map<BirthPlace, Summary> perBP;
	acqDebugLock();
	{
		for (auto i = population->begin(); i != population->end(); i++) {
			BirthPlace bp(i->second.file, i->second.line);
			auto c = perBP.find(bp);
			if (c == perBP.end()) { perBP.insert(std::make_pair(bp, Summary())); c = perBP.find(bp); }
			c->second.count++;
			c->second.size += i->second.size;
		}
	}
	relDebugLock();
	std::cerr << "Census taken! populationCount:" << populationCount << std::endl;
	std::cout << "Census taken! populationCount:" << populationCount << std::endl;
	for (auto i = perBP.begin(); i != perBP.end(); i++) {
		std::cerr << "    " << i->first.file << ":" << i->first.line << " " << i->second.count << " " << i->second.size/1000. << "KB" << std::endl;
		std::cout << "    " << i->first.file << ":" << i->first.line << " " << i->second.count << " " << i->second.size/1000. << "KB" << std::endl;
	}

	for (auto cdp = censusDetailsProviders; cdp; cdp = cdp->next()) {
		cdp->showPopulation(std::cerr);
		cdp->showPopulation(std::cout);
	}
#endif
}

static void censusWkr(void* p, bool isVector, bool died, size_t size = 0, const char* file = nullptr, size_t line = 0) {
    if (died && !p) return;
    
    if (died)
#pragma omp atomic
        populationCount--;
    else
#pragma omp atomic
        populationCount++;
    
#if defined(DEBUG_LEAKS)

#include "./util_heap_undefs.h"
	if (!population) 
		population = new Population;
#include "./util_heap_defs.h"

	bool shouldShowPopulation = false;
	acqDebugLock();
	{
        // Note: Ptrs can be freed on different threads than the creator
        //
        static size_t watchTick = 0;
        static void*  watchP    = nullptr;
        static size_t tick      = 0;
        
        auto i = population->find(p);
        if (i == population->end()) {
            if (died) {
                std::cerr << "Undocumented alien died" << std::endl;
            } else {
                auto ii = population->insert(i, std::make_pair(p,CensusCountable(isVector, size, file, int(line), ++tick)));
                if (watchTick == tick) {
                    std::cerr << "p:" << (void*)p << " being made at tick:" << tick << " " << ii->second.file << ":" << ii->second.line << std::endl;
                    watchP = p;
                }
            }
        } else {
            if (died) {
				if (isVector != i->second.isVector) {
					std::cerr << "new [] / delete [] mismatch, made at " << i->second.file << ":" << i->second.line << std::endl;
					exitAbnormally(__FILE__,__LINE__);
				}
                population->erase(i);
                if (watchP == p) {
                    std::cerr << "watchP:" << (void*)p << " died" << std::endl;
                    watchP = nullptr;
                }
            } else {
                std::cerr << "Identify theft p:" << (void*)p << " previously made at tick:" << i->second.tick << " " << i->second.file << ":" << i->second.line << std::endl;
                EXIT_ABNORMALLY;
            }
        }
        
        if (!died && populationCount > populationShowTrigger) {
            if (populationShowTrigger < 1000000) populationShowTrigger *= 10; else populationShowTrigger += 1000000;
			shouldShowPopulation = true;
        }
    }
	relDebugLock();
	if (shouldShowPopulation) showPopulation();
#endif
}

void noteAllocatedWkr(void* p, size_t size, bool isVector, const char* file, int line) {
	censusWkr(p, isVector, false, size, file, line);
}

void noteFreedWkr(void* p, bool isVector) {
	censusWkr(p, isVector, true);
}


void* aMallocWkr(size_t size, size_t alignment, const char* file, int line) {
#include "./resmap_util_heap_undefs.h"
	if (alignment > 64) {
		std::cerr << "Can't align more than 64" << std::endl;
		EXIT_ABNORMALLY;
	}

//#define USE_GUARDS
#ifndef USE_GUARDS
	auto p = _wrapped_aligned_malloc(size, alignment);
	noteAllocatedWkr(p, size, false, file, line);
	return p;
#else
	// Put in a 64 byte guard at the beginning and the end
	// and put some info in it
	auto pc = (char*)_wrapped_aligned_malloc(size + 128, 64);
	auto ps = (size_t*)pc;
	ps[ 0] = 77;
	ps[ 1] = size;
	pc[63] = 37;		// not ps
	noteAllocatedWkr(pc+64, size, false, file, line);
	pc[64+size] = 37;
	return pc+64;
#endif

#include "./resmap_util_heap_defs.h"
}


void  aFreeWkr(void* p) {
#include "./resmap_util_heap_undefs.h"
#ifndef USE_GUARDS
	noteFreedWkr(p, false);
	_wrapped_aligned_free(p);
#else
	// Remove the 64 byte guard at the beginning
	auto pc = ((char*)p) - 64;
	auto ps = (size_t*)pc;
	// Check the info we put there is still ok
	assert(ps[ 0] == 77);
	size_t size = ps[1];
	assert(pc[63]      == 37);
	assert(pc[64+size] == 37);
	noteFreedWkr(pc+64, false);
	_wrapped_aligned_free(pc);
#endif
#include "./resmap_util_heap_defs.h"
}

namespace Heap {

    static double* malloc2DDoubles(const char* name, size_t sizeOfItem, size_t s) {
        assert(s % sizeOfItem == 0);
        auto p = mallocDoubles(s);
        return p;
    }
    static double* mallocZeroed2DDoubles(const char* name, size_t sizeOfItem, size_t s) {
        assert(s % sizeOfItem == 0);
        auto p = mallocZeroedDoubles(s);
        return p;
    }
    
    double* allocDoublesForImages(const char* name, size_t nr_images, size_t s) { return malloc2DDoubles(name, nr_images, s); }
    double* allocDoublesForClasses(const char* name, size_t nr_classes, size_t s) { return malloc2DDoubles(name, nr_classes, s); }
    double* allocZeroedDoublesForClasses(const char* name, size_t nr_classes, size_t s) { return mallocZeroed2DDoubles(name, nr_classes, s); }
};
