/***************************************************************************
 *
 * Authors: "Yongbei(Glow) Ma,Jiayi (Timmy) Wu, Youdong (Jack) Mao"
 * Dana-Farber Cancer Institute, Harvard Medical School and Peking University
 * Authors: "Brett, Bevin"
 * Intel Corporation
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

#include "memory.h"

static size_t populationCount;


//#define DEBUG_LEAKS
#if defined(DEBUG_LEAKS)

static omp_lock_t* _debug_lock;
static void acqDebugLock() {
#include "./util_heap_undefs.h"
	if (!_debug_lock) omp_init_lock(_debug_lock = new omp_lock_t);
#include "./util_heap_defs.h"
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
#include "./util_heap_undefs.h"
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

#include "./util_heap_defs.h"
}


void  aFreeWkr(void* p) {
#include "./util_heap_undefs.h"
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
#include "./util_heap_defs.h"
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



// L2 cache modelling is used to try to understand our L2 miss rates in production runs without distort performance
//
namespace L2CacheModel {

	class SeqAcc {
	public:
		SeqAcc* _next;
#define SEP			
#define ELT(T,N)	T N;
		L2_ACCESS_ELTS
#undef ELT
#undef SEP
	};

	struct PerThread {
		size_t  _active;
		SeqAcc* _oldest;
		SeqAcc* _newest;
		PerThread() : _active(0), _oldest(nullptr), _newest(nullptr) {}
	};

	std::vector<PerThread> perThreads(omp_get_max_threads());

	class CHiP {
	public:
		static const size_t cacheLineSizeLog2	=  6;	// 64 Byte cache lines
		static const size_t cacheLineSize		= (size_t(1) << cacheLineSizeLog2);
		static const size_t cacheLineOffsetMask	= cacheLineSize - 1;
		static const size_t maxCacheSizeLog2	= 		// KNL 1MB max interesting cache size,			20
												  		// but shared by two cores						19
													17;	// but BDW and SKL L2 caches are only 256KB		17
		static const size_t usableCacheSize     = (size_t(1) << (maxCacheSizeLog2-1));
														// Only try to half-fill the available cache 
		static const size_t depthLog2			=  4;	// allows for an efficient implementation
		static const size_t widthLog2			= maxCacheSizeLog2 - cacheLineSizeLog2 - depthLog2;
		static const size_t width				= 1 << widthLog2;
		static const size_t depth				= 1 << depthLog2;

		__declspec(align(64))
		size_t cols    [width*depth];
		size_t rowCount[depth];
		size_t misses;
		CHiP() {
			for (size_t i = 0; i < width*depth; i++) cols    [i] = 0;
			for (size_t r = 0; r <       depth; r++) rowCount[r] = 0;
			misses = 0;
		}
		size_t total() {
			size_t total(misses);
			for (size_t r = 0; r < depth; r++) {
				total += rowCount[r];
			}
			return total;
		}
		float hitRate() {
			size_t sizeSoFar(0);
			size_t total(0);
			size_t r = 0;
			for (; r < depth; r++) {
				total += rowCount[r];
				sizeSoFar += width*cacheLineSize;
				if (sizeSoFar > usableCacheSize) break;
			}
			auto hits = total;
			for (; r < depth; r++) {
				total += rowCount[r];
			}
			total += misses;
			return hits / float(std::max(size_t(1),total));
		}
		float missPct() {
			return 0.1f * int(1000.0f*hitRate());
		}
		bool goodHitRate() {
			return missPct() < 0.5;
		}
		void showSummary(std::ostream & os, Life* life) {
			const size_t total = this->total();
			auto missPct = this->missPct();
			os << "L2 misses/total " << misses << "/" << total
				<< "(" << ((total == 0) ? 0 : missPct) << "%)";
			if (life) {
				os << " during " << life->path();
			}
			os << " tid:" << omp_get_thread_num() << std::endl;
		}
		void showDetails(std::ostream & os, Life* life) {
			const size_t total = this->total();
			if (total == 0) return;
			size_t size(0), sum(0);
			for (size_t r = 0; r < depth; r++) {
				size += width*cacheLineSize; sum += rowCount[r];
				auto pct = 0.1f * int(1000.0 * sum / total);
				os << "    " << std::setw(6) << size << " " << pct << std::endl;
				if (pct >= 99.9) break;
			}
		}
 		void access(size_t addr, size_t count, size_t amount, size_t stride) {
			if (amount == stride) {
				while (count > 0 && !(count & 1) && (amount < cacheLineSize)) {
					amount *= 2;
					stride *= 2;
					count  /= 2;
				}
			}
			while (count--) {
				size_t remains = amount;
				while (remains) {
					size_t o     = addr & cacheLineOffsetMask;
					size_t eaten = std::min(amount, cacheLineSize - o);
					access(addr >> cacheLineSizeLog2);
					addr    += eaten;
					remains -= eaten;
				}
				addr += stride - amount;
			}
		}
		void access(size_t line) {
			size_t c = line & (width-1);
			size_t* col = &cols[c * depth];
			size_t insert = line;
			for (size_t r = 0; r < depth; r++) {
				size_t prev = col[r];
				col[r] = insert;
				if (prev == line) {
					rowCount[r]++;
					return;
				}
				insert = prev;
			}
			misses++;
		}
	};

#ifdef L2_CACHE_MODELING
	void seqAccWkr(
#define SEP			,
#define ELT(T,N)	T N
		L2_ACCESS_ELTS
#undef ELT
#undef SEP
	) {
		auto & pt = perThreads[omp_get_thread_num()];
		auto &_active = pt._active;
		auto &_oldest = pt._oldest;
		auto &_newest = pt._newest;

		if (!_active) return;
		auto p = sNew(SeqAcc);
		p->_next = nullptr;
#define SEP			
#define ELT(T,N)	p->N = N;
		L2_ACCESS_ELTS
#undef ELT
#undef SEP
		if (_newest) _newest->_next = p; else _oldest = p;
		_newest = p;

		// BEVIN
		// std::cerr << "seqAccWkr list is now " << std::endl;
		// for (auto p = _oldest; !!p; p = p->_next) {
		// 	std::cerr << "   " << p->nameBase << "[" << p->nameIndex << "]" << std::endl;
		// }
	}
#endif

	class Interval::Impl {
	public:
		size_t  _depth;
		SeqAcc* _beginNewest;
		Impl() : _depth(0), _beginNewest(nullptr) {}
	};

	void Interval::begin() {
#ifdef L2_CACHE_MODELING
		if (!!_impl) {
			_impl->_depth++;
			return;
		}
		auto & pt = perThreads[omp_get_thread_num()];
		auto &_active = pt._active;
		auto &_oldest = pt._oldest;
		auto &_newest = pt._newest;
		_impl = sNew(Impl);
		_impl->_beginNewest = _newest;
		_active++;
#endif
	}

	void Interval::endWkr() {
		if (--_impl->_depth) return;
		sDelete(_impl);
		auto & pt = perThreads[omp_get_thread_num()];
		auto &_active = pt._active;
		auto &_oldest = pt._oldest;
		auto &_newest = pt._newest;
		--_active;
		if (_active > 0) return;
		while (_oldest) {
			auto next = _oldest->_next;
			sDelete(_oldest);
			_oldest = next;
		}
		_newest = nullptr;
	}

	size_t Interval::accessesRecordedCount() {
#ifndef L2_CACHE_MODELING
		return 1;
#else
		if (!_impl) return 0;
		auto & pt = perThreads[omp_get_thread_num()];
		auto &_active = pt._active;
		auto &_oldest = pt._oldest;
		auto &_newest = pt._newest;
		size_t sum(0);
		for (auto p = !!_impl->_beginNewest ? _impl->_beginNewest->_next : _oldest; p; p = p->_next) {
			sum++;
	}
		return sum;
#endif

	}

	size_t Interval::cacheProbeCount() {
#ifndef L2_CACHE_MODELING
		return 1;
#else
		if (!_impl) return 0;
		auto & pt = perThreads[omp_get_thread_num()];
		auto &_active = pt._active;
		auto &_oldest = pt._oldest;
		auto &_newest = pt._newest;
		size_t sum(0);
		for (auto p = !!_impl->_beginNewest ? _impl->_beginNewest->_next : _oldest; p; p = p->_next) {
			sum += (p->count*p->amount);
		}
		return sum / CHiP::cacheLineSize;		// not right but good enough
#endif
	}

	bool Interval::showL2CacheIntervalUnlocked(Life* life) {
		       showL2CacheInterval(std::cerr, life);
		return showL2CacheInterval(std::cout, life);
	}

	float Interval::hitRate() {
#ifndef L2_CACHE_MODELING
		return 1.0;
#else
		if (!_impl) return 1.0;
		auto & pt = perThreads[omp_get_thread_num()];
		auto &_active = pt._active;
		auto &_oldest = pt._oldest;
		auto &_newest = pt._newest;
		CHiP chip;
		for (auto p = !!_impl->_beginNewest ? _impl->_beginNewest->_next : _oldest; p; p = p->_next) {
			chip.access(size_t(p->start), p->count, p->amount, p->stride);
		}
		return chip.hitRate();
#endif
	}

	bool Interval::showL2CacheInterval(std::ostream & os, Life* life) {
#ifndef L2_CACHE_MODELING
		return true;
#else
		if (!_impl) return true;

		auto & pt = perThreads[omp_get_thread_num()];
		auto &_active = pt._active;
		auto &_oldest = pt._oldest;
		auto &_newest = pt._newest;
		CHiP chip;
		for (auto p = !!_impl->_beginNewest ? _impl->_beginNewest->_next : _oldest; p; p = p->_next) {
			chip.access(size_t(p->start), p->count, p->amount, p->stride);
		}

		if (chip.goodHitRate()) return true;

		chip.showSummary(os, life);

		// eliminate duplicates
		//
		struct PerSite {
			float  missPct;
			size_t duplicates;
			PerSite() : missPct(0.0), duplicates(0) {}
		};

		auto const path = life->path();
		auto const missPct = chip.missPct();

		static std::map<const std::string, PerSite> perSitePerFile[2];
		auto& perSite = perSitePerFile[(&os == &std::cout) ? 0 : 1];

		auto p = perSite.find(path);
		if (p == perSite.end()) {
			perSite.insert({path,PerSite()}); p = perSite.find(path);
		}
		auto & prev = p->second;

		if (prev.missPct >= missPct) {
			prev.duplicates++;
			return false;
		}
		if (prev.duplicates > 1) os << "...prev detailed cache hit behavior bettered or duplicated " << prev.duplicates - 1 << " times" << std::endl;
		prev.missPct = missPct;
		prev.duplicates = 1;

		// not a duplicate
		//
		chip.showDetails(os, life);	// doesn't show anything if good hit rate

		if (true) {
			os << "Traffic causing this" << std::endl;
			struct Data { SeqAcc* seqAcc; size_t count; Data(SeqAcc* seqAcc, size_t count) : seqAcc(seqAcc), count(count) {} };
			std::map<void*, Data> variables;
			for (auto p = !!_impl->_beginNewest ? _impl->_beginNewest->_next : _oldest; p; p = p->_next) {
				auto i = variables.find(p->start);
				if (i != variables.end()) i->second.count++; else variables.insert(std::make_pair(p->start, Data(p,1)));
			}
			for (auto & i : variables) {
				auto seqAcc = i.second.seqAcc;
				os << "    " << i.first << " " << seqAcc->nameBase << " " << i.second.count 
					<< " uses accessing " << seqAcc->count*seqAcc->amount << " bytes" 
					<< std::endl;
			}
		}
		return false;
#endif
	}
}
