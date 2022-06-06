/***************************************************************************
 *
 * Dana-Farber Cancer Institute, Harvard Medical School and Peking University
 * Intel® Parallel Computing Center for Structural Biology
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
#ifndef _L2CACHEMODEL_H_
#define _L2CACHEMODEL_H_

#include "checker.h"

#include "../src/resmap/resmap_util.h"

//#define L2_CACHE_MODELING

// L2 cache modelling is used to try to understand our L2 miss rates in production runs without distort performance
//
namespace L2CacheModel {
    
#define L2_ACCESS_ELTS \
ELT(const char* , nameBase	) SEP \
ELT(int			, nameIndex	) SEP \
ELT(void*		, start		) SEP \
ELT(size_t		, amount	) SEP \
ELT(size_t		, count		) SEP \
ELT(size_t		, stride	) \
// end of macro
    
    class Interval {
        class Impl;
        Impl* _impl;
    public:
        Interval() : _impl(nullptr) {}
        ~Interval() { if (_impl) end(); }
        void begin();
        void end() { if (_impl) endWkr(); }
        size_t accessesRecordedCount();
        size_t cacheProbeCount();
        float hitRate();
        bool showL2CacheIntervalLocked(Life* life) {
            bool result;
#pragma omp critical
            result = showL2CacheIntervalUnlocked(life);
            return result;
        }
        bool showL2CacheIntervalUnlocked(Life* life);
        // returns true if a good hit rate
    private:
        void endWkr();
        bool showL2CacheInterval(std::ostream & os, Life* life);
    };
    

#ifndef L2_CACHE_MODELING
    static
#endif
    void seqAccWkr(
#define SEP			,
#define ELT(T,N)	T N
                   L2_ACCESS_ELTS
#undef ELT
#undef SEP
    )
#ifdef L2_CACHE_MODELING
    ;
#else
    {}
#endif
    
    template <typename Elt>
    void seqAcc(const char* nameBase, int nameIndex, Elt* start, size_t count = 1, size_t stride = 1) {
        seqAccWkr(nameBase, nameIndex, (void*)start, sizeof(Elt), count, sizeof(Elt));
    }
};

#endif