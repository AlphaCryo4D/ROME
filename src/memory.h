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

#ifndef MEMORY_H_
#define MEMORY_H_

#include "./util.h"

#ifdef _WIN32 // Bevin is using windows for alloc data
#define _wrapped_aligned_malloc(s,a) _aligned_malloc(s,a)
#define _wrapped_aligned_free(p) _aligned_free(p)
#else
#define _wrapped_aligned_malloc(s,a) _mm_malloc(s,a)
#define _wrapped_aligned_free(p) _mm_free(p)
#endif

 // Heap support
//
template<typename T>
T* mallocCacheAlignedTemplate(size_t len, const char* file, size_t line) {
	auto size = sizeof(T)*len;
    return (T*)aMallocWkr(size, 64, file, int(line));
}
#define mallocCacheAligned(T,L) mallocCacheAlignedTemplate<T>(L, __FILE__, __LINE__)


namespace Heap {
    template <class ScalarType>
    ScalarType* allocScalars(size_t s, const char* file, size_t line) {
        return mallocCacheAlignedTemplate<ScalarType>(s, file, line);
    }
    template <class ScalarType>
    ScalarType* allocInitedScalars(size_t s, ScalarType v, const char* file, size_t line) {
        auto p = allocScalars<ScalarType>(s, file, line);
        for (size_t i = 0; i < s; i++) p[i] = v;
        return p;
    }
    template <class ScalarType>
    ScalarType* allocZeroedScalars(size_t s, const char* file, size_t line) {
        return allocInitedScalars<ScalarType>(s, ScalarType(0), file, line);
    }
    template <class ScalarType>
    void freeScalars(ScalarType* & p) {
        aFree(p); p = NULL;
    }
    static auto allocChars = allocScalars < char > ;
    static auto allocFloats = allocScalars < float > ;
    static auto allocDoubles = allocScalars < double > ;
    static auto allocZeroedDoubles = allocZeroedScalars < double > ;
    static auto freeDoubles = freeScalars < double > ;
    static auto freeChars = freeScalars < char > ;
    static auto freeFloats = freeScalars < float > ;
    static auto freeInts = freeScalars < int > ;
    double* allocDoublesForImages(const char* name, size_t nr_images, size_t s);
    double* allocDoublesForClasses(const char* name, size_t nr_classes, size_t s);
    double* allocZeroedDoublesForClasses(const char* name, size_t nr_classes, size_t s);
    void traceDoubles(double* ptr, const char* name, size_t sizeOfItem, size_t s);
    void traceFloats(float * ptr, const char* name, size_t sizeOfItem, size_t s);
    void untrace(void  * ptr);
};

#define mallocDoubles(len)			Heap::allocDoubles      (len,		__FILE__, __LINE__)
#define mallocFloats(len)			Heap::allocFloats       (len,		__FILE__, __LINE__)
#define mallocValDoubles(len,val)   Heap::allocInitedDoubles(len, val,	__FILE__, __LINE__)
#define mallocValFloats(len,val)    Heap::allocInitedFloats (len, val,	__FILE__, __LINE__)
#define mallocZeroedDoubles(len)	Heap::allocZeroedDoubles(len,		__FILE__, __LINE__)
#define mallocZeroedFloats(len)   	Heap::allocZeroedFloats (len,		__FILE__, __LINE__)


#ifdef USEMCDRAM
// high bandwidth memory
#include  <hbwmalloc.h>   // hbwmalloc interface
#endif
//
static void *
allocMemory(size_t __size, size_t __align,bool __inHBM, const char* file, int line)
{
    void *__mallocedMemory;
#ifdef USEMCDRAM
    if (__inHBM) int ret = hbw_posix_memalign((void**) &__mallocedMemory, 64, __size);
    else __mallocedMemory = (void*)aMallocWkr(__size, __align, file, line);
#else
    __mallocedMemory = (void*)aMallocWkr(__size, __align, file, line);
#endif
    return __mallocedMemory;
}
//
static void
freeMemory(void *__p,bool __inHBM = false)
{
#ifdef USEMCDRAM
    if (__inHBM) hbw_free(__p);
    else aFree(__p);
#else
    aFree(__p);
#endif
}


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

//#define L2_CACHE_MODELING
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


#endif /* defined(MEMORY_H_) */
