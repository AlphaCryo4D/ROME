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

#ifndef UTIL_H_
#define UTIL_H_

// Include all the external header files that are needed in header files here
// so that the precede the #include "./util_defs.h" calls
// in case they have new, delete, ... uses that we can not replace
//
// If an external header file is included into a .cpp, it should be done before the first of our header files is included
//
#include <assert.h>

#include <mkl.h>
#include <mkl_dfti.h>
#include <omp.h>
#ifndef __APPLE__
#include <malloc.h>
#endif
#include <memory>
#include <cassert>
#include <cmath>
#include <cstdlib>
#include <cstring>

// denormal number control
#include <pmmintrin.h>
#include <xmmintrin.h>

#include <exception>

#include <algorithm>
#include <atomic>
#include <complex>
#include <initializer_list>
#include <limits>
#include <map>
#include <random>
#include <string>
#include <tuple>
#include <vector>

#include <iostream>
#include <iomanip>
#include <sstream>
#include <fstream>

//--------------- header files for resmap
#include <array>
#include <functional>
#include <type_traits>
#include <numeric>
#include <iterator>
#include <list>
#include <set>
#include <regex>
#include <utility>
#include <cfloat>
#include <ios>
#include <thread>
//----------------------

#ifndef _CYGWIN
#ifdef __APPLE__
#include <limits.h>
#elif defined(_WIN32)
 //do nothing
#else
#include <values.h>
#endif
#endif


#if defined(_WIN32)
#include <Winsock2.h>
#include <Windows.h>
#undef min //winodws min/max conflicts with std::min/std::max
#undef max
#pragma comment(lib, "Ws2_32.lib")
#include <direct.h>
#else
#include <pthread.h>
#include <unistd.h>
#endif

// ntohl(),for swaping binary data order
#ifdef _WIN32
#else
#include <arpa/inet.h>
#endif


typedef double Microseconds;

#ifdef _WIN32
#include <time.h>
#undef GROUP_NAME
int gettimeofday(struct timeval *tp, void *tzp);
#elif __MACH__
#include <sys/time.h>
#define CLOCK_REALTIME 0
#define CLOCK_MONOTONIC 0
// clock_gettime is not implemented on OSX,accuracy is 1 microsecond
int clock_gettime(int /*clk_id*/, struct timespec* t);
#else
#include <time.h>
#include <unistd.h>
#include <sys/time.h>
#endif


#ifdef _mm_malloc
#undef _mm_malloc
#undef _mm_free
#endif

#undef mkl_alloc
#undef mkl_free
#define mkl_alloc use_aMalloc_instead_of_mkl_alloc
#define mkl_free  use_aFree_instead_of_mkl_free

#ifdef __APPLE__ 
typedef __SIZE_TYPE__ size_t;
#endif


// Unfortunately Windows and Linux take a different approach to macros in pragmas
// This hides the difference
//
#ifdef _WIN32
#define PragmaInMacro(unquoted_text) __pragma(unquoted_text)
#define NOINLINE __declspec(noinline)
#else
#define PragmaInMacro(unquoted_text) _Pragma(#unquoted_text)
#define NOINLINE __attribute__((noinline))
#endif

// Sets denormal results from floating-point calculations to zero,FTZ (flush-to-zero)
// Treats denormal values used as input to floating-point instructions as zero,DAZ (denormals-are-zero)
#define DISABLE_DENORMALS \
_MM_SET_DENORMALS_ZERO_MODE(_MM_DENORMALS_ZERO_ON); \
_MM_SET_FLUSH_ZERO_MODE(_MM_FLUSH_ZERO_ON);

#define ENABLE_DENORMALS \
_MM_SET_DENORMALS_ZERO_MODE(_MM_DENORMALS_ZERO_OFF); \
_MM_SET_FLUSH_ZERO_MODE(_MM_FLUSH_ZERO_OFF);

//===============================================================================================================================
// Some useful additions
//
#define EXIT_ABNORMALLY exitAbnormally(__FILE__,__LINE__)
void exitAbnormally(const char* file, int line);
	//
	// Calling exit(....) directly is hard to debug, so instead this should be used

static bool isVectorAligned(void const* p) { return (size_t(p) % 64 == 0); }

template <class T>
T square(T v) {
	return v*v;
}

template <class T>
int divRoundUp(T lhs, T rhs) { return (lhs + rhs - 1) / rhs; }

template <class XT, class YT>
YT interpolate(XT x, XT xLo, YT yLo, XT xHi, YT yHi) {
	return YT(double(x - xLo)/double(xHi - xLo)*double(yHi) + double(xHi - x)/double(xHi - xLo)*double(yLo));
}

struct DoublePair {
	double lo;
	double hi;
	DoublePair() : lo(0), hi(0) {}
	DoublePair(double lo, double hi) : lo(lo), hi(hi) {}
	DoublePair(std::string);
	std::string string() const;
	bool operator==(DoublePair const & rhs) const { return lo == rhs.lo && hi == rhs.hi; }
	bool operator!=(DoublePair const & rhs) const { return !(*this == rhs); }
	double mid() const { return (lo + hi) / 2.0; }
	static DoublePair midFrac(double mid, double frac) { auto diff = std::abs(frac*mid); return DoublePair(mid - diff, mid + diff); }
};


class DoSomePerIter {
public:
	int count, limit, currIter;
	DoSomePerIter() : currIter(-1) {}
	template <typename Callable>
	void note(const bool doIt, int iter, Callable callable)
	{
		if (doIt && omp_get_thread_num() == 0)
#pragma omp critical
		{
			if (currIter != iter) { currIter = iter; count = 0; limit = 1; }
			if (count++ >= limit) {
				if (limit < 4) limit++; else limit *= 2;
				callable(count);
			}
		}
	}
};


template <class T>
struct InitZeroSingleton {
	T v;
	InitZeroSingleton() : v(0) {}
};


class PerformanceCounter {
public:
	PerformanceCounter();
	virtual ~PerformanceCounter();
	virtual void show(std::ostream & os) = 0;
	virtual void reset() = 0;
	void showAndReset(std::ostream & os) { show(os); reset(); }
	static void showAll(std::ostream & os, bool andReset = false);
private:
	PerformanceCounter* _prev;
	PerformanceCounter* _next;
};

class IntPerformanceCounter : public PerformanceCounter {
	const char* const name;
public:
	IntPerformanceCounter(const char* const name) : name(name) {}
	InitZeroSingleton<std::atomic<int>> count;
	virtual void show(std::ostream & os);
	virtual void reset();
};

inline int omp_get_max_threads_or_one_if_no_openmp() {
#if defined(_OPENMP)
    int num_threads = omp_get_max_threads();
#else
    int num_threads = 1;
#endif
    return num_threads;
}
//===============================================================================================================================
// By using this as a member or as a base class, the code will get an error message
// if an attempt is made to copy the object.  
//
// Use when copying doesn't make sense, because its use causes obscure bugs to become compile-time errors.
//
#include "resmap_util_nocopy.h"

//===============================================================================================================================
// Parallelism control
//
int initial_omp_get_max_threads();
#define omp_get_max_threads initial_omp_get_max_threads
void setGoSerial  (std::string label);
void setGoParallel(std::string label);
void maybeGoSerial(const char* label);
class Lock : NoCopy {
	char const * volatile _heldByFile;	// so a debugger can examine them
	int			 volatile _heldByLine;
	omp_lock_t			  _lock;
public:
	Lock() : _heldByFile(nullptr), _heldByLine(-1) { omp_init_lock(&_lock);  }
	~Lock() {
		if (_heldByFile) EXIT_ABNORMALLY;
		omp_destroy_lock(&_lock);
	}
	void acquire(char const * file, int line) {
		omp_set_lock(&_lock);
		_heldByFile = file;
		_heldByLine = line;
	}
	void release() {
		_heldByLine = -1;
		_heldByFile = nullptr;
		omp_unset_lock(&_lock);
	}
};
#if defined(_WIN32)
static void yield(bool andSleep) { if (andSleep) Sleep(1); else SwitchToThread(); }
#else
static void yield(bool andSleep) { if (andSleep) sleep(1); else sched_yield(); }
#endif


class ScopedAcquire {
public:
	Lock & lock;
	ScopedAcquire(Lock & lock, char const * file, int line) : lock(lock) { lock.acquire(file, line); }
	~ScopedAcquire() { lock.release(); }
};

//===============================================================================================================================
// Heap allocation and census support
//		It is easy to leak memory, and hard to keep track of how much is allocated from where in the code
//		The code notes all its allocations and frees to track at least the memory we allocate

void noteAllocatedWkr(void* p, size_t s, bool isVector, const char* file, int line);
void noteFreedWkr    (void* p,           bool isVector);

template <typename T> 
T* allocated(T* p, size_t s, bool isVector, const char* file, int line) {
	noteAllocatedWkr(p, s, isVector, file, line);
	return p;
}
template<typename T>
T* allocated(T* p, const char* file, size_t line) {
	noteAllocatedWkr(p, sizeof(*p), false, file, line);
	return p;
}


// Census support
//
class CensusDetailsProvider : public NoCopy {
public:
	CensusDetailsProvider();
	virtual ~CensusDetailsProvider();
	virtual void showPopulation(std::ostream & os) = 0;
	CensusDetailsProvider* next() const { return _next; }
private:
	CensusDetailsProvider* _next;
};


size_t census();
void showPopulation();
	// memory.cpp


//		To stop accidental uses of new and delete and other mechanisms, those keywords and functions are #defined to give an error
//		This causes a problem with external header files that use them, so most external header files are included above
//		before the macros are defined - alternatively the #include could be bracketed with util_heap_(un)defs.h inclusions

#include "./resmap_util_heap_defs.h"
	// This defines a lot of the likely ways of allocating heap storage as macros saying not to use them

#include "./resmap_util_heap_undefs.h"
	// This undoes those defines.  These two includes can bracket a region where one of those ways should be allowed
	// like this one, where we are defining the templates that allow a syntax that doesn't need the macros

// The following convenience mechanisms avoid most of the need for these #includes

// Bytes allocator and deallocator
//
void* aMallocWkr(size_t size, size_t alignment, const char* file, int line);
void  aFreeWkr(void* p);

#define aMalloc(S,A) aMallocWkr(S,A,__FILE__,__LINE__)					  //	  E.g.		  aMalloc(sizeof(double)*100,64)


template<typename T>
void aFree(T* volatile &p) {
	if (!p) return;
	aFreeWkr(p); 
	p = nullptr; 
}
template <typename T> 
void aFreeConst(T* const p) {
	if (!p) return;
	aFreeWkr(p);
}


// Singleton allocator and deallocator - no constructor parameters
//	If you need parameters, use the following style
//																			//		E.g.		#include "./util_heap_undefs.h"
//																			//					sNewA(Class,(arg1,arg2))
//																			//					#include "./util_heap_defs.h"
//	It is best to put in a static member function of the class
//
template <typename T> 
T* sNewTemplate(const char* file, int line) {
	return allocated(new T, sizeof(T), false, file, line);
}
#define sNew(T) sNewTemplate<T>(__FILE__,__LINE__)							//		E.g.		sNew(Class)
#define sNewA(T,A) allocated(new T A, sizeof(T), false, __FILE__, __LINE__) //		E.g.		sNewA(Class,(arg1,arg2))
template <typename T> 
void sDelete(T* volatile &p) {
	if (!p) return;
	noteFreedWkr((void*)p, false);
	delete p; 
	p = nullptr; 
}
template <typename T> 
void sDeleteConst(T* const p) {
	if (!p) return;
	noteFreedWkr((void*)p, false);
	delete p; 
}
#define born(p) allocated(p,__FILE__,__LINE__)
template <typename T> 
T* died(T* p) {
	if (p) noteFreedWkr((void*)p, false);
	return p;
}

// Vector allocator and deallocator - no constructor parameters
//
template <typename T, typename CT>
T* vNewTemplate(CT c, const char* file, int line) {
	return allocated(new T[c], c*sizeof(T), true, file, line);
}
#define vNew(T,C) vNewTemplate<T>(C, __FILE__, __LINE__)					//		E.g.		vNew(char,500)

template <typename T> 
void vDelete(T* volatile &p) {
	if (!p) return;
	noteFreedWkr((void*)p, true);
	delete [] p;
	p = NULL; 
}
template <typename T> 
void vFreeConst(T* const p) {
	if (!p) return;
	noteFreedWkr((void*)p, true);
	delete [] p;
}


// By default the above must be used
//
#include "./resmap_util_heap_defs.h"



//===============================================================================================================================
// Convenient vector operations
//
template <class _Elt>
class VectorOfNoCopyable {
public:
	typedef _Elt Elt;
	VectorOfNoCopyable(size_t size) : _size(size), _elts(vNew(Elt,size)) {}
	~VectorOfNoCopyable() { vDelete(_elts); }
	size_t size() const { return _size; }
	Elt       & operator[](size_t index)       { return _elts[check(index)]; }
	Elt const & operator[](size_t index) const { return _elts[check(index)]; }
private:
	size_t _size;
	Elt*   _elts;
	size_t check(size_t index) const {
		assert(index < _size);
		return index;
	}
};


template<typename T>
class Checksum {
	std::vector<T> copy;
	int size;
	int capacity;
	T x[3];
	double hash;
public:
	Checksum() : capacity(0), size(0), copy(), hash(0.0) {
		for (int i = 0; i < 3; i++) x[i] = 0.0; 
	}
	~Checksum() { fini(); }
	void fini() {
		copy.resize(0);
		capacity = size = 0;
	}
	void init(T* vec, int length) {
		if (length > capacity) {
			fini();
			capacity = length * 1.5;
			copy.resize(capacity);
		}
		x[0] = vec[length * 1 / 5];
		x[1] = vec[length * 3 / 5];
		x[2] = vec[length * 4 / 5];
		hash = 0.0;
		for (int i = 0; i < length; i++) {
			hash += std::abs(vec[i]+1.0);
			copy[i] = vec[i];
		}
		size = length;
	}
	bool operator==(Checksum const & rhs) {
		return 
			x[0] == rhs.x[0] &&
			x[1] == rhs.x[1] &&
			x[2] == rhs.x[2] &&
			hash == rhs.hash;
	}
	void print(std::ostream & os) const {
		os << "Checksum for " << std::endl;
		for (int i = 0; i < size; ) {
			os << i << ":" << copy[i] << std::endl;
			if (i < 10  || size < i +  10) i+= 1; else
			if (i < 100 || size < i + 100) i+=10; else
										   i+=100;
		}
	}
	void printDiff(std::ostream & os, Checksum const & rhs) const {
		os << "Checksum for " << std::endl;
		int count(0);
		for (int i = 0; i < size; i++) {
			if (copy[i] == rhs.copy[i]) continue;
			os << i << ":" << copy[i] << " != " << rhs.copy[i] << std::endl;
			if (count++ == 30) i = std::max(i,size - 30);
		}
	}
};


class MeanAndStdDeviation {
	size_t _n;
	double _sum, _sumOfSquares;
	bool   _resultComputed;
	float  _mean, _stdDev;
public:
	MeanAndStdDeviation() { init(); }
	void init() { _n = 0; _sum = _sumOfSquares = 0.0; _resultComputed = false; }
	float mean  () { compute(); return _mean; }
	float stdDev() { compute(); return _stdDev; }
	void insert(float f) {
		_resultComputed = false;
		_n++;
		_sum += f;						// double adjusts to a certain extent for the numerical instability of these sums
		_sumOfSquares += square(f);
	}
private:
	void compute();
};


template<typename T>
static T sumvec(T *vec, int length)
{
	T sum = 0;
	for (int i = 0; i < length; i++) {
		sum += fabs(vec[i]);
	}
	return sum;
}

template<typename T>
static T sumvec(T *real_vec, T *imag_vec, int length)
{
    T sum = 0;
    for (int i = 0; i < length; i++) {
        sum += fabs(real_vec[i]);
        sum += fabs(imag_vec[i]);
    }
    return sum;
}

template<typename T>
static void printvec(std::ostream& out, T* vec, int length)
{
	for (int i = 0; i < std::min(100, length); i++)
		out << vec[i] << " ";
	out << std::endl;
};


// Some convenient 2D operations
//
template <class Coeff>
class Raster {
public:
	typedef size_t V;
	Raster(V w, V h) : _width(w), _height(h) {}
	virtual ~Raster() {}
	const V _width;
	const V _height;
	virtual Coeff operator()(V w, V h) = 0;
};
typedef Raster<float> FloatRaster;



// Much faster than a vector<char> or string
//
struct FastCharBuf {
	FastCharBuf(size_t initialCapacity)
		: capacity(std::max(size_t(64), initialCapacity)),
		ptr(vNew(char,capacity)),
		size(0)
	{}
	~FastCharBuf() { vDelete(ptr); }

	size_t capacity;
	char*  ptr;
	size_t size;

	void push_back(char c) {
		if (size == capacity) grow();
		ptr[size++] = c;
	}

	void grow() {
		capacity += capacity / 2;
		char* p = vNew(char, capacity);
		memcpy(p, ptr, size);
		vDelete(ptr);
		ptr = p;
	}
};


// Floating point comparisons for equality are too risky in many cases
// so this provides a slightly less code-generation-sensitive test
//
template <typename T>
bool nearEnoughTemplate(T const & tentative, T const & known, bool info = false) {
	if (info) { std::cerr << "nearEnoughTemplate<T> " << tentative << " != " << known << std::endl; }
	return tentative == known;
}
template <> bool nearEnoughTemplate(double const & tentative, double const & known, bool info);
template <> bool nearEnoughTemplate(float const & tentative, float const & known, bool info);


//===============================================================================================================================
// Some convenience io support
//
class ifstreamCheckingExistence : public std::ifstream {
public:
	ifstreamCheckingExistence(const char* fileName);
	static ifstreamCheckingExistence* make(std::string const & fileName) {
		return ifstreamCheckingExistence::make(fileName.c_str());
	}
	static ifstreamCheckingExistence* make(const char* fileName) {
#include "./resmap_util_heap_undefs.h"
		return new ifstreamCheckingExistence(fileName);
#include "./resmap_util_heap_defs.h"
	}
};

class ofstreamCheckingCreated : public std::ofstream {
public:
	ofstreamCheckingCreated(const char* fileName);
	static ofstreamCheckingCreated* make(std::string const & fileName) {
		return ofstreamCheckingCreated::make(fileName.c_str());
	}
	static ofstreamCheckingCreated* make(const char* fileName) {
#include "./resmap_util_heap_undefs.h"
		return new ofstreamCheckingCreated(fileName);
#include "./resmap_util_heap_defs.h"
	}
};


static std::string currentDirectory() {
	char buffer[FILENAME_MAX];
	if (!
#if defined(_WIN32)
		_getcwd
#else
		getcwd
#endif
		(buffer, FILENAME_MAX)) {
		return "<currentDirectory() failed>";
	}
	return buffer;
}


#endif

