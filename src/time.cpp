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

#include "util.h"		// used for building precompiled headers on Windows

#include "time.h"

#ifdef _WIN32

	static int gettimeofday(struct timeval *tp, void *tzp)
	{
	    time_t clock;
	    struct tm tm;
	    SYSTEMTIME wtm;
	    GetLocalTime(&wtm);
	    tm.tm_year     = wtm.wYear - 1900;
	    tm.tm_mon      = wtm.wMonth - 1;
	    tm.tm_mday     = wtm.wDay;
	    tm.tm_hour     = wtm.wHour;
	    tm.tm_min      = wtm.wMinute;
	    tm.tm_sec      = wtm.wSecond;
	    tm. tm_isdst   = -1;
	    clock = mktime(&tm);
	    tp->tv_sec = long(clock);
	    tp->tv_usec = wtm.wMilliseconds * 1000;
	    return 0;
	}

	void sleepInSecs(unsigned int seconds) { Sleep(DWORD(seconds)*1000); } 
#else
	void sleepInSecs(unsigned int seconds) { sleep(seconds); } 
#endif

#if __MACH__
// clock_gettime is not implemented on OSX,accuracy is 1 microsecond
int clock_gettime(int /*clk_id*/, struct timespec* t) {
    struct timeval now;
    int rv = gettimeofday(&now, NULL);
    if (rv) return rv;
    t->tv_sec  = now.tv_sec;
    t->tv_nsec = now.tv_usec * 1000;
    return 0;
}
#endif


double dtime(){
	static bool baseDtimeSet = false;
	static double baseDtime;
    double tseconds=0.0;
    struct timeval mytime;
    gettimeofday(&mytime,(struct timezone*)0);
    tseconds=(double)(mytime.tv_sec+mytime.tv_usec*1.0e-6);
	if (!baseDtimeSet) { baseDtimeSet = true; baseDtime = tseconds; }
    return tseconds - baseDtime;
}

// The constructor is not called on OXS if define here
// static AccurateTimer microsecondsTimer;

static AccurateTimer microsecondsTimer;
Microseconds timeInMicroseconds() {
	return double(microsecondsTimer.sinceInited()) / microsecondsTimer.countPerMicrosecond();
}


// Will move out once working on Linux

AccurateTimer::AccurateTimer() {
    _countPerMicrosecond = staticCountPerMicrosecond();
    init();
}

#ifndef _WIN32
void AccurateTimer::init() {
	    _inited = cgt();
	}
	AccurateTimer::Cycles AccurateTimer::sinceInited() {
	    struct timespec current = cgt();
	    return (current.tv_sec-_inited.tv_sec)*1000000000+(current.tv_nsec-_inited.tv_nsec);
	}
	double AccurateTimer::staticCountPerMicrosecond() {
	    return 1e3;
	}
#else
	void AccurateTimer::init() {
    	_inited = qpc();
	}

	AccurateTimer::Cycles AccurateTimer::sinceInited() {
    	return qpc() - _inited;

	}

	double AccurateTimer::staticCountPerMicrosecond() {
		#ifdef USE_QPC
    		LARGE_INTEGER hertz;
    		QueryPerformanceFrequency(&hertz);
    		return double(hertz.QuadPart) / 1000000.0;
		#else
    		return 1e3;
		#endif
	}


	#ifndef USE_QPC
		AccurateTimer::Cycles AccurateTimer::qpc() {
			Cycles tcycles;
			#pragma omp critical
		    {
				double tseconds = dtime();									// in seconds, accurate to milliseconds

				static bool  laterCall = false;
				static double baseTseconds;
		    	if (laterCall) {
		    		tseconds    -= baseTseconds;
		    	} else {
		    		laterCall    = true;
		    		baseTseconds = tseconds;								// start the tseconds at 0 so the 1e-9 below is significant
		    		tseconds     = 0.0;
		    	}

				tcycles = Cycles(tseconds*1e9);								// assume a gigahertz clock, see staticCountPerMicrosecond

				static Cycles minNextTcycles = 0;
				if (minNextTcycles >= tcycles) tcycles = minNextTcycles;
				minNextTcycles = tcycles + 1000;							// smallest interval is 1000 cycles
			}
			return tcycles;
		}
	#endif





#endif


static const int topCapacity = 10;
static int topSize = 1;
int topIter(-1);
double topElapsed[topCapacity+1];		// initializes to 0.0's, support writing one beyond "end"
int    topLine   [topCapacity+1];
bool   topBad    [topCapacity+1];

void Tuning_Scope_Para::init() {
	startTime = dtime();
}

void Tuning_Scope_Para::fini() {
	const double elapsed = dtime() - startTime;
	const bool   enough = (tripCount/5 >= 60);
	if (topIter != iter) { 
		if (topIter >= 0) {
			std::cerr << "TUNING_SCOPE_PARA topElapsed:";
			for (int i = 0; i < topSize; i++) std::cerr << topElapsed[i] << "@" << topLine[i] << (topBad[i]?"_BAD " : " ");
			std::cerr << std::endl;
		}
		topIter = iter; 
		topSize = 0; 
	}
	// Look for a matching one to merge with
	int i = topSize; 
	for (int j = topSize-1; j >= 0; j--) {
		if (topLine   [j] != line   ) continue;
		if (topBad    [j] != !enough) continue;
		// Matches
		if (topElapsed[j] > elapsed ) return;		
			// it already is the largest.  If should be using the smallest, the replacement below is wrong.
		// Replace it by using it as the empty slot
		i = j;
		break;
	}
	for (; i > 0; i--) {
		double next = topElapsed[i-1];
		if (next >= elapsed) break;
		topElapsed[i] = next;
		topLine   [i] = topLine[i-1];
		topBad    [i] = topBad [i-1];
	}
	topElapsed[i] = elapsed;
	topLine   [i] = line;
	topBad    [i] = !enough;
	if (topSize < topCapacity) topSize++;
	if (i == topCapacity) return;
	if (enough) return;
	std::cerr << "TUNING_SCOPE_PARA iter:" << iter << " " << name << "@" << line << " tripCount:" << tripCount
		<< " elapsed:" << elapsed
		<< " on " << omp_get_max_threads() << " threads. "
		<< (!enough ? " <<<<<<<<<<<<<< NOT ENOUGH FOR 60 CORES" : "")
		<< std::endl;
}

#if defined(TUNING)

namespace Tuning {

	std::ostream* ftrace_os = NULL;
	AccurateTimer* ftrace_accurateTimer;

#ifdef assert
#undef assert
#endif

#define ASSERT(C) assert((C), #C)
	static void assert(bool c, const char* text) {
		if (c) return;
		std::cout << "assert(" << text << ") failed" << std::endl;
		std::cerr << "assert(" << text << ") failed" << std::endl;
		if (!!ftrace_os) *ftrace_os << "assert(" << text << ") failed" << std::endl;
		*(int*)(-1) = 0;
	}
	void setFTraceFilename(const char* name) {
		ftrace_os = new std::ofstream(name);
		(*ftrace_os) << "#           TASK-PID   CPU#      TIMESTAMP  FUNCTION" << std::endl;
		(*ftrace_os) << "#              | |       |          |         |" << std::endl;
		ftrace_accurateTimer = new AccurateTimer();
	}               //  "                -0001  [000] 123456.123456: %s: %s%s%n"

	typedef unsigned __int64 U64;
	void digits(U64 number, int numberOfDigits) {
		if (numberOfDigits > 1) digits(number / 10, numberOfDigits - 1);
		(*ftrace_os) << char('0' + number % 10);
	}

	class ScopeBase::FTraceLine {
	private:
		~FTraceLine() { *(int*)(-1) = 0; }	// do not call!
	protected:
		friend class FTraceLinePool;
		friend class FTraceLines;
		FTraceLine() : prev(NULL), next(NULL) {}

		void release() {
			ASSERT(!bOrEOrNULL);
			ASSERT(!next);
			ASSERT(!prev);
		}
	public:
		FTraceLine *prev, *next;
		U64  nanoTime;
		int  tid;
		bool isB;
		FTraceLine* bOrEOrNULL;
		const char* file; int line; const char* funcEtc;
		void disconnectBorE() {
			bOrEOrNULL->bOrEOrNULL = NULL;
			bOrEOrNULL = NULL;
		}
		static void elide(
			FTraceLine* oldPairE,
			FTraceLine* newPairB) {
			auto oldPairB = oldPairE->bOrEOrNULL;
			auto newPairE = newPairB->bOrEOrNULL;
			oldPairB->bOrEOrNULL = newPairE;
			newPairE->bOrEOrNULL = oldPairB;
			oldPairE->bOrEOrNULL = NULL;
			newPairB->bOrEOrNULL = NULL;
		}
		void init(
			FTraceLine* bOrNULL,
			const char* file, int line, const char* funcEtc) {
			this->isB = (!bOrNULL);
			this->bOrEOrNULL = bOrNULL;
			this->file = file;
			this->line = line;
			this->funcEtc = funcEtc;
			this->next = NULL;
			this->prev = NULL;
			if (bOrNULL) {
				ASSERT(!bOrNULL->bOrEOrNULL);
				bOrNULL->bOrEOrNULL = this;
			}
		}
		void write() {
			auto & os = *ftrace_os;
			U64 seconds = nanoTime / 1000000000LL;
			U64 fraction = nanoTime % 1000000000LL / 1000;
#pragma omp critical
			{
				os << "                -";
				digits(tid, 4);      // tid that wrote the trace record
				os << "  [";
				digits(tid, 3);
				os << "] ";
				digits(seconds, 6);
				os << ".";
				digits(fraction, 6);
				os << ": tracing_mark_write: "
					<< (isB ? "B|" : "E|")
					<< tid          // tid
					<< "|";
				if (isB) {
					os << funcEtc << "@" << line;
				} else if (bOrEOrNULL) {
					os << bOrEOrNULL->funcEtc << "@" << bOrEOrNULL->line;
				}
				os << std::endl;
			}
		}
	};

	class FTraceLinePool {
	public:
		FTraceLinePool() : available(omp_get_max_threads()) {}
		ScopeBase::FTraceLine* acquire() {
			auto& ptr = available[omp_get_thread_num()];
			auto result = ptr;
			if (!ptr) {
				static const int count = 10000;
				auto results = new ScopeBase::FTraceLine[count];	// Make more than one to cut down on malloc locking
				for (int i = 0; i < count; i++) {
					results[i].next = ptr;
					ptr = &results[i];
				}
			}
			result = ptr;
			ptr = result->next;
			return result;
		}
		void release(ScopeBase::FTraceLine* ftl) {
			ftl->release();
			auto& ptr = available[omp_get_thread_num()];
			ftl->next = ptr;
			ptr = ftl;
		}
	private:
		std::vector<ScopeBase::FTraceLine*> available;
	} fTraceLinePool;

	struct FTraceLines : private omp_lock_t {
		FTraceLines() {
			omp_init_lock(this);
			_head = NULL;
			_tail = NULL;
		}
		~FTraceLines() {
			flush();
		}
		void flush() {
			while (auto ftl = removeHead()) {
				ftl->write();
				ASSERT(ftl->bOrEOrNULL);
				if (!ftl->isB) {
					ftl->disconnectBorE();
					if (ftl->bOrEOrNULL) fTraceLinePool.release(ftl->bOrEOrNULL);
					fTraceLinePool.release(ftl);
				}
			}
		}
		ScopeBase::FTraceLine* removeHead() {
			acquire();
			auto ftl = _head;
			if (!!ftl) removeNoLock(ftl);
			release();
			return ftl;
		}
		void insert(ScopeBase::FTraceLine * ftl) {
			acquire();
			ftl->nanoTime =
				U64(
				double(ftrace_accurateTimer->sinceInited()) /
				(ftrace_accurateTimer->countPerMicrosecond()*1e-3));
			ftl->tid = omp_get_thread_num();
			ftl->prev = _tail;
			ftl->next = NULL;
			if (_tail) _tail->next = ftl; else _head = ftl;
			_tail = ftl;
			release();
		}
		void remove(ScopeBase::FTraceLine * ftl) {
			acquire();
			removeNoLock(ftl);
			release();
		}
	private:
		ScopeBase::FTraceLine *_head, *_tail;
		void acquire() { omp_set_lock(this); }
		void release() { omp_unset_lock(this); }
		void removeNoLock(ScopeBase::FTraceLine * ftl) {
			if (ftl->prev) ftl->prev->next = ftl->next; else _head = ftl->next;
			if (ftl->next) ftl->next->prev = ftl->prev; else _tail = ftl->prev;
			ftl->next = ftl->prev = NULL;
		}
	} ftraceLines;

	// TBD - buffer these until a non-critical time before spending all the time writing them
	//
	static bool suppressed = false;
	ScopeBase::ScopeBase(const char* file, int line, const char* funcEtc)
		: file(file), line(line), ftraceLineB(NULL), tracing(!suppressed && !!ftrace_os && omp_get_thread_num()<6) {		// BEVIN 6
		ftraceLineB = BorE(true, file, line, funcEtc);
	}
	ScopeBase::~ScopeBase() {
	}
	ScopeBase::FTraceLine* ScopeBase::BorE(bool makeB, const char* file, int line, const char* funcEtc) {
		if (!tracing) return NULL;
		ASSERT(makeB ^ !!ftraceLineB);
		auto ftl = fTraceLinePool.acquire();
		ftl->init(makeB ? NULL : ftraceLineB, file, line, funcEtc);
		ftraceLines.insert(ftl);
		return ftl;
	}

	ScopeStep::ScopeStep(const char* file, int line, const char* funcEtc)
		: ScopeBase(file, line, funcEtc) {
	}
	ScopeStep::~ScopeStep() {
		BorE(false, file, line, "");
	}

	ScopeParallel::ScopeParallel(
		ScopeStep& scopeStep, const char* file, int line, const char* funcEtc)
		:
		ScopeBase(file, line, funcEtc) {
		omp_init_lock(&lock);
	}
	ScopeParallel::~ScopeParallel() {
		BorE(false, file, line, "");
		for (auto i = threadToIterationsRingBuffer.begin();
			i != threadToIterationsRingBuffer.end();
			i++) {
			delete i->second;
		}
		threadToIterationsRingBuffer.clear();
		omp_destroy_lock(&lock);
	}
	void ScopeParallel::iteration(FTraceLine* b, FTraceLine* e) {
		if (!b || !e) return;
		omp_set_lock(&lock);
		IterationsRingBuffer* iterationsRingBuffer = NULL;
		{
			auto tid = omp_get_thread_num();
			auto i = threadToIterationsRingBuffer.find(tid);
			if (i != threadToIterationsRingBuffer.end()) {
				iterationsRingBuffer = i->second;
			} else {
				iterationsRingBuffer = new IterationsRingBuffer;
				threadToIterationsRingBuffer.insert(std::make_pair(tid, iterationsRingBuffer));
			}
		}
		iterationsRingBuffer->iteration(b, e);
		omp_unset_lock(&lock);
	}
	void ScopeParallel::IterationsRingBuffer::iteration(FTraceLine* b, FTraceLine* e) {

		// There are a LOT of these, so we only report the first and last few
		// by having a shared queue in the ScopeParallel and removing middle items from the queue.
		// This is done by keeping the first few, then putting the remainder into a ring buffer,
		// removing the oldest from both the ring buffer and the queue when the space is needed
		// for a later one.  The dropped iterations are represented by a single combined one,
		// using the tail iteration for this purpose.

		// Keep the first ones
		if (initialIterationCount < iterationsCapacity) {
			initialIterationCount++;
			return;
		}

		// Make space to store these if the buffer is full
		// by combining the last oldest two entries
		if (iterationsHead == iterationsTail) {
			auto & oldest = v[iterationsTail];
			iterationsTail = (iterationsTail + 1) % iterationsCapacity;
			auto & merge = v[iterationsTail];
			if (!!oldest.b) {       // the ringBuffer starts with an empty element in it

				if (false) {
					// just remove the oldest
					oldest.b->disconnectBorE();
					ftraceLines.remove(oldest.b); fTraceLinePool.release(oldest.b);
					ftraceLines.remove(oldest.e); fTraceLinePool.release(oldest.e);
				} else {
					// merge
					FTraceLine::elide(oldest.e, merge.b);
					ftraceLines.remove(oldest.e); fTraceLinePool.release(oldest.e);
					ftraceLines.remove(merge.b);  fTraceLinePool.release(merge.b);
					merge.b = oldest.b;
					merge.b->funcEtc = "repeated";
				}
				oldest.b = NULL;
				oldest.e = NULL;
			}
		}

		auto & head = v[iterationsHead];
		iterationsHead = (iterationsHead + 1) % iterationsCapacity;
		head.b = b;
		head.e = e;
	}

	// These are when the code spreads the iterations across the threads
	//
	ScopeForBody::ScopeForBody(ScopeParallel& scopeParallel, const char* file, int line, const char* funcEtc) :
		ScopeBase(file, line, funcEtc), scopeParallel(scopeParallel) {
	}

	ScopeForBody::~ScopeForBody() {
		scopeParallel.iteration(ftraceLineB, BorE(false, file, line, ""));
	};

	static int ftrace_debug_count = 0;
	void flush() {
		ftraceLines.flush();
		if (!ftrace_accurateTimer) return;
		auto elapsed = ftrace_accurateTimer->sinceInited()/ftrace_accurateTimer->countPerMicrosecond()/1e6;
		static double nextTransition = 15.0;	// collect these first seconds		// BEVIN
		static double nextSkip		 = 15.0;	// then skip this amount			// BEVIN
		if (elapsed < nextTransition) {
			if (ftrace_debug_count++ < 10) std::cerr << "Tuning::flush elapsed=" << elapsed << " < nextTransition:" << nextTransition << std::endl;
			return;
		}
		std::cerr << "Tuning::flush elapsed=" << elapsed << " >= nextTransition:" << nextTransition << (suppressed?"was suppressed":"") << std::endl;
		if (suppressed) {
			suppressed = false;
			nextTransition = elapsed + 15.0;						// collect another 15 seconds
		} else {
			suppressed = true;
			nextTransition = elapsed + nextSkip; nextSkip *= 2.0;	// skip this amount, twice as much next time
		}
	}
}

#endif
