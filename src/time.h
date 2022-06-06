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

#ifndef TIME_H_
#define TIME_H_

#include "./util.h"



// Slow inaccurate timing


void sleepInSecs(unsigned int seconds);
double dtime();

Microseconds timeInMicroseconds();

struct TimePoint {
    int line;
    const char* function;
    const char* file;
    double startTime = 0;
    double endTime = 0;
    
    void init(const char* _file,const char* _func,int _line){
        line = _line;
        file = _file;
        function = _func;
        tick();
    }
    
    void tick(){
        if(startTime == 0)
            endTime = startTime = timeInMicroseconds();
        else
            endTime = timeInMicroseconds();
    }
    
};

#if 0 // Yongbei

#define TIMEPOINT_INIT std::map<int,int> timePointIndex;std::vector<TimePoint> timePoints;int pre_timePoint = 0;

#define TIMEPOINT \
if(timePointIndex.find(__LINE__) != timePointIndex.end()) { timePoints[timePointIndex[__LINE__]].tick(); } \
else { TimePoint tp ;tp.init(__FILE__,__FUNCTION__,__LINE__);timePoints.push_back(tp);timePointIndex[__LINE__]=timePoints.size()-1; } \
std::cout<<"[ "<<timePoints[pre_timePoint].line<<" ---> "<<__LINE__<<" : "<<(timePoints[timePointIndex[__LINE__]].endTime-timePoints[pre_timePoint].endTime)<<" ] "; \
pre_timePoint = timePointIndex[__LINE__];

#define TIMEPOINT_FINA std::cout<<std::endl<<"line from --->  to  , time_costs(microseconds) "<<std::endl; \
for(int i = 0;i < timePoints.size()-1;i++) \
std::cout<<std::setw(6)<<" "<<timePoints[i].line<<" ---> "<<timePoints[i+1].line<<" ,  "<<(timePoints[i+1].endTime-timePoints[i].startTime)<<std::endl;

#else

#define TIMEPOINT_INIT
#define TIMEPOINT
#define TIMEPOINT_FINA

#endif

// Fast precise timing


class AccurateTimer {
public:
#ifdef _WIN32
	#define USE_QPC
    typedef __int64 Cycles;
#else
    typedef long long Cycles;
#endif
    
    AccurateTimer();
    void init();
    Cycles sinceInited();
	double countPerMicrosecond() { return _countPerMicrosecond; }    
    static double staticCountPerMicrosecond();
    
private:
    double        _countPerMicrosecond;

#ifdef _WIN32
    Cycles        _inited;
    
	Cycles qpc()
	#ifdef USE_QPC
		{   LARGE_INTEGER now; QueryPerformanceCounter(&now); return now.QuadPart; }
	#else
		;
	#endif

#else
    struct timespec        _inited;
    
    struct timespec cgt() {
        struct timespec curr_time = {0, 0};
        clock_gettime(CLOCK_MONOTONIC, &curr_time);
        return curr_time;
    }
#endif
};

// #define TUNING

#if !defined(TUNING)
	class Tuning_Scope_Para { 
	public:
		static const bool tuning = false;
		const bool doing;
		const char* const name;
		const int line;
		const int iter;
		const int tripCount;
		double startTime;
		Tuning_Scope_Para(const char* name, int line, int iter, long int tripCount) 
		  : doing(Tuning_Scope_Para::tuning && omp_get_thread_num()==0), 
			name(name), line(line), iter(iter), tripCount(tripCount) 
		{ 
			if (doing) init(); 
		}
		~Tuning_Scope_Para() {
			if (doing) fini();
		}
	private:
		void init();
		void fini();
	};
	#define TUNING_FLUSH
	#define TUNING_SCOPE_STEP(NAME)
	#define TUNING_SCOPE_STEP_BEGIN(NAME)
	#define TUNING_SCOPE_STEP_END
	#define TUNING_SCOPE_PARA_BEGIN(NAME,TRIPCOUNT)		//	Tuning_Scope_Para tuning_Scope_Para_##NAME(#NAME,__LINE__,iter,TRIPCOUNT);
	#define TUNING_SCOPE_PARA_END
	#define TUNING_SCOPE_ITER_BEGIN(NAME)
	#define TUNING_SCOPE_ITER_END
#else
	#define TUNING_FLUSH					Tuning::flush(); 
    #define TUNING_SCOPE_STR(NAME) #NAME
    #define TUNING_SCOPE_FUNCETC(NAME,SUFFIX) TUNING_SCOPE_STR(NAME##_##SUFFIX)
	#define TUNING_SCOPE_STEP(NAME)          Tuning::ScopeStep NAME##_tuningScopeStep(__FILE__, __LINE__, TUNING_SCOPE_STR(NAME));
	#define TUNING_SCOPE_STEP_BEGIN(NAME)   {Tuning::ScopeStep NAME##_tuningScopeStep(__FILE__, __LINE__, TUNING_SCOPE_STR(NAME));
	#define TUNING_SCOPE_STEP_END }
	#define TUNING_SCOPE_PARA_BEGIN(NAME,TRIPCOUNT)   {Tuning::ScopeParallel NAME##_tuningScopeParallel(NAME##_tuningScopeStep, __FILE__, __LINE__, TUNING_SCOPE_FUNCETC(NAME,para));
	#define TUNING_SCOPE_PARA_END			}
	#define TUNING_SCOPE_ITER_BEGIN(NAME)   {Tuning::ScopeForBody  NAME##_tuningScopeForBody(NAME##_tuningScopeParallel, __FILE__, __LINE__, TUNING_SCOPE_FUNCETC(NAME,iter));
	#define TUNING_SCOPE_ITER_END           }
	namespace Tuning {
		void flush();
		void setFTraceFilename(const char* name);
		class ScopeBase {
		public:
			const bool tracing;
			const char* const file; const int line;
			ScopeBase(const char* file, int line, const char* funcEtc);
			~ScopeBase();
			struct FTraceLine;
		protected:
			FTraceLine* BorE(bool makeB, const char* file, int line, const char* funcEtc);
			FTraceLine* ftraceLineB;
		};
		class ScopeStep : public ScopeBase {
		public:
			ScopeStep(const char* file, int line, const char* funcEtc);
			~ScopeStep();
		};
		class ScopeParallel : public ScopeBase {
		public:
			ScopeParallel(ScopeStep& scopeStep, const char* file, int line, const char* funcEtc);
			~ScopeParallel();
			void iteration(FTraceLine* b, FTraceLine* e);
		private:
			struct Iteration {
				FTraceLine* b;
				FTraceLine* e;
				Iteration() : b(NULL), e(NULL) {}
			};
			static const size_t iterationsCapacity = 10;
			omp_lock_t lock;
			struct IterationsRingBuffer {
				Iteration v[iterationsCapacity];
				size_t initialIterationCount;
				size_t iterationsHead, iterationsTail;
				IterationsRingBuffer() : initialIterationCount(0), iterationsTail(0), iterationsHead(1) {}
				void iteration(FTraceLine* b, FTraceLine* e);
			};
			std::map<int,IterationsRingBuffer*> threadToIterationsRingBuffer;
		};																						// end of macro
		class ScopeForBody : public ScopeBase {
			ScopeParallel& scopeParallel;
		public:
			ScopeForBody(ScopeParallel& scopeParallel, const char* file, int line, const char* funcEtc);
			~ScopeForBody();
		};																						// end of macro
	}
#endif

#endif
