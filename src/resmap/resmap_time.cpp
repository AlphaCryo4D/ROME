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

#include "resmap_time.h"

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

#else // elif not defined _WIN32

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
#endif //end if node define USE_QPC

#endif
