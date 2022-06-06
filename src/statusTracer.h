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

#ifndef STATUSTRACER_H_
#define STATUSTRACER_H_

#include "util.h"

class StatusTracer {
protected:
	static double* mallocDoubleAligned(size_t size){
		return (double*)aMalloc(size*sizeof(double), 64);
	}
    static float* mallocFloatAligned(size_t size){
        return (float*)aMalloc(size*sizeof(float), 64);
    }
	static int* mallocIntAligned(size_t size){
		return (int*)aMalloc(size*sizeof(int), 64);
	}
	template<typename T>
	static void freeAligned(T* ptr){
		aFree(ptr);
	}

public:
	StatusTracer()  { }
	~StatusTracer() { clear(); }

    void clear(){ doublePtr.resize(0); intPtr.resize(0); floatPtr.resize(0);}

	void appendDoublePtr(double* data, int len, std::string what){
		doublePtr.push_back(DoublePtr({ data, len, what }));
	}
    void appendFloatPtr(float* data, int len, std::string what, bool compareFloatWithDouble = false){
        floatPtr.push_back(FloatPtr({data, len, what, compareFloatWithDouble}));
    }
	void appendIntPtr(int* data, int len, std::string what){
		intPtr.push_back(IntPtr({ data, len, what }));
	}

	void backupStatus(std::string status_fn);

	void recoveryStatus(std::string status_fn);
		// read the tracing data to _backup.model file

	void checkStatus(std::string guidance_fn);

private:
	struct DoublePtr { double* data; int len; std::string what; };
	std::vector<DoublePtr> doublePtr;
    struct FloatPtr {float* data;int len;std::string what;bool compareFloatWithDouble;};
    std::vector<FloatPtr> floatPtr;
	struct IntPtr { int* data; int len; std::string what; };
	std::vector<IntPtr> intPtr;

};


#endif
