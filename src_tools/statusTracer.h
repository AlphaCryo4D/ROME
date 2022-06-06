/***************************************************************************
 *
 * Dana-Farber Cancer Institute, Harvard Medical School and Peking University
 * Intel® Parallel Computing Center for Structural Biology
 *
 * Authors: "Yong Bei Ma(galowma@gmail.com)"
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

#include "../src/resmap/resmap_util.h"
#include "../src/resmap/resmap_mpi.h"
#include "../src/resmap/resmap_error.h"

#define append_double_ptr(d,l,w) appendDoublePtr(d,l,w)
#define append_float_ptr(d,l,w) appendFloatPtr(d,l,w)
#define append_int_ptr(d,l,w) appendIntPtr(d,l,w)
#define append_bool_ptr(d,l,w) appendBoolPtr(d,l,w)

template<typename T>
static bool checkvec(std::ostream& out, T* vec1, T* vec2, int length){
    size_t count = 0;
    out.precision(12);
    T max_vec1,max_vec2,max_index;
    T max_error = -1;
    for (int i = 0; i < length; i++)
    {
        T error = fabs(vec1[i] - vec2[i]);
        if (error > max_error) {
            max_vec1 = vec1[i];
            max_vec2 = vec2[i];
            max_index = i;
            max_error = error;
        }
        if (error > 1e-8) count++;
    }
    /*MASTERNODE*/ if (count > 0) {
        out << "#### maximum error(index,value1,value2,error,precision) : (" << max_index << "," << max_vec1 << "," << max_vec2 << ","
        	<<max_error<<","<<max_error/max_vec1<<")";
        out << " (error_count,data_length,percentage) : (" << count << "," << length << "," << (double)count / length*100. << "%) ";
    }
    if (count > 0) return false;
    else return true;
};

template<typename T1,typename T2>
static bool checkDiffTypeVec(std::ostream& out, T1* vec1, T2* vec2, int length){
    size_t count = 0;
    out.precision(12);
    double max_vec1,max_vec2,max_index;
    double max_error = -1;
    for (int i = 0; i < length; i++)
    {
        double error = fabs(vec1[i] - vec2[i]);
        if (error > max_error) {
            max_vec1 = vec1[i];
            max_vec2 = vec2[i];
            max_index = i;
            max_error = error;
        }
        if (error > 1e-8) count++;
    }
    /*MASTERNODE*/ if (count > 0) {
        out << "#### maximum error(index,value1,value2,error) : (" << max_index << "," << max_vec1 << "," << max_vec2 << ","<<max_error<<") ";
        out << "(error_count,data_length,percentage) : (" << count << "," << length << "," << (double)count / length*100. << "%) ";
    }
    if (count > 0) return false;
    else return true;
};

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
    static bool* mallocBoolAligned(size_t size){
        return (bool*)aMalloc(size*sizeof(bool), 64);
    }
	template<typename T>
	static void freeAligned(T* ptr){
		aFree(ptr);
	}

public:
	StatusTracer()  { }
	~StatusTracer() { clear(); }

    void clear(){ doublePtr.resize(0); intPtr.resize(0); floatPtr.resize(0);boolPtr.clear();}

	void appendDoublePtr(double* data, int len, std::string what){
		doublePtr.push_back(DoublePtr({ data, len, what }));
	}
    void appendFloatPtr(float* data, int len, std::string what, bool compareFloatWithDouble = false){
        floatPtr.push_back(FloatPtr({data, len, what, compareFloatWithDouble}));
    }
	void appendIntPtr(int* data, int len, std::string what){
		intPtr.push_back(IntPtr({ data, len, what }));
	}
    void appendBoolPtr(bool* data, int len, std::string what){
        boolPtr.push_back(BoolPtr({ data, len, what }));
    }
	void backupStatus(std::string status_fn);

	void recoveryStatus(std::string status_fn);
		// read the tracing data to _backup.model file

	void checkStatus(std::string guidance_fn,bool do_recovery = true);

private:
	struct DoublePtr { double* data; int len; std::string what; };
	std::vector<DoublePtr> doublePtr;
    struct FloatPtr {float* data;int len;std::string what;bool compareFloatWithDouble;};
    std::vector<FloatPtr> floatPtr;
	struct IntPtr { int* data; int len; std::string what; };
	std::vector<IntPtr> intPtr;
    struct BoolPtr { bool* data; int len; std::string what; };
    std::vector<BoolPtr> boolPtr;
};


#endif
