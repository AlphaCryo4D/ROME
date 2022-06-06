/***************************************************************************
 *
 * Intel® Parallel Computing Center for Structural Biology
 * Principal Investigator : Youdong (Jack) Mao (Youdong_Mao@dfci.harvard.edu)
 * Dana-Farber Cancer Institute, Harvard Medical School and Peking University
 *
 * Authors: "Yong Bei Ma(galowma@gmail.com) Bevin R Brett(bevin_brett@hotmail.com)"
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

#ifndef ERROR_H_
#define ERROR_H_

#include "./resmap_util.h"

// Nans and infinities can propagate a long way before causing damage
// Support for instrumenting the code to catch them soon after creation
//
//#define CHECK_FOR_NANS_ETC
static double checkDoubleValue(double value, const char * name) {
#ifdef CHECK_FOR_NANS_ETC
    if (std::abs(value) > 1e30) {
        static int count;
        if (count++ < 10) {
            std::cerr << name << " absurdly high!  value:" << value << std::endl;
        }
    }
#endif
    return value;
}
static double checkIrefsValue  (double value) { return checkDoubleValue(value, "IrefsValue"); }
static double checkFrefPadValue(double value) { return checkDoubleValue(value, "FrefPadsValue"); }
static double checkFrefRotValue(double value) { return checkDoubleValue(value, "FrefRotValue"); }
static double checkFrefCtfValue(double value) { return checkDoubleValue(value, "FrefCtfValue"); }

static double checkNonNegative(double x) { assert(x >= 0.0); return x; };

// check nan and infinite value
template <typename T>
T checkNotNanOrInf(T v, int line) {
#ifdef CHECK_FOR_NANS_ETC
    if (std::isnan(v) || std::isinf(v)) {
        std::cerr << "Nan or inf " << v << " found at line " << line << std::endl;
        std::cout << "Nan or inf " << v << " found at line " << line << std::endl;
        assert(false);
        EXIT_ABNORMALLY;
    }
#endif
    return v;
}
#define CHECK_NOT_NAN_OR_INF(F) checkNotNanOrInf(F, __LINE__)

// check indefinite value
static double checkNotIndef(double x, const char* expr, const char* file, int line) {
#ifdef CHECK_FOR_NANS_ETC
    if (x*x < 0) {
        std::cerr << file << ":" << line << " checkNotIndef() found an indefinite " << expr << ":" << x << std::endl;
        *(int*)(-1) = 0;
    }
#endif
    return x;
};
#define CHECK_NOT_IND(F) checkNotIndef(F, #F, __FILE__, __LINE__)

// ----------------------------------------------------------------

template <class T>
void assertAllEltsZero(size_t len, const T* ptr) {
#ifdef NDEBUG
    for (size_t i = 0; i < len; i+=len/15) {		// BEVIN
        if (T(0) == ptr[i]) continue;
        std::cerr << "assertAllEltsZero failed" << std::endl;
        bad_exit(1);
    }
#else
    for (size_t i = 0; i < len; i++) {
        assert(T(0) == ptr[i]);
    }
#endif
}

// ----------------------------------------------------------------

#ifndef ERROR_REPORT
#define ERROR_REPORT(what) { rome_error::assertFalse(what,__FILE__,__LINE__); }
#endif

#ifndef ERROR_CHECK
#define ERROR_CHECK(condition,what) if(condition) { ERROR_REPORT(what); }
#endif

//
class rome_error { // : public std::exception
public:
    explicit rome_error( const char*        _what_arg, const char* _file, int _line);
    explicit rome_error( const std::string& _what_arg, const char* _file, int _line);
    ~rome_error(){}
    const char* what() const /* noexcept */ { return what_arg; }
	static void assertFalse(const char*        _what_arg, const char* _file, int _line);
	static void assertFalse(const std::string& _what_arg, const char* _file, int _line);
private:
    const char* what_arg;
    int			line;
    const char* file;
};

#endif
