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

#ifndef CHECKER_H_
#define CHECKER_H_

#include "checkerbase.h"

namespace Checker {
    //
    class Benchmark;
    static const int maxStudiedTids = 256;
    
    // SECTION TWO: BEING ABLE TO COLLECT THE DATA TO BE PROCESSED BY SECTION ONE
    
    void flush();
    // causes any buffering to be flushed to disk files
    // Call this infrequently but often enough
    // SECTION THREE: The ftrace file gets a timing log of a subset of the calls of a subset of threads
    //
    void setFtraceFile(std::string const & checker_ftrace);
    void unsetFtraceFile();
}

//===============================================================================================================================
// Scopes and the lives of a Scope
//
// The execution of a thread enters and exits scopes over time, the entering and exiting of a scope is a Life
// These lives are used to
//		provide stack traces during production code execution
//		measure execution time, including parallel v. serial execution
//
class Scope_Flags {
public:
    Scope_Flags() : v(0) {}
#define SCOPE_FLAGS \
    ELT(None,      0) SEP		/* The following are or'ed into this						*/ \
    ELT(Wait,	   1) SEP		/* The scope is the thread idling waiting for something     */ \
    ELT(Loop,	   2) SEP		/* The scope is a loop                                      */ \
    ELT(Iter,	   4) SEP		/* The scope is a loop body                                 */ \
    ELT(ChkBal,    8)			/* Check if the loop is load balanced						*/ \
    // end of macro
#define ELT(N,V)				\
    static Scope_Flags N()			 { return Scope_Flags(V); }	\
    bool        contains##N() { return !!(v & V);	  }
#define SEP
    SCOPE_FLAGS
#undef SEP
#undef ELT
    std::string toString() const;
    Scope_Flags operator|(Scope_Flags const rhs) const { return Scope_Flags(v | rhs.v); }
    Scope_Flags operator&(Scope_Flags const rhs) const { return Scope_Flags(v&rhs.v); }
private:
    int v;
    Scope_Flags(int v) : v(v) {}
};


class Scope : public NoCopy {
public:
    Scope(const char * name, Scope_Flags flags)
    : _name(name), _flags(flags), _file(nullptr), _line(0),
    _noteLoadBalanceCount(0), _worstLoadBalance(1.0)
    {}
    // Objects are constructed somewhere that no locking is required
    Scope& init(const char * file, int line) {
        // init only requires locking once
        if (file != _file || line != _line) initWkr(file, line);
        return *this;
    }
    const char* name() const { return _name; }
    const char* file() const { return _file; }
    int			line() const { return _line; }
    Scope_Flags flags() const { return _flags; }
    void noteLoadBalance(float loadBalance, Microseconds elapsed);
private:
    void initWkr(const char * file, int line);
    const char* const _name;
    Scope_Flags const _flags;
    const char*       _file;
    int			      _line;
    // internal implementation details
    class Impl; class Impl* _impl;
    size_t			  _noteLoadBalanceCount;
    float			  _worstLoadBalance;
};

class Life : public NoCopy {
public:
    Life(Scope& scope, Life* parent = NULL, Checker::Benchmark* benchmark = nullptr);	// parent is useful when known, needed across a omp_para boundary
    virtual ~Life();
    std::string path() const;
protected:
    int const					tid; // must be first so init'ed early
    Microseconds				iterUsed[Checker::maxStudiedTids];
    Scope&						scope;
    Life*				const   parent;
    Checker::Benchmark*	const   _benchmark;
    std::string  mutable		_path;
    Microseconds				startTime;
    // internal implementation details
    class Impl; class Impl* _impl;
    bool debug();
    void iterStarted (int tid, Microseconds now);
    void iterFinished(int tid, Microseconds used);
};

//#define DO_CHECKER_SCOP
#ifdef DO_CHECKER_SCOP
#define MAJOR_SCOPE2(NAME, LIFE)																	\
    static Scope checker_scope_##_LINE_(#NAME, Scope_Flags::None());								\
    Life  LIFE(checker_scope_##_LINE_.init(__FILE__, __LINE__));									\
// end of macro
#define MAJOR_SCOPE(NAME)																			\
    MAJOR_SCOPE2(NAME, checker_scope_##_LINE_##_live)												\
// end of macro
#else
#define MAJOR_SCOPE2(NAME, LIFE)																	
// end of macro
#define MAJOR_SCOPE(NAME)																		
// end of macro
#endif

// SECTION FOUR: MACROS TO MAKE IT EASY TO USE OR NOT USE THE ABOVE
//
//#define DO_CHECKER
#ifdef DO_CHECKER /* if do checker */

#define CHECKER_PARALLEL_COUNTED_FOR1(TYPE1, NAME1, LO1, HI1)													\
{																												\
	static Scope checker_loop(#NAME1, Scope_Flags::Loop()|Scope_Flags::ChkBal());								\
	static Scope checker_iter(#NAME1, Scope_Flags::Iter()|Scope_Flags::ChkBal());								\
	Life         checker_loop_live(checker_loop.init(__FILE__, __LINE__));										\
	auto & checker_iter_scope = checker_iter.init(__FILE__, __LINE__);											\
	PragmaInMacro(omp parallel for schedule(dynamic))															\
	for (TYPE1 NAME1 = LO1; NAME1 < HI1; NAME1++) {																\
	Life  checker_iter_live(checker_iter_scope, &checker_loop_live);											\
// end of macro

#define CHECKER_PARALLEL_COUNTED_FOR2(TYPE1, NAME1, LO1, HI1, TYPE2, NAME2, LO2, HI2)							\
{																												\
	static Scope checker_loop(#NAME1, Scope_Flags::Loop()|Scope_Flags::ChkBal());								\
	static Scope checker_iter(#NAME1, Scope_Flags::Iter()|Scope_Flags::ChkBal());								\
	Life         checker_loop_live(checker_loop.init(__FILE__, __LINE__));										\
	auto & checker_iter_scope = checker_iter.init(__FILE__, __LINE__);											\
	PragmaInMacro(omp parallel for collapse(2) schedule(dynamic))												\
	for (TYPE1 NAME1 = LO1; NAME1 < HI1; NAME1++)																\
	for (TYPE2 NAME2 = LO2; NAME2 < HI2; NAME2++) {																\
	Life  checker_iter_live(checker_iter_scope, &checker_loop_live);											\
// end of macro

#define CHECKER_PARALLEL_COUNTED_FOR3(TYPE1, NAME1, LO1, HI1, TYPE2, NAME2, LO2, HI2, TYPE3, NAME3, LO3, HI3) 	\
{																												\
    static Scope checker_loop(#NAME1, Scope_Flags::Loop()|Scope_Flags::ChkBal());								\
    static Scope checker_iter(#NAME1, Scope_Flags::Iter()|Scope_Flags::ChkBal());								\
    Life         checker_loop_live(checker_loop.init(__FILE__, __LINE__));										\
    auto & checker_iter_scope = checker_iter.init(__FILE__, __LINE__);											\
    PragmaInMacro(omp parallel for collapse(3) schedule(dynamic))												\
    for (TYPE1 NAME1 = LO1; NAME1 < HI1; NAME1++)																\
    for (TYPE2 NAME2 = LO2; NAME2 < HI2; NAME2++)																\
    for (TYPE3 NAME3 = LO3; NAME3 < HI3; NAME3++) {																\
    Life  checker_iter_live(checker_iter_scope, &checker_loop_live);											\
// end of macro

#define CHECKER_PARALLEL_COUNTED_FOR4(TYPE1, NAME1, LO1, HI1, TYPE2, NAME2, LO2, HI2, TYPE3, NAME3, LO3, HI3, TYPE4, NAME4, LO4, HI4) \
{																												\
    static Scope checker_loop(#NAME1, Scope_Flags::Loop()|Scope_Flags::ChkBal());								\
    static Scope checker_iter(#NAME1, Scope_Flags::Iter()|Scope_Flags::ChkBal());								\
    Life         checker_loop_live(checker_loop.init(__FILE__, __LINE__));										\
    auto & checker_iter_scope = checker_iter.init(__FILE__, __LINE__);											\
    PragmaInMacro(omp parallel for collapse(4) schedule(dynamic))												\
    for (TYPE1 NAME1 = LO1; NAME1 < HI1; NAME1++)																\
    for (TYPE2 NAME2 = LO2; NAME2 < HI2; NAME2++)																\
    for (TYPE3 NAME3 = LO3; NAME3 < HI3; NAME3++)																\
    for (TYPE4 NAME4 = LO4; NAME4 < HI4; NAME4++) {																\
    Life  checker_iter_live(checker_iter_scope, &checker_loop_live);											\
// end of macro

#define CHECKER_PARALLEL_COUNTED_FOR5(TYPE1, NAME1, LO1, HI1, TYPE2, NAME2, LO2, HI2, TYPE3, NAME3, LO3, HI3, TYPE4, NAME4, LO4, HI4, TYPE5, NAME5, LO5, HI5) \
{																												\
    static Scope checker_loop(#NAME1, Scope_Flags::Loop()|Scope_Flags::ChkBal());								\
    static Scope checker_iter(#NAME1, Scope_Flags::Iter()|Scope_Flags::ChkBal());								\
    Life         checker_loop_live(checker_loop.init(__FILE__, __LINE__));										\
    auto & checker_iter_scope = checker_iter.init(__FILE__, __LINE__);											\
    PragmaInMacro(omp parallel for collapse(5) schedule(dynamic))												\
    for (TYPE1 NAME1 = LO1; NAME1 < HI1; NAME1++)																\
    for (TYPE2 NAME2 = LO2; NAME2 < HI2; NAME2++)																\
    for (TYPE3 NAME3 = LO3; NAME3 < HI3; NAME3++)																\
    for (TYPE4 NAME4 = LO4; NAME4 < HI4; NAME4++)																\
    for (TYPE5 NAME5 = LO5; NAME5 < HI5; NAME5++) {																\
    Life  checker_iter_live(checker_iter_scope, &checker_loop_live);											\
// end of macro

#define CHECKER_PARALLEL_FOR_END																				\
	}}																											\
// end of macro

#else /* if not do checker */

#define CHECKER_PARALLEL_COUNTED_FOR1(TYPE1, NAME1, LO1, HI1)													\
	CHECKER_PARALLEL_COUNTED_FOR1_PLUS(TYPE1, NAME1, LO1, HI1, 1)												\
// end of macro

#define CHECKER_PARALLEL_COUNTED_FOR1_PLUS(TYPE1, NAME1, LO1, HI1, P1)											\
{																												\
	CheckerForBounds<TYPE1,int,int,int> checkerForBounds(__FILE__,__LINE__,LO1,HI1,P1);							\
	PragmaInMacro(omp parallel for schedule(dynamic))															\
	for (TYPE1 NAME1 = checkerForBounds.lo1; NAME1 < checkerForBounds.hi1; NAME1+=checkerForBounds.p1)			\
// end of macro

#define CHECKER_PARALLEL_COUNTED_FOR2(TYPE1, NAME1, LO1, HI1, TYPE2, NAME2, LO2, HI2)							\
	CHECKER_PARALLEL_COUNTED_FOR2_PLUS(TYPE1, NAME1, LO1, HI1, 1, TYPE2, NAME2, LO2, HI2, 1)					\
// end of macro

#define CHECKER_PARALLEL_COUNTED_FOR2_PLUS(TYPE1, NAME1, LO1, HI1, P1, TYPE2, NAME2, LO2, HI2, P2)				\
{																												\
	CheckerForBounds<TYPE1,TYPE2,int,int> checkerForBounds(__FILE__,__LINE__,LO1,HI1,P1,LO2,HI2,P2);			\
	PragmaInMacro(omp parallel for collapse(2) schedule(dynamic))												\
	for (TYPE1 NAME1 = checkerForBounds.lo1; NAME1 < checkerForBounds.hi1; NAME1+=checkerForBounds.p1)			\
	for (TYPE2 NAME2 = checkerForBounds.lo2; NAME2 < checkerForBounds.hi2; NAME2+=checkerForBounds.p2)			\
// end of macro

#define CHECKER_PARALLEL_COUNTED_FOR3(TYPE1, NAME1, LO1, HI1, TYPE2, NAME2, LO2, HI2, TYPE3, NAME3, LO3, HI3)	\
	CHECKER_PARALLEL_COUNTED_FOR3_PLUS(TYPE1, NAME1, LO1, HI1, 1, TYPE2, NAME2, LO2, HI2, 1, TYPE3, NAME3, LO3, HI3, 1) \
// end of macro

#define CHECKER_PARALLEL_COUNTED_FOR3_PLUS(TYPE1, NAME1, LO1, HI1, P1, TYPE2, NAME2, LO2, HI2, P2, TYPE3, NAME3, LO3, HI3, P3) \
{																											\
	CheckerForBounds<TYPE1,TYPE2,TYPE3,int> checkerForBounds(__FILE__,__LINE__,LO1,HI1,P1,LO2,HI2,P2,LO3,HI3,P3);	\
	PragmaInMacro(omp parallel for collapse(3) schedule(dynamic))											\
	for (TYPE1 NAME1 = checkerForBounds.lo1; NAME1 < checkerForBounds.hi1; NAME1+=checkerForBounds.p1)		\
	for (TYPE2 NAME2 = checkerForBounds.lo2; NAME2 < checkerForBounds.hi2; NAME2+=checkerForBounds.p2)		\
	for (TYPE3 NAME3 = checkerForBounds.lo3; NAME3 < checkerForBounds.hi3; NAME3+=checkerForBounds.p3)		\
// end of macro

#define CHECKER_PARALLEL_COUNTED_FOR4(TYPE1, NAME1, LO1, HI1, TYPE2, NAME2, LO2, HI2, TYPE3, NAME3, LO3, HI3, TYPE4, NAME4, LO4, HI4) \
	CHECKER_PARALLEL_COUNTED_FOR4_PLUS(TYPE1, NAME1, LO1, HI1, 1, TYPE2, NAME2, LO2, HI2, 1, TYPE3, NAME3, LO3, HI3, 1, TYPE4, NAME4, LO4, HI4, 1) \
// end of macro

#define CHECKER_PARALLEL_COUNTED_FOR4_PLUS(TYPE1, NAME1, LO1, HI1, P1, TYPE2, NAME2, LO2, HI2, P2, TYPE3, NAME3, LO3, HI3, P3, TYPE4, NAME4, LO4, HI4, P4) \
{																											\
	CheckerForBounds<TYPE1,TYPE2,TYPE3,TYPE4> checkerForBounds(__FILE__,__LINE__,LO1,HI1,P1,LO2,HI2,P2,LO3,HI3,P3,LO4,HI4,P4);	\
	PragmaInMacro(omp parallel for collapse(4) schedule(dynamic))											\
	for (TYPE1 NAME1 = checkerForBounds.lo1; NAME1 < checkerForBounds.hi1; NAME1+=checkerForBounds.p1)		\
	for (TYPE2 NAME2 = checkerForBounds.lo2; NAME2 < checkerForBounds.hi2; NAME2+=checkerForBounds.p2)		\
	for (TYPE3 NAME3 = checkerForBounds.lo3; NAME3 < checkerForBounds.hi3; NAME3+=checkerForBounds.p3)		\
	for (TYPE4 NAME4 = checkerForBounds.lo4; NAME4 < checkerForBounds.hi4; NAME4+=checkerForBounds.p4)		\
// end of macro

#define CHECKER_PARALLEL_COUNTED_FOR5(TYPE1, NAME1, LO1, HI1, TYPE2, NAME2, LO2, HI2, TYPE3, NAME3, LO3, HI3, TYPE4, NAME4, LO4, HI4, TYPE5, NAME5, LO5, HI5) \
	CHECKER_PARALLEL_COUNTED_FOR5_PLUS(TYPE1, NAME1, LO1, HI1, 1, TYPE2, NAME2, LO2, HI2, 1, TYPE3, NAME3, LO3, HI3, 1, TYPE4, NAME4, LO4, HI4, 1, TYPE5, NAME5, LO5, HI5, 1) \
// end of macro

#define CHECKER_PARALLEL_COUNTED_FOR5_PLUS(TYPE1, NAME1, LO1, HI1, P1, TYPE2, NAME2, LO2, HI2, P2, TYPE3, NAME3, LO3, HI3, P3, TYPE4, NAME4, LO4, HI4, P4, TYPE5, NAME5, LO5, HI5, P5) \
{																												\
	CheckerForBounds<TYPE1,TYPE2,TYPE3,TYPE4,TYPE5> checkerForBounds(__FILE__,__LINE__,LO1,HI1,P1,LO2,HI2,P2,LO3,HI3,P3,LO4,HI4,P4,LO5,HI5,P5);	\
	PragmaInMacro(omp parallel for collapse(5) schedule(dynamic))											\
	for (TYPE1 NAME1 = checkerForBounds.lo1; NAME1 < checkerForBounds.hi1; NAME1+=checkerForBounds.p1)		\
	for (TYPE2 NAME2 = checkerForBounds.lo2; NAME2 < checkerForBounds.hi2; NAME2+=checkerForBounds.p2)		\
	for (TYPE3 NAME3 = checkerForBounds.lo3; NAME3 < checkerForBounds.hi3; NAME3+=checkerForBounds.p3)		\
	for (TYPE4 NAME4 = checkerForBounds.lo4; NAME4 < checkerForBounds.hi4; NAME4+=checkerForBounds.p4)		\
	for (TYPE5 NAME5 = checkerForBounds.lo5; NAME5 < checkerForBounds.hi5; NAME5+=checkerForBounds.p5)		\
// end of macro

#define CHECKER_PARALLEL_FOR_END																				\
	}																											\
// end of macro

#define CHECKER_PARALLEL_ITER


#endif /* end of DO_CHECKER */


#endif