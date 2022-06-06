/***************************************************************************
 *
 * Authors: "Yongbei(Glow) Ma,Jiayi (Timmy) Wu, Youdong (Jack) Mao"
 *		Dana-Farber Cancer Institute, Harvard Medical School and Peking University
 *		"Brett, Bevin" Intel Corporation
 *		"Brett, Bevin" After retiring from Intel Corporation
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

// Overview:  See docs/checker.doc


#ifndef CHECKER_H_
#define CHECKER_H_

#include "./time.h"


namespace Checker {

	void configureFTraceSkipping(
		double collectSecs		= 15.0,
		double firstSkipSecs	= 15.0);
		//
		// FTrace files would be too large if they collected everything
		// so they skip over portions of the call subtree

	// SECTION ONE: BEING ABLE TO WRITE, READ, AND COMPARE DATA
	//
	class Instance;
	class Instance_Impl;

	template <typename T>
	class DimHdrs : public NoCopy {
	public:
		struct Hdr { 
		protected:
			friend class DimHdrs;
			friend class Instance_Impl;
			size_t    _index; 
			Instance* _instance; 
			Hdr(size_t index, Instance* instance) : _index(index), _instance(instance) {} 
		public:
			Hdr() : _index(0), _instance(NULL) {}
			size_t    index()    const { return _index;      }
			Instance* instance() const { return _instance;   }
			bool operator()()          { return !!_instance; }
		};

		typedef std::map<std::string, size_t>	StringToHeaderIndex;
		typedef std::vector<std::string>		HeaderIndexToString;

		StringToHeaderIndex const & stringToHeaderIndex() const { return _stringToHeaderIndex; }
		HeaderIndexToString const & headerIndexToString() const { return _headerIndexToString; }

		size_t size() const { return _headerIndexToString.size(); }

		void clear() {
			_stringToHeaderIndex.clear();
			_headerIndexToString.resize(0);
		}

	protected:
		friend class Instance;
		Hdr addIfMissing(std::string header, Instance* instance) {
			auto i = _stringToHeaderIndex.find(header);
			if (i != _stringToHeaderIndex.end()) return Hdr(i->second, instance);
			auto hdr = Hdr(_headerIndexToString.size(),instance);
			_headerIndexToString.push_back(header);
			_stringToHeaderIndex.insert(std::make_pair(header, hdr.index()));
			return hdr;
		}

	private:
		StringToHeaderIndex _stringToHeaderIndex;
		HeaderIndexToString _headerIndexToString;
	};

	class Row; typedef DimHdrs<Row> RowHdrs; typedef RowHdrs::Hdr RowHdr;
	class Col; typedef DimHdrs<Col> ColHdrs; typedef ColHdrs::Hdr ColHdr;
	typedef RowHdrs::Hdr RowHdr;
	typedef ColHdrs::Hdr ColHdr;

	class Instance : NoCopy {
	public:
		friend class Instance_Impl;

		Instance();
		~Instance();

		// Set the input and output streams
		// An output can be used as an input to a later instance
		//
		void read (std::istream & input);
		void write(std::ostream & output);

		RowHdrs const & rowHdrs() const;
		ColHdrs const & colHdrs() const;

		RowHdr addRowIfMissing(std::string const & header);
		ColHdr addColIfMissing(std::string const & header);

		ColHdr setDefaultColumnHeader(ColHdr const & colHdr);	// returns the previous default ColHdr
		ColHdr setDefaultColumnHeader(std::string const & header) { return setDefaultColumnHeader(addColIfMissing(header)); }
 
		// The cells can be either numbers or strings
		// If numbers, there is an acceptable range of values that will match
		//
		std::string getStr(RowHdr const & r, ColHdr const & c);				            std::string getStr(RowHdr const & r);
		void        set   (RowHdr const & r, ColHdr const & c, std::string const & to); void        set   (RowHdr const & r, std::string const & to);
		DoublePair  getDbl(RowHdr const & r, ColHdr const & c);	                        DoublePair  getDbl(RowHdr const & r);
		void        set   (RowHdr const & r, ColHdr const & c, DoublePair const & to);	void        set   (RowHdr const & r, DoublePair const & to);

	private:
		Instance_Impl* impl;
	};


	// The benchmark is the current instance.
	// If there is an input, it is read to initialize the instance
	//		and those values are used to check the validity of the updates
	// If there is an output, that output is what the write calls create
	//
	class Benchmark : public Instance {
	public:
		Benchmark();
		~Benchmark();
		void setFiles(std::string const & inputFilename, std::string const & outputFilename);
		bool hasAFile() const { return _hasAFile; }
		bool hasIFile() const { return inputFilename .size() > 0; }
		bool hasOFile() const { return outputFilename.size() > 0; }
		void write();
		void writeIfEnoughTimeElapsed();
		void write(Microseconds now);
		bool writeIfEnoughTimeElapsed(Microseconds now);	// true if written

		static Benchmark* defaultBenchmark();
		static void setDefaultBenchmark(Benchmark* to);

	private:
		Microseconds gapUntilNextWrite;
		Microseconds nextWriteTime;
		bool		 _hasAFile;
		std::string  inputFilename;
		std::string  outputFilename;
	};


	// SECTION TWO: BEING ABLE TO COLLECT THE DATA TO BE PROCESSED BY SECTION ONE

	void flush();
		// causes any buffering to be flushed to disk files
		// Call this infrequently but often enough


	// SECTION THREE: The ftrace file gets a timing log of a subset of the calls of a subset of threads
	//
	void setFtraceFile(std::string const & checker_ftrace);
	void unsetFtraceFile();

}


// SECTION FOUR: MACROS TO MAKE IT EASY TO USE OR NOT USE THE ABOVE 
//
#ifdef DO_CHECKER /* if do checker */

#define CHECKER_PARALLEL_COUNTED_FOR1(TYPE1, NAME1, LO1, HI1)														\
	{																												\
		static Scope checker_loop(#NAME1, Scope_Flags::Loop()|Scope_Flags::ChkBal());	\
		static Scope checker_iter(#NAME1, Scope_Flags::Iter()|Scope_Flags::ChkBal());	\
		Life         checker_loop_live(checker_loop.init(__FILE__, __LINE__));								\
		auto & checker_iter_scope = checker_iter.init(__FILE__, __LINE__);											\
		PragmaInMacro(omp parallel for schedule(dynamic))																	\
		for (TYPE1 NAME1 = LO1; NAME1 < HI1; NAME1++) {																\
			Life  checker_iter_live(checker_iter_scope, &checker_loop_live);								\
	// end of macro

#define CHECKER_PARALLEL_COUNTED_FOR2(TYPE1, NAME1, LO1, HI1, TYPE2, NAME2, LO2, HI2)								\
	{																												\
		static Scope checker_loop(#NAME1, Scope_Flags::Loop()|Scope_Flags::ChkBal());	\
		static Scope checker_iter(#NAME1, Scope_Flags::Iter()|Scope_Flags::ChkBal());	\
		Life         checker_loop_live(checker_loop.init(__FILE__, __LINE__));								\
		auto & checker_iter_scope = checker_iter.init(__FILE__, __LINE__);											\
		PragmaInMacro(omp parallel for collapse(2) schedule(dynamic))														\
		for (TYPE1 NAME1 = LO1; NAME1 < HI1; NAME1++)																\
		for (TYPE2 NAME2 = LO2; NAME2 < HI2; NAME2++) {																\
			Life  checker_iter_live(checker_iter_scope, &checker_loop_live);								\
	// end of macro

#define CHECKER_PARALLEL_COUNTED_FOR3(TYPE1, NAME1, LO1, HI1, TYPE2, NAME2, LO2, HI2, TYPE3, NAME3, LO3, HI3) \
	{																												\
		static Scope checker_loop(#NAME1, Scope_Flags::Loop()|Scope_Flags::ChkBal());	\
		static Scope checker_iter(#NAME1, Scope_Flags::Iter()|Scope_Flags::ChkBal());	\
		Life         checker_loop_live(checker_loop.init(__FILE__, __LINE__));								\
		auto & checker_iter_scope = checker_iter.init(__FILE__, __LINE__);											\
		PragmaInMacro(omp parallel for collapse(3) schedule(dynamic))														\
		for (TYPE1 NAME1 = LO1; NAME1 < HI1; NAME1++)																\
		for (TYPE2 NAME2 = LO2; NAME2 < HI2; NAME2++)																\
		for (TYPE3 NAME3 = LO3; NAME3 < HI3; NAME3++) {																\
			Life  checker_iter_live(checker_iter_scope, &checker_loop_live);								\
	// end of macro

#define CHECKER_PARALLEL_COUNTED_FOR4(TYPE1, NAME1, LO1, HI1, TYPE2, NAME2, LO2, HI2, TYPE3, NAME3, LO3, HI3, TYPE4, NAME4, LO4, HI4) \
	{																												\
		static Scope checker_loop(#NAME1, Scope_Flags::Loop()|Scope_Flags::ChkBal());	\
		static Scope checker_iter(#NAME1, Scope_Flags::Iter()|Scope_Flags::ChkBal());	\
		Life         checker_loop_live(checker_loop.init(__FILE__, __LINE__));								\
		auto & checker_iter_scope = checker_iter.init(__FILE__, __LINE__);											\
		PragmaInMacro(omp parallel for collapse(4) schedule(dynamic))														\
		for (TYPE1 NAME1 = LO1; NAME1 < HI1; NAME1++)																\
		for (TYPE2 NAME2 = LO2; NAME2 < HI2; NAME2++)																\
		for (TYPE3 NAME3 = LO3; NAME3 < HI3; NAME3++)																\
		for (TYPE4 NAME4 = LO4; NAME4 < HI4; NAME4++) {																\
			Life  checker_iter_live(checker_iter_scope, &checker_loop_live);								\
	// end of macro

#define CHECKER_PARALLEL_FOR_END																					\
	}}																												\
	// end of macro

#else /* if not do checker */

#define CHECKER_PARALLEL_COUNTED_FOR1(TYPE1, NAME1, LO1, HI1)													\
		CHECKER_PARALLEL_COUNTED_FOR1_PLUS(TYPE1, NAME1, LO1, HI1, 1)											\
    // end of macro

#define CHECKER_PARALLEL_COUNTED_FOR1_PLUS(TYPE1, NAME1, LO1, HI1, P1)											\
	{																											\
		CheckerForBounds<TYPE1,int,int,int> checkerForBounds(__FILE__,__LINE__,LO1,HI1,P1);						\
	PragmaInMacro(omp parallel for schedule(dynamic))														\
		for (TYPE1 NAME1 = checkerForBounds.lo1; NAME1 < checkerForBounds.hi1; NAME1+=checkerForBounds.p1)		\
    // end of macro

#define CHECKER_PARALLEL_COUNTED_FOR2(TYPE1, NAME1, LO1, HI1, TYPE2, NAME2, LO2, HI2)							\
		CHECKER_PARALLEL_COUNTED_FOR2_PLUS(TYPE1, NAME1, LO1, HI1, 1, TYPE2, NAME2, LO2, HI2, 1)				\
    // end of macro

#define CHECKER_PARALLEL_COUNTED_FOR2_PLUS(TYPE1, NAME1, LO1, HI1, P1, TYPE2, NAME2, LO2, HI2, P2)				\
	{																											\
		CheckerForBounds<TYPE1,TYPE2,int,int> checkerForBounds(__FILE__,__LINE__,LO1,HI1,P1,LO2,HI2,P2);		\
		PragmaInMacro(omp parallel for collapse(2) schedule(dynamic))											\
		for (TYPE1 NAME1 = checkerForBounds.lo1; NAME1 < checkerForBounds.hi1; NAME1+=checkerForBounds.p1)		\
		for (TYPE2 NAME2 = checkerForBounds.lo2; NAME2 < checkerForBounds.hi2; NAME2+=checkerForBounds.p2)		\
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

#define CHECKER_PARALLEL_FOR_END																				\
    }																											\
    // end of macro

#define CHECKER_PARALLEL_ITER

template <class T1, class T2, class T3, class T4>
struct CheckerForBounds {
	const char* const file;
	int const line;
	T1 const lo1, hi1, p1;
	T2 const lo2, hi2, p2;
	T3 const lo3, hi3, p3;
	T4 const lo4, hi4, p4;
	CheckerForBounds(
		const char* const file,
		int const line,
		T1 const lo1,     T1 const hi1,     T1 const p1,
		T2 const lo2 = 0, T2 const hi2 = 1, T2 const p2 = 1,
		T3 const lo3 = 0, T3 const hi3 = 1, T3 const p3 = 1,
		T4 const lo4 = 0, T4 const hi4 = 1, T4 const p4 = 1
		) 
	  :
		file(file), line(line),
		lo1(lo1), hi1(hi1), p1(p1),
		lo2(lo2), hi2(hi2), p2(p2),
		lo3(lo3), hi3(hi3), p3(p3),
		lo4(lo4), hi4(hi4), p4(p4)
	{
		size_t iterations =
			(lo1 >= hi1 || lo2 >= hi2 || lo3 >= hi3 || lo4 >= hi4)
			? 0
			: ( divRoundUp(hi1 - lo1, p1) *
				divRoundUp(hi2 - lo2, p2) *
				divRoundUp(hi3 - lo3, p3) *
				divRoundUp(hi4 - lo4, p4)
			  );
        bool debug_flag = false;
		if (debug_flag && iterations < omp_get_max_threads()) {
			auto msg = [&](std::ostream & os) {
				MASTERNODE os << "Load balance problem - only iterations:" << iterations << " at " << file << ":" << line << std::endl;
			};
			msg(std::cout);
			msg(std::cerr);
		}
	}
};


#endif /* end of DO_CHECKER */


#endif /* end of CHECKER_H_ */
