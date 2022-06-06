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
// Overview:  See docs/checker.doc

#ifndef CHECKERBASE_H_
#define CHECKERBASE_H_

#include "../src/resmap/resmap_time.h"

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

    // Row and Col header
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
	//typedef RowHdrs::Hdr RowHdr; // redundant declaration
	//typedef ColHdrs::Hdr ColHdr;

    
    // instance is table , example :
    // ------------------------------------
    // |      | col1 | col2 | col3 |
    // | row1 |  3.1 |  "t" |      |
    // | row2 |      |      |      |
    // ------------------------------------
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
}

// Check the load balance problem
template <class T1, class T2, class T3, class T4, class T5 = int>
struct CheckerForBounds {
    const char* const file;
    int const line;
    T1 const lo1, hi1, p1;
    T2 const lo2, hi2, p2;
    T3 const lo3, hi3, p3;
    T4 const lo4, hi4, p4;
    T5 const lo5, hi5, p5;
    CheckerForBounds(
         const char* const file,
         int const line,
         T1 const lo1,     T1 const hi1,     T1 const p1,
         T2 const lo2 = 0, T2 const hi2 = 1, T2 const p2 = 1,
         T3 const lo3 = 0, T3 const hi3 = 1, T3 const p3 = 1,
         T4 const lo4 = 0, T4 const hi4 = 1, T4 const p4 = 1,
         T5 const lo5 = 0, T5 const hi5 = 1, T5 const p5 = 1
         )
    :
    file(file), line(line),
    lo1(lo1), hi1(hi1), p1(p1),
    lo2(lo2), hi2(hi2), p2(p2),
    lo3(lo3), hi3(hi3), p3(p3),
    lo4(lo4), hi4(hi4), p4(p4),
    lo5(lo5), hi5(hi5), p5(p5)
    {
        size_t iterations =
        (lo1 >= hi1 || lo2 >= hi2 || lo3 >= hi3 || lo4 >= hi4 || lo5 >= hi5)
        ? 0
        : ( divRoundUp(hi1 - lo1, p1) *
           divRoundUp(hi2 - lo2, p2) *
           divRoundUp(hi3 - lo3, p3) *
           divRoundUp(hi4 - lo4, p4) *
           divRoundUp(hi5 - lo5, p5)
           );
        bool debug_flag = false;
        if (debug_flag && iterations < omp_get_max_threads()) {
            auto msg = [&](std::ostream & os) {
                os << "Load balance problem - only iterations:" << iterations << " at " << file << ":" << line << std::endl;
            };
            msg(std::cout);
            msg(std::cerr);
        }
    }
};


#endif /* end of CHECKER_H_ */
