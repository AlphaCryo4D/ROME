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

#include "util.h"		// used for building precompiled headers on Windows

#include "./checker.h"

#include <algorithm>
#include <omp.h>

namespace Checker {

	class Instance_Impl : NoCopy {
	public:
		ColHdr  defaultColHdr;
		RowHdrs rowHdrs;
		ColHdrs colHdrs;

		std::vector<std::vector<std::string>> grid;

		void clear() {
			colHdrs.clear();
			rowHdrs.clear();
			grid.resize(0);
		}

		std::string& operator()(size_t r, size_t c) {
			if (grid.size() <= r) grid.resize(rowHdrs.size());
			auto& row = grid[r];
			if (row.size() <= c) row.resize(colHdrs.size());
			return row[c];
		}

		std::string& operator()(RowHdr const & r, ColHdr const & c) {
			assert(r.instance()->impl == this);
			assert(c.instance()->impl == this);
			return (*this)(r.index(), c.index());
		}
	};

	Instance::Instance() : impl(sNew(Instance_Impl)) {}
	Instance::~Instance() { sDelete(impl); }

	// Read and write
	// An output file can be used as an input file to a later run
	//
	void Instance::read(std::istream & input) {

		impl->clear();

		std::string line;
		FastCharBuf buffer(1024);

		size_t row(0);
		size_t col(0);

		auto append = [&]() {
			std::string bufstr(buffer.ptr, buffer.ptr + buffer.size);
			if (row == 0) {
				if (col > 0) {
					auto h = impl->colHdrs.addIfMissing(bufstr, this);
					assert(h.index() == col - 1);
				}
			}
			else {
				if (col == 0) {
					auto h = impl->rowHdrs.addIfMissing(bufstr, this);
					assert(h.index() == row - 1);
				}
				else {
					(*impl)(row - 1, col - 1) = bufstr;
				}
			}
			col++;
		};

		while (std::getline(input, line)) {
			auto   p = line.c_str();
			buffer.size = 0;
			while (auto c = *p++) {
				if (c == ',') {
					append();
					buffer.size = 0;
				}
				else if (c != '\\') {
					buffer.push_back(c);
				}
				else {
					if (*p) buffer.push_back(*p++);
				}
			}
			append();
			col = 0; row++;
		}
	}

	void Instance::write(std::ostream & output) {
		for (size_t col = 0; col < impl->colHdrs.size(); col++)
			output << "," << impl->colHdrs.headerIndexToString()[col];
		output << std::endl;
		for (size_t row = 0; row < impl->grid.size(); row++) {
			output << impl->rowHdrs.headerIndexToString()[row];
			for (size_t col = 0; col < impl->colHdrs.size(); col++)
				output << "," << impl->grid[row][col];
			output << std::endl;
		}
	}

	// The files are .csv formatted files
	// The rows are named in the first column
	// The columns are the various runs of the collection - typically each column is a specific platform and s/w-version combination
	//		eg: KNLbU_n4p4t8-Rel_2016Jun19a 
	//				might mean a KNL rev b running Unix using 4 mpi nodes, each with 4 processes, with 4 threads/process
	//				running the Release code checked out early on 2016 Jun 19
	// This code does not attribute any meanings to the column headers
	//
	RowHdr Instance::addRowIfMissing(std::string const & header) {
		return impl->rowHdrs.addIfMissing(header, this);
	}
	ColHdr Instance::addColIfMissing(std::string const & header) {
		return impl->colHdrs.addIfMissing(header, this);
	}

	ColHdr Instance::setDefaultColumnHeader(ColHdr const & colHdr) {
		auto old = impl->defaultColHdr;
		impl->defaultColHdr = colHdr;
		return old;
	}


	// The cells can be either numbers or strings
	// If numbers, there is an acceptable range of values that will match
	//
	std::string Instance::getStr(RowHdr const & r, ColHdr const & c) {
		return (*impl)(r, c);
	}
	void		Instance::set(RowHdr const & r, ColHdr  const & c, std::string const & to) {
		(*impl)(r, c) = to;
	}

	DoublePair  Instance::getDbl(RowHdr const & r, ColHdr const & c) { return DoublePair(getStr(r, c)); }
	void        Instance::set(RowHdr const & r, ColHdr const & c, DoublePair const & to) {
		auto expectedStr = getStr(r, c);
		if (expectedStr.size() > 0) {
			auto expected = DoublePair(expectedStr);
			if (to.mid() < expected.lo) std::cerr << impl->colHdrs.headerIndexToString()[c.index()] << " less than expected: " << to.mid() << " < " << expected.mid() << std::endl;
			if (to.mid() > expected.hi) std::cerr << impl->colHdrs.headerIndexToString()[c.index()] << " more than expected: " << to.mid() << " > " << expected.mid() << std::endl;
		}
		(*impl)(r, c) = to.string();
	}

	std::string Instance::getStr(RowHdr const & r) { return getStr(r, impl->defaultColHdr); }
	void		Instance::set(RowHdr const & r, std::string const & to) { set(r, impl->defaultColHdr, to); }
	DoublePair  Instance::getDbl(RowHdr const & r) { return getDbl(r, impl->defaultColHdr); }
	void        Instance::set(RowHdr const & r, DoublePair const & to) { set(r, impl->defaultColHdr, to); }

	static Benchmark* _defaultBenchmark = nullptr;
	Benchmark* Benchmark::defaultBenchmark() { return _defaultBenchmark; }
	void       Benchmark::setDefaultBenchmark(Benchmark* to) { _defaultBenchmark = to; }

	Benchmark::Benchmark() : _hasAFile(false), gapUntilNextWrite(1e6), nextWriteTime(0.0) {
	}

	Benchmark::~Benchmark() {
		write();
	}

	void Benchmark::write(Microseconds now) {
		if (outputFilename.size() == 0) return;
		ofstreamCheckingCreated ofstream(outputFilename.c_str());
		Instance::write(ofstream);
		gapUntilNextWrite = (std::min)(60e6, gapUntilNextWrite * 2);
		nextWriteTime = now + gapUntilNextWrite;
		std::cerr << "Benchmark::write() wrote " << outputFilename
			<< " at " << size_t(now*1e-6) << " secs"
			<< ", next will be after " << size_t(nextWriteTime*1e-6) << " secs"
			<< std::endl;
	}
	void Benchmark::write() {
		if (outputFilename.size() == 0) return;
		write(timeInMicroseconds());
	}

	bool Benchmark::writeIfEnoughTimeElapsed(Microseconds now) {
		if (now <= nextWriteTime) return false;
		write(now);
		return true;
	}

	void Benchmark::writeIfEnoughTimeElapsed() {
		if (outputFilename.size() == 0) return;
		writeIfEnoughTimeElapsed(timeInMicroseconds());
	}

	void Benchmark::setFiles(std::string const & inputFilename, std::string const & outputFilename) {
		this->inputFilename = inputFilename;
		this->outputFilename = outputFilename;
		_hasAFile = (inputFilename.size() > 0 || outputFilename.size() > 0);
	}


	// Ftrace files
	//	These would get very large if all the call tree was kept, especially if kept for all threads.
	//  For any node in the tree, we would like to understand its first few, some middle, and its last few children.
	//
	typedef unsigned __int64 Nanosecs;

	static const int hiTidStudied = (std::min)(omp_get_max_threads(), maxStudiedTids);

	static double ftraceCollectSecs = 15.0;
	static double ftraceFirstSkipSecs = 15.0;
	void configureFTraceSkipping(
		double collectSecs,
		double firstSkipSecs) {
		ftraceCollectSecs = collectSecs;
		ftraceFirstSkipSecs = firstSkipSecs;
	}

	class Ftrace {
		static bool const debug = false;

		AccurateTimer*	accurateTimer;
		bool			suppressed;
		int				flushCount;

		static void digits(std::ostream & os, __int64 number, int numberOfDigits) {
			if (numberOfDigits > 1) digits(os, number / 10, numberOfDigits - 1);
			os << char('0' + number % 10);
		}

	public:
		Ftrace()
			: accurateTimer(NULL),
			suppressed(false),
			flushCount(0),
			perTids(hiTidStudied)
		{
			for (int tid = 0; tid < hiTidStudied; tid++) perTids[tid].init(this, tid);
		}

		~Ftrace() {
			std::cerr << "BEVIN " << __FILE__ << "~Ftrace()" << std::endl;
		}

		Nanosecs nanoTimeSinceEra() {
			if (!accurateTimer) return 0;
			return
				Nanosecs(
					double(accurateTimer->sinceInited()) /
					(accurateTimer->countPerMicrosecond()*1e-3));
		}

		bool isBeingFtraced(int tid) {
			return (tid < hiTidStudied) && !suppressed && !!accurateTimer;
		}

		bool isBeingFtraced() { return isBeingFtraced(omp_get_thread_num()); }

		bool isTracing() const { return !!accurateTimer; }

		void setFtraceFilename(std::string const & checker_ftrace) {
			accurateTimer = sNew(AccurateTimer);
			for (size_t tid = 0; tid < hiTidStudied; tid++) perTids[tid].setFtraceFilename(checker_ftrace);
		}

		void unsetFtraceFilename() {
			flush();
			for (size_t tid = 0; tid < hiTidStudied; tid++) perTids[tid].unsetFtraceFilename();
			sDelete(accurateTimer);
		}

		class PerTid {
			Ftrace*			parent;
			int				tid;
			std::ostream*	ftrace_os;
		public:
			PerTid() : parent(NULL), ftrace_os(NULL) {}

			class Line;
			class LinePool;
			class Lines;
			class Region;

			class RegionBase {
			protected:
				Ftrace*		const ftrace;
				const char* const file;
				int			const line;
				const char* const funcEtc;
				int			const tid;
				bool        const tracing;
				bool		const debug;
				Line*		firstTentativeLine;
			public:
				RegionBase(Ftrace* const ftrace, const char* file, int line, const char* funcEtc, int tid)
					: ftrace(ftrace),
					file(file),
					line(line),
					funcEtc(funcEtc),
					tid(tid),
					firstTentativeLine(NULL),
					debug(false && (tid == 1)),
					tracing(ftrace->isBeingFtraced(tid)) {
					assert(!!file);
					assert(line > 0);
					assert(!!funcEtc);
				}

				bool shouldWriteRegion(Nanosecs nanoTime) {
					if (!tracing) return false;
					auto regionDuration = nanoTime - firstTentativeLine->nanoTime;
					return (regionDuration >= nanoTime / 1000);
				}

				void definitelyWrite() {
					assert(tracing);
					firstTentativeLine = NULL;
				}
			};

			class Line : public NoCopy {
				Ftrace*		 parent;
				int			_tid;
				RegionBase* _inTentativeWrittenLiveRegion;
			public:
				Line() : parent(NULL), newer(NULL) { }
				~Line() { *(int*)(-1) = 0; }	// do not call!

				void init(Ftrace* const parent, int tid) { this->parent = parent; this->_tid = tid; }

				int tid() { return _tid; }

				Line *newer, *older;
				Nanosecs nanoTime;
				bool     isB;
				const char* file; int line; const char* funcEtc;
				void init(
					bool isB,
					const char* file, int line, const char* funcEtc) {
					assert(!!file);
					assert(line > 0);
					assert(!!funcEtc);
					this->isB = isB;
					this->_inTentativeWrittenLiveRegion = NULL;
					this->file = file;
					this->line = line;
					this->funcEtc = funcEtc;
					this->newer = NULL;
					this->older = NULL;
				}
				void setInTentativeWrittenLiveRegion(RegionBase* to) {
					_inTentativeWrittenLiveRegion = to;
				}
				RegionBase* isInTentativeWrittenLiveRegion() {
					return _inTentativeWrittenLiveRegion;
				}
				void definitelyWrite() {
					_inTentativeWrittenLiveRegion->definitelyWrite();
					_inTentativeWrittenLiveRegion = NULL;
				}
				void write(std::ostream & os) {
					unsigned __int64 seconds = nanoTime / 1000000000LL;
					unsigned __int64 fraction = nanoTime % 1000000000LL / 1000;
#pragma omp critical
					{
						os << "                -";
						digits(os, _tid, 4);      // tid that wrote the trace record
						os << "  [";
						digits(os, _tid, 3);
						os << "] ";
						digits(os, seconds, 6);
						os << ".";
						digits(os, fraction, 6);
						os << ": tracing_mark_write: "
							<< (isB ? "B|" : "E|")
							<< _tid
							<< "|"
							<< funcEtc << "@" << line
							<< std::endl;
					}
				}
			};

			void init(Ftrace* const parent, int tid) {
				this->parent = parent;
				this->tid = tid;
				linePool.init(this);
				lines.init(this);
			}

			void setFtraceFilename(std::string const & checker_ftrace) {
				char suffix[3]; suffix[0] = '_'; suffix[1] = ('0' + tid); suffix[2] = 0;
				ftrace_os = ofstreamCheckingCreated::make(checker_ftrace + suffix + ".txt");
				std::cerr << "Creating " << checker_ftrace + suffix << (ftrace_os->fail() ? " - failed" : "") << std::endl;
				(*ftrace_os) << "#           TASK-PID   CPU#      TIMESTAMP  FUNCTION" << std::endl;
				(*ftrace_os) << "#              | |       |          |         |" << std::endl;
			}               //  "                -0001  [000] 123456.123456: %s: %s%s%n"

			void unsetFtraceFilename() {
				auto old = ftrace_os;
				ftrace_os = nullptr;
				sDelete(old);
			}

			class LinePool {
				PerTid* parent;
			public:
				LinePool() : available(nullptr) {}
				void init(PerTid* const parent) { this->parent = parent; }
				Line* remove() {
					auto result = available;
					if (!result) {
						static const int count = 10000;
						auto results = vNew(Line, count);			// Make more than one to cut down on malloc locking
						for (int i = 0; i < count; i++) {
							results[i].init(parent->parent, parent->tid);
							results[i].newer = available;
							available = &results[i];
						}
						result = available;
					}
					available = result->newer;
					return result;
				}
				void insert(Line* oldest, Line* newest) {
					assert(parent->tid == oldest->tid());
					newest->newer = available;
					available = oldest;
				}
			private:
				Line* available;
			} linePool;

			class Lines {
			public:
				PerTid* parent;
				bool    flushPending;

				Lines() : parent(NULL), flushPending(false) {
					_oldest = NULL;
					_newest = NULL;
				}
				void init(PerTid* parent) {
					this->parent = parent;
				}
				~Lines() {
					if (parent) flush(parent->parent->nanoTimeSinceEra());
				}

				void flush(Nanosecs nanoTime) {
					if (!_oldest || parent->tid != omp_get_thread_num()) {
						flushPending = true;
						return;
					}
					flushPending = false;
					Line* firstWrittenIfAny = _oldest;
					Line* lastWritten = NULL;
					for (auto p = _oldest; !!p; p = p->newer) {
						auto r = p->isInTentativeWrittenLiveRegion();
						if (!!r) {
							if (r->shouldWriteRegion(nanoTime)) {
								r->definitelyWrite();
							}
							else {
								break;
							}
						}
						p->write(*parent->ftrace_os);
						lastWritten = p;
					}
					if (!lastWritten) {
						// std::cerr << "DEBUG" << std::endl;	// BEVIN
						return;
					}
					_oldest = lastWritten->newer;
					if (!_oldest) {
						_newest = NULL;
					}
					parent->linePool.insert(firstWrittenIfAny, lastWritten);
				}
				void flushIfPending(Nanosecs nanoTime) { if (flushPending) flush(nanoTime); }

				void insert(Line * ftl, Nanosecs nanoTime, int tid) {
					flushIfPending(nanoTime);
					ftl->nanoTime = nanoTime;
					ftl->newer = NULL;
					ftl->older = _newest;
					if (_newest) _newest->newer = ftl; else _oldest = ftl;
					_newest = ftl;
				}
				void removeLines(Line* oldest, Nanosecs nanoTime) {
					Line* newest = this->_newest;
					assert(!!newest);
					auto prev = oldest->older;
					auto& prevNewer = prev ? prev->newer : _oldest;
					prevNewer = NULL;
					_newest = prev;
					parent->linePool.insert(oldest, newest);
					flushIfPending(nanoTime);
				}
			private:
				Line *_oldest, *_newest;
			} lines;
		};
		std::vector<PerTid> perTids;

#define REGION_TYPES ELT(Step) SEP ELT(Loop) SEP ELT(Iter)

#define ELT(NAME) class Region##NAME;
#define SEP
		REGION_TYPES
#undef SEP
#undef ELT

			class Region : public PerTid::RegionBase {
			public:
				Region(Ftrace* const ftrace, const char* file, int line, const char* funcEtc, int tid)
					: RegionBase(ftrace, file, line, funcEtc, tid)
				{
					if (!tracing) return;
					firstTentativeLine = BorE(true);
					firstTentativeLine->setInTentativeWrittenLiveRegion(this);
				}

				virtual ~Region() {
					if (!tracing) return;
					BorE(false);
					if (!!firstTentativeLine) {
						firstTentativeLine->setInTentativeWrittenLiveRegion(NULL);
					}
				}

				PerTid::Line* BorE(bool makeB) {
					if (debug) std::cerr << "Checker::Region::BorE tid:" << tid << (makeB ? " begin" : " end") << (!!tracing ? "" : " not tracing") << std::endl;

					auto nanoTime = ftrace->nanoTimeSinceEra();

					auto& perTid = ftrace->perTids[tid];
					if (!makeB && firstTentativeLine && !shouldWriteRegion(nanoTime)) {
						if (debug) std::cerr << "Checker::Region::BorE tid:" << tid << " discarding region" << std::endl;
						perTid.lines.removeLines(firstTentativeLine, nanoTime);
						return NULL;
					}

					auto ftl = perTid.linePool.remove();
					ftl->init(makeB, file, line, funcEtc);
					perTid.lines.insert(ftl, nanoTime, tid);
					if (debug) std::cerr << "Checker::Region::BorE tid:" << tid << " inserting line" << std::endl;

					return ftl;
				}

#define ELT(NAME)																				\
			Region##NAME* toRegion##NAME() { return this?this->toRegion##NAME##Wkr():NULL; }	\
			virtual Region##NAME* toRegion##NAME##Wkr() { return NULL; }						\
			// end of macro
#define SEP
				REGION_TYPES
#undef SEP
#undef ELT
		};

		class RegionLoop : public Region {
		public:
			RegionLoop(
				Ftrace* const ftrace,
				const char* file, int line, const char* funcEtc, int tid)
				: Region(ftrace, file, line, funcEtc, tid) {
			}

			~RegionLoop() {
			}
			static RegionLoop* make(Ftrace* const ftrace,
				const char* file, int line, const char* funcEtc, int tid) {
#include "./util_heap_undefs.h"
				return sNewA(RegionLoop,(ftrace, file, line, funcEtc, tid));
#include "./util_heap_defs.h"
			}
		};


		// These are when the code spreads the iterations across the threads
		//
		class RegionIter : public Region {
			RegionLoop& regionLoop;
		public:
			RegionIter(Ftrace* const ftrace, RegionLoop& regionLoop, const char* file, int line, const char* funcEtc, int tid)
				: Region(ftrace, file, line, funcEtc, tid), regionLoop(regionLoop) {
			}
			~RegionIter() {
			}
			static RegionIter* make(Ftrace* const ftrace, RegionLoop& regionLoop, const char* file, int line, const char* funcEtc, int tid) {
#include "./util_heap_undefs.h"
				return sNewA(RegionIter,(ftrace, regionLoop, file, line, funcEtc, tid));
#include "./util_heap_defs.h"
			}
		};

		class RegionOtherScope : public Region {
		public:
			RegionOtherScope(Ftrace* const ftrace, const char* file, int line, const char* funcEtc, int tid)
				: Region(ftrace, file, line, funcEtc, tid) {
			}
			~RegionOtherScope() {
			}
			static RegionOtherScope* make(Ftrace* const ftrace, const char* file, int line, const char* funcEtc, int tid) {
#include "./util_heap_undefs.h"
				return sNewA(RegionOtherScope,(ftrace, file, line, funcEtc, tid));
#include "./util_heap_defs.h"
			}
		};

		void flush(int tid) {
			if (isBeingFtraced(tid)) perTids[tid].lines.flush(nanoTimeSinceEra());
		}

		void flush() {
			auto nanoTime = nanoTimeSinceEra();
			for (size_t tid = 0; tid < hiTidStudied; tid++) perTids[tid].lines.flush(nanoTime);
			if (!accurateTimer) return;
			auto elapsed = accurateTimer->sinceInited() / accurateTimer->countPerMicrosecond() / 1e6;
			static double nextTransition = ftraceCollectSecs;	// collect these first seconds
			static double nextSkip = ftraceFirstSkipSecs;	// then skip this amount
			if (elapsed < nextTransition) {
				if (debug && (flushCount++ < 10)) std::cerr << "Tuning::flush elapsed=" << elapsed << " < nextTransition:" << nextTransition << std::endl;
				return;
			}
			if (debug) std::cerr << "Tuning::flush elapsed=" << elapsed << " >= nextTransition:" << nextTransition << (suppressed ? "was suppressed" : "") << std::endl;
			if (suppressed) {
				suppressed = false;
				nextTransition = elapsed + ftraceCollectSecs;						// collect another 15 seconds
			}
			else {
				suppressed = true;
				nextTransition = elapsed + nextSkip; nextSkip *= 2.0;	// skip this amount, twice as much next time
			}
		}
	};
	Ftrace* ftrace() {
		static Ftrace* ptr;
		if (!ptr) ptr = sNew(Ftrace);
		return ptr;
	}

	void setFtraceFile(std::string const & checker_ftrace) {
		ftrace()->setFtraceFilename(checker_ftrace);
	}

	void unsetFtraceFile() {
		ftrace()->unsetFtraceFilename();
	}

	void flush() {
		ftrace()->flush();
	}
}

using namespace Checker;

// Scopes
//	There are several planned uses for these
//
//	Some of them will have their elapsed times recorded in the files 
//	so that they can be compared across runs easily
//	  
//	Many (all?) of them will have their start and finish times recorded in the 
//	performance traces that are used for diagnosing where the time is going
//	  
//	Some of them will be used for diagnosing mpi load imbalances 
//
class Scope::Impl : public NoCopy {
};

static const char* persist(std::string const & s) {
	auto r = vNew(char, s.size() + 1);
	strcpy(r, s.c_str());
	return r;
}

void Scope::noteLoadBalance(float loadBalance, Microseconds elapsed) {
	// if (_worstLoadBalance < loadBalance) return;
	_worstLoadBalance = (std::min)(_worstLoadBalance, loadBalance);
	std::cerr << std::endl
		<< "Scope::noteLoadBalance"
		<< this->file() << ":" << this->line() << " " << this->name()
		<< " elapsed:" << elapsed
		<< " used cycles:" << int(loadBalance*100.0) << "%"
		<< std::endl;
}

void Scope::initWkr(const char * file, int line) {
	// Should only be called once per Scope, 
	// The Scopes are statically allocated so should not bottlenck
#pragma omp critical
	if (_line > 0) {
		assert(_line == line);
		assert(_file == file);
	}
	else {
		assert(line > 0);
		_line = line;
		_file = file;
	}
}

bool Life::debug() {
	return false;
}

class Life::Impl : public NoCopy {
public:
	Ftrace::Region* ftraceRegionBase;
	Impl(Life* parent) : ftraceRegionBase(NULL) {}
	~Impl() {
		sDelete(ftraceRegionBase);
	}
	static Impl* make(Life* parent) {
#include "./util_heap_undefs.h"
		return sNewA(Impl,(parent));
#include "./util_heap_defs.h"
	}
};

static std::vector<Life*> * volatile topLivesPtr;
std::vector<Life*>& topLives() {
	if (!topLivesPtr) {
		auto p = sNew(std::vector<Life*>);
		p->resize(hiTidStudied);
		for (auto i : *p) i = nullptr;
		topLivesPtr = p;
	}
	return *topLivesPtr;
}


Life::Life(Scope& scope, Life* parent, Benchmark* benchmark)
	: tid(omp_get_thread_num()),
	parent(parent ? parent : topLives()[tid]),
	scope(scope),
	_impl(NULL),
	_benchmark(benchmark)
{
	if (scope.flags().containsChkBal()) {
		startTime = timeInMicroseconds();
		if (tid < hiTidStudied && parent && scope.flags().containsIter()) {
			parent->iterStarted(tid, startTime);
		}
		if (scope.flags().containsLoop()) {
			for (int t = 0; t < hiTidStudied; t++) {
				iterUsed[t] = 0;
			}
		}
	}

	auto ftracePtr = ftrace();

	if (tid != 0) {
		if (!ftracePtr->isBeingFtraced(tid)) {
			return;
		}
	}

	auto debug = this->debug();
	if (debug) std::cerr << "ScopeLive::ScopeLive " << scope.name();
	_impl = Impl::make(this);
	{
		Ftrace::Region* r = NULL;
		if (!ftracePtr->isTracing()) {
			// nothing to do
		} else if (scope.flags().containsIter()) {
			auto impl = parent->_impl;
			auto frb = impl->ftraceRegionBase;
			auto loop = frb->toRegionLoop();
			r = Ftrace::RegionIter::make(ftracePtr, *loop, scope.file(), scope.line(), scope.name(), tid);
		} else if (scope.flags().containsLoop()) {
			r = Ftrace::RegionLoop::make(ftracePtr, scope.file(), scope.line(), scope.name(), tid);
		} else {
			auto fileX = scope.file();
			if (!fileX)
				scope.file();
			r = Ftrace::RegionOtherScope::make(ftracePtr, scope.file(), scope.line(), scope.name(), tid);
		}
		_impl->ftraceRegionBase = r;
	}
	topLives()[tid] = this;
	if (tid == 0) {
		if (_benchmark) {
			if (debug) std::cerr << " is a benchmark";
			if (_benchmark->hasAFile()) {
				if (debug) std::cerr << " that is being timed" << std::endl;
				startTime = timeInMicroseconds();
			}
		}
		if (debug) std::cerr << std::endl;
	}
}

Life::~Life() {
	auto debug = this->debug();
	if (debug) std::cerr << "ScopeLive::~ScopeLive " << scope.name();

	if (scope.flags().containsChkBal()) {
		auto now = timeInMicroseconds();
		if (parent && scope.flags().containsIter()) {
			parent->iterFinished(tid, now - startTime);
		}
		if (tid < hiTidStudied && scope.flags().containsLoop()) {
			auto elapsed = now - startTime;
			Microseconds used = 0;
			for (int t = 0; t < hiTidStudied; t++) {
				used += iterUsed[t];
			}
			auto loadBalance = float(used) / (float(elapsed) * omp_get_max_threads());
			scope.noteLoadBalance(loadBalance, elapsed);
		}
	}

	bool ftrace_flush(false);
	if (tid != 0) {
		ftrace()->flush(tid);
	}
	else {
		if (_benchmark) {
			ftrace_flush = true;
			auto & bm = *_benchmark;
			if (bm.hasAFile()) {
				auto now = timeInMicroseconds();
				auto elapsed = now - startTime;
				bm.set(
					bm.addRowIfMissing(path()),
					DoublePair::midFrac(elapsed, 0.05));
				auto wrote = bm.writeIfEnoughTimeElapsed(now);
				if (debug) std::cerr << " after writing csv";
			}
		}
		if (debug) std::cerr << std::endl;
		topLives()[tid] = parent;
	}
	if (!!_impl) sDelete(_impl);
	if (ftrace_flush) ftrace()->flush();
}

std::string Life::path() const {
	if (_path.size() == 0) {
		if (parent) { _path += parent->path(); _path += "->"; }
		_path += scope.name();
	}
	return _path;
}

std::string Scope_Flags::toString() const {
	std::string result = "{";
	const char* prefix = "";
	auto append = [&](const char* name) { result += prefix; result += name; prefix = ","; };
#define ELT(N,V) if (v & V) append(#V);
#define SEP
	SCOPE_FLAGS
#undef SEP
#undef ELT
		result += "}";
	return result;
}

void Life::iterStarted(int tid, Microseconds now) {
	assert(scope.flags().containsLoop());
}

void Life::iterFinished(int tid, Microseconds used) {
	assert(scope.flags().containsLoop());
	assert(tid < hiTidStudied);
	iterUsed[tid] += used;
}

