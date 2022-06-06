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

#include "./checkerbase.h"

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

    // benchmark
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
}
