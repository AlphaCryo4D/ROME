/***************************************************************************
 *
 * Dana-Farber Cancer Institute, Harvard Medical School and Peking University
 * Intel® Parallel Computing Center for Structural Biology
 *
 * Authors: "Bevin R Brett(bevin_brett@hotmail.com) 2017-feb-18"
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
#ifndef PEAKFINDER_H_
#define PEAKFINDER_H_

#include "../src/resmap/resmap_util.h"

//===============================================================================================================================
// PeakFinder finds the local maximum in a multi-dimensional sampled space of noisy measurements
//
class PeakFinder : public NoCopy {
public:
	struct Bound { int begin; int end; };
	typedef std::vector<Bound>	Bounds;
	typedef std::vector<int>	Coordinates;

	PeakFinder(Bounds const & bounds);
	~PeakFinder();

	bool hasFoundPeak() const;
	float volume() const;
	size_t numberAssignmentsReported() const;

	typedef int ResultIndex;
	static const ResultIndex noAssignment = -1;
	class Assignment {
	public:
		~Assignment() {}
		virtual void coords(ResultIndex resultIndex, Coordinates const & coordinates) = 0;
	};

	ResultIndex askForAssignment(Assignment & assignment);
	void assignmentMustBeRedone(ResultIndex resultIndex);
	void reportResult(ResultIndex resultIndex, float result);

	Coordinates const & bestCoordinates() const;
	float bestResult() const;

private:
	struct Internals;
	Internals* const internals;
};

#endif

