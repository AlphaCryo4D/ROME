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
//#include "util.h"		// used for building precompiled headers on Windows

#include "peakFinder.h"

struct PeakFinder::Internals {
	static const int	numberOfTries = 1;
	Lock				lock;
	Bounds const &		bounds;
	std::vector<int>	divisionBegin;
	std::vector<int>	divisionStep;
	std::vector<int>	divisionEnd;
	Coordinates			nextCoords;
	int					nextCount;
	float				nextResult;
	Coordinates			bestCoords;
	float				bestResult;
	bool				done;
	size_t				numberAssignmentsGiven;
	size_t				numberAssignmentsReported;
	std::vector<Coordinates> assignedCoordinates;
	std::vector<ResultIndex> assignedCoordinatesAvailable;
	std::vector<ResultIndex> assignedCoordinatesToRedo;

	MeanAndStdDeviation	stats;

	Internals(Bounds const & bounds)
	  : bounds           (bounds),
		divisionBegin    (bounds.size()),
		divisionStep	 (bounds.size()),
		divisionEnd		 (bounds.size()),
		nextCoords       (bounds.size()),
		nextCount (numberOfTries),
		bestResult(-std::numeric_limits<float>::max()),
		done(false),
		numberAssignmentsGiven(0),
		numberAssignmentsReported(0),
		assignedCoordinates(omp_get_max_threads())
	{
		for (int i = 0; i < bounds.size(); i++) {
			divisionBegin[i] = bounds[i].begin;
			divisionEnd  [i] = bounds[i].end;
			divisionStep [i] = 1;
			while (divisionStep[i] < (divisionEnd[i] - divisionBegin[i]) / 2) divisionStep[i] *= 2;
			nextCoords   [i] = bounds[i].begin;
		}
		for (int i = 0; i < assignedCoordinates.size(); i++)
			assignedCoordinatesAvailable.push_back(i);
	}

	void assignmentMustBeRedone(ResultIndex resultIndex) {
		assignedCoordinatesToRedo.push_back(resultIndex);
	}
	void reportResult(int resultIndex, float result) {
		numberAssignmentsReported++;
		stats.insert(result);
		if (bestResult < result) {
			bestResult = result;
			bestCoords = assignedCoordinates[resultIndex];
		}
		assignedCoordinatesAvailable.push_back(resultIndex);
	}

	int askForAssignment(Assignment & assignment) {
		if (!assignedCoordinatesToRedo.empty()) {
			auto ri = assignedCoordinatesToRedo.back();
			assignedCoordinatesToRedo.pop_back();
			assignment.coords(ri, assignedCoordinates[ri]);
			return ri;
		}
		if (done || assignedCoordinatesAvailable.empty()) return noAssignment;

		auto giveAssignment = [&]()->int {
			numberAssignmentsGiven++;
			auto ri = assignedCoordinatesAvailable.back();
			assignedCoordinatesAvailable.pop_back();
			assignedCoordinates[ri] = nextCoords;
			assignment.coords(ri, nextCoords);
			return ri;
		};

		// Repeat as often as desired
		if (nextCount > 0) {
			nextCount--;
			return giveAssignment();
		}
		
		// Bump the odometer
		nextCount = numberOfTries;
		for (int i = 0; i < bounds.size(); i++) {
			if ((nextCoords[i]+=divisionStep[i]) < divisionEnd[i]) {
				return giveAssignment();
			}
		FoundNext:
			nextCoords[i] = divisionBegin[i];
		}

		// Until all the results are in,  can not decide what to do when odometer overflows
		if (numberAssignmentsGiven < numberAssignmentsReported) {
			return noAssignment;
		}

		// Have received all the values at this resolution
		//
		bool foundAPeak = (bestResult - stats.mean()) > 0.5*stats.stdDev();
		std::cerr 
			<< "bestResult:" << bestResult
			<< " stats.mean():" << stats.mean()
			<< " stats.stdDev():" << stats.stdDev()
			<< (foundAPeak ? " foundAPeak" : " plateau")
			<< std::endl;

		if (foundAPeak) {
			// Bound the search
			for (int i = 0; i < bounds.size(); i++) {
				divisionBegin[i] = std::max(bounds[i].begin, bestCoords[i] - divisionStep[i]    );
				divisionEnd  [i] = std::min(bounds[i].end,   bestCoords[i] + divisionStep[i] + 1);
				nextCoords   [i] = divisionBegin[i];
			}
		}

		// Search at a higher resolution
		bool zoomed(false);
		for (int i = 0; i < bounds.size(); i++) {
			if (divisionStep[i] > 1) {
				divisionStep[i] /= 2;
				zoomed = true;
			}
		}

		if (!zoomed) {
			done = true;
			return noAssignment;
		}

		if (1)
		#pragma omp critical 
		{
			std::cerr << "Next divisions are " << std::endl;
			std::cout << "Next divisions are " << std::endl;
			for (int i = 0; i < divisionStep.size(); i++) {
				std::cerr << divisionBegin[i] << ".." << divisionEnd[i] << " step " << divisionStep[i] << std::endl;
				std::cout << divisionBegin[i] << ".." << divisionEnd[i] << " step " << divisionStep[i] << std::endl;
			}
		}

		stats.init();

		return giveAssignment();
	}

	static Internals* make(Bounds const & bounds) {
#include "../src/resmap/resmap_util_heap_undefs.h"
		return sNewA(Internals, (bounds));
#include "../src/resmap/resmap_util_heap_defs.h"
	}
};

PeakFinder::PeakFinder(Bounds const & bounds) : internals(Internals::make(bounds)) {
}

PeakFinder::~PeakFinder() {
	sDeleteConst(internals);
}

float PeakFinder::volume() const {
	float v = 1.0;
	for (auto & b : internals->bounds) v *= float(double(b.end) - double(b.begin));
	return v;
}

size_t PeakFinder::numberAssignmentsReported() const {
	return internals->numberAssignmentsReported;
}

PeakFinder::Coordinates const & PeakFinder::bestCoordinates() const {
	return internals->bestCoords;
}

float PeakFinder::bestResult() const {
	return internals->bestResult;
}

bool PeakFinder::hasFoundPeak() const {
	ScopedAcquire acquire(internals->lock, __FILE__, __LINE__);
	return internals->done;
}

PeakFinder::ResultIndex PeakFinder::askForAssignment(Assignment & assignment) {
	ScopedAcquire acquire(internals->lock, __FILE__, __LINE__);
	return internals->askForAssignment(assignment);
}

void PeakFinder::assignmentMustBeRedone(ResultIndex resultIndex) {
	ScopedAcquire acquire(internals->lock, __FILE__, __LINE__);
	internals->assignmentMustBeRedone(resultIndex);
}
void PeakFinder::reportResult(PeakFinder::ResultIndex resultIndex, float result) {
	ScopedAcquire acquire(internals->lock, __FILE__, __LINE__);
	internals->reportResult(resultIndex, result);
}
