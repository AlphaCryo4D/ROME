/***************************************************************************
 *
 * Dana-Farber Cancer Institute, Harvard Medical School and Peking University
 * Intel® Parallel Computing Center for Structural Biology
 *
 * Authors: "Bevin R Brett(bevin_brett@hotmail.com) 2012-09-21"
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
#pragma once

// This code is included into the MAP#DOptimizer_original and MAP#DOptimizer_new namespaces
// to make their tandem algorithms.  The function called by run() will differ depending
// on which one it is included in...

#define ALGORITHM(PHASE,pHASE) 																				\
class PHASE##Algorithm : public TandemExecution::Algorithm {												\
public:																										\
	PHASE##Algorithm(TandemExecution::AlgorithmsP algorithms) : TandemExecution::Algorithm(algorithms) {	\
	}																										\
																											\
	virtual void run() {																					\
		pHASE();																							\
	}																										\
};																											\
																											\
TandemExecution::AlgorithmP pHASE##Algorithm(TandemExecution::AlgorithmsP algorithms) {						\
	static PHASE##Algorithm* p;																				\
	if (!p) p = born(new PHASE##Algorithm(algorithms));														\
	return p;																								\
}																											\
// end of macro

#include "../src/resmap/resmap_util_heap_undefs.h"
ALGORITHM(Prepare,prepare)
ALGORITHM(Iterate,iterate)
#include "../src/resmap/resmap_util_heap_defs.h"
