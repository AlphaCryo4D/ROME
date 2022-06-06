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

#include "./util_heap_undefs.h"
ALGORITHM(Prepare,prepare)
ALGORITHM(Iterate,iterate)
#include "./util_heap_defs.h"
