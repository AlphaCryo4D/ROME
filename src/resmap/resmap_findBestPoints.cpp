/***************************************************************************
*
* Authors: Bevin Brett
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

#include "resmap_util.h"		// used for building precompiled headers on Windows

#include "./resmap_findBestPoints.h"
#include "resmap_time.h"



// Support comparing several different implementations easily
//
struct FindBestPoints::Impl {
	virtual ~Impl() {}

	virtual void insert(double w) = 0;
	virtual void prepareToPop() = 0;
	virtual double pop() = 0;

	int size;
	int weight_count;
	double target_weight;
	double sum_weight;

	Impl() : size(0), weight_count(0), target_weight(0.0), sum_weight(0.0) {}
	void setWeightCount(int weight_count, double target_weight, double sum_weight) {
		this->size			= 0;
		this->weight_count  = weight_count;
		this->target_weight = target_weight;
		this->sum_weight    = sum_weight;
		setWeightCountWkr();
	}

	virtual double * get_sorted_weights_ptr_if_avail() { return nullptr; }
	virtual void setWeightCountWkr() { }
};


// Implementation that sorts all the weights
//
struct FindBestPoints_Impl_Sort : public FindBestPoints::Impl {

	std::vector<double> sorted_weight;

	void setWeightCountWkr() {
		if (sorted_weight.size() < weight_count) {
			sorted_weight.resize(weight_count);
		}
	}

	double * get_sorted_weights_ptr_if_avail() { 
		return &sorted_weight[0]; 
	}

	void insert(double w) {
		sorted_weight[size] = w;
		size++;
	}

	void prepareToPop() {
		int info = LAPACKE_dlasrt('I', weight_count, &sorted_weight[0]);
		ERROR_CHECK(info<0, "the "+std::to_string((long long)-info)+"-th argument had an illegal value");
	}

	double pop() {
		return sorted_weight[--size];
	}
};



// Implementation that buckets the weights into logarithmically smaller buckets
// stopping when sure they won't be useful
//
struct FindBestPoints_Impl_Buckets : public FindBestPoints::Impl {

	struct Bucket {
		double hiWeight;
		int    capacity;
		int    weightIndex;
		int    size;
		double sumWeightsSoFar;
	};

	size_t              bucketsSize;
	std::vector<Bucket> buckets;
	std::vector<double> weights;

	FindBestPoints_Impl_Buckets() : buckets(32) {}
	~FindBestPoints_Impl_Buckets() {}

	virtual void setWeightCountWkr() {
		// Build the buckets
		int    weightIndex = 0;
		int    capacity    = 2;
		double hiWeight    = this->sum_weight;
		for (bucketsSize = 1; ; bucketsSize++) {
			if (buckets.size() == bucketsSize-1) {
				buckets.resize(2*bucketsSize);
			}
			auto& bucket = buckets[bucketsSize-1];
			bucket.hiWeight        = hiWeight;
			bucket.capacity        = capacity + 1 + capacity/128;	// allow for rounding errors
			bucket.weightIndex     = weightIndex;
			bucket.size            = 0;
			bucket.sumWeightsSoFar = 0.0;
			weightIndex += bucket.capacity;
			if (capacity >= weight_count) break;
			capacity *= 2;
			hiWeight *= 0.5;
		}
		if (weights.size() < weightIndex) {
			weights.resize(weightIndex);
		}
	}

	void insert(double w) {
		assert(w >= 0.0);
		// search for the bucket this weight belongs in, adjusting all the sum weights
		// and trimming buckets that this weight proves not needed
		int bucketsIndex = bucketsSize-1;
		for (; bucketsIndex >= 0; bucketsIndex--) {
			auto& bucket = buckets[bucketsIndex];
			bucket.sumWeightsSoFar += w;	
				// this is going into this or one of the lower-indexed (higher weight) buckets
			if (bucket.sumWeightsSoFar >= this->target_weight) {
				// discard the higher-indexed buckets because they are not needed
				if (bucketsSize > bucketsIndex+1) {
					if (false) {
						static int count = 0;
						if (count < 50 && count++ < 50) {
							std::cerr << "Shrunk bucketsSize from bucketsSize" << bucketsSize << " to " << bucketsIndex+1 << std::endl;
						}
					}
					bucketsSize = bucketsIndex+1;
				}
			}
			//
			// The buckets are in descending  0:[8..4] 1:[4..2] 2:[2..1] ...
			// so stop searching towards the lower indices (higher weights) when w is less than this bucket's hiWeight
			if (w <= bucket.hiWeight) break;
		}
		// found the bucket
		assert(bucketsIndex >= 0);	// bucket[0] has sum_weight as its low weight
		auto& bucket = buckets[bucketsIndex];
		weights[bucket.weightIndex + bucket.size++] = w;
	}

	void partiallySort(int begin, int end, double sumUptoBegin) {
		// Now that values are partially sorted into buckets, with each bucket holding values that are smaller than those in the lower indexed buckets
		// The sum of the values in the buckets is known, and we know that the sum of the values in all but the highest indexed bucket is less than or equal to the target weight
		// Now partially sort the values in the last bucket into a high value set and a low value set, 
		//		where the sum of the lower indexed bucketed values plus the sum of the high set is less than or equal to the target weight
		//		and adding the smallest value from the low set will push it over the target weight
		// This is done using a variant of quick sort that only sorts the partitions whose values are not certainly in the high set or the low set
		//
		if (begin + 1 >= end) return;
		double highSetSum = 0.0;
		// choose the partitioning weight
		auto w = (weights[begin] + weights[(begin + end)/2] + weights[end-1]) / 3.0;
		// partition
		int lo = begin, hi = end;
		for (;;) {
			// skip all the lo's that are in the right place
			while (lo < hi && weights[lo  ] <= w) {
				//highSetSum += weights[lo];						// sadly the high values are sorted into the lo indices...
				lo++;											// elements in [begin, lo) are less than or equal to w
			}
			// skip all the hi's that are in the right place	// either lo == hi or weights[lo] > w
            while (lo < hi && weights[hi-1] >  w) {highSetSum += weights[hi-1];hi--;}			// elements in [hi, end) are greater than w
			if (lo == hi) break;    							// either lo == hi or weights[hi - 1] <= w, which means hi - 1 can not be lo because weights[lo] > w 
			std::swap(weights[lo], weights[hi-1]);
		}
		// only partiallySort the partition that contains the element needed to exceed the target_weight, the others don't need further ordering
		if (sumUptoBegin + highSetSum > target_weight) partiallySort(lo, 	end,  sumUptoBegin);
		else                                           partiallySort(begin, lo, sumUptoBegin + highSetSum);
	}

	void prepareToPop() {
		// only really have to sort the last bucket since the values in the others are definitely needed
		// sort all the buckets to getting the same pops as the other algorithms
		for (int bucketsIndex = 
//#ifdef NDEBUG
			bucketsSize-1;
//#else
//			0; 
//#endif
			bucketsIndex < bucketsSize; bucketsIndex++) {
			auto& bucket = buckets[bucketsIndex];
			if (bucket.size > 0) {
				auto doSort = [&]() {
					//std::sort(&weights[bucket.weightIndex + 0], &weights[bucket.weightIndex + bucket.size]);
					partiallySort(bucket.weightIndex, bucket.weightIndex + bucket.size, bucketsIndex == 0 ? 0.0 : buckets[bucketsIndex-1].sumWeightsSoFar);
				};
				if (true) {
					doSort();
				} else {
					#pragma omp critical
					std::cerr << __FILE__ << ":" << __LINE__ << " prepareToPop sorting " << bucket.size << " items" << std::endl;
					AccurateTimer at;
					doSort();
					if (bucket.size > 10000) {
						std::cerr << "    Sorting " << bucket.size << " items took "
							<< std::setprecision(2) << at.sinceInited()/at.countPerMicrosecond() * 1e-6 << " seconds"
							<< std::endl;
					}
				}
			}
		}
		popBucketsIndex = 0;
		popWeightIndex = 0;
	}

	int popBucketsIndex;
	int popWeightIndex;

	double pop() {
		for (;popBucketsIndex < bucketsSize; popBucketsIndex++) {
			auto& bucket = buckets[popBucketsIndex];
			if (popWeightIndex < bucket.size) {
				return weights[bucket.weightIndex + bucket.size - ++popWeightIndex];
			}
			popWeightIndex = 0;
		}
		assert(false);	// asked for too many
		return 0.0;
	}
};


// Implementations that heap the weights
//
static void heapInsert(std::vector<double>& heap, int& size, double w) {
	//	0			i
	//	1	2		2i+1  2i+2
	//	3 4 5 6
	//
	int i = size;
	while (i > 0) {
		auto parentIndex = (i-1)/2;
		auto& parent = heap[parentIndex];
		if (w <= parent) break;
		heap[i] = parent;
		i = parentIndex;
	}
	heap[i] = w;
	size++;
}

struct FindBestPoints_Impl_HeapBase : public FindBestPoints::Impl {
	FindBestPoints_Impl_HeapBase() {}
	~FindBestPoints_Impl_HeapBase() {}

	std::vector<double> heap_weight;

	virtual void setWeightCountWkr() {
		if (heap_weight.size() < weight_count) {
			heap_weight.resize(weight_count);
		}
	}

	void insert(double w) {
		heapInsert(heap_weight, size, w);
	}

	void check_heap_weight() {
		for (auto i = 1; i < size; i++) {
			auto parent = (i-1)/2;
			assert(heap_weight[i] < heap_weight[parent]);
		}
	}
};

struct FindBestPoints_Impl_Heap : public FindBestPoints_Impl_HeapBase {

	void prepareToPop() {}

	double pop() {
		size--;
		auto result = heap_weight[0];

		// extracting the root is done by replacing it by the greatest child and extracting that child from its subtree
		int i = 0;
		for (;;) {
			// promote the greater child
			auto greaterChildIndex = i*2+1;		if (greaterChildIndex >= weight_count) break;
			auto greaterChild      = heap_weight[greaterChildIndex];
			auto otherChildIndex   = greaterChildIndex+1;
			if (otherChildIndex < weight_count && heap_weight[otherChildIndex] > greaterChild) { greaterChildIndex = otherChildIndex; greaterChild = heap_weight[otherChildIndex]; }
			heap_weight[i] = greaterChild;
			i = greaterChildIndex;
			if (greaterChild == 0.0) break;		// This effectively deletes the deeper parts of the heap, speeding up many extracts case
		}
		heap_weight[i] = 0.0;					// This effectively deletes the deeper parts of the heap, speeding up many extracts case

		if (0) check_heap_weight();

		return result;
	}

};


// Implementation that heaps all the weights
// but the runs a frontier through the heap
// Never seems to win - sigh
//
struct FindBestPoints_Impl_HeapFrontier : public FindBestPoints_Impl_HeapBase {
	FindBestPoints_Impl_HeapFrontier() {}
	~FindBestPoints_Impl_HeapFrontier() {}

	std::vector<int>    heap_frontier;
	int					heap_frontier_count;

	virtual void setWeightCountWkr() {
		FindBestPoints_Impl_HeapBase::setWeightCountWkr();
		if (heap_frontier.size() < weight_count) {
			heap_frontier.resize(weight_count);
		}
	}

	void prepareToPop() {
		check_heap_weight();
		heap_frontier_count = 0;
		if (weight_count > 0) {
			heap_frontier_count = 1;
			heap_frontier[0] = 0;		// The root is the frontier
		}
		check_heap_frontier();
	}

	double pop() {
		size--;
		assert(heap_frontier_count > 0);
		auto biggestInFrontier = heap_frontier[0];

		// replace the frontier root with its first child, and push it down to the right place
		// if there is no first child, use the last in the frontier
		//
		auto toBeInserted = biggestInFrontier*2 + 1;
		if (toBeInserted >= weight_count) {
			if (--heap_frontier_count == 0) return heap_weight[biggestInFrontier];	// the frontier is now empty
			toBeInserted = heap_frontier[heap_frontier_count];
		}
		auto toBeInsertedWeight = heap_weight[toBeInserted];

		int i = 0;
		for (;;) {
			// find the greater child
			auto greaterChildIndex  = i*2+1;		if (greaterChildIndex >= heap_frontier_count) break;
			auto greaterChild       = heap_frontier[greaterChildIndex];
			auto greaterChildWeight = heap_weight[greaterChild];
			auto otherChildIndex    = i*2+2;
			if (otherChildIndex < heap_frontier_count) {
				auto otherChild = heap_frontier[otherChildIndex];
				auto otherChildWeight = heap_weight[otherChild];
				if (otherChildWeight > greaterChildWeight) {
					greaterChildIndex  = otherChildIndex;
					greaterChild       = otherChild;
					greaterChildWeight = otherChildWeight;
				}
			}
			// promote the greater child if the toBeInserted can't go here 
			if (greaterChildWeight < toBeInsertedWeight) break;
			heap_frontier[i] = greaterChild;
			i = greaterChildIndex;
		}
		heap_frontier[i] = toBeInserted;

		if (0) check_heap_frontier();

		toBeInserted = biggestInFrontier*2 + 2;
		if (toBeInserted < weight_count) {
			toBeInsertedWeight = heap_weight[toBeInserted];

			auto i = heap_frontier_count++;
			while (i > 0) {
				auto parentIndex = (i-1)/2;
				auto parent = heap_frontier[parentIndex];
				auto parentWeight = heap_weight[parent];
				if (toBeInsertedWeight <= parentWeight) break;
				heap_frontier[i] = parent;
				i = parentIndex;
			}
			heap_frontier[i] = toBeInserted;
		}

		if (0) check_heap_frontier();

		// return the one that was removed
		return heap_weight[biggestInFrontier];
	}

	void check_heap_frontier() {
		for (auto i = 1; i < heap_frontier_count; i++) {
			auto parent = (i-1)/2;
			assert(heap_weight[heap_frontier[i]] < heap_weight[heap_frontier[parent]]);
		}
	}

};


// Connect the external view to the implementation
//
static const int ImplToUse_begin = 0;
static const int ImplToUse_end = 4;
static int implToUse      = 1;				// Try the buckets approach
FindBestPoints::Impl* makeImpl() {
	switch (implToUse) {
	case 0: return sNew(FindBestPoints_Impl_Sort);
	case 1: return sNew(FindBestPoints_Impl_Buckets);
	case 2: return sNew(FindBestPoints_Impl_Heap);
	case 3: return sNew(FindBestPoints_Impl_HeapFrontier);
	}
	return nullptr;
}

FindBestPoints::FindBestPoints() : _impl(makeImpl()) {}

FindBestPoints::~FindBestPoints() {
	sDelete(_impl); _impl = nullptr;
}

void FindBestPoints::setWeightCount(int to, double target_weight, double sum_weight) {
	_impl->setWeightCount(to, target_weight, sum_weight);
}

int FindBestPoints::get_weight_count() const { 
	return _impl->weight_count; 
}

double * FindBestPoints::get_sorted_weights_ptr_if_avail() { 
	return _impl->get_sorted_weights_ptr_if_avail(); 
}

void FindBestPoints::insert(double w) {
	_impl->insert(w);
}

void FindBestPoints::prepareToPop() {
	_impl->prepareToPop();

}

double FindBestPoints::pop() {
	return _impl->pop();
}


// Test it
//
#include "./time.h"

void findBestPointsUnitTestUnitTest(FindBestPoints & findBestPoints, double& rightAnswerOrNegative) {

	static const int weight_count = 0xfffff;
	static const int weight_xor   = 0x5e5e5;

	double sum_weight = 0.0;
	double target_weight = -1.0;
	for (int pass = 0; pass < 2; pass++) {
		if (pass == 1) findBestPoints.setWeightCount(weight_count, target_weight, sum_weight);
		for (int i = 0; i < weight_count; i++) {
			auto x = double(i ^ weight_xor);
			double w = std::abs(1+x+x*x+x*x*x+x*x*x*x);
			if (pass == 0)	sum_weight += w;
			else			findBestPoints.insert(w);
		}
		if (pass == 0) target_weight = sum_weight * 0.1;
	}
	findBestPoints.prepareToPop();
	double soFar = 0.0;
	for (int i = 0; i < weight_count; i++) {
		auto result = findBestPoints.pop();
		//	when put ascending numbers in as the w		assert(result == double(limit - i));
		soFar += result;
		if (soFar >= target_weight) {
			if (rightAnswerOrNegative < 0.0) {
				rightAnswerOrNegative = result;
			} else {
				std::cerr << "Had to pop element [" << i << "] to exceed target_weight" << std::endl;
				if (rightAnswerOrNegative != result) {
					std::cerr << "Different answer, rightAnswerOrNegative:" << rightAnswerOrNegative << " result:" << result << std::endl;
				}
			}
			break;
		}
	}
}

void findBestPointsUnitTestUnitTest() {

	std::cerr << __FILE__ << ":" << __LINE__ 
		<< "FindBestPointsUnitTest" << std::endl;

	double rightAnswerOrNegative(-1.0);

	if (false) {	// used for correctness debugging when the following loop detects an error
		implToUse = 3;
		FindBestPoints findBestPoints;
		findBestPointsUnitTestUnitTest(findBestPoints, rightAnswerOrNegative);
	}

	for (implToUse = ImplToUse_begin; implToUse < ImplToUse_end; implToUse++) {
		FindBestPoints findBestPoints;
		double bestElapsed = std::numeric_limits<double>::max();
		for (int trial = 0; trial < 1; trial++) {
			AccurateTimer accurateTimer;
			findBestPointsUnitTestUnitTest(findBestPoints, rightAnswerOrNegative);
			auto elapsed = accurateTimer.sinceInited() / accurateTimer.countPerMicrosecond();
			bestElapsed = std::min(bestElapsed, elapsed);
		}
		std::cerr << "Algorithm " << implToUse << ": elapsed " << bestElapsed << std::endl;
	}

	EXIT_ABNORMALLY;
}

#if 0
class FindBestPointsUnitTestUnitTest {
public:
	FindBestPointsUnitTestUnitTest() {
		findBestPointsUnitTestUnitTest();
	}
} unitTest_for_FindBestPointsUnitTestUnitTest;
#endif