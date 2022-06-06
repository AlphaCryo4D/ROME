/***************************************************************************
 *
 * Authors: "Yongbei(Glow) Ma,Jiayi (Timmy) Wu, Youdong (Jack) Mao"
 * Dana-Farber Cancer Institute, Harvard Medical School and Peking University
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

#include "resmap_exp_model.h"

IntPerformanceCounter FourierShellTranslation::transformFinePerformanceCounter("transformFine");
IntPerformanceCounter FourierShellTranslation::transformCoarsePerformanceCounter("transformCoarse");


 // rearrange the image array to resolution increase array....
void FourierShellTranslation::setCoarseImageBoundary(const int* _Mresol_coarse, int _size)
{
	changeUid();	// helps the Shifter (below) detect that the deferring it is doing is wrong
					// because one of the inputs to the deferred or redone calculation has changed

	int outSideResCount = 0;
	oldFsize2Coarse = _size; newFsize2Coarse = 0;
	for (int i = 0; i < oldFsize2Coarse; i++) {
		if (_Mresol_coarse[i] <= 0) outSideResCount++;
		else newFsize2Coarse++;
	}
	assert(outSideResCount + newFsize2Coarse == oldFsize2Coarse);
	nIndexCoarse = (int*)aMalloc(sizeof(int)*newFsize2Coarse, 64);
	// std::cout<<" "<<newFsize2Coarse<<" "<<outSideResCount<<" "<<oldFsize2Coarse<<std::endl;
	int index = 0;
	for (int i = 0; i < oldFsize2Coarse; i++) {
		if (_Mresol_coarse[i]>0) {
			nIndexCoarse[index] = i;
			index++;
		}
	}
	assert(index == newFsize2Coarse);
	assert(index + outSideResCount == oldFsize2Coarse);
}

void FourierShellTranslation::setNormCorrectionBoundary(const int* _Mresol_fine, int _size)
{
	changeUid();

	min_iRes = (std::numeric_limits<int>::max)();
	max_iRes = (std::numeric_limits<int>::min)();
	int outSideResCount = 0;
	oldFsize2Fine = _size; newFsize2Fine = 0;
	for (int i = 0; i < oldFsize2Fine; i++) {
		if (_Mresol_fine[i] == -1) {
			outSideResCount++;
		}
		else {
			newFsize2Fine++;
			if (_Mresol_fine[i]<min_iRes) min_iRes = _Mresol_fine[i];
			if (_Mresol_fine[i]>max_iRes) max_iRes = _Mresol_fine[i];
		}
	}
	assert(outSideResCount + newFsize2Fine == oldFsize2Fine);
	//
	nIndexFine = (int*)aMalloc(sizeof(int)*newFsize2Fine, 64);
	iResFine = (int*)aMalloc(sizeof(int)*newFsize2Fine, 64);
	// std::cout<<min_iRes<<" "<<max_iRes<<" "<<newFineFsize2<<" "<<outSideResCount<<" "<<oldFineFsize2<<std::endl;
	int index = 0;
	for (int ires = min_iRes; ires <= max_iRes; ires++) {
		for (int i = 0; i < oldFsize2Fine; i++) {
			if (_Mresol_fine[i] == ires) {
				nIndexFine[index] = i;
				iResFine[index] = ires;
				index++;
			}
		}
	}
	assert(index == newFsize2Fine);
	assert(index + outSideResCount == oldFsize2Fine);
	//
	normCorrLo = 0;
	normCorrHi = newFsize2Fine;
}

void FourierShellTranslation::setScaleCorrectionBoundary(const int* _Mresol_fine,const VectorOfArray2d<FDOUBLE>& data_vs_prior_class)
{
	changeUid();

	int nr_class = data_vs_prior_class.size();
	scaleCorrThresholdRes.resize(nr_class);
	scaleCorrLo.resize(nr_class);
	scaleCorrHi.resize(nr_class);
    scaleCorrExcludingIndex.resize(nr_class,nullptr);
	for (int iclass = 0; iclass < nr_class; iclass++)
	{
		scaleCorrThresholdRes[iclass] = -1;
		int thresholdResCount = 0;

		for (int ires = min_iRes; ires <= max_iRes; ires++) {
			if (data_vs_prior_class[iclass][ires] > 3.) {
				scaleCorrThresholdRes[iclass] = ires;
				thresholdResCount++;
			}
		}
		assert(thresholdResCount>0);

        auto print_data_vs_prior_class = [&]() {
            std::cout << "data_vs_prior_class [ ";
            for (int ires = min_iRes; ires <= max_iRes; ires++)
                std::cout << data_vs_prior_class[iclass][ires] << ", ";
            std::cout << " ]"<<std::endl;
        };
        
        if (thresholdResCount == 0) {
            MASTERNODE std::cout << "WARNING : " << __FILE__ << " " << __LINE__ << " data_va_prior has no larger than 3!!!" << std::endl;
            MASTERNODE print_data_vs_prior_class();
        }
        bool is_continuous = true;
        for (int ires = min_iRes; ires <= scaleCorrThresholdRes[iclass]; ires++) {
            if (data_vs_prior_class[iclass][ires] <= 3.) {
                is_continuous = false;
            }
        }
        
        if (!is_continuous)
        {
            MASTERNODE std::cout << "WARNING : " << __FILE__ << " " << __LINE__ << " data_va_prior is not continuously decreasing!!!" << std::endl;
            MASTERNODE print_data_vs_prior_class();
            //std::cout << "WARNING : " << " RELION excluding <=3 points between highest~3,but ROME not!!! " <<std::endl;
            MASTERNODE std::cout << "WARNING : " << " you can remove '-scale' option to prevent scale correction!!!" <<std::endl;
            MASTERNODE std::cout << "WARNING : " << " If you are doing autorefine,ignore this warning!!!"<<std::endl;
        }
        
        if (!is_continuous && nr_class==1)
        {
            // TODO : Solution
            // 1) add excluding index if data_vs_prior is not monotone decreasing
            // scaleCorrExcludingIndex[iclass] = (int*)aMalloc(sizeof(int)*newFsize2Fine, 64);
            // TODO
            // 2) can reorder the nIndexFine if nr_classes==1(do autorefine)
            int index = 0;
            std::vector<int> ires_temp;
            auto addIres = [&](int& ires) {
                for (int i = 0; i < oldFsize2Fine; i++)
                    if (_Mresol_fine[i] == ires) {
                        nIndexFine[index] = i;iResFine[index] = ires;index++;
                    }
            };
            for (int ires = min_iRes; ires <= max_iRes; ires++) {
                if (data_vs_prior_class[iclass][ires] > 3.) addIres(ires);
                else ires_temp.push_back(ires);
            }
            for (auto& ires : ires_temp) {
                assert(data_vs_prior_class[iclass][ires] <= 3.);
                addIres(ires);
            }
            assert(index == newFsize2Fine);
        }
		//
		scaleCorrLo[iclass] = 0;
        bool has_found = false;
		for (int i = 0; i < newFsize2Fine; i++) {
            // reorder the ires,so only pick ires < scaleCorrThresholdRes
            if (nr_class==1 && !is_continuous && has_found && iResFine[i] < scaleCorrThresholdRes[iclass]) break;
            // need to include all ires < 3 between (ires=0)~(ires=last_ires>3)
			if (iResFine[i] <= scaleCorrThresholdRes[iclass]) {
				scaleCorrHi[iclass] = i + 1;
                if (iResFine[i] == scaleCorrThresholdRes[iclass]) has_found = true;
			}
		}
	}
}

void FourierShellTranslation::clear()
{
	changeUid();

	if (nIndexFine) { aFree(nIndexFine); nIndexFine = nullptr; }
	if (iResFine) { aFree(iResFine); iResFine = nullptr; }
	if (nIndexCoarse) { aFree(nIndexCoarse); nIndexCoarse = nullptr; }
	newFsize2Fine = 0; oldFsize2Fine = 0;
	newFsize2Coarse = 0; oldFsize2Coarse = 0;
	min_iRes = max_iRes = -1;
	normCorrLo = normCorrHi = -1;
	scaleCorrThresholdRes.resize(0);
	scaleCorrLo.resize(0);
	scaleCorrHi.resize(0);
    //
    for(auto & v : scaleCorrExcludingIndex){if(v) aFree(v);v=nullptr;}
    scaleCorrExcludingIndex.resize(0);
	for (int thread = 0; thread < TempDataThread.size(); thread++)
		aFree(TempDataThread[thread]);
	TempDataThread.clear();
}

void doTransformNow(
	const ParticleModel::State  & desiredState,
	ParticleModel::Shifter* shifter) {

	const int tid = omp_get_thread_num();

	auto & currentState = shifter->currentState;

	currentState.fourierShellTrans = desiredState.fourierShellTrans;
	currentState.fourierShellTrans_uid = desiredState.fourierShellTrans_uid;
	currentState.transformToCoarse = desiredState.transformToCoarse;
	currentState.exp_current_Fsize2 = desiredState.exp_current_Fsize2;

	auto & fourierShellTrans = *currentState.fourierShellTrans;
	auto exp_current_Fsize2 = currentState.exp_current_Fsize2;

#if defined(DEBUG_SHIFTER_POOL)
	if (desiredState.debugTransformed) 
#pragma omp critical
	{
		// Check that the inputs are the same
		//
		Checksum<int>    fourierShellTrans_checksum_coarse, fourierShellTrans_checksum_fine;
		fourierShellTrans.checksum(fourierShellTrans_checksum_coarse, fourierShellTrans_checksum_fine);
		if (!(currentState.transformToCoarse
			? fourierShellTrans_checksum_coarse == desiredState.fourierShellTrans_checksum_coarse
			: fourierShellTrans_checksum_fine == desiredState.fourierShellTrans_checksum_fine)) {
			EXIT_ABNORMALLY;
		}

		Checksum<double> aTable_checksum_before;
		aTable_checksum_before.init(shifter->aTable_wptr(), shifter->tableSize());
		if (aTable_checksum_before == desiredState.aTable_checksum_before) {
			std::cerr << tid << ": Got the right inputs for the debugTransformed" << std::endl;
			std::cout << tid << ": Got the right inputs for the debugTransformed" << std::endl;
		} else {
			auto p = [&](std::ostream& os) {
				os << tid << ": Got the WRONG inputs for the debugTransformed" << std::endl;
				os << "aTable_checksum_before" << std::endl;
				aTable_checksum_before.print(os);
				os << "desiredState.aTable_checksum_before" << std::endl;
				desiredState.aTable_checksum_before.print(os);
				os << "    DIFF   current : desired" << std::endl;
				aTable_checksum_before.printDiff(os, desiredState.aTable_checksum_before);
				os << std::endl;
			};
			p(std::cerr);
			p(std::cout);
	    }
	}
#endif

	if (currentState.transformToCoarse) {
		fourierShellTrans.transformCoarse(tid, shifter->aTable_wptr(), exp_current_Fsize2);
		fourierShellTrans.transformCoarse(tid, shifter->bTable_wptr(), exp_current_Fsize2);
	}
	else {
		fourierShellTrans.transformFine(tid, shifter->aTable_wptr(), exp_current_Fsize2);
		fourierShellTrans.transformFine(tid, shifter->bTable_wptr(), exp_current_Fsize2);
	}
	shifter->currentState.status = ParticleModel::State::transformDone;

#if defined(DEBUG_SHIFTER_POOL)
	if (desiredState.debugTransformed) 
#pragma omp critical
	{
		Checksum<double> aTable_checksum_after;
		aTable_checksum_after.init(shifter->aTable_wptr(), shifter->tableSize());

		if (aTable_checksum_after == desiredState.aTable_checksum_after) {
			std::cerr << tid << ": aTable_checksum_after == desiredState.aTable_checksum_after" << std::endl;
			std::cout << tid << ": aTable_checksum_after == desiredState.aTable_checksum_after" << std::endl;
		} else {
			auto p = [&](std::ostream& os) {
				os << tid << ": Got the WRONG outputs from the debugTransformed" << std::endl;
				os << "aTable_checksum_after" << std::endl;
				aTable_checksum_after.print(os);
				os << "desiredState.aTable_checksum_after" << std::endl;
				desiredState.aTable_checksum_after.print(os);
				os << "    DIFF   current : desired" << std::endl;
				aTable_checksum_after.printDiff(os, desiredState.aTable_checksum_after);
				os << std::endl;
			};
			p(std::cerr);
			p(std::cout);
		}
	}
#endif
}

static IntPerformanceCounter shifterPoolAcquirePeak         ("shifterPoolPeak");
static IntPerformanceCounter shifterPoolAcquireRightReusable("shifterPoolAcquireRightReusable");
static IntPerformanceCounter shifterPoolAcquireWrongReusable("shifterPoolAcquireWrongReusable");

class ParticleModel::ShifterPool : public NoCopy {
	// Idea:	Reuse the cache memory needed for the shifters by recomputing them rather than saving them
	//			Keep only as many as will fit in the cache, along with the other things that are needed
	//			Of course, all the ones in use must be kept
	//			Each thread has its own free list
	//			For now, requires the acquiring and releasing thread be the same (checked)
	//
	class PerThread : public NoCopy {
		static const size_t desired = 1024*1024;		// BEVIN TBD set the right number
		ShifterPool* parent;
		int			 made;
		int			 inUseCount;
		Shifter*	 reusable;
		std::vector< ParticleModel::Shifter* > v;
	public:
		PerThread() : parent(nullptr), made(0), inUseCount(0), reusable(nullptr) { }
		~PerThread() {
			if (inUseCount > 0) EXIT_ABNORMALLY;
			while (auto p = reusable) {
				reusable = p->freeNext;
				sDelete(p);
			}
		}
		void setParent(ShifterPool* parent) { 
			this->parent = parent;  
			resize();
		}
		Shifter* acquire(int i) {
			Shifter* shifter = v[i];
			if (shifter) {
				if (shifter->acquiredCount == 0) {
					shifterPoolAcquireRightReusable.count.v++;
					// remove from the list
					if (reusable == shifter) reusable = shifter->freeNext;
					if (!!shifter->freeNext) shifter->freeNext->freePrev = shifter->freePrev;
					if (!!shifter->freePrev) shifter->freePrev->freeNext = shifter->freeNext;
					shifter->freeNext = shifter->freePrev = nullptr;
				} else {
					// assert not in the list
					assert(reusable != shifter);
					assert(!shifter->freeNext && !shifter->freePrev);
				}
			} else {
				shifter = (made >= desired) ? reusable : nullptr;
				if (shifter) {
					shifterPoolAcquireWrongReusable.count.v++;
					// The indexes are revisited smallest first, so reuse the latest freed
					assert(shifter->acquiredCount == 0);
					v[shifter->index] = nullptr;
					reusable = shifter->freeNext;
					if (reusable) reusable->freePrev = nullptr;
					shifter->freeNext = shifter->freePrev = nullptr;
				} else {
					made++;
					shifter = sNew(Shifter);
					while (made > shifterPoolAcquirePeak.count.v) shifterPoolAcquirePeak.count.v = made;
				}
				shifter->index = i;
				v[i] = shifter;
			}
			if (shifter->acquiredCount++ == 0) inUseCount++;
			return shifter;
		}
		void release(Shifter* shifter) {
			auto i = shifter->index;
			assert(shifter == v[i]);
			if (shifter->acquiredCount-- == 0)
				EXIT_ABNORMALLY;
			if (shifter->acquiredCount > 0) return;
			inUseCount--;
			shifter->freeNext = reusable;
			assert(!shifter->freePrev);
			if (reusable) reusable->freePrev = shifter;
			reusable = shifter;
		}
		void resize() {
			if (inUseCount > 0) EXIT_ABNORMALLY;
			v.resize(parent->size);
			for (auto & sp : v) if (sp != nullptr) EXIT_ABNORMALLY;
		}
	};

	VectorOfNoCopyable<PerThread> forEachThread;
	int size;
public:
	ShifterPool() : forEachThread(omp_get_max_threads()), size(0) {
		for (size_t i = 0; i < forEachThread.size(); i++) forEachThread[i].setParent(this);
	}

	void resize(int size) {		// size is the number of i, not the size of the pool
		this->size = size;
		for (size_t i = 0; i < forEachThread.size(); i++) {
			forEachThread[i].resize();
		}
	}

	void release(Shifter* shifter) {
		forEachThread[omp_get_thread_num()].release(shifter);
	}

	Shifter* acquire(int i, State const & desiredState) {
		auto shifter = forEachThread[omp_get_thread_num()].acquire(i);

		auto & currentState = shifter->currentState;

#if defined(DEBUG_SHIFTER_POOL)
		bool interesting = false;

		if (desiredState.debugTransformed) {
			static std::atomic<int> count;
			if (count++ < 10) interesting = true;
		}
		if (interesting) {
			std::cerr << "ParticleModel::ShifterPool::acquire the debugTransformed" << std::endl;
			std::cout << "ParticleModel::ShifterPool::acquire the debugTransformed" << std::endl;
		}
#endif

		if (currentState == desiredState)  {
#if defined(DEBUG_SHIFTER_POOL)
			if (interesting) {
				std::cerr << "no work needed to acquire debugTransformed index:" << i << std::endl;
				std::cout << "no work needed to acquire debugTransformed index:" << i << std::endl;
			}
#endif
			return shifter;
		}


		if (currentState.status == State::transformDone) {
			// Beware - moving this AFTER the next "if" gets the wrong answers - TBD why
			// Because the wrong contents are now in the table
			shifter->currentState.status = State::initTableDone;
		}

		if (currentState.status > desiredState.status) currentState.status = desiredState.status;

		assert(desiredState.status >= State::initTableDone);
		if (currentState.status < State::initTableDone ||
			currentState.fine_Fsize2 != desiredState.fine_Fsize2) {

			currentState.fine_Fsize2  = desiredState.fine_Fsize2;
			shifter->setNeededTableCapacity(currentState.fine_Fsize2);
			currentState.status = State::initTableDone;
		}

		assert(desiredState.status >= State::initDone);
		if (currentState.status < State::initDone
			|| currentState.current_size !=	desiredState.current_size
			|| currentState.shiftx       !=	desiredState.shiftx      
			|| currentState.shifty       !=	desiredState.shifty      
			|| currentState.ori_size     != desiredState.ori_size) {

			assert(desiredState.status >= State::initDone);

			currentState.current_size = desiredState.current_size;
			currentState.shiftx       = desiredState.shiftx;
			currentState.shifty       = desiredState.shifty;
			currentState.ori_size     = desiredState.ori_size;
			shifter->init(currentState.current_size, currentState.shiftx, currentState.shifty, currentState.ori_size);
			shifter->currentState.status = State::initDone;
		}

		if (desiredState.status >= State::transformDone
			&& (currentState.status < State::transformDone
				|| currentState.fourierShellTrans != desiredState.fourierShellTrans
				|| currentState.fourierShellTrans_uid != desiredState.fourierShellTrans_uid
				|| currentState.transformToCoarse != desiredState.transformToCoarse
				|| currentState.exp_current_Fsize2 != desiredState.exp_current_Fsize2)) {

            doTransformNow(desiredState, shifter);
		}

		assert(shifter->currentState == desiredState);
		return shifter;
	}

};

ParticleModel::ShiftImageAssistor::ShiftImageAssistor(ParticleModel & parent) 
  : max_dim0(0), max_dim1(0), dim0(0), dim1(0), desiredStates(nullptr),
	shifterPool(sNew(ShifterPool))
{
}

ParticleModel::ShiftImageAssistor::~ShiftImageAssistor() {
	vDelete(desiredStates);
	sDeleteConst(shifterPool);
}

void ParticleModel::ShiftImageAssistor::setMaxDims(int max_exp_nr_trans, int max_exp_nr_over_trans) {
	max_dim0 = max_exp_nr_trans;
	max_dim1 = max_exp_nr_over_trans;
	if (!!desiredStates) {
		if (max_exp_nr_trans != 0 || max_exp_nr_over_trans != 0) EXIT_ABNORMALLY;
		vDelete(desiredStates);
	} else {
		desiredStates = vNew(State, max_dim0*max_dim1);
	}
	shifterPool->resize(max_dim0*max_dim1);	// wrong, but will fix later
}

void ParticleModel::ShiftImageAssistor::initTable(int fine_Fsize2) {
	if (omp_get_thread_num() != 0) EXIT_ABNORMALLY;
	if (omp_get_num_threads() != 1) EXIT_ABNORMALLY;
	for (int i = 0; i < max_dim0*max_dim1; i++) {
		auto & state = desiredStates[i];
		state.status = State::initTableDone;
		state.fine_Fsize2 = fine_Fsize2;
	}
}

void ParticleModel::ShiftImageAssistor::setCurrDims(int exp_nr_trans, int exp_nr_over_trans) {
	dim0 = exp_nr_trans;
	dim1 = exp_nr_over_trans;
	assert(dim0 <= max_dim0);
	assert(dim1 <= max_dim1);
}
void ParticleModel::ShiftImageAssistor::setInitArgs(int itrans, int iover_trans, int itrans_over_trans, int current_size, double shiftx, double shifty, int ori_size) {
	auto i = index(itrans, iover_trans, itrans_over_trans);
	auto & desiredState = desiredStates[i];
	if (omp_get_thread_num() != 0) EXIT_ABNORMALLY;
	if (omp_get_num_threads() != 1) EXIT_ABNORMALLY;
#if defined(DEBUG_SHIFTER_POOL)
	if (desiredState.debugTransformed)
		std::cerr << "setInitArgs on one that has been reverted" << std::endl;
#endif
	desiredState.status = State::initDone;
	desiredState.current_size	= current_size; 
	desiredState.shiftx			= shiftx; 
	desiredState.shifty			= shifty; 
	desiredState.ori_size		= ori_size;
}

static IntPerformanceCounter ShiftImageAssistor_transformCounter("ShiftImageAssistor_transformCounter");

void ParticleModel::ShiftImageAssistor::transform(int itrans, int iover_trans, int itrans_over_trans,

	const int exp_current_Fsize2,
	FourierShellTranslation const & fourierShellTrans,
	bool isCoarse) {

	ShiftImageAssistor_transformCounter.count.v++;

	auto i = index(itrans, iover_trans, itrans_over_trans);
	auto & desiredState = desiredStates[i];

	assert(desiredState.status >= State::initDone);
	auto shifter = acquireShifterW(itrans, iover_trans, itrans_over_trans);
	assert(shifter->currentState.status >= State::initDone);

#if defined(DEBUG_SHIFTER_POOL)
	// The debug version compares the checksums of some of them to see if this recomputation is different
	//
#pragma omp critical
	{
		static int count;
		static int nextInteresting = 1;
		if (++count == nextInteresting) {
			desiredState.debugTransformed = true;
			nextInteresting *= 2;
		}
	}

	if (desiredState.debugTransformed) {

		// Before changing the desired state
		// get the acquire the current state, checksum it
		// do the transform now to a copy, and checksum that result also
		//
		std::cerr << "Checksumming a " 
			<< ((desiredState.status > State::initDone) ? "already transformed" : "initDone")
			<< " one" << std::endl;
		std::cout << "Checksumming a "
			<< ((desiredState.status > State::initDone) ? "already transformed" : "initDone")
			<< " one" << std::endl;

		const int  tid		   = omp_get_thread_num();
		const int  tableSize   = shifter->tableSize();
		const auto aTable_wptr = shifter->aTable_wptr();

		// checksum the inputs
		//
		desiredState.aTable_checksum_before.init(aTable_wptr, tableSize);
		fourierShellTrans.checksum(
			desiredState.fourierShellTrans_checksum_coarse,
			desiredState.fourierShellTrans_checksum_fine);

		std::vector<double> copy(tableSize);
		for (int i = 0; i < tableSize; i++) copy[i] = aTable_wptr[i];

		if (isCoarse) {
			fourierShellTrans.transformCoarse(tid, copy.data(), exp_current_Fsize2);
		} else {
			fourierShellTrans.transformFine(tid, copy.data(), exp_current_Fsize2);
		}

		// checksum the expected output
		//
		desiredState.aTable_checksum_after.init(copy.data(), tableSize);
	}
#endif

	desiredState.status = State::transformDone;
	desiredState.fourierShellTrans     = &fourierShellTrans;
	desiredState.fourierShellTrans_uid = fourierShellTrans.uid();
	desiredState.transformToCoarse     = isCoarse;
	desiredState.exp_current_Fsize2    = exp_current_Fsize2;

	// By changing the currentState back to initTableDone,
	// the init step will be done again in the acquire, 
	// and the transform step also
	//
	// The debug version compares the checksums to see if this recomputation is different
	//
	shifter->currentState.status = State::initTableDone;

#if defined(DEBUG_SHIFTER_POOL)
	if (desiredState.debugTransformed) {
		std::cerr << "transform does: shifter->currentState.status = State::initTableDone" << std::endl;
		std::cout << "transform does: shifter->currentState.status = State::initTableDone" << std::endl;
	}
#endif
	releaseShifter(shifter);
}

ParticleModel::Shifter const * ParticleModel::ShiftImageAssistor::acquireShifter(int itrans, int iover_trans, int itrans_over_trans) {
	return acquireShifterW(itrans, iover_trans, itrans_over_trans);
}

ParticleModel::Shifter * ParticleModel::ShiftImageAssistor::acquireShifterW(int itrans, int iover_trans, int itrans_over_trans) {	// makes it easy to find the writers
	auto i = index(itrans, iover_trans, itrans_over_trans);
	auto & desiredState = desiredStates[i];
	
	ParticleModel::Shifter* shifter = shifter = shifterPool->acquire(i, desiredState);
	assert(shifter->currentState.status >= State::initDone);

	return shifter;
}

void ParticleModel::ShiftImageAssistor::releaseShifter(ParticleModel::Shifter const * & shifterConst) {
	auto shifter = const_cast<Shifter*>(shifterConst);
	shifterConst = nullptr;
	releaseShifter(shifter);
}

void ParticleModel::ShiftImageAssistor::releaseShifter(ParticleModel::Shifter * & shifter) {
	shifterPool->release(shifter);
	shifter = nullptr;
}


//
void ParticleModel::initialize(int _ori_size,double _pixel_size,double _particle_diameter,int _width_mask_edge,
                               double _sigma2_fudge,int _random_seed,bool _do_norm_correction,bool _do_zero_mask,
                               bool _do_shifts_onthefly,int _nr_threads,DataStream* _global_data_stream)
{
    ori_size = _ori_size;ori_size2 = _ori_size*ori_size;
    ori_Fsize = (_ori_size/2+1);ori_Fsize2 = (_ori_size/2+1)*_ori_size;
    pixel_size = _pixel_size;
    particle_diameter = _particle_diameter;
    width_mask_edge = _width_mask_edge;
    do_norm_correction = _do_norm_correction;
    do_zero_mask = _do_zero_mask;
    random_seed = _random_seed;
    sigma2_fudge = _sigma2_fudge;
    // Prepare thread data
    nr_threads = _nr_threads;//omp_get_max_threads();
    threadFimages_real	.init(nr_threads, ori_size2);	threadFimages_real	.fill_with_first_touch(0.);
    threadFimages_imag	.init(nr_threads, ori_size2);	threadFimages_imag	.fill_with_first_touch(0.);
    threadImages		.init(nr_threads, ori_size2);	threadImages		.fill_with_first_touch(0.);
    threadFFTtransformer.resize(nr_threads);
    threadCTFer			.resize(nr_threads);
    //
#if defined(FLOAT_PRECISION)
    for(auto &ffter : threadFFTtransformer) ffter = FFTWFTransformer::make(ori_size,ori_size);
#else
    for(auto &ffter : threadFFTtransformer) ffter = FFTWTransformer::make(ori_size,ori_size);
#endif
    for(auto &ctfer : threadCTFer) ctfer = sNew(CTF);
    //
    do_shifts_onthefly = _do_shifts_onthefly;
    global_data_stream = _global_data_stream;
}

//
void ParticleModel::finalize()
{
    threadFimages_real	.fini();
    threadFimages_imag	.fini();
    threadImages		.fini();
    exp_metadata		.fini();
    
    for (auto &ffter : threadFFTtransformer) sDelete(ffter);
    for (auto &ctfer : threadCTFer		   ) sDelete(ctfer);
    threadFFTtransformer.resize(0);
    threadCTFer			.resize(0);
}

// Set up the Images Transformer
void ParticleModel::setup(int _nr_pool,int _current_size,int _coarse_size,int exp_nr_trans/* = 0*/,int exp_nr_over_trans/* = 0*/)
{
    if (exp_nr_trans==0) assert(exp_nr_over_trans==0 && do_shifts_onthefly==false);
    if (exp_nr_trans!=0) assert(exp_nr_over_trans!=0 && do_shifts_onthefly==true);
    exp_nr_images = _nr_pool;
    fine_size = _current_size;fine_size2 = _current_size*_current_size;
    fine_Fsize = (_current_size/2+1);fine_Fsize2 = (_current_size/2+1)*_current_size;
    coarse_size = _coarse_size;coarse_size2 = _coarse_size*_coarse_size;
    coarse_Fsize = (_coarse_size/2+1);coarse_Fsize2 = (_coarse_size/2+1)*_coarse_size;
    
    // Prepare the images data
    // TODO : since each thread will get all images on KNL,fix this first touch
    exp_images	.init(exp_nr_images, ori_size2);	exp_images	.fill_with_first_touch(0.);
    Fctfs		.init(exp_nr_images, fine_Fsize2);	Fctfs		.fill_with_first_touch(0.);
    FctfsOne	.init(exp_nr_images, fine_Fsize2);	FctfsOne	.fill_with_first_touch(1.);
    exp_metadata.init(exp_nr_images);
    
    // Prepare the fourier images data
    Fimages_mask_fine_real	.init(exp_nr_images, fine_Fsize2);	Fimages_mask_fine_real	.fill_with_first_touch(0.);
    Fimages_mask_fine_imag	.init(exp_nr_images, fine_Fsize2);	Fimages_mask_fine_imag	.fill_with_first_touch(0.);
    Fimages_mask_coarse_real.init(exp_nr_images, coarse_Fsize2);Fimages_mask_coarse_real.fill_with_first_touch(0.);
    Fimages_mask_coarse_imag.init(exp_nr_images, coarse_Fsize2);Fimages_mask_coarse_imag.fill_with_first_touch(0.);
    Fimages_nomask_real		.init(exp_nr_images, fine_Fsize2);	Fimages_nomask_real		.fill_with_first_touch(0.);
    Fimages_nomask_imag		.init(exp_nr_images, fine_Fsize2);	Fimages_nomask_imag		.fill_with_first_touch(0.);

	shiftAssistorsSetMaxDims(exp_nr_trans,exp_nr_over_trans);
	shiftAssistorsInitTable(fine_Fsize2);
}

void ParticleModel::destroy()
{
    exp_images				.fini();
    Fctfs					.fini();
    FctfsOne				.fini();
    Fimages_mask_fine_real	.fini();
    Fimages_mask_fine_imag	.fini();
    Fimages_mask_coarse_real.fini();
    Fimages_mask_coarse_imag.fini();
    Fimages_nomask_real		.fini();
    Fimages_nomask_imag		.fini();
    //
    shiftImageAssistor  .setMaxDims(0,0);
#if defined(DOUBLE_TRANSLATION) || defined(TRIPLE_TRANSLATION)
    largeShiftedABTable	.setMaxDims(0,0);
    largeXShiftedABTable.setMaxDims(0,0);
    largeYShiftedABTable.setMaxDims(0,0);
    smallShiftedABTable	.setMaxDims(0,0);
#endif
}

//
void NOINLINE ParticleModel::getLargeShiftedMaskImageOneTile(int iimage,FDOUBLE* Fimgs_shifted_mask_real,FDOUBLE* Fimgs_shifted_mask_imag,
                         					 		int n_start,int n_end,int itrans)
{
    assert(iimage<exp_nr_images);
    // Shift through phase-shifts in the Fourier transform
    // Note that the shift search range is centered around (exp_old_xoff, exp_old_yoff)
	auto shifter = largeShiftedABTable.acquireShifter(itrans,0,itrans);
	shifter->transformOneTile(n_start, n_end,
                                                 Fimages_mask_fine_real[iimage].wptrAll(),
                                                 Fimages_mask_fine_imag[iimage].wptrAll(),
                                                 Fimgs_shifted_mask_real, Fimgs_shifted_mask_imag);
	largeShiftedABTable.releaseShifter(shifter);
}

void NOINLINE ParticleModel::getLargeShiftedNomaskImageOneTile(int iimage,FDOUBLE* Fimgs_shifted_nomask_real,FDOUBLE* Fimgs_shifted_nomask_imag,
                                                      int n_start,int n_end,int itrans)
{
    assert(iimage<exp_nr_images);
    // Shift through phase-shifts in the Fourier transform
    // Note that the shift search range is centered around (exp_old_xoff, exp_old_yoff)
	auto shifter = largeShiftedABTable.acquireShifter(itrans, 0, itrans);
	shifter->transformOneTile(n_start, n_end,
                                                 Fimages_nomask_real[iimage].wptrAll(),
                                                 Fimages_nomask_imag[iimage].wptrAll(),
                                                 Fimgs_shifted_nomask_real, Fimgs_shifted_nomask_imag);
	largeShiftedABTable.releaseShifter(shifter);
}

void  NOINLINE ParticleModel::shiftOneFimageByXYabTable(const FDOUBLE* Fimgs_in_real,const FDOUBLE* Fimgs_in_imag,
                                       			  	 FDOUBLE* Fimgs_out_real,FDOUBLE* Fimgs_out_imag,
                                       			  	 int n_start,int n_end,int itrans,SamplingGrid& samplingGrid)
{
    auto shiftx = samplingGrid.exp_trans_xyshift[itrans].first;
    auto shifty = samplingGrid.exp_trans_xyshift[itrans].second;
    ShiftImageInFourierTransformNew<FDOUBLE, FDOUBLE> transform;
    if (shiftx.status!=0 && shifty.status!=0) {
        int itransx = samplingGrid.exp_positive_shift_index[shiftx.shift];
        int itransy = samplingGrid.exp_positive_shift_index[shifty.shift];
		auto xShifter = largeXShiftedABTable.acquireShifter(itransx, 0, itransx);
		auto yShifter = largeYShiftedABTable.acquireShifter(itransy, 0, itransy);
		transform.transformOneTileByShiftXY(n_start, n_end, Fimgs_in_real, Fimgs_in_imag, Fimgs_out_real, Fimgs_out_imag,
                                            xShifter->aTable_rptr(), 
											xShifter->bTable_rptr(),
                                            xShifter->tableSize(), shiftx.status,
                                            yShifter->aTable_rptr(), 
											yShifter->bTable_rptr(),
                                            yShifter->tableSize(), shifty.status);
		largeYShiftedABTable.releaseShifter(yShifter);
		largeXShiftedABTable.releaseShifter(xShifter);
    }
    else if(shiftx.status!=0 && shifty.status==0){
        int itransx = samplingGrid.exp_positive_shift_index[shiftx.shift];
		auto xShifter = largeXShiftedABTable.acquireShifter(itransx, 0, itransx);
		transform.transformOneTileByShiftX(n_start, n_end, Fimgs_in_real, Fimgs_in_imag, Fimgs_out_real, Fimgs_out_imag,
											xShifter->aTable_rptr(), 
											xShifter->bTable_rptr(),
											xShifter->tableSize(), shiftx.status);
		largeXShiftedABTable.releaseShifter(xShifter);
	}
    else if(shiftx.status==0 && shifty.status!=0){
        int itransy = samplingGrid.exp_positive_shift_index[shifty.shift];
		auto yShifter = largeYShiftedABTable.acquireShifter(itransy, 0, itransy);
		transform.transformOneTileByShiftY(n_start, n_end, Fimgs_in_real, Fimgs_in_imag, Fimgs_out_real, Fimgs_out_imag,
											yShifter->aTable_rptr(),
											yShifter->bTable_rptr(),
											yShifter->tableSize(), shifty.status);
		largeYShiftedABTable.releaseShifter(yShifter);
	}
    else{// just copy
        transform.transformOneTileByCopy(n_start, n_end, Fimgs_in_real, Fimgs_in_imag, Fimgs_out_real, Fimgs_out_imag);
    }
}

void ParticleModel::getLargeShiftedMaskImageDecompOneTile(int iimage,FDOUBLE* Fimgs_shifted_mask_real,FDOUBLE* Fimgs_shifted_mask_imag,
                                                   		  int n_start,int n_end,int itrans,SamplingGrid& samplingGrid)
{
    assert(iimage<exp_nr_images);
    // Shift through phase-shifts in the Fourier transform
    // Note that the shift search range is centered around (exp_old_xoff, exp_old_yoff)
    shiftOneFimageByXYabTable(Fimages_mask_fine_real[iimage].rptrAll(), Fimages_mask_fine_imag[iimage].rptrAll(),
                              Fimgs_shifted_mask_real, Fimgs_shifted_mask_imag, n_start, n_end, itrans, samplingGrid);
}

void ParticleModel::getLargeShiftedNomaskImageDecompOneTile(int iimage,FDOUBLE* Fimgs_shifted_nomask_real,FDOUBLE* Fimgs_shifted_nomask_imag,
                                                            int n_start,int n_end,int itrans,SamplingGrid& samplingGrid)
{
    assert(iimage<exp_nr_images);
    // Shift through phase-shifts in the Fourier transform
    // Note that the shift search range is centered around (exp_old_xoff, exp_old_yoff)
    shiftOneFimageByXYabTable(Fimages_nomask_real[iimage].rptrAll(), Fimages_nomask_imag[iimage].rptrAll(),
                              Fimgs_shifted_nomask_real, Fimgs_shifted_nomask_imag, n_start, n_end, itrans, samplingGrid);
}

void ParticleModel::getShiftedMaskImageOneTileCoarse(int iimage,FDOUBLE* Fimgs_shifted_mask_real,FDOUBLE* Fimgs_shifted_mask_imag,
                                                   	 int n_start,int n_end, int itrans, int iover_trans, int itrans_over_trans)
{
    assert(iimage<exp_nr_images);
    // Shift through phase-shifts in the Fourier transform
    // Note that the shift search range is centered around (exp_old_xoff, exp_old_yoff)
	auto shifter = shiftImageAssistor.acquireShifter(itrans, iover_trans, itrans_over_trans);
	shifter->transformOneTile(	n_start, n_end,
                      	Fimages_mask_coarse_real[iimage].wptr(fine_Fsize),
                      	Fimages_mask_coarse_imag[iimage].wptr(fine_Fsize),
                      	Fimgs_shifted_mask_real, Fimgs_shifted_mask_imag);
	shiftImageAssistor.releaseShifter(shifter);
}

void ParticleModel::getShiftedMaskImageOneTileFine(int iimage,FDOUBLE* Fimgs_shifted_mask_real,FDOUBLE* Fimgs_shifted_mask_imag,
                                                   int n_start,int n_end, int itrans, int iover_trans, int itrans_over_trans)
{
    assert(iimage<exp_nr_images);
    // Shift through phase-shifts in the Fourier transform
    // Note that the shift search range is centered around (exp_old_xoff, exp_old_yoff)
	auto shifter = shiftImageAssistor.acquireShifter(itrans, iover_trans, itrans_over_trans);
    shifter->transformOneTile(	n_start, n_end,
                        Fimages_mask_fine_real[iimage].wptr(fine_Fsize),
                        Fimages_mask_fine_imag[iimage].wptr(fine_Fsize),
                        Fimgs_shifted_mask_real, Fimgs_shifted_mask_imag);
	shiftImageAssistor.releaseShifter(shifter);
}

void ParticleModel::getShiftedNomaskImageOneTile(int iimage,FDOUBLE* Fimgs_shifted_nomask_real,FDOUBLE* Fimgs_shifted_nomask_imag,
                                          		 int n_start,int n_end, int itrans, int iover_trans, int itrans_over_trans)
{
    assert(iimage<exp_nr_images);
    // Shift through phase-shifts in the Fourier transform
    // Note that the shift search range is centered around (exp_old_xoff, exp_old_yoff)
	auto shifter = shiftImageAssistor.acquireShifter(itrans, iover_trans, itrans_over_trans);
	shifter->transformOneTile(n_start, n_end,
                        Fimages_nomask_real[iimage].wptr(fine_Fsize2),
                        Fimages_nomask_imag[iimage].wptr(fine_Fsize2),
                        Fimgs_shifted_nomask_real, Fimgs_shifted_nomask_imag);
	shiftImageAssistor.releaseShifter(shifter);
}

//
void ParticleModel::testDoubleOrTripleTranslation(SamplingGrid& samplingGrid,int exp_current_size)
{
    int exp_nr_trans = samplingGrid.exp_nr_trans();
    int exp_nr_over_trans = samplingGrid.exp_nr_over_trans();
    int exp_current_Fsize2 = exp_current_size*(exp_current_size/2+1);
    FDOUBLE* fout_single_real = (FDOUBLE*)aMalloc(sizeof(FDOUBLE)*fine_Fsize2,64);
    FDOUBLE* fout_single_imag = (FDOUBLE*)aMalloc(sizeof(FDOUBLE)*fine_Fsize2,64);
    FDOUBLE* fout_double_real1 = (FDOUBLE*)aMalloc(sizeof(FDOUBLE)*fine_Fsize2,64);
    FDOUBLE* fout_double_imag1 = (FDOUBLE*)aMalloc(sizeof(FDOUBLE)*fine_Fsize2,64);
    FDOUBLE* fout_double_real2 = (FDOUBLE*)aMalloc(sizeof(FDOUBLE)*fine_Fsize2,64);
    FDOUBLE* fout_double_imag2 = (FDOUBLE*)aMalloc(sizeof(FDOUBLE)*fine_Fsize2,64);
    int diff_count = 0;
    // TODO test shift
    for (int itrans = 0; itrans < exp_nr_trans; itrans++) {
        FDOUBLE shiftx1,shifty1,shiftxOver,shiftyOver;
        samplingGrid.getShiftxy(shiftx1, shifty1, shiftxOver, shiftyOver, itrans, 0);
        auto shiftx2 = samplingGrid.exp_trans_xyshift[itrans].first;
        auto shifty2 = samplingGrid.exp_trans_xyshift[itrans].second;
        assert(shiftx1==shiftx2.status*shiftx2.shift);
        assert(shifty1==shifty2.status*shifty2.shift);
    }
    // TODO
    for (int iimage = 0; iimage < exp_nr_images; iimage++)
    {
        for (int itrans = 0; itrans < exp_nr_trans; itrans++)
        {
//#define TEST_DOUBLE_TRANSLATION
#ifdef TEST_DOUBLE_TRANSLATION
            // double translation 1
            largeShiftedABTable[itrans]
            .transformOneTile(0, exp_current_Fsize2,
                              Fimages_mask_fine_real[iimage].wptrAll(),
                              Fimages_mask_fine_imag[iimage].wptrAll(),
                              fout_double_real1, fout_double_imag1);
#else
            shiftOneFimageByXYabTable(Fimages_mask_fine_real[iimage].wptrAll(), Fimages_mask_fine_imag[iimage].wptrAll(),
                                      fout_double_real1, fout_double_imag1, 0, exp_current_Fsize2, itrans, samplingGrid);
#endif
            //
            for (int iover_trans = 0; iover_trans < exp_nr_over_trans; iover_trans++)
            {
                
                // single translation
				auto shifter = shiftImageAssistor.acquireShifter(itrans, iover_trans, itrans*exp_nr_over_trans+iover_trans);
				shifter->transformOneTile(0, exp_current_Fsize2,
                                   Fimages_mask_fine_real[iimage].wptrAll(),
                                   Fimages_mask_fine_imag[iimage].wptrAll(),
                                   fout_single_real, fout_single_imag);
				shiftImageAssistor.releaseShifter(shifter);
                // double or triple translation 2
				auto smallShifter = smallShiftedABTable.acquireShifter(0, iover_trans, iover_trans);
				smallShifter->transformOneTile(0, exp_current_Fsize2,
                                  fout_double_real1,fout_double_imag1,
                                  fout_double_real2, fout_double_imag2);
				smallShiftedABTable.releaseShifter(smallShifter);

                // compare
                double accuracy = 1e-2;
                for (int i = 0; i < exp_current_Fsize2; i++) {
                    if (fabs(fout_double_real2[i]-fout_single_real[i])/fabs(fout_single_real[i]) > accuracy||
                        fabs(fout_double_imag2[i]-fout_single_imag[i])/fabs(fout_single_imag[i]) > accuracy) {
                        std::cout<<"larger error..."<<itrans<<","<<iover_trans<<","<<fout_double_real2[i]<<","<<fout_single_real[i]<<","
                        		 <<fout_double_imag2[i]<<","<<fout_single_imag[i]<<std::endl;
                        diff_count++;
                    }
                }
            }
        }
    }
    std::cout	<<"####### diff_count percentage = "<<( (double)diff_count /( (double)exp_current_Fsize2*exp_nr_images*exp_nr_trans*exp_nr_over_trans ) )
    			<<" total = "<<(exp_nr_images*exp_nr_trans*exp_nr_over_trans)<<std::endl;
    aFree(fout_single_real);aFree(fout_single_imag);
    aFree(fout_double_real1);aFree(fout_double_imag1);
    aFree(fout_double_real2);aFree(fout_double_imag2);
}

// TODO,use optimized one
void ParticleModel::preShiftedImagesCtfsAndInvSigma2s(Aligned3dArray<FDOUBLE>& exp_Fimgs_shifted_real,
                                                      Aligned3dArray<FDOUBLE>& exp_Fimgs_shifted_imag,
                                                      VectorOfArray2d<FDOUBLE>& exp_local_Minvsigma2s,
                                                      VectorOfArray1d<FDOUBLE>& exp_local_sqrtXi2,
                                                      VectorOfArray2d<FDOUBLE>& exp_local_Fctfs,
                                                      int exp_current_size,bool do_cc,bool do_coarse_search,
                                                      SamplingGrid& samplingGrid,MLModel& mlModel,
                                                      VectorOfInt& Mresol_coarse,VectorOfInt& Mresol_fine)
{
#ifdef DATA_STREAM
    global_data_stream->foutDouble(0, "##################start_precalculateShiftedImagesCtfsAndInvSigma2s#####################", __FILE__, __LINE__);
#endif
    
    const int exp_current_Fsize2 = exp_current_size*(exp_current_size/2+1);
    // wrong because the coarse_size maybe be equal to fine_size
    // bool do_coarse_search = (exp_current_size == coarse_size);
	const int exp_nr_trans       = samplingGrid.exp_nr_trans();
	const int exp_nr_over_trans  = samplingGrid.exp_nr_over_trans();
    // prepare all AB matrix for shift image on-fly
    if (do_shifts_onthefly)
    {
        if (do_coarse_search)
        {
            // This is very fast and not locked		#pragma omp parallel for
            for (int itrans = 0; itrans < exp_nr_trans; itrans++)
            {
                FDOUBLE shiftx,shifty;
                samplingGrid.getShiftxy(shiftx, shifty, itrans, 0);
				shiftImageAssistor.setInitArgs(itrans, 0, itrans*exp_nr_over_trans + 0,
					exp_current_size, shiftx, shifty, ori_size);
            }
        }
        else
        {
			// This is very fast and not locked		#pragma omp parallel for collapse(2)
            for (int itrans = 0; itrans < exp_nr_trans; itrans++){
                for (int iover_trans = 0; iover_trans < exp_nr_over_trans; iover_trans++)
                {
                    FDOUBLE shiftx,shifty;
                    samplingGrid.getShiftxy(shiftx, shifty, itrans, iover_trans);
                    shiftImageAssistor.setInitArgs(itrans, iover_trans, itrans*exp_nr_over_trans + iover_trans,
						exp_current_size, shiftx, shifty, ori_size);
                }
            }
#if defined(DOUBLE_TRANSLATION)
            // double translation
            for (int itrans = 0; itrans < exp_nr_trans; itrans++) {
                FDOUBLE shiftx,shifty,shiftxOver,shiftyOver;
                samplingGrid.getShiftxy(shiftx, shifty, shiftxOver, shiftyOver, itrans, 0);
                largeShiftedABTable.setInitArgs(itrans,0,itrans, exp_current_size, shiftx, shifty, ori_size);
            }
#endif
#if defined(DOUBLE_TRANSLATION) || defined(TRIPLE_TRANSLATION)
            for (int iover_trans = 0; iover_trans < exp_nr_over_trans; iover_trans++) {
                FDOUBLE shiftx,shifty,shiftxOver,shiftyOver;
                samplingGrid.getShiftxy(shiftx, shifty, shiftxOver, shiftyOver, 0, iover_trans);
                smallShiftedABTable.setInitArgs(0,iover_trans, iover_trans, exp_current_size, shiftxOver, shiftyOver, ori_size);
            }
#endif
#if defined(TRIPLE_TRANSLATION)
            // triple translation
            auto & exp_positive_shift_index = samplingGrid.exp_positive_shift_index;
            for (auto & shift_index : exp_positive_shift_index) {
                double shift = shift_index.first;assert(shift>0);
                int index = shift_index.second;
                largeXShiftedABTable.setInitArgs(index,0,index, exp_current_size, shift, 0, ori_size);
                largeYShiftedABTable.setInitArgs(index,0,index, exp_current_size, 0, shift, ori_size);
            }
#endif
        }
    }
    //
#pragma omp parallel for
    for (int iimage = 0; iimage < exp_nr_images; iimage++)
    {
        bool get_shifted_image = false;
#ifdef DATA_STREAM
        global_data_stream->foutInt(iimage, "precalculateShiftedImagesCtfsAndInvSigma2s()_start------", __FILE__, __LINE__);
        get_shifted_image = true;
#endif
        int tid = omp_get_thread_num();
        
        int igroup = exp_metadata[iimage].GROUP_NO-1;

        auto tempFimage_real = threadFimages_real[tid].wptr(ori_size2);
        auto tempFimage_imag = threadFimages_imag[tid].wptr(ori_size2);
        
        if (!do_shifts_onthefly || get_shifted_image || do_cc)
        {
            // Downsize Fimg and Fctf to exp_current_image_size, also initialise Fref and Fimg_shift to the right size
            // In the second pass of the adaptive approach this will have no effect,
            // since then exp_current_image_size will be the same as the size of exp_Fctfs
            windowFourierTransform(Fimages_mask_fine_real[iimage].wptr(fine_Fsize2),
                                   Fimages_mask_fine_imag[iimage].wptr(fine_Fsize2), fine_size,
                                   tempFimage_real, tempFimage_imag, exp_current_size);
        }
        
#ifdef DATA_STREAM
        global_data_stream->foutDouble(tempFimage_real, exp_current_Fsize2, "precalculateShiftedImagesCtfsAndInvSigma2s()_Fimg_real", __FILE__, __LINE__);
        global_data_stream->foutDouble(tempFimage_imag, exp_current_Fsize2, "precalculateShiftedImagesCtfsAndInvSigma2s()_Fimg_imag", __FILE__, __LINE__);
#endif
        
        if (do_cc)
        {
            double sumxi2 = 0.;
            for (int n = 0; n < exp_current_Fsize2; n++) {
                sumxi2 += tempFimage_real[n]*tempFimage_real[n] \
                        + tempFimage_imag[n]*tempFimage_imag[n];
            }
            // Normalised cross-correlation coefficient: divide by power of reference (power of image is a constant)
            exp_local_sqrtXi2.wptrAll()[iimage] = sqrt(sumxi2);
        }
        
        // always need to do even coarse_size=fine_size to avoid transformFourierArrayToShellIncreasedCoarse .. Fine bug....
        windowTransform(Fctfs[iimage].wptr(fine_Fsize2), fine_size, exp_local_Fctfs[iimage].wptr(exp_current_Fsize2), exp_current_size);
        
#ifdef DATA_STREAM
        global_data_stream->foutDouble(exp_local_Fctfs[iimage].wptr(exp_current_Fsize2), exp_current_Fsize2, "precalculateShiftedImagesCtfsAndInvSigma2s()_exp_local_Fctfs", __FILE__, __LINE__);
        global_data_stream->check();global_data_stream->flush();
#endif
      
        // only need if not do-shift-on-fly or do-data-stream comparison
        if (!do_shifts_onthefly || get_shifted_image)
        {
            for (int itrans = 0; itrans < exp_nr_trans; itrans++)
            {
                for (int iover_trans = 0; iover_trans < exp_nr_over_trans; iover_trans++)
                {
                    // Shift through phase-shifts in the Fourier transform
                    // Note that the shift search range is centered around (exp_old_xoff, exp_old_yoff)
                    FDOUBLE shiftx,shifty;
                    samplingGrid.getShiftxy(shiftx, shifty, itrans, iover_trans);
                    
                    int itrans_over_trans = itrans*exp_nr_over_trans+iover_trans;
                    auto exp_Fimgs_shifted_real_aux = exp_Fimgs_shifted_real.wptr(iimage, itrans_over_trans, exp_current_Fsize2);
                    auto exp_Fimgs_shifted_imag_aux = exp_Fimgs_shifted_imag.wptr(iimage, itrans_over_trans, exp_current_Fsize2);
                    
                    shiftImageInFourierTransform(tempFimage_real, tempFimage_imag,
                                                 exp_Fimgs_shifted_real_aux,exp_Fimgs_shifted_imag_aux,
                                                 exp_current_size,shiftx,shifty,ori_size);
                }
            }
        }

        auto Minvsigma2 = exp_local_Minvsigma2s[iimage].wptr(exp_current_Fsize2);
        auto myMresol = do_coarse_search ? Mresol_coarse.rptrAll() : Mresol_fine.rptrAll();
        auto sigma2_noise_igroup = mlModel.sigma2_noise[igroup].wptr(ori_Fsize);
        // With group_id and relevant size of Fimg, calculate inverse of sigma^2 for relevant parts of Mresol
        for (int n = 0; n < exp_current_Fsize2; n++)
        {
            int ires = myMresol[n];
            // Exclude origin (ires==0) from the Probability-calculation
            // This way we are invariant to additive factors
            if (ires > 0){
                Minvsigma2[n] = CHECK_NOT_NAN_OR_INF(1. / (sigma2_fudge * sigma2_noise_igroup[ires]));
            }
            else
                Minvsigma2[n] = 0;
        }
#ifdef DATA_STREAM
        global_data_stream->foutDouble(Minvsigma2, exp_current_Fsize2, "precalculateShiftedImagesCtfsAndInvSigma2s()_exp_local_Minvsigma2s", __FILE__, __LINE__);
        global_data_stream->foutDouble(exp_Fimgs_shifted_real.wptr(iimage, 0, exp_current_Fsize2), exp_current_Fsize2, "precalculateShiftedImagesCtfsAndInvSigma2s()_exp_Fimgs_shifted_real1", __FILE__, __LINE__);
        global_data_stream->foutDouble(exp_Fimgs_shifted_imag.wptr(iimage, 0, exp_current_Fsize2), exp_current_Fsize2, "precalculateShiftedImagesCtfsAndInvSigma2s()_exp_Fimgs_shifted_imagn", __FILE__, __LINE__);
        global_data_stream->foutDouble(exp_Fimgs_shifted_real.wptr(iimage, exp_nr_trans*exp_nr_over_trans-1, exp_current_Fsize2), exp_current_Fsize2, "precalculateShiftedImagesCtfsAndInvSigma2s()_exp_Fimgs_shifted_realN", __FILE__, __LINE__);
        global_data_stream->foutDouble(exp_Fimgs_shifted_imag.wptr(iimage, exp_nr_trans*exp_nr_over_trans-1, exp_current_Fsize2), exp_current_Fsize2, "precalculateShiftedImagesCtfsAndInvSigma2s()_exp_Fimgs_shifted_imagN", __FILE__, __LINE__);
        global_data_stream->check();global_data_stream->flush();
#endif
    }
}


void ParticleModel::pregetShiftedImagesNomask(Aligned3dArray<FDOUBLE>& exp_nomask_Fimgs_shifted_real,
                                              Aligned3dArray<FDOUBLE>& exp_nomask_Fimgs_shifted_imag,
                                              int exp_current_size,SamplingGrid& samplingGrid)
{
    assert(exp_current_size==fine_size);
    int exp_current_Fsize2 = exp_current_size*(exp_current_size/2+1);
    int exp_nr_trans = samplingGrid.exp_nr_trans();
    int exp_nr_over_trans = samplingGrid.exp_nr_over_trans();
	
    // NOTE : in map3d_optimizer_new do-shift-on-fly,sot this is not called
	CHECKER_PARALLEL_COUNTED_FOR3(
		int, iimage, 0, exp_nr_images,
		int, itrans, 0, exp_nr_trans,
		int, iover_trans, 0, exp_nr_over_trans)
// #pragma omp parallel for collapse(3)
//     for (int iimage = 0; iimage < exp_nr_images; iimage++)
//     {
//         for (int itrans = 0; itrans < exp_nr_trans; itrans++)
//         {
//             for (int iover_trans = 0; iover_trans < exp_nr_over_trans; iover_trans++)
            {
                int itrans_over_trans = itrans*exp_nr_over_trans + iover_trans;
                
                int tid = omp_get_thread_num();
                
                // auto tempFimage_real = threadFimages_real[tid];
                // auto tempFimage_imag = threadFimages_imag[tid];
                
                // exp_current_size is same as current_size
                // windowFourierTransform(Fimages_nomask_real[iimage], Fimages_nomask_imag[iimage], fine_size,
                //                        tempFimage_real, tempFimage_imag, exp_current_size);
                
                // Shift through phase-shifts in the Fourier transform
                // Note that the shift search range is centered around (exp_old_xoff, exp_old_yoff)
                FDOUBLE shiftx,shifty;
                samplingGrid.getShiftxy(shiftx, shifty, itrans, iover_trans);
                
                // NOTE : using exp_Fimgs_shifted as exp_Fimg_shifted_nomask
                auto exp_Fimgs_shifted_nomask_real_aux = exp_nomask_Fimgs_shifted_real.wptr(iimage, itrans_over_trans, exp_current_Fsize2);
                auto exp_Fimgs_shifted_nomask_imag_aux = exp_nomask_Fimgs_shifted_imag.wptr(iimage, itrans_over_trans, exp_current_Fsize2);
                
                //shiftImageInFourierTransform(Fimg_nomask_aux, exp_Fimgs_shifted_nomask_aux,exp_current_size,shiftx,shifty,ori_size);
                shiftImageInFourierTransform(Fimages_nomask_real[iimage].wptr(fine_Fsize2),
                                             Fimages_nomask_imag[iimage].wptr(fine_Fsize2),
                                             exp_Fimgs_shifted_nomask_real_aux,exp_Fimgs_shifted_nomask_imag_aux,
                                             exp_current_size,shiftx,shifty,ori_size);
//            }
//        }
    } CHECKER_PARALLEL_FOR_END
}
    
void ParticleModel::unitTest()
{
    // Vector2D test;
}
