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

#include "l2cachemodel.h"
// L2 cache modelling is used to try to understand our L2 miss rates in production runs without distort performance
//
namespace L2CacheModel {
    
    class SeqAcc {
    public:
        SeqAcc* _next;
#define SEP
#define ELT(T,N)	T N;
        L2_ACCESS_ELTS
#undef ELT
#undef SEP
    };
    
    struct PerThread {
        size_t  _active;
        SeqAcc* _oldest;
        SeqAcc* _newest;
        PerThread() : _active(0), _oldest(nullptr), _newest(nullptr) {}
    };
    
    std::vector<PerThread> perThreads(omp_get_max_threads());
    
    class CHiP {
    public:
        static const size_t cacheLineSizeLog2	=  6;	// 64 Byte cache lines
        static const size_t cacheLineSize		= (size_t(1) << cacheLineSizeLog2);
        static const size_t cacheLineOffsetMask	= cacheLineSize - 1;
        static const size_t maxCacheSizeLog2	= 		// KNL 1MB max interesting cache size,			20
        // but shared by two cores						19
        17;	// but BDW and SKL L2 caches are only 256KB		17
        static const size_t usableCacheSize     = (size_t(1) << (maxCacheSizeLog2-1));
        // Only try to half-fill the available cache
        static const size_t depthLog2			=  4;	// allows for an efficient implementation
        static const size_t widthLog2			= maxCacheSizeLog2 - cacheLineSizeLog2 - depthLog2;
        static const size_t width				= 1 << widthLog2;
        static const size_t depth				= 1 << depthLog2;
        
        __declspec(align(64))
        size_t cols    [width*depth];
        size_t rowCount[depth];
        size_t misses;
        CHiP() {
            for (size_t i = 0; i < width*depth; i++) cols    [i] = 0;
            for (size_t r = 0; r <       depth; r++) rowCount[r] = 0;
            misses = 0;
        }
        size_t total() {
            size_t total(misses);
            for (size_t r = 0; r < depth; r++) {
                total += rowCount[r];
            }
            return total;
        }
        float hitRate() {
            size_t sizeSoFar(0);
            size_t total(0);
            size_t r = 0;
            for (; r < depth; r++) {
                total += rowCount[r];
                sizeSoFar += width*cacheLineSize;
                if (sizeSoFar > usableCacheSize) break;
            }
            auto hits = total;
            for (; r < depth; r++) {
                total += rowCount[r];
            }
            total += misses;
            return hits / float(std::max(size_t(1),total));
        }
        float missPct() {
            return 0.1f * int(1000.0f*hitRate());
        }
        bool goodHitRate() {
            return missPct() < 0.5;
        }
        void showSummary(std::ostream & os, Life* life) {
            const size_t total = this->total();
            auto missPct = this->missPct();
            os << "L2 misses/total " << misses << "/" << total
            << "(" << ((total == 0) ? 0 : missPct) << "%)";
            if (life) {
                os << " during " << life->path();
            }
            os << " tid:" << omp_get_thread_num() << std::endl;
        }
        void showDetails(std::ostream & os, Life* life) {
            const size_t total = this->total();
            if (total == 0) return;
            size_t size(0), sum(0);
            for (size_t r = 0; r < depth; r++) {
                size += width*cacheLineSize; sum += rowCount[r];
                auto pct = 0.1f * int(1000.0 * sum / total);
                os << "    " << std::setw(6) << size << " " << pct << std::endl;
                if (pct >= 99.9) break;
            }
        }
        void access(size_t addr, size_t count, size_t amount, size_t stride) {
            if (amount == stride) {
                while (count > 0 && !(count & 1) && (amount < cacheLineSize)) {
                    amount *= 2;
                    stride *= 2;
                    count  /= 2;
                }
            }
            while (count--) {
                size_t remains = amount;
                while (remains) {
                    size_t o     = addr & cacheLineOffsetMask;
                    size_t eaten = std::min(amount, cacheLineSize - o);
                    access(addr >> cacheLineSizeLog2);
                    addr    += eaten;
                    remains -= eaten;
                }
                addr += stride - amount;
            }
        }
        void access(size_t line) {
            size_t c = line & (width-1);
            size_t* col = &cols[c * depth];
            size_t insert = line;
            for (size_t r = 0; r < depth; r++) {
                size_t prev = col[r];
                col[r] = insert;
                if (prev == line) {
                    rowCount[r]++;
                    return;
                }
                insert = prev;
            }
            misses++;
        }
    };
    
#ifdef L2_CACHE_MODELING
    void seqAccWkr(
#define SEP			,
#define ELT(T,N)	T N
                   L2_ACCESS_ELTS
#undef ELT
#undef SEP
    ) {
        auto & pt = perThreads[omp_get_thread_num()];
        auto &_active = pt._active;
        auto &_oldest = pt._oldest;
        auto &_newest = pt._newest;
        
        if (!_active) return;
        auto p = sNew(SeqAcc);
        p->_next = nullptr;
#define SEP
#define ELT(T,N)	p->N = N;
        L2_ACCESS_ELTS
#undef ELT
#undef SEP
        if (_newest) _newest->_next = p; else _oldest = p;
        _newest = p;
        
        // BEVIN
        // std::cerr << "seqAccWkr list is now " << std::endl;
        // for (auto p = _oldest; !!p; p = p->_next) {
        // 	std::cerr << "   " << p->nameBase << "[" << p->nameIndex << "]" << std::endl;
        // }
    }
#endif
    
    class Interval::Impl {
    public:
        size_t  _depth;
        SeqAcc* _beginNewest;
        Impl() : _depth(0), _beginNewest(nullptr) {}
    };
    
    void Interval::begin() {
#ifdef L2_CACHE_MODELING
        if (!!_impl) {
            _impl->_depth++;
            return;
        }
        auto & pt = perThreads[omp_get_thread_num()];
        auto &_active = pt._active;
        auto &_oldest = pt._oldest;
        auto &_newest = pt._newest;
        _impl = sNew(Impl);
        _impl->_beginNewest = _newest;
        _active++;
#endif
    }
    
    void Interval::endWkr() {
        if (--_impl->_depth) return;
        sDelete(_impl);
        auto & pt = perThreads[omp_get_thread_num()];
        auto &_active = pt._active;
        auto &_oldest = pt._oldest;
        auto &_newest = pt._newest;
        --_active;
        if (_active > 0) return;
        while (_oldest) {
            auto next = _oldest->_next;
            sDelete(_oldest);
            _oldest = next;
        }
        _newest = nullptr;
    }
    
    size_t Interval::accessesRecordedCount() {
#ifndef L2_CACHE_MODELING
        return 1;
#else
        if (!_impl) return 0;
        auto & pt = perThreads[omp_get_thread_num()];
        auto &_active = pt._active;
        auto &_oldest = pt._oldest;
        auto &_newest = pt._newest;
        size_t sum(0);
        for (auto p = !!_impl->_beginNewest ? _impl->_beginNewest->_next : _oldest; p; p = p->_next) {
            sum++;
        }
        return sum;
#endif
        
    }
    
    size_t Interval::cacheProbeCount() {
#ifndef L2_CACHE_MODELING
        return 1;
#else
        if (!_impl) return 0;
        auto & pt = perThreads[omp_get_thread_num()];
        auto &_active = pt._active;
        auto &_oldest = pt._oldest;
        auto &_newest = pt._newest;
        size_t sum(0);
        for (auto p = !!_impl->_beginNewest ? _impl->_beginNewest->_next : _oldest; p; p = p->_next) {
            sum += (p->count*p->amount);
        }
        return sum / CHiP::cacheLineSize;		// not right but good enough
#endif
    }
    
    bool Interval::showL2CacheIntervalUnlocked(Life* life) {
        showL2CacheInterval(std::cerr, life);
        return showL2CacheInterval(std::cout, life);
    }
    
    float Interval::hitRate() {
#ifndef L2_CACHE_MODELING
        return 1.0;
#else
        if (!_impl) return 1.0;
        auto & pt = perThreads[omp_get_thread_num()];
        auto &_active = pt._active;
        auto &_oldest = pt._oldest;
        auto &_newest = pt._newest;
        CHiP chip;
        for (auto p = !!_impl->_beginNewest ? _impl->_beginNewest->_next : _oldest; p; p = p->_next) {
            chip.access(size_t(p->start), p->count, p->amount, p->stride);
        }
        return chip.hitRate();
#endif
    }
    
    bool Interval::showL2CacheInterval(std::ostream & os, Life* life) {
#ifndef L2_CACHE_MODELING
        return true;
#else
        if (!_impl) return true;
        
        auto & pt = perThreads[omp_get_thread_num()];
        auto &_active = pt._active;
        auto &_oldest = pt._oldest;
        auto &_newest = pt._newest;
        CHiP chip;
        for (auto p = !!_impl->_beginNewest ? _impl->_beginNewest->_next : _oldest; p; p = p->_next) {
            chip.access(size_t(p->start), p->count, p->amount, p->stride);
        }
        
        if (chip.goodHitRate()) return true;
        
        chip.showSummary(os, life);
        
        // eliminate duplicates
        //
        struct PerSite {
            float  missPct;
            size_t duplicates;
            PerSite() : missPct(0.0), duplicates(0) {}
        };
        
        auto const path = life->path();
        auto const missPct = chip.missPct();
        
        static std::map<const std::string, PerSite> perSitePerFile[2];
        auto& perSite = perSitePerFile[(&os == &std::cout) ? 0 : 1];
        
        auto p = perSite.find(path);
        if (p == perSite.end()) {
            perSite.insert({path,PerSite()}); p = perSite.find(path);
        }
        auto & prev = p->second;
        
        if (prev.missPct >= missPct) {
            prev.duplicates++;
            return false;
        }
        if (prev.duplicates > 1) os << "...prev detailed cache hit behavior bettered or duplicated " << prev.duplicates - 1 << " times" << std::endl;
        prev.missPct = missPct;
        prev.duplicates = 1;
        
        // not a duplicate
        //
        chip.showDetails(os, life);	// doesn't show anything if good hit rate
        
        if (true) {
            os << "Traffic causing this" << std::endl;
            struct Data { SeqAcc* seqAcc; size_t count; Data(SeqAcc* seqAcc, size_t count) : seqAcc(seqAcc), count(count) {} };
            std::map<void*, Data> variables;
            for (auto p = !!_impl->_beginNewest ? _impl->_beginNewest->_next : _oldest; p; p = p->_next) {
                auto i = variables.find(p->start);
                if (i != variables.end()) i->second.count++; else variables.insert(std::make_pair(p->start, Data(p,1)));
            }
            for (auto & i : variables) {
                auto seqAcc = i.second.seqAcc;
                os << "    " << i.first << " " << seqAcc->nameBase << " " << i.second.count 
                << " uses accessing " << seqAcc->count*seqAcc->amount << " bytes" 
                << std::endl;
            }
        }
        return false;
#endif
    }
}
