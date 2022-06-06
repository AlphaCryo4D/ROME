/***************************************************************************
*
* Authors: "Yongbei(Glow) Ma,Jiayi (Timmy) Wu, Youdong (Jack) Mao"
* Dana-Farber Cancer Institute, Harvard Medical School and Peking University
* "Bevin Brett"
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

#include "./resmap_exp_mweight.h"
	// must follow the system includes

#ifndef DONT_INCLUDE_SAMPLING
#include "./resmap_sampling.h"
#endif

#ifdef Mweight_DEBUGGING

DebuggingMweights::DebuggingMweights() : _state("ctr"), _count(1) {
	sgnt_iimage = -1;
	sgnt_iclass = -1;
	sgnt_idir = -1;
	sgnt_ipsi = -1;
	sgnt_iover_rot = -1;
	sgnt_itrans = -1;
	sgnt_iover_trans = -1;
	sgnt_index = -1;
}

bool DebuggingMweights::setState(const char* to, bool interesting) {
	if (to == _state && !interesting) { 
		_count++; 
		if (_count == 10) std::cerr << "..." << std::endl;
		return _count < 10; 
	}
	std::cerr << "DebuggingMweights::setState ";
	if (to != _state) std::cerr << "changing from " << _state << " to ";
	std::cerr << to;
	if (_count > 1) std::cerr << " after " << _count;
	std::cerr << std::endl;
	_state = to;
	_count = 1;
	return true;
}

std::ostream& showDimensions(int iimage, int iclass, int idir, int ipsi) {
	std::cerr 
		<<  "image:" << iimage 
		<< " class:" << iclass
		<< " dir:" << idir
		<< " psi:" << ipsi
		<< " ";
	return std::cerr;
}

std::ostream& DebuggingMweights::showDimensions(MweightsForSomeSpinsAndSlides* in) {
	return ::showDimensions(in->_iimage, in->_iclass, in->_idir, in->_ipsi);
}

void DebuggingMweights::clear(int iimage) {
	if (setState("clear")) std::cerr << "image:" << iimage << " cleared" << std::endl;
}

void DebuggingMweights::insert(MweightsForSomeSpinsAndSlides* in, size_t index, PackedOverRotTransOverTrans key, double to) {
	if (setState("insert")) showDimensions(in) << "index:" << index << " set to " << to << std::endl;
}

void DebuggingMweights::setValue(MweightsForSomeSpinsAndSlides* in, size_t index, double to) {
	if (setState("setValue")) showDimensions(in) << "index:" << index << " set to " << to << "  count:" << _count << std::endl;
}

void DebuggingMweights::setSignificant(
	int iimage, int iclass, int idir, int ipsi, int iover_rot, int itrans, int iover_trans, int index, bool value) {
	if (!setState("setSignificant")) return;
	::showDimensions(iimage, iclass, idir, ipsi) 
		<< " iover_rot:" << iover_rot
		<< " itrans:" << itrans
		<< " iover_trans:" << iover_trans
		<< " index:" << index
		<< " wrote:" << value
		<< std::endl;
	sgnt_iimage = iimage;
	sgnt_iclass = iclass;
	sgnt_idir = idir;
	sgnt_ipsi = ipsi;
	sgnt_iover_rot = iover_rot;
	sgnt_itrans = itrans;
	sgnt_iover_trans = iover_trans;
	if (value) sgnt_index = index;
}

void DebuggingMweights::getSignificant(
	int iimage, int iclass, int idir, int ipsi, int iover_rot, int itrans, int iover_trans, int index, bool value) {

	bool interesting =
		true
#define P(N) ((N == -1) || (N == sgnt_##N))
		&& P(iimage)
		&& P(iclass)
		&& P(idir)
		&& P(ipsi)
		&& P(iover_rot)
		&& P(itrans)
		&& P(iover_trans);
	interesting |=
		P(index);
#undef P
	if (!setState("getSignificant", interesting)) return;
	::showDimensions(iimage, iclass, idir, ipsi)
		<< " iover_rot:" << iover_rot
		<< " itrans:" << itrans
		<< " iover_trans:" << iover_trans
		<< " index:" << index;
	if (index == sgnt_index)
		std::cerr << " ******** read the significant index";
	std::cerr << " read:" << value << std::endl;
}


DebuggingMweights debuggingMweights;

#endif

std::string to_string(MweightsUpperIndex const & mui) {
	std::string s;
#ifdef NDEBUG
	s = "MweightsUpperIndex_unknown";
#else
	s = "{";
	s+= "iimage:"     + std::to_string((long long)mui.iimage())
		+ ", iclass:" + std::to_string((long long)mui.iclass())
		+ ", idir:"	  + std::to_string((long long)mui.idir  ())
		+ ", ipsi:"	  + std::to_string((long long)mui.ipsi  ())
		+ "}";
#endif
	return s;
}


MweightsForSomeSpinsAndSlides::MweightsForSomeSpinsAndSlides()
  : _capacity(0),
#ifdef  MweightsForSomeSpinsAndSlides_USE_VECTOR
	_keys(0),
	_data(0),
#else
	_keys(nullptr),
	_data(nullptr),
#endif
#ifdef MweightsForSomeSpinsAndSlides_INITED_CHECKS
	_initedSize(0),
#ifdef  MweightsForSomeSpinsAndSlides_USE_VECTOR
	_inited(0),
#else
	_inited(nullptr),
#endif
#endif
	_size    (0)
{
}

void MweightsForSomeSpinsAndSlides::setParent(Exp_Mweight_new* parent, MweightsUpperIndex const & mui) {
	mui_connect(parent, mui);
	auto capacity = parent->nr_over_rot * parent->nr_trans * parent->nr_over_trans;
	if (_capacity < capacity) {
		_capacity = capacity;
#ifdef  MweightsForSomeSpinsAndSlides_USE_VECTOR
		_keys.resize(_capacity);
		_data.resize(_capacity);
#else
		vDelete(_keys); _keys = vNew(PackedOverRotTransOverTrans, _capacity);
		vDelete(_data); _data = vNew(double, _capacity);
#endif
#ifdef MweightsForSomeSpinsAndSlides_INITED_CHECKS
		_initedSize = (_capacity + 8 * sizeof(size_t) - 1) / (8 * sizeof(size_t));
#ifdef  MweightsForSomeSpinsAndSlides_USE_VECTOR
		_inited(_initedSize),
#else
		vDelete(_inited); _inited = vNew(size_t,_initedSize);
#endif
#endif
	}
	_size = 0;
	uninit();
	if (interesting()) std::cerr << "uid:" << uid() << " having its parent set, _capacity:" << _capacity << std::endl;
}

MweightsForSomeSpinsAndSlides::~MweightsForSomeSpinsAndSlides() {
#ifdef  MweightsForSomeSpinsAndSlides_USE_VECTOR
	_data.resize(0);
	_keys.resize(0);
#ifdef MweightsForSomeSpinsAndSlides_INITED_CHECKS
	_inited.resize(0);
#endif
#else
	vFreeConst(_data); 
	vFreeConst(_keys); 
#ifdef MweightsForSomeSpinsAndSlides_INITED_CHECKS
	delete[] _inited;
#endif
#endif
}

size_t MweightsForSomeSpinsAndSlides::bytesUsed() const {
	size_t sum = sizeof(*this) + _capacity*(sizeof(PackedOverRotTransOverTrans) + sizeof(double));
	return sum;
}

void MweightsForSomeSpinsAndSlides::uninit() {
#ifdef MweightsForSomeSpinsAndSlides_INITED_CHECKS
	for (size_t i = 0; i <_initedSize; i++) _inited[i] = 0;
#endif
}

void MweightsForSomeSpinsAndSlides::clear() {
	if (interesting()) {
		std::cerr << "MweightsForSomeSpinsAndSlides::clear() interesting(). " << std::endl;
	}
	uninit();
	_size = 0;
}

void MweightsForSomeSpinsAndSlides::zero() {
	if (interesting()) {
		std::cerr << "MweightsForSomeSpinsAndSlides::zero() interesting(). " << std::endl;
	}
	for (size_t i = 0; i < _size; i++) _data[i] = 0.0;
#ifdef MweightsForSomeSpinsAndSlides_INITED_CHECKS
	for (size_t i = 0; i <_initedSize; i++) _inited[i] = -1;
#endif
}

void MweightsForSomeSpinsAndSlides::zeroSomeMakeSureRestUnused() {
	uninit();
}

void MweightsForSomeSpinsAndSlides::insert(
	PackedOverRotTransOverTrans key,
	double to) {
	if (interesting()) {
		int iover_rot, itrans, iover_trans;
		key.decode(iover_rot, itrans, iover_trans);
		std::cerr << "insertion " << _size+1 << " into interesting() uid:" << uid()
			<< " mui:"			<< to_string(mui())
			<< " iover_rot:"	<< iover_rot
			<< " itrans:"		<< itrans
			<< " iover_trans:"	<< iover_trans
			<< std::endl;
	}
	if (_size >= _capacity) {
		std::cerr << "MweightsForSomeSpinsAndSlides::insert no room in "
			<< " mui:"		 << to_string(mui())
			<< " uid:"       << uid()
			<< " _size:"     << _size
			<< " _capacity:" << _capacity 
			<< std::endl;
		assert(false);
	}
	assert(key.inited());
	_keys[_size] = key;
	_data[_size] = to;
	inited(_size);
#ifdef Mweight_DEBUGGING
	debuggingMweights.insert(this, _size, key, to);
#endif
	_size++;
}


namespace Exp_MWeight_Namespace {

    // thread data Heap
	template <class _Node>
	class NodePoolTemplate : public CensusDetailsProvider {
		typedef _Node Node;
		__declspec(align(64)) // so each gets its own cache lines on a core to avoid false sharing
		Node*		_head;
		const char* _name;
		int			_poolIndex;
		size_t		_inserts, _removes;
		size_t		_poolSize;
	public:
		NodePoolTemplate() : _name("anon"), _poolIndex(-1), _head(nullptr), _inserts(0), _removes(0), _poolSize(0) { }
		static Node* remove()           { 
			auto& p = pool(); auto & _head = p._head; 
			Node* n = nullptr; 
			if (_head) { n = _head; _head = n->_next; p._removes++; p._poolSize--; }
			assert(!n || n->mui_iimage() == -1);
			if (n && n->interesting()) std::cerr << "_uid:" << n->uid() << " removed from pool" << std::endl;
			return n;
		}
		static void  insert(Node*& n)	{ 
			assert(n->mui_iimage() == -1);
			assert(n->childCount() == 0);
			if (n && n->interesting()) std::cerr << "_uid:" << n->uid() << " inserted into pool" << std::endl;
			auto& p = pool(); auto & _head = p._head; 
			n->_next = _head; _head = n; n = nullptr; 
			p._inserts++; p._poolSize++;
		}
		virtual void showPopulation(std::ostream & os) {
			
			return;	//	Only needed during debugging and tuning

			bool canShow = (omp_get_num_threads() == 1);
			static bool prevCanShow = true;
			if (!canShow) {
				if (prevCanShow) os << "Can't show exp_mweight pool " << _name << "[" << _poolIndex << "] because threads running" << std::endl;
			} else {
				size_t count = 0;
				for (auto n = _head; n; n = n->_next) count++;
				if (count + _inserts + _removes > 0) {
					os << "Exp_Mweight pool:" << _name << "[" << _poolIndex << "]" << ", Population:" << count << ", inserts:" << _inserts << ", removes:" << _removes << std::endl;
				}
			}
			prevCanShow = canShow;
		}
		static void setNameEtc(const char* name) {
			for (size_t poolIndex = 0; poolIndex < maxPools; poolIndex++) pool(poolIndex).setNameEtc(name,poolIndex);
		}
		void setNameEtc(const char* name, int poolIndex) {
			_name		= name;
			_poolIndex  = poolIndex;
		}
        size_t inline getPoolSize() {
            return _poolSize;
        }
		const char* getNameEtc() {
			return _name;
		}
        static size_t show(bool debug_flag) {
            size_t total_count = 0;
            if(debug_flag) std::cout<< pool(0).getNameEtc() << ",Each thread's population  : ";
            for (size_t poolIndex = 0; poolIndex < omp_get_max_threads(); poolIndex++){
                auto count = pool(poolIndex).getPoolSize();
                total_count += count;
                if(debug_flag) std::cout<<" "<<count<<" ";
            }
            if(debug_flag) std::cout<<" Total : "<<total_count<<std::endl;
            return total_count;
        }
		static bool hasNode() {
			return pool().getPoolSize() > 0;
		}
		// clear NodeForMweightsForSomeSpinsAndSlides pool if possible
		//void clear(double percentage) {
		//	assert(_name == "NodeForMweightsForSomeSpinsAndSlides");
		//	if (!_head) return;
		//	size_t size = get_count()*percentage;
		//	size_t i = 0;
		//	auto& pre = _head;
		//	//
		//	while(pre->_next) {
		//		auto next = pre->_next;
		//		sDelete(pre); pre = nullptr;
		//		_head = next; pre = _head;
		//		if (++i == size) break;
		//	}
		//}
  //      static void clearPool(double percentage) {
		//	#pragma omp parallel for
		//	for (size_t poolIndex = 0; poolIndex < omp_get_max_threads(); poolIndex++) {
		//		pool(poolIndex).clear(percentage);
		//	}
		//}
	private:
		static const size_t maxPools = 256;
		static NodePoolTemplate& pool() {
			static int _nextPoolIndex = 0;
			auto poolIndex = 
				omp_get_num_threads() > 1 
				? omp_get_thread_num()												// so don't need locks
				: ( _nextPoolIndex = (_nextPoolIndex+1) % omp_get_max_threads() );	// so spread evenly across the pools
			return pool(poolIndex);
		}
		static NodePoolTemplate& pool(int poolIndex) {
			if (poolIndex >= maxPools) {
				std::cerr << __FILE__ << ":" << __LINE__ << " too many threads, poolIndex:" << poolIndex << ", maxPools:" << maxPools << std::endl;
				EXIT_ABNORMALLY;
			}
			static NodePoolTemplate p[maxPools];
			return p[poolIndex];
		}
	};

    // Fix a bug : 0x0000000003c0d28c in _INTERNAL_22_______src_kmp_lock_cpp_a132da10::__kmp_unset_indirect_lock(unsigned int*, int) ()
    // centos 7.3
    static Lock							_exp_lock;
    
	template <class _NextLevel,int _sizes_size>
	class NodeTemplate : public ForMweightsUpperIndex {
	public:
        // initialized by template parameter to avoid Intel complier crash on Mac
        // static const size_t	   _sizes_size = _NextLevel::_sizes_size + 1;
	private:
		typedef _NextLevel NextLevel;
		typedef NextLevel* NextLevelPtr;
		size_t							_sizes[_sizes_size];
#ifdef  MweightsForSomeSpinsAndSlides_USE_VECTOR
		std::vector<NextLevelPtr>		_data;
		bool			 volatile		_dataSized;
#else
		NextLevelPtr volatile* volatile	_data;
#endif
        // use global  Lock instead
//		Lock							_lock;

	public:
		// NOTE : 
		// Bevin : should not be used, should be returned to the pool
		// YongBei : need to delete if it run out os memory
		~NodeTemplate() { /*assert(!"should not be used, should be returned to the pool");*/ }

		NodeTemplate*	_next;	// used for the pool

		static NodeTemplate* make(Exp_Mweight_new* parent, MweightsUpperIndex const & mui, size_t sizes_size, size_t* sizes) {
			auto n = NodePoolTemplate<NodeTemplate>::remove();
			if (!n) n = sNew(NodeTemplate);
			n->setParent(parent, mui, sizes_size, sizes);
			assert(n->mui().matches(mui));
			return n;
		}
		NodeTemplate() 
			: 
#ifdef  MweightsForSomeSpinsAndSlides_USE_VECTOR
			  _data(0), _dataSized(false)
#else
			  _data(nullptr)
#endif
		{
			for (size_t i = 0; i < _sizes_size; i++) _sizes[i] = 0;
		}
		void setParent(Exp_Mweight_new* parent, MweightsUpperIndex const & mui, size_t sizes_size, size_t* sizes) {
			mui_connect(parent, mui, true);
			assert(_sizes_size == sizes_size);
#ifdef  MweightsForSomeSpinsAndSlides_USE_VECTOR
			_data.resize(0);
			_dataSized = false;
#else
			if (!!_data && _sizes[0] < sizes[0]) {
				auto data = _data;
				_data = nullptr;
				for (size_t i = 0; i < _sizes[0]; i++) assert(!data[i]);
				vDelete(data);
			}
#endif
			for (size_t i = 0; i < _sizes_size; i++) _sizes[i] = sizes[i];
		}
		bool hasData() const {
#ifdef  MweightsForSomeSpinsAndSlides_USE_VECTOR
			return _dataSized;
#else
			return _data != nullptr;
#endif

		}
		size_t childCount() const {
			size_t sum = 0;
			if (hasData()) {
				for (size_t i = 0; i < _sizes[0]; i++) if (_data[i]) sum++;
			}
			return sum;
		}
		size_t bytesUsed() const {
			size_t sum = sizeof(*this);
			if (hasData()) {
				for (size_t i = 0; i < _sizes[0]; i++) if (_data[i]) sum += _data[i]->bytesUsed();
			}
			return sum;
		}
		void uninit() {
			if (!hasData()) return;
			for (size_t i = 0; i < _sizes[0]; i++) if (_data[i]) _data[i]->uninit();
		}
		void clear() {
			if (!hasData()) return;
			for (size_t i = 0; i < _sizes[0]; i++) if (_data[i]) _data[i]->clear();
		}
		void zero() {
			if (!hasData()) return;
			for (size_t i = 0; i < _sizes[0]; i++) if (_data[i]) _data[i]->zero();
		}
		void fini() {
			if (!hasData()) return;
			deleteAllExistingValues_noLockNeeded();
#ifdef  MweightsForSomeSpinsAndSlides_USE_VECTOR
			_dataSized = false;
			_data.resize(0);
#else
			vDelete(_data);
#endif
		}

		void deleteExistingValue(size_t i, bool interesting) {
            if (!_data[i]) return;
			_exp_lock.acquire(__FILE__,__LINE__);													// keep the lock for as short a time as possible
			auto n = _data[i];
			_data[i] = nullptr;
			_exp_lock.release();
			n->deleteAllExistingValues_noLockNeeded();
			n->mui_disconnect();
			NodePoolTemplate<NextLevel>::insert(n);
		}

		void deleteAllExistingValues_noLockNeeded() {
			if (!_data) return;
			for (size_t i = 0; i < _sizes[0]; i++) {
				auto n = _data[i];
				if (!n) continue;
				_data[i] = nullptr;
				n->deleteAllExistingValues_noLockNeeded();
				n->mui_disconnect();
				NodePoolTemplate<NextLevel>::insert(n);
			}
			assert(childCount() == 0);
		}

		NextLevel* value(size_t i, MweightsUpperIndex const & mui, bool createIfNullptr = true) {
			// Create the _data if necessary
			if (i >= _sizes[0] || !hasData()) {
				assert(_sizes[0]);
				if (!hasData()) {
					if (!createIfNullptr) return nullptr;
					_exp_lock.acquire(__FILE__, __LINE__);
					if (!hasData()) {
#ifdef  MweightsForSomeSpinsAndSlides_USE_VECTOR
						_data.resize(_sizes[0]);
						for (size_t i = 0; i < _sizes[0]; i++) _data[i] = nullptr;
						_dataSized = true;
#else
						auto data = vNew(NextLevelPtr, _sizes[0]);
						for (size_t i = 0; i < _sizes[0]; i++) data[i] = nullptr;
						_data = data;
#endif
					}
					_exp_lock.release();
				}
				assert(i < _sizes[0]);
			}
			// Create the _data[i] if necessary
			auto howMade = "initial _data[i]";		// debugging only
			auto result = _data[i];
			if (!result) {
				if (!createIfNullptr) return nullptr;
				_exp_lock.acquire(__FILE__, __LINE__);
				howMade = "locked _data[i]";
				if (!(result = _data[i])) {
					howMade = "NextLevel::make";
					_data[i] = result = NextLevel::make(parent(), mui, _sizes_size - 1, &_sizes[1]);
				}
				_exp_lock.release();
			}
			assert(result->mui().matches(mui));
			// Done
			return result;
		}
	};

	class NodeForMweightsForSomeSpinsAndSlides : public MweightsForSomeSpinsAndSlides {
	public:
		// NOTE : 
		// Bevin : should not be used, should be returned to the pool !!!
		// YongBei : need to delete if it run out of memory
		~NodeForMweightsForSomeSpinsAndSlides() {};
		NodeForMweightsForSomeSpinsAndSlides* _next;	// used for the pool

		static const size_t	   _sizes_size = 0;

		NodeForMweightsForSomeSpinsAndSlides() 
		  : _next(nullptr) {
		}
		void setParent(Exp_Mweight_new* parent, MweightsUpperIndex const & mui, size_t sizes_size, size_t* sizes) {
			MweightsForSomeSpinsAndSlides::setParent(parent, mui);
			assert(sizes_size == _sizes_size);
			assert(sizes_size == 0);
		}
		size_t childCount() const {
			return size();
		}
		void deleteAllExistingValues_noLockNeeded() { 
			clear();
			assert(childCount() == 0);
		}

		static NodeForMweightsForSomeSpinsAndSlides* make(Exp_Mweight_new* parent, MweightsUpperIndex const & mui, size_t sizes_size, size_t* sizes) {
			auto n = NodePoolTemplate<NodeForMweightsForSomeSpinsAndSlides>::remove();
			if (!n) n = sNew(NodeForMweightsForSomeSpinsAndSlides);
			n->setParent(parent, mui, sizes_size, sizes);
			assert(n->mui().matches(mui));
			return n;
		}
	};

    typedef NodeTemplate<NodeForMweightsForSomeSpinsAndSlides,1>	ForEachPsi;
    typedef NodeTemplate<ForEachPsi,2>								ForEachDir;
    typedef NodeTemplate<ForEachDir,3>								ForEachClass;
    typedef NodeTemplate<ForEachClass,4>							ForEachImage;
};

using namespace Exp_MWeight_Namespace;


class Exp_Mweight_new::SparseData : public ForEachImage {
public:
	SparseData* _next;
	SparseData() : _next(nullptr) {
		NodePoolTemplate<NodeForMweightsForSomeSpinsAndSlides>::setNameEtc("NodeForMweightsForSomeSpinsAndSlides");
		NodePoolTemplate<ForEachPsi>						  ::setNameEtc("ForEachPsi");
		NodePoolTemplate<ForEachDir>						  ::setNameEtc("ForEachDir");
		NodePoolTemplate<ForEachClass>						  ::setNameEtc("ForEachClass");
		NodePoolTemplate<SparseData>						  ::setNameEtc("SparseData");
	}
	static SparseData* make(Exp_Mweight_new* parent, size_t sizes_size, size_t* sizes) {
		auto p = NodePoolTemplate<SparseData>::remove();
		if (!p) p = sNew(SparseData);
		p->setParent(parent, MweightsUpperIndex(), sizes_size, sizes);
		return p;
	}
};


Exp_Mweight_new::Exp_Mweight_new() :
	_sparseData(nullptr),
	nr_images (_sparseDataDimensions[0]),	// fan out fast for fewest lock collisions, also to implement clear_image easily
	nr_classes(_sparseDataDimensions[1]),	// also try to create the fewest leaves
	nr_dir    (_sparseDataDimensions[2]),
	nr_psi    (_sparseDataDimensions[3]),
	nr_over_rot	 (_nr_over_rot  ),
	nr_trans	 (_nr_trans     ),
	nr_over_trans(_nr_over_trans),

	current_fn("NULL"),
	os_file	  (),
	max_nr_significant_ihidden(0),
	max_density(0),
	max_size   (0)
{
	for (size_t i = 0; i < numberOfSparseDimensions; i++) _sparseDataDimensions[i] = 0;
}

Exp_Mweight_new::~Exp_Mweight_new() {
	fini();
	NodePoolTemplate<SparseData>::insert(_sparseData);
}

void Exp_Mweight_new::init(
	int nr_images,
	int nr_classes,
	int nr_dir,
	int nr_psi,
	// The following size the MweightsForSomeSpinsAndSlides
	// int elementCapacity
	int nr_over_rot,
	int nr_trans,
	int nr_over_trans)
{
	if (true
		&& this->nr_images     == nr_images
		&& this->nr_classes    == nr_classes
		&& this->nr_dir        == nr_dir
		&& this->nr_psi        == nr_psi
		&& this->nr_over_rot   == nr_over_rot
		&& this->nr_trans      == nr_trans
		&& this->nr_over_trans == nr_over_trans) {

		if (_sparseData) _sparseData->uninit();

		return;
	}

	fini();

	static_assert(numberOfSparseDimensions == 4, "assuming 4 dims");
	_sparseDataDimensions[0] = nr_images;		assert(this->nr_images     == nr_images);
	_sparseDataDimensions[1] = nr_classes;		assert(this->nr_classes    == nr_classes);
	_sparseDataDimensions[2] = nr_dir;			assert(this->nr_dir        == nr_dir);
	_sparseDataDimensions[3] = nr_psi;			assert(this->nr_psi        == nr_psi);
	_nr_over_rot             = nr_over_rot;		assert(this->nr_over_rot   == nr_over_rot);
	_nr_trans                = nr_trans;		assert(this->nr_trans      == nr_trans);
	_nr_over_trans           = nr_over_trans;	assert(this->nr_over_trans == nr_over_trans);

	_sparseData = SparseData::make(this, numberOfSparseDimensions, &_sparseDataDimensions[0]);
}

void Exp_Mweight_new::fini() {
	if (_sparseData) _sparseData->fini();
}

void Exp_Mweight_new::clear(int iimage) {
#ifdef Mweight_DEBUGGING
	debuggingMweights.clear(iimage);
#endif
	auto deeper = _sparseData->value(iimage, MweightsUpperIndex(iimage), false);
	if (deeper) {
		// deeper->clear();
		_sparseData->deleteExistingValue(iimage, false);
	}
}

void Exp_Mweight_new::clearAllImages() {
	//#pragma omp parallel for
	for (int iimage = 0; iimage < nr_images; iimage++) {
		clear(iimage);
	}
}

void Exp_Mweight_new::statisticHeap() {
	// Actually it not really clear the data,it only push the data to 'thread Heap', and pull the data back when needed.
	// When 'thread Heap' is empty , it need to create new data
	// Because the imbalace of image distribution,it always need to create new data for some empty 'thread Heap'
	// this will cause the program need more memory!
	// NodePoolTemplate<NodeForMweightsForSomeSpinsAndSlides>::clearPool(0.5);
	std::cout << __FILE__ << " " << __LINE__ << " statisticHeap() " << std::endl;
	std::cout << NodePoolTemplate<NodeForMweightsForSomeSpinsAndSlides>::show(1) << " , " ;
	std::cout << NodePoolTemplate<ForEachPsi>::show(1) << " , ";
	std::cout << NodePoolTemplate<ForEachDir>::show(1) << " , ";
	std::cout << NodePoolTemplate<ForEachClass>::show(1) << " , ";
	std::cout << NodePoolTemplate<SparseData>::show(1) << std::endl;
}

bool Exp_Mweight_new::hasNodeForMweightsForSomeSpinsAndSlides() {
	return NodePoolTemplate<NodeForMweightsForSomeSpinsAndSlides>::hasNode();
}
//  void Exp_Mweight_new::zero(int iimage) {
//  	if (auto deeper = _sparseData->value(iimage, false)) deeper->zero();
//  }
//  
//  void Exp_Mweight_new::zeroAllImages() {
//  	for (int iimage = 0; iimage < nr_images; iimage++) {
//  		zero(iimage);
//  	}
//  }

MweightsForSomeSpinsAndSlides* Exp_Mweight_new::mweightsForSomeSpinsAndSlidesOrNullptr(int iimage, int iclass, int idir, int ipsi, bool forceNotNullPtr) {
	if (!_sparseData) EXIT_ABNORMALLY;
	auto a = _sparseData     ->value(iimage, MweightsUpperIndex(iimage),					 forceNotNullPtr);
	auto b = !a ? nullptr : a->value(iclass, MweightsUpperIndex(iimage, iclass),			 forceNotNullPtr);
	auto c = !b ? nullptr : b->value(idir,   MweightsUpperIndex(iimage, iclass, idir),		 forceNotNullPtr);
	auto d = !c ? nullptr : c->value(ipsi,   MweightsUpperIndex(iimage, iclass, idir, ipsi), forceNotNullPtr);
	if (!d) return nullptr;
#ifdef Mweight_DEBUGGING
	d->v.setDimensions(iimage, iclass, idir, ipsi);
#endif
	return d;
}

int Exp_Mweight_new::sum_sizes_non_zero(
	int iimage,
	int class_begin, int class_end,
	int dir_begin, int dir_end,
	int psi_begin, int psi_end) {

	int sum = 0;

	auto a = _sparseData->value(iimage, MweightsUpperIndex(iimage), false);
	if (!a) return sum;

	for (auto iclass = class_begin; iclass < class_end; iclass++) {
		auto b = a->value(iclass, MweightsUpperIndex(iimage, iclass), false);
		if (!b) continue;
		for (int idir = dir_begin; idir < dir_end; idir++) {
			auto c = b->value(idir, MweightsUpperIndex(iimage, iclass, idir), false);
			if (!c) continue;
			for (int ipsi = psi_begin; ipsi < psi_end; ipsi++) {
				auto d = c->value(ipsi, MweightsUpperIndex(iimage, iclass, idir, ipsi), false);
				if (d) sum += d->sizeNonZero();
			}
		}
	}

	return sum;
}

void Exp_Mweight_new::analysis(std::string output_fn, std::string note)
{
	double per_vector_size = 3.5;
	auto resetStatic = [&]() {
		max_nr_significant_ihidden = 0;
		max_density = 0;
		max_size = 0;
	};
	// write the head
	auto writeHead = [&]() {
		os_file << std::setw(20) << " " << std::setw(20) << "nr_signif_ihidden"
			<< std::setw(20) << "density" << std::setw(20) << "memory_size(GB)" << std::endl;
	};
	// write the tail
	auto writeTail = [&]() {
		os_file << std::setw(20) << "max above : " << std::setw(20) << max_nr_significant_ihidden
			<< std::setw(20) << max_density
			<< std::setw(20) << max_size << std::endl;
		os_file << std::setw(20) << " ---------- " << std::setw(20) << " --------------- "
			<< std::setw(20) << " ----------- "
			<< std::setw(20) << " ----------- " << std::endl;
	};
	//
	if (current_fn == "NULL") {
		writeTail();
		current_fn = output_fn;
		os_file.open((current_fn + ".txt").c_str(), std::ios::out);
		writeHead();
		resetStatic();
	}
	else if (current_fn != output_fn) {
		writeTail();
		os_file.close();
		current_fn = output_fn;
		os_file.open((current_fn + ".txt").c_str(), std::ios::out);
		writeHead();
		resetStatic();
	}
	// write out the info
	std::atomic<size_t> nr_significant_ihidden = 0;
	if (_sparseData) {
		#pragma omp parallel for
		for (int iimage = 0; iimage < nr_images; iimage++) {
			size_t nr_significant_ihidden_partial = 0;
			auto a = _sparseData->value(iimage, MweightsUpperIndex(iimage), false);
			if (!a) continue;
			for (int iclass = 0; iclass < nr_classes; iclass++) {
				auto b = a->value(iclass, MweightsUpperIndex(iimage, iclass), false);
				if (!b) continue;
				for (int idir = 0; idir < nr_dir; idir++) {
					auto c = b->value(idir, MweightsUpperIndex(iimage, iclass, idir), false);
					if (!c) continue;
					for (int ipsi = 0; ipsi < nr_psi; ipsi++) {
						auto d = c->value(ipsi, MweightsUpperIndex(iimage, iclass, idir, ipsi), false);
						if (!d) continue;
						nr_significant_ihidden_partial += d->size();
					}
				}
			}
			nr_significant_ihidden += nr_significant_ihidden_partial;
		}
	}

	double density = nr_significant_ihidden / (double(nr_images)*nr_classes*nr_dir*nr_psi*nr_over_rot*nr_trans*nr_over_trans);
	double sizeInGB = (_sparseData ? _sparseData->bytesUsed() : 0) / 1024. / 1024. / 1024.;
	os_file << std::setw(20) << note << std::setw(20) << nr_significant_ihidden
		<< std::setw(20) << density
		<< std::setw(20) << sizeInGB << "GB" << std::endl;
	// compare with the max
	if (nr_significant_ihidden > max_nr_significant_ihidden) {
		max_nr_significant_ihidden = nr_significant_ihidden;
		max_density = density;
		max_size = sizeInGB;
	}
}

#ifndef DONT_INCLUDE_SAMPLING
void Exp_Mweight_new::printSamplingInfo() {
	std::cout << "//" << std::setw(15) << "healpix_order," << std::setw(15) << "oversampling," << std::setw(10) << "nr_images," << std::setw(10) << "nr_class," \
		<< std::setw(15) << "nr_direction," << std::setw(15) << "nr_rotation," << std::setw(15) << "nr_over_rot" << std::setw(15) << "nr_trans," \
		<< std::setw(15) << "nr_over_trans," << std::setw(15) << "angle_precision," << std::setw(20) << "size," << std::setw(18) << "size(GB)" << std::endl;
	for (int i = 1; i < 11; i++) {
		HealpixSampler sampling3d;
		sampling3d.initialize(2, 10, -1, i);
		int adaptive_oversampling = 1;
		size_t nr_images = 8;
		size_t nr_classes = 8;
		size_t nr_dir = sampling3d.NrDir();
		size_t nr_rot = sampling3d.NrPsi();
		size_t nr_trans = sampling3d.NrTrans();
		size_t nr_over_rot = sampling3d.oversamplingFactorOrientations(adaptive_oversampling);
		size_t nr_over_trans = sampling3d.oversamplingFactorTranslations(adaptive_oversampling);
		size_t size = nr_images*nr_classes*nr_dir*nr_rot*nr_trans*nr_over_rot*nr_over_trans;
		double size_gb = size*8. / 1024. / 1024. / 1024.;
		std::cout << "//" << std::setw(14) << i << "," << std::setw(14) << adaptive_oversampling << "," << std::setw(9) << nr_images << "," << std::setw(9) << nr_classes << "," \
			<< std::setw(14) << nr_dir << "," << std::setw(14) << nr_rot << "," << std::setw(14) << nr_over_rot << "," << std::setw(14) << nr_trans << "," \
			<< std::setw(14) << nr_over_trans << "," << std::setw(14) << 360. / sampling3d.NrPsi(adaptive_oversampling) << "," << std::setw(20) << size << "," \
			<< std::setw(14) << size_gb << "(GB)" << std::endl;
	}
}
#endif

//------------ Exp_Mcoarse_Rot_significant_base
Exp_Mcoarse_Rot_significant_base::Exp_Mcoarse_Rot_significant_base(){}
Exp_Mcoarse_Rot_significant_base::~Exp_Mcoarse_Rot_significant_base() {
}
void Exp_Mcoarse_Rot_significant_base::init(int _nr_classes,int _nr_dir,int _nr_psi,int _nr_trans,int _nr_images){
    nr_classes = _nr_classes;nr_dir = _nr_dir;nr_psi = _nr_psi;
    nr_images = _nr_images;nr_trans = _nr_trans;
}
void Exp_Mcoarse_Rot_significant_base::fini(){
    nr_classes = 0;nr_dir = 0;nr_psi = 0;
    nr_images = 0;nr_trans = 0;
}

//------------ Exp_Mcoarse_Rot_significant_full  ------------------------------- //
Exp_Mcoarse_Rot_significant_full::Exp_Mcoarse_Rot_significant_full(){}
Exp_Mcoarse_Rot_significant_full::~Exp_Mcoarse_Rot_significant_full() {
}
void Exp_Mcoarse_Rot_significant_full::init(int _nr_classes,int _nr_dir,int _nr_psi,int _nr_trans,int _nr_images){
    Exp_Mcoarse_Rot_significant_base::init(_nr_classes, _nr_dir, _nr_psi, _nr_trans, _nr_images);
    exp_Mcoarse_significant_full.init(nr_classes*nr_dir*nr_psi, nr_trans*nr_images);
    exp_Rot_significant_full.init(nr_classes*nr_dir*nr_psi, 1);
}
void Exp_Mcoarse_Rot_significant_full::fini(){
    Exp_Mcoarse_Rot_significant_base::fini();
    exp_Mcoarse_significant_full.fini();
    exp_Rot_significant_full.fini();
}
void Exp_Mcoarse_Rot_significant_full::resetMcoarseSignificant()
{
	exp_Mcoarse_significant_full.fill_with_first_touch(false);
}
void Exp_Mcoarse_Rot_significant_full::allocMcoarseSignificant(int iclass,int idir,int ipsi)
{
    // donot free memory last time
}
const char* Exp_Mcoarse_Rot_significant_full::isMcoarseSignificantRptrAll(int iclass,int idir,int ipsi)
{
    const char* significant_rptr =  exp_Mcoarse_significant_full[(iclass*nr_dir + idir)*nr_psi + ipsi].rptrAll();
    return significant_rptr;
}
void Exp_Mcoarse_Rot_significant_full::assignMcoarseSignificant(int iclass,int idir,int ipsi,std::vector<int>& itrans_image,char value)
{
    for (auto &i : itrans_image) exp_Mcoarse_significant_full[(iclass*nr_dir + idir)*nr_psi + ipsi].wptrAll()[i] = value;
}
void Exp_Mcoarse_Rot_significant_full::resetRotSignificant()
{
    exp_Rot_significant_full.fill_with_first_touch(false);
}
void Exp_Mcoarse_Rot_significant_full::allocAndAssignRotSignificant(int iclass,int idir,int ipsi,char value)
{
    auto set_full =[&](int iclass,int idir,int ipsi,bool value)
    {
        exp_Rot_significant_full[(iclass*nr_dir+idir)*nr_psi+ipsi].wptrAll()[0] = value;
    };
    auto exp_Mcoarse_significant_rptr = isMcoarseSignificantRptrAll(iclass, idir, ipsi);
    // set exp_nr_images = nr_images,in real case,exp_nr_images < nr_images!!!
    int exp_nr_images = nr_images;
    for (int i = 0; i < nr_trans*exp_nr_images; i++)
    {
        if (exp_Mcoarse_significant_rptr[i]){
            set_full(iclass,idir,ipsi,value);
            return;
        }
    }
}
bool Exp_Mcoarse_Rot_significant_full::isRotSignificant(int iclass,int idir,int ipsi)
{
    return exp_Rot_significant_full[(iclass*nr_dir+idir)*nr_psi+ipsi].wptrAll()[0];
}

//--------------  Exp_Mcoarse_Rot_significant_sparse  ------------------------------- //
Exp_Mcoarse_Rot_significant_sparse::Exp_Mcoarse_Rot_significant_sparse(){}
Exp_Mcoarse_Rot_significant_sparse::~Exp_Mcoarse_Rot_significant_sparse() {
    assert(exp_Mcoarse_significant_sparse.size()==0);
    assert(exp_Rot_significant_sparse.size()==0);
}
void Exp_Mcoarse_Rot_significant_sparse::init(int _nr_classes,int _nr_dir,int _nr_psi,int _nr_trans,int _nr_images){
    Exp_Mcoarse_Rot_significant_base::init(_nr_classes, _nr_dir, _nr_psi, _nr_trans, _nr_images);
    {/*do nothing*/}
}
void Exp_Mcoarse_Rot_significant_sparse::fini(){
    Exp_Mcoarse_Rot_significant_base::fini();
    for (auto & v : exp_Mcoarse_significant_sparse) aFree(v.second);
    exp_Mcoarse_significant_sparse.clear();
    exp_Rot_significant_sparse.clear();
}
void Exp_Mcoarse_Rot_significant_sparse::resetMcoarseSignificant()
{
    for (auto & v : exp_Mcoarse_significant_sparse) aFree(v.second);
    exp_Mcoarse_significant_sparse.clear();
}
void Exp_Mcoarse_Rot_significant_sparse::allocMcoarseSignificant(int iclass,int idir,int ipsi)
{
    auto orientation = std::make_tuple(iclass,idir,ipsi);
    auto it = exp_Mcoarse_significant_sparse.find(orientation);
    ERROR_CHECK(it != exp_Mcoarse_significant_sparse.end(), "exp_Mcoarse_significant already alloc.")
    char* significant_itrans_images = (char*)aMalloc(nr_trans*nr_images, 64);
    auto newit = exp_Mcoarse_significant_sparse.insert(it, std::make_pair(orientation,significant_itrans_images));
    for (int i = 0; i < nr_trans*nr_images; i++)
        newit->second[i] = false;
}
const char* Exp_Mcoarse_Rot_significant_sparse::isMcoarseSignificantRptrAll(int iclass,int idir,int ipsi)
{
    const char* significant_rptr = nullptr;
    auto orientation = std::make_tuple(iclass,idir,ipsi);
    auto it = exp_Mcoarse_significant_sparse.find(orientation);
    if (it != exp_Mcoarse_significant_sparse.end()) {
        significant_rptr = &it->second[0];
    }
    return significant_rptr;
}
void Exp_Mcoarse_Rot_significant_sparse::assignMcoarseSignificant(int iclass,int idir,int ipsi,std::vector<int>& itrans_image,char value)
{
    auto orientation = std::make_tuple(iclass,idir,ipsi);
    auto it = exp_Mcoarse_significant_sparse.find(orientation);
    ERROR_CHECK(it == exp_Mcoarse_significant_sparse.end(), "cound found exp_Mcoarse_significant.")
    for (auto& i : itrans_image) it->second[i] = value;
}
void Exp_Mcoarse_Rot_significant_sparse::resetRotSignificant()
{
    exp_Rot_significant_sparse.clear();
}
void Exp_Mcoarse_Rot_significant_sparse::allocAndAssignRotSignificant(int iclass,int idir,int ipsi,char value)
{
    auto set_sparse = [&](int iclass,int idir,int ipsi,bool value)
    {
        auto orientation = std::make_tuple(iclass,idir,ipsi);
        auto it = exp_Rot_significant_sparse.find(orientation);
        assert(it==exp_Rot_significant_sparse.end());
        exp_Rot_significant_sparse.insert(it, std::make_pair(orientation,value));
    };
    auto exp_Mcoarse_significant_rptr = isMcoarseSignificantRptrAll(iclass, idir, ipsi);
    if (exp_Mcoarse_significant_rptr != nullptr){
        // set exp_nr_images = nr_images,in real case,exp_nr_images < nr_images!!!
        int exp_nr_images = nr_images;
        for (int i = 0; i < nr_trans*exp_nr_images; i++) {
            if (exp_Mcoarse_significant_rptr[i]){
                set_sparse(iclass,idir,ipsi,value);
                return;
            }
        }
    }
}
bool Exp_Mcoarse_Rot_significant_sparse::isRotSignificant(int iclass,int idir,int ipsi)
{
    return !( exp_Rot_significant_sparse.end() == exp_Rot_significant_sparse.find(std::make_tuple(iclass,idir,ipsi)) );
}



// --------------------------------- Exp_Mcoarse_Rot_significant_new ------------------------------- //
Exp_Mcoarse_Rot_significant_new::Exp_Mcoarse_Rot_significant_new(){}
Exp_Mcoarse_Rot_significant_new::~Exp_Mcoarse_Rot_significant_new(){}
void Exp_Mcoarse_Rot_significant_new::init(int _nr_classes,int _nr_dir,int _nr_psi,int _nr_trans,int _nr_images,int do_local_searching)
{
    use_full_data = !do_local_searching;
    use_sparse_data = do_local_searching;
    if (use_full_data) {
        exp_Mcoarse_Rot_significant_full = sNew(Exp_Mcoarse_Rot_significant_full);
        exp_Mcoarse_Rot_significant_full->init(_nr_classes, _nr_dir, _nr_psi, _nr_trans, _nr_images);
    }
    if (use_sparse_data) {
        exp_Mcoarse_Rot_significant_sparse = sNew(Exp_Mcoarse_Rot_significant_sparse);
        exp_Mcoarse_Rot_significant_sparse->init(_nr_classes, _nr_dir, _nr_psi, _nr_trans, _nr_images);
    }
}
void Exp_Mcoarse_Rot_significant_new::fini()
{
    if (use_full_data){
        exp_Mcoarse_Rot_significant_full->fini();
        sDelete(exp_Mcoarse_Rot_significant_full);
    }
    if (use_sparse_data){
        exp_Mcoarse_Rot_significant_sparse->fini();
        sDelete(exp_Mcoarse_Rot_significant_sparse);
    }
}
// exp_Mcoarse_significant
void Exp_Mcoarse_Rot_significant_new::resetMcoarseSignificant()
{
    if (use_full_data) exp_Mcoarse_Rot_significant_full->resetMcoarseSignificant();
    if (use_sparse_data) exp_Mcoarse_Rot_significant_sparse->resetMcoarseSignificant();
}
void Exp_Mcoarse_Rot_significant_new::allocMcoarseSignificant(int iclass,int idir,int ipsi)
{
    if (use_full_data) exp_Mcoarse_Rot_significant_full->allocMcoarseSignificant(iclass, idir, ipsi);
    if (use_sparse_data) exp_Mcoarse_Rot_significant_sparse->allocMcoarseSignificant(iclass, idir, ipsi);
}
const char* Exp_Mcoarse_Rot_significant_new::isMcoarseSignificantRptrAll(int iclass,int idir,int ipsi)
{
    assert(use_full_data || use_sparse_data);
	const char* ptr = nullptr;
    if (use_full_data) ptr = exp_Mcoarse_Rot_significant_full->isMcoarseSignificantRptrAll(iclass, idir, ipsi);
    if (use_sparse_data) ptr = exp_Mcoarse_Rot_significant_sparse->isMcoarseSignificantRptrAll(iclass, idir, ipsi);
	return ptr;
}
void Exp_Mcoarse_Rot_significant_new::assignMcoarseSignificant(int iclass,int idir,int ipsi,std::vector<int>& itrans_image,char value)
{
    assert(use_full_data || use_sparse_data);
    if (use_full_data) exp_Mcoarse_Rot_significant_full->assignMcoarseSignificant(iclass, idir, ipsi, itrans_image, value);
    if (use_sparse_data) exp_Mcoarse_Rot_significant_sparse->assignMcoarseSignificant(iclass, idir, ipsi, itrans_image, value);
}
// exp_Rot_significant
void Exp_Mcoarse_Rot_significant_new::resetRotSignificant()
{
    if (use_full_data) exp_Mcoarse_Rot_significant_full->resetRotSignificant();
    if (use_sparse_data) exp_Mcoarse_Rot_significant_sparse->resetRotSignificant();
}
void Exp_Mcoarse_Rot_significant_new::allocAndAssignRotSignificant(int iclass,int idir,int ipsi,char value)
{
    if (use_full_data) exp_Mcoarse_Rot_significant_full->allocAndAssignRotSignificant(iclass, idir, ipsi, value);
    if (use_sparse_data) exp_Mcoarse_Rot_significant_sparse->allocAndAssignRotSignificant(iclass, idir, ipsi, value);
    if (debug) {
        assert(use_full_data && use_sparse_data);
#pragma omp critical
        {
            ERROR_CHECK(exp_Mcoarse_Rot_significant_full->isRotSignificant(iclass, idir, ipsi) !=
                        exp_Mcoarse_Rot_significant_sparse->isRotSignificant(iclass, idir, ipsi), "isRotSignificant Not Equal.");
        }
    }
}
bool Exp_Mcoarse_Rot_significant_new::isRotSignificant(int iclass,int idir,int ipsi)
{
    if (debug) {
        assert(use_full_data && use_sparse_data);
#pragma omp critical
        {
            ERROR_CHECK(exp_Mcoarse_Rot_significant_full->isRotSignificant(iclass, idir, ipsi) !=
                        exp_Mcoarse_Rot_significant_sparse->isRotSignificant(iclass, idir, ipsi), "isRotSignificant Not Equal.");
        }
    }
	bool isRotSignificantOrNot = false;
	if (use_full_data) isRotSignificantOrNot = exp_Mcoarse_Rot_significant_full->isRotSignificant(iclass, idir, ipsi);
	if (use_sparse_data) isRotSignificantOrNot = exp_Mcoarse_Rot_significant_sparse->isRotSignificant(iclass, idir, ipsi);
	return isRotSignificantOrNot;
}

//void Exp_Mcoarse_Rot_significant_new::compare()
//{
//    assert(use_full_data && use_sparse_data);
//    // compare Mcoarse
//    for (int iimage = 0; iimage < nr_images; iimage++) {
//        for (int iclass = 0; iclass < nr_classes; iclass++) {
//            for (int idir = 0; idir < nr_dir; idir++) {
//                for (int ipsi = 0; ipsi < nr_psi; ipsi++) {
//                    //
//                    isMcoarseSignificantRptrAll()
//                }
//            }
//        }
//    }
//    // compare Rot
//    for (int iclass = 0; iclass < nr_classes; iclass++) {
//        for (int idir = 0; idir < nr_dir; idir++) {
//            for (int ipsi = 0; ipsi < nr_psi; ipsi++) {
//                //
//                
//            }
//        }
//    }
//}
