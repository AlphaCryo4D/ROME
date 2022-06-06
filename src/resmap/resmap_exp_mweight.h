#pragma once

// #define Mweight_DEBUGGING


#include "./resmap_util.h"
#include "./resmap_error.h"
#include "./resmap_string.h"
#include "./resmap_mpi.h"
#include "./resmap_time.h"
#include "./resmap_array_vector.h"

// this matrix for exp_Mweight(nr_images*nr_classes*nr_directions*nr_psi*nr_over_rot*nr_translations*nr_over_trans)
// this is the matrix size for different oversampling and healpix-order:
// healpix_order,  oversampling,nr_images, nr_class,  nr_direction,           nr_rotation,    nr_over_rot      nr_trans, nr_over_trans,angle_precision,               size,          size(GB)
//             1,             1,        8,        8,            48(12*4)  ,            12,             8,            81,             4,            15,            95551488,      0.711914(GB)
//             2,             1,        8,        8,           192(12*4^2),            24,             8,            81,             4,           7.5,           764411904,       5.69531(GB)
//             3,             1,        8,        8,           768(12*4^3),            48,             8,            81,             4,          3.75,          6115295232,       45.5625(GB)
//             4,             1,        8,        8,          3072(12*4^4),            96,             8,            81,             4,         1.875,         48922361856,         364.5(GB)
//             5,             1,        8,        8,         12288(12*4^5),           192,             8,            81,             4,        0.9375,        391378894848,          2916(GB)
//             6,             1,        8,        8,         49152(12*4^6),           384,             8,            81,             4,       0.46875,       3131031158784,         23328(GB)
//             7,             1,        8,        8,        196608(12*4^7),           768,             8,            81,             4,      0.234375,      25048249270272,        186624(GB)
//             8,             1,        8,        8,        786432(12*4^8),          1536,             8,            81,             4,      0.117188,     200385994162176,   1.49299e+06(GB)
//             9,             1,        8,        8,       3145728(12*4^9),          3072,             8,            81,             4,     0.0585938,    1603087953297408,   1.19439e+07(GB)
//            10,             1,        8,        8,      12582912(12*4^10),         6144,             8,            81,             4,     0.0292969,   12824703626379264,   9.55515e+07(GB)
// NOTE : TILE size for direction can be 1,4,16,16*16,....for 'NEST healpix scheme'
// also need consider the opposite direction will share the same memory access for the 3D volume

class Exp_Mweight_new;

class PackedOverRotTransOverTrans {
	static const int over_rot_end   = 1<<9;		// see HealpixSampler::oversamplingFactorOrientations (order <= 3)
	static const int over_trans_end = 1<<6;		// see HealpixSampler::oversamplingFactorTranslations (order <= 3)
	static const int trans_end      = 2147483647/(over_rot_end*over_trans_end);
	int v;
public:
	PackedOverRotTransOverTrans() 
#ifndef NDEBUG
		: v(-1) 
#endif
	{ 
		assert(!inited()); 
	}
	PackedOverRotTransOverTrans(
		int iover_rot,
		int itrans,
		int iover_trans) {
		assert(0 <= iover_rot   && iover_rot   < over_rot_end);
		assert(0 <= itrans      && itrans      < trans_end);
		assert(0 <= iover_trans && iover_trans < over_trans_end);
		v = (itrans*over_rot_end + iover_rot)*over_trans_end + iover_trans;
		assert(inited());
	}
#ifndef NDEBUG
	bool inited() { return v >= 0; }
#endif
	bool operator==(PackedOverRotTransOverTrans rhs) const { return v == rhs.v; }
	bool operator!=(PackedOverRotTransOverTrans rhs) const { return v != rhs.v; }
	void decode(
		int & iover_rot,
		int & itrans,
		int & iover_trans
	) {
		assert(inited());
		auto r = v;
		iover_trans = r % over_trans_end; r /= over_trans_end;
		iover_rot   = r % over_rot_end;   r /= over_rot_end;
		itrans      = r;
	}
	int itrans() {
		auto r = v / (over_trans_end*over_rot_end);
		return r;
	}
};

#ifdef Mweight_DEBUGGING
class MweightsForSomeSpinsAndSlides;

class DebuggingMweights {
	static const bool debuggingMweights = true;
	const char* _state;
	size_t		_count;
	bool setState(const char* to, bool interesting = false);
	std::ostream& showDimensions(MweightsForSomeSpinsAndSlides*);
public:
	DebuggingMweights();
	void clear(int iimage);
	void insert(MweightsForSomeSpinsAndSlides* in, size_t index, PackedOverRotTransOverTrans key, double to);
	void setValue(MweightsForSomeSpinsAndSlides* in, size_t index, double to);
	void setSignificant(
		int iimage, int iclass, int idir, int ipsi, int iover_rot, int itrans, int iover_trans, int index, bool value);
	void getSignificant(
		int iimage, int iclass, int idir, int ipsi, int iover_rot, int itrans, int iover_trans, int index, bool value);
private:
	int sgnt_iimage, sgnt_iclass, sgnt_idir, sgnt_ipsi, sgnt_iover_rot, sgnt_itrans, sgnt_iover_trans, sgnt_index;
};

extern DebuggingMweights debuggingMweights;
#endif

class MweightsUpperIndex {
#ifdef NDEBUG
public:
	MweightsUpperIndex() {}
	explicit MweightsUpperIndex(int iimage, int iclass = -1, int idir = -1, int ipsi = -1) {
	}
#else
	int _iimage, _iclass, _idir, _ipsi;
public:
	MweightsUpperIndex() : _iimage(-1), _iclass(-1), _idir(-1), _ipsi(-1) {}
	explicit MweightsUpperIndex(int iimage, int iclass = -1, int idir = -1, int ipsi = -1) : _iimage(iimage), _iclass(iclass), _idir(idir), _ipsi(ipsi) {
		assert(valid(true));
	}
	int  iimage() const { return _iimage; }
	int  iclass() const { return _iclass; }
	int  idir  () const { return _idir;   }
	int  ipsi  () const { return _ipsi;   }
	bool valid(bool allowPartial) const {
		bool result			= true;
		bool allowNegative  = allowPartial;
		if (ipsi  () < 0) result &= allowNegative; else allowNegative = false;
		if (idir  () < 0) result &= allowNegative; else allowNegative = false;
		if (iclass() < 0) result &= allowNegative; else allowNegative = false;
		if (iimage() < 0) result &= allowNegative; else allowNegative = false;
		return result;
	}
	bool matches(MweightsUpperIndex const & pattern) const {
		return 
			match(ipsi  (), pattern.ipsi  ()) &&
			match(idir  (), pattern.idir  ()) &&
			match(iclass(), pattern.iclass()) &&
			match(iimage(), pattern.iimage());
	}
#endif
private:
	static bool match(int m, int n) {
		return (m < 0) || (n < 0) || (m == n);
	}
};

std::string to_string(MweightsUpperIndex const & mui);

class ForMweightsUpperIndex : public NoCopy {
	Exp_Mweight_new*   _parent;
	MweightsUpperIndex _mui;
public:
	ForMweightsUpperIndex() : _parent(nullptr) {}
	void mui_disconnect() { 
		assert(!!_parent);
		_parent = nullptr;
		_mui = MweightsUpperIndex(); 
	}
	void mui_connect(Exp_Mweight_new* parent, MweightsUpperIndex const & mui, bool allowPartial = false) {
		assert(!_parent);
		assert(parent);
		assert(mui.valid(allowPartial));
		_parent = parent;
		_mui    = mui;
	}
	Exp_Mweight_new* parent() const	 { return _parent;       }
	MweightsUpperIndex const & mui() { return _mui; }
#ifndef NDEBUG
	int  mui_iimage() const			 { return _mui.iimage(); }
	int  mui_iclass() const			 { return _mui.iclass(); }
	int  mui_idir  () const			 { return _mui.idir();   }
	int  mui_ipsi  () const			 { return _mui.ipsi(); }
#endif
};


class MweightsForSomeSpinsAndSlides : public ForMweightsUpperIndex {
private:
	Exp_Mweight_new*					_parent;
	size_t 								_capacity;

//#define MweightsForSomeSpinsAndSlides_USE_VECTOR
#ifdef  MweightsForSomeSpinsAndSlides_USE_VECTOR
	std::vector<PackedOverRotTransOverTrans> _keys;
	std::vector<double>                      _data;
#else
	PackedOverRotTransOverTrans* 		_keys;
	double*								_data;
#endif

//#define MweightsForSomeSpinsAndSlides_INITED_CHECKS
#ifdef MweightsForSomeSpinsAndSlides_INITED_CHECKS
	size_t const						_initedSize;
#ifdef  MweightsForSomeSpinsAndSlides_USE_VECTOR
	std::vector<size_t>					_inited;
#else
	size_t*	const						_inited;
#endif
#endif

	size_t								_size;

#ifdef Mweight_DEBUGGING
public:
	int _iimage, _iclass, _idir, _ipsi;
	void setDimensions(int iimage, int iclass, int idir, int ipsi) {
		_iimage = iimage; _iclass = iclass; _idir = idir; _ipsi = ipsi;
	}
#endif

public:
	MweightsForSomeSpinsAndSlides();
	void setParent(Exp_Mweight_new* parent, MweightsUpperIndex const & mui);
	~MweightsForSomeSpinsAndSlides();

	size_t bytesUsed() const;
	void uninit();

	void clear();
		// removes all elements
	void zero();
		// zeros all elements

	void zeroSomeMakeSureRestUnused();
	// Zero is called for the significant elts
	// and I suspect the others are never used...
	// This helps verify that assumption
	void zero(
		int iover_rot,
		int itrans,
		int iover_trans) {
		PackedOverRotTransOverTrans key(iover_rot, itrans, iover_trans);
		for (size_t index = 0; index < _size; index++) {
			if (_keys[index] != key) continue;
			setValue(index, 0.0);
			return;
		}
		assert(false);
	}
	void zeroMatchingTrans(int itrans) {
		for (size_t index = 0; index < _size; index++) {
			if (_keys[index].itrans() != itrans) continue;
			setValue(index, 0.0);
		}
        assert(false);// ??
	}
	void insert(
		PackedOverRotTransOverTrans key,
		double to);
	void insert(
		int iover_rot,
		int itrans,
		int iover_trans,
		double to) {
		insert(PackedOverRotTransOverTrans(iover_rot, itrans, iover_trans), to);
	}
	size_t capacity() const { return _capacity; }
	size_t size    () const { return _size;     }
    size_t sizeNonZero() const {
        size_t __size = 0;
        for (int index = 0; index < _size; index++) {
            if (_data[index] != 0) {
                assert(_data[index] > 0);__size++;
            }
        }
        return __size;
    }
	PackedOverRotTransOverTrans key(size_t index) const {
		assert(index < _size);
		return _keys[validatedIndexR(index)];
	}
	double value(size_t index) const {
		return _data[validatedIndexR(index)];
	}
	void setValue(size_t index, double to) {
#ifdef Mweight_DEBUGGING
		debuggingMweights.setValue(this, index, to);
#endif
		_data[validatedIndexW(index)] = to;
	}
	void key(size_t index, int & iover_rot, int & itrans, int & iover_trans) const {
		key(index).decode(iover_rot, itrans, iover_trans);
	}
private:
	void inited(size_t index) {
#ifdef MweightsForSomeSpinsAndSlides_INITED_CHECKS
		auto i = index / (8 * sizeof(size_t));
		auto j = index % (8 * sizeof(size_t));
		assert(i < _initedSize);
		_inited[i] |= size_t(1) << j;
#endif
	}
	bool isInited(size_t index) const {
#ifndef MweightsForSomeSpinsAndSlides_INITED_CHECKS
		return true;
#else
		auto i = index / (8 * sizeof(size_t));
		auto j = index % (8 * sizeof(size_t));
		assert(i < _initedSize);
		return (_inited[i] & (size_t(1) << j));
#endif
	}
	size_t validatedIndexW(size_t index) {
		assert(index < _size);
#ifdef MweightsForSomeSpinsAndSlides_INITED_CHECKS
		inited(index);
#endif
		return index;
	}
	size_t validatedIndexR(size_t index) const {
		assert(index < _size);
#ifdef MweightsForSomeSpinsAndSlides_INITED_CHECKS
		assert(isInited(index));
#endif
		return index;
	}
};

class Exp_Mweight_new : public NoCopy {
	class SparseData;
	static const size_t numberOfSparseDimensions = 4;
	size_t          _sparseDataDimensions[numberOfSparseDimensions];
	int             _nr_over_rot;
	int             _nr_trans;
	int             _nr_over_trans;
	SparseData *    _sparseData;

	std::string		current_fn;
	std::ofstream	os_file;
	double			max_nr_significant_ihidden;
	double			max_density;
	double			max_size;

public:
	size_t const & nr_images;
	size_t const & nr_classes;
	size_t const & nr_dir;
	size_t const & nr_psi;
	int    const & nr_over_rot;
	int    const & nr_trans;
	int    const & nr_over_trans;

	Exp_Mweight_new();
	~Exp_Mweight_new();
	void init(
		int nr_images,
		int nr_classes,
		int nr_dir,
		int nr_psi,
		// The following size the MweightsForSomeSpinsAndSlides
		// int elementCapacity
		int nr_over_rot,
		int nr_trans,
		int nr_over_trans);
	void fini();
	
	void clear(int iimage);		// removes all elements
	void clearAllImages();
	void statisticHeap();
	bool hasNodeForMweightsForSomeSpinsAndSlides();
	//void zero(int iimage);	// zeros   all elements
	//void zeroAllImages();

	MweightsForSomeSpinsAndSlides* mweightsForSomeSpinsAndSlidesOrNullptr(int iimage, int iclass, int idir, int ipsi, bool forceNotNullPtr = false);
	MweightsForSomeSpinsAndSlides& mweightsForSomeSpinsAndSlides         (int iimage, int iclass, int idir, int ipsi) {
		return *mweightsForSomeSpinsAndSlidesOrNullptr(iimage, iclass, idir, ipsi, true);
	}

	int sum_sizes_non_zero(
		int iimage,
		int class_begin, int class_end,
		int dir_begin,   int dir_end,
		int psi_begin,   int psi_end);

public:

	void analysis(std::string output_fn, std::string note);
		// do some analysis for exp_Mweight....(the compare matrix)

#ifndef DONT_INCLUDE_SAMPLING
	void printSamplingInfo();
#endif
};

// memorize all significant points which found on coarse searching setp
class Exp_Mcoarse_Rot_significant_base : public NoCopy{
protected:
    int nr_classes;
    int nr_dir;
    int nr_psi;
    int nr_images;
    int nr_trans;
public:
    int get_nr_classes(){return nr_classes;}
    int get_nr_dir(){return nr_dir;}
    int get_nr_psi(){return nr_psi;}
    int get_nr_images(){return nr_images;}
    int get_nr_trans(){return nr_trans;}
public:
    Exp_Mcoarse_Rot_significant_base();
    ~Exp_Mcoarse_Rot_significant_base();
    void init(int _nr_classes,int _exp_nr_dir,int _exp_nr_psi,int exp_nr_trans,int nr_images);
    void fini();
    // exp_Mcoarse_significant
    virtual void resetMcoarseSignificant() = 0;
    // alloc memory space
    virtual void allocMcoarseSignificant(int iclass,int idir,int ipsi) = 0;
    // assign
    virtual void assignMcoarseSignificant(int iclass,int idir,int ipsi,std::vector<int>& itrans_image,char value) = 0;
    virtual const char* isMcoarseSignificantRptrAll(int iclass,int idir,int ipsi) = 0;
    // exp_Rot_significant
    virtual void resetRotSignificant() = 0;
    virtual void allocAndAssignRotSignificant(int iclass,int idir,int ipsi,char value) = 0;
    virtual bool isRotSignificant(int iclass,int idir,int ipsi) = 0;
};

class Exp_Mcoarse_Rot_significant_full : public Exp_Mcoarse_Rot_significant_base{
private:
    // nr_classes*exp_nr_dir*exp_nr_psi  exp_nr_trans*nr_images
    VectorOfArray2d< char > exp_Mcoarse_significant_full;
    VectorOfArray2d< char > exp_Rot_significant_full;
public:
    Exp_Mcoarse_Rot_significant_full();
    ~Exp_Mcoarse_Rot_significant_full();
    void init(int _nr_classes,int _exp_nr_dir,int _exp_nr_psi,int exp_nr_trans,int nr_images);
    void fini();
    // exp_Mcoarse_significant
    void resetMcoarseSignificant();
    void allocMcoarseSignificant(int iclass,int idir,int ipsi);
    void assignMcoarseSignificant(int iclass,int idir,int ipsi,std::vector<int>& itrans_image,char value);
    const char* isMcoarseSignificantRptrAll(int iclass,int idir,int ipsi);
    // exp_Rot_significant
    void resetRotSignificant();
    void allocAndAssignRotSignificant(int iclass,int idir,int ipsi,char value);
    bool isRotSignificant(int iclass,int idir,int ipsi);
};

class Exp_Mcoarse_Rot_significant_sparse : public Exp_Mcoarse_Rot_significant_base{
private:
    // nr_classes*exp_nr_dir*exp_nr_psi  exp_nr_trans*nr_images
    std::map< std::tuple<int, int, int>, char* > exp_Mcoarse_significant_sparse;// each image's exp_Mcoarse_significant
    std::map < std::tuple<int, int, int>, bool > exp_Rot_significant_sparse;
public:
    Exp_Mcoarse_Rot_significant_sparse();
    ~Exp_Mcoarse_Rot_significant_sparse();
    void init(int _nr_classes,int _exp_nr_dir,int _exp_nr_psi,int exp_nr_trans,int nr_images);
    void fini();
    // exp_Mcoarse_significant
    void resetMcoarseSignificant();
    void allocMcoarseSignificant(int iclass,int idir,int ipsi);
    void assignMcoarseSignificant(int iclass,int idir,int ipsi,std::vector<int>& itrans_image,char value);
    const char* isMcoarseSignificantRptrAll(int iclass,int idir,int ipsi);
    // exp_Rot_significant
    void resetRotSignificant();
    void allocAndAssignRotSignificant(int iclass,int idir,int ipsi,char value);
    bool isRotSignificant(int iclass,int idir,int ipsi);
};

class Exp_Mcoarse_Rot_significant_new : public NoCopy {
private:
    bool debug = false;//true;
    Exp_Mcoarse_Rot_significant_full* exp_Mcoarse_Rot_significant_full;
    Exp_Mcoarse_Rot_significant_sparse* exp_Mcoarse_Rot_significant_sparse;
public:
    bool use_full_data;
    bool use_sparse_data;
    Exp_Mcoarse_Rot_significant_new();
    ~Exp_Mcoarse_Rot_significant_new();
    void init(int _nr_classes,int _nr_dir,int _nr_psi,int _nr_trans,int _nr_images,int do_local_searching);
    void fini();
    // exp_Mcoarse_significant
    void resetMcoarseSignificant();
    void allocMcoarseSignificant(int iclass,int idir,int ipsi);
    void assignMcoarseSignificant(int iclass,int idir,int ipsi,std::vector<int>& itrans_image,char value);
    const char* isMcoarseSignificantRptrAll(int iclass,int idir,int ipsi);
    char* isMcoarseSignificantWptrAll(int iclass,int idir,int ipsi);
    // exp_Rot_significant
    void resetRotSignificant();
    void allocAndAssignRotSignificant(int iclass,int idir,int ipsi,char value);
    bool isRotSignificant(int iclass,int idir,int ipsi);
};
