/***************************************************************************
 *
 * Intel® Parallel Computing Center for Structural Biology
 * Principal Investigator : Youdong (Jack) Mao (Youdong_Mao@dfci.harvard.edu)
 * Dana-Farber Cancer Institute, Harvard Medical School and Peking University
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

#ifndef ARRAY_VECTOR_H_
#define ARRAY_VECTOR_H_

#include "./resmap_util.h"
#include "./resmap_memory.h"
#include "./resmap_macros.h"

template <class T>
void fill(std::vector<T> & v, T with) {
    for (auto i = v.begin(); i != v.end(); i++) i->fill(with);
}

template <class T>
void zero(std::vector<T> & v) {
    for (auto i = v.begin(); i != v.end(); i++) i->zero();
}

template <class T>
void copy(T* lhs_p, int lhs_len, T const * rhs_p, int rhs_len) {
    assert(lhs_len == rhs_len);
    for (int i = 0; i < lhs_len; i++) {
        lhs_p[i] = rhs_p[i];
    }
}

template <class T1, class T2>
static void copyConverted(T1* lhs_p, int lhs_len, T2 const * rhs_p, int rhs_len) {
    assert(lhs_len == rhs_len);
    for (int i = 0; i < lhs_len; i++) {
        lhs_p[i] = T1(rhs_p[i]);
    }
}
// Other Utilities
//
template<typename T>
T sumvec(T const *vec,int length){
    if (length == 0) return 0;
    auto sum = fabs(vec[0]);				// can sum a 'const T*' vector
    for (int i = 1; i < length; i++) {
        sum += fabs(vec[i]);
    }
    return sum;
}

template<typename T>
void zerovec(T & vec) {
    for (auto i = 0; i != vec.size(); i++)
        vec[i].zero();
}

template<typename T1, typename T2>
void fillvec(T1 & vec, T2 const & fill) {
    for (auto i = 0; i < vec.size(); i++)
        vec[i] = fill;
}

// The code uses lots of vectors of scalars, and vectors of vectors of scalars
// Support efficient operations on them
//
template <typename T>
class VectorOfSingle {
protected:
    bool   _written;
    size_t _capacity;
    size_t _size;
    T*     _ptr;
    VectorOfSingle()                           : _written(false), _size(0), _capacity(0), _ptr(nullptr) { }
    VectorOfSingle(int size)                   : _written(false), _size(0), _capacity(0), _ptr(nullptr) { }
    VectorOfSingle(VectorOfSingle const & rhs) : _written(false), _size(0), _capacity(0), _ptr(nullptr) { init(rhs.size()); *this = rhs; }
public:
    typedef T Elt;
    
    VectorOfSingle(std::vector<T> & v) : _capacity(0), _size(v.size()), _ptr(v.data()) { }
    // This hack is used to make the comparison code easier
    
    ~VectorOfSingle() {
        assert(_capacity == 0);
        assert(_size     == 0);
        assert(_ptr      == nullptr);
    }
    
    void fini() {
        _written = false;
    }
    
    // Read, Modify, Write
    // But until fixed most Write uses may be reading existing contents
    //
    bool written() const { return _written; }
    size_t size()  const { return _size;    }
    
    const Elt* rptr(int lenToRead) const { assert(lenToRead <= size());  assert(_written); return _ptr; }
    Elt*       mptr(int lenToWrite)      { assert(lenToWrite <= size()); assert(_written); return _ptr; }
    Elt*       wptr(int lenToWrite)      { assert(lenToWrite <= size()); _written = true;  return _ptr; }
    
    const Elt* rptrAll()           const {                               assert(_written); return _ptr; }
    Elt*       mptrAll()                 {                               assert(_written); return _ptr; }
    Elt*       wptrAll()                 {                               _written = true;  return _ptr; }
    
    T  operator[](int i) const { assert(0 <= i && i < _size); assert(_written); return _ptr[i]; }
    T& operator[](int i)	   { assert(0 <= i && i < _size); _written = true;  return _ptr[i]; }
    
    void operator=(VectorOfSingle const & rhs) {
        auto rp = rhs.rptrAll();
        assert(_size == rhs.size());
#pragma vector aligned
        for (int i = 0; i < _size; i++)
            _ptr[i] = rp[i];
        _written = true;
    }
};

template <typename T>
class VectorOfScalar : public VectorOfSingle<T> {
public:
    typedef T Elt;
    VectorOfScalar() : VectorOfSingle<T>() { }
    VectorOfScalar(int size) : VectorOfSingle<T>(size) { init(size); }
    VectorOfScalar(VectorOfScalar const & rhs) { init(rhs.size()); *this = rhs; }
    
    ~VectorOfScalar() {
        if (this->_capacity > 0) { Heap::freeScalars<Elt>(this->_ptr); this->_capacity = 0; }
        this->_size = 0;
        this->_ptr = nullptr;
    }
    
    void init(size_t size) {											// not resize because the values are not kept
        if (size == 0 || size > this->_capacity) {
            if (this->_capacity > 0) Heap::freeScalars(this->_ptr);		// == 0 when ctr with a vector<T>
            this->_size = this->_capacity = 0;
        }
        if (this->_capacity < size) {
            this->_capacity = size;
            this->_ptr = Heap::allocScalars<Elt>(this->_capacity, __FILE__, __LINE__);
        }
        this->_size = size;
        this->_written = false;												// scalars do not self-init
    }
    
    void fill(double with) {
#pragma vector aligned
        for (int i = 0; i < this->_size; i++) this->_ptr[i] = with;
        this->_written = true;
    }
    void zero() {
        fill(0.0);
    }
    void operator=(VectorOfScalar const & rhs) {
        auto rp = rhs.rptrAll();
        assert(this->_size == rhs.size());
#pragma vector aligned
        for (int i = 0; i < this->_size; i++)
            this->_ptr[i] = rp[i];
        this->_written = true;
    }
    void operator+=(VectorOfScalar const & rhs) {
        auto rp = rhs.rptrAll();
        assert(this->_size == rhs.size());
#pragma vector aligned
        for (int i = 0; i < this->_size; i++)
            this->_ptr[i] += rp[i];
        this->_written = true;
    }
    void operator*=(VectorOfScalar const & rhs) {
        auto rp = rhs.rptrAll();
        assert(this->_size == rhs.size());
#pragma vector aligned
        for (int i = 0; i < this->_size; i++)
            this->_ptr[i] *= rp[i];
        this->_written = true;
    }
    void operator/=(Elt rhs) {
#pragma vector aligned
        for (int i = 0; i < this->_size; i++)
            this->_ptr[i] /= rhs;
        this->_written = true;
    }
};

template <typename T>
class VectorOfStruct : public VectorOfSingle<T>, NoCopy {
public:
    typedef T Elt;
    VectorOfStruct()         : VectorOfSingle<T>()     { }
    VectorOfStruct(int size) : VectorOfSingle<T>(size) { init(size); }
    
    ~VectorOfStruct() {
        if (this->_capacity > 0) {	// == 0 when ctr with a vector<T>
			vDelete(this->_ptr); this->_capacity = 0; 
		}
        this->_size = 0;
        this->_ptr  = nullptr;
    }
    
    void init(size_t size) {											// not resize because the values are not kept
        if (size == 0 || size > this->_capacity) {
            if (this->_capacity > 0) vDelete(this->_ptr);					// == 0 when ctr with a vector<T>
            this->_size = this->_capacity = 0;
        }
        if (this->_capacity < size) {
            this->_capacity = size;
            this->_ptr = vNew(Elt,this->_capacity);
        }
        this->_size = size;
        this->_written = true;												// assumes structs are constructed
    }
};

template <typename T>
class VectorOfVector : public std::vector<T>, NoCopy {
public:
    VectorOfVector() {}
    VectorOfVector(int size) { init(size); }
    void init(size_t size) {
        resize(0);		// deletes previous contents
        resize(size);
    }
    void fill(double with) {
        for (auto i = this->begin(); i != this->end(); i++) (*i).fill(with);
    }
    void zero() {
        fill(0.0);
    }
private:
    void resize(size_t size) {
        size_t initialSize = this->size();
        std::vector<T>::resize(size);
        for (auto i = initialSize; i < this->size(); i++) (*this)[i].init();
    }
};

// warning:  vector<bool> is special cased by C++ so don't use it
//
typedef VectorOfScalar<char>   VectorOfChar;
typedef VectorOfScalar<int>    VectorOfInt;
typedef VectorOfScalar<double> VectorOfDouble;
typedef VectorOfScalar<float>  VectorOfFloat;

//
template<typename T>
class VectorOfArray1d : public VectorOfScalar<T>,NoCopy {
public:
    VectorOfArray1d(){} ~VectorOfArray1d(){fini();}
    void init(int size){VectorOfScalar<T>::init(size);}
    void fini(){VectorOfScalar<T>::init(0);}
    void fill(T with){VectorOfScalar<T>::fill(with);}
private: // undefine this slow operation,4X slow when call operator[]
    T  operator[](int i) const { assert(0 <= i && i < this->_size); assert(this->_written); return this->_ptr[i]; }
    T& operator[](int i)	   { assert(0 <= i && i < this->_size); this->_written = true;  return this->_ptr[i]; }
};

template<typename T>
class VectorOfArray2d : public std::vector< VectorOfScalar<T> >,NoCopy {
public:
    VectorOfArray2d(){} ~VectorOfArray2d(){fini();}
    void init(int size,int length){
        std::vector< VectorOfScalar<T> >::resize(size);
        for(int i = 0;i < size;i++) (*this)[i].init(length);
    }
    void fini(){std::vector< VectorOfScalar<T> >::resize(0);}
    void fill(T with){for(int i = 0;i < (*this).size();i++) (*this)[i].fill(with);}
    void fill_with_first_touch(double with){
#pragma omp parallel for
        for(int i = 0;i < (*this).size();i++) {
            (*this)[i].fill(with);
        }
    }
};

template<typename T>
class VectorOfArray3d : public std::vector< std::vector< VectorOfScalar<T> > >,NoCopy {
public:
    VectorOfArray3d(){} ~VectorOfArray3d(){fini();}
    void init(int nr_array2d,int size,int length){
		std::vector< std::vector< VectorOfScalar<T> > >::resize(nr_array2d);
        for (int n = 0 ; n < nr_array2d; n++) {
            (*this)[n].resize(size);
            for(int i = 0;i < size;i++) (*this)[n][i].init(length);
        }
    }
    void fini(){std::vector< std::vector< VectorOfScalar<T> > >::resize(0);}
    void fill(T with){
        for(int n = 0;n < (*this).size();n++)
            for (int i = 0; i < (*this)[n].size(); i++)
                (*this)[n][i].fill(with);
    }
    void fill_with_first_touch(double with){
#pragma omp parallel for
        for(int n = 0;n < (*this).size();n++){
            for (int i = 0; i < (*this)[n].size(); i++)
                (*this)[n][i].fill(with);
        }
    }
};

template<typename T>
class GenricImage : public VectorOfScalar<T>,NoCopy
{
public:
    GenricImage(){}
    ~GenricImage(){fini();}
    void init(int size1,int size2){assert(size1==size2);_width=size1;VectorOfScalar<T>::init(size1*size2);}
    void fini(){VectorOfScalar<T>::init(0);}
    void fill(T with){VectorOfScalar<T>::fill(with);}
    int width(){return _width;}
private:
    int _width;
};

#if defined(FLOAT_PRECISION)
typedef VectorOfFloat Image;
typedef VectorOfFloat VectorOfFDOUBLE;
#else
typedef VectorOfDouble Image;
typedef VectorOfDouble VectorOfFDOUBLE;
#endif

typedef std::vector<Image> Images;

#endif /* defined(ARRAY_VECTOR_H_) */
