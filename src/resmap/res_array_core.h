/***************************************************************************
 *
 * Intel® Parallel Computing Center for Structural Biology
 * Principal Investigator : Youdong (Jack) Mao (Youdong_Mao@dfci.harvard.edu)
 * Dana-Farber Cancer Institute, Harvard Medical School and Peking University
 *
 * Authors: "Jian Wang(wj_hust08@hust.edu.cn)"
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

#include "resmap_util.h"
#include "resmap_error.h"
#include "resmap_mpi.h"
#include "resmap_string.h"

#include "res_fft.h"
#include "res_traits.h"

#ifdef _WIN32
#undef max
#undef min
#endif

#define JN_ENABLE(_cond) typename std::enable_if<_cond, int>::type = 1
#define DIE(...) do {std::cerr << strMerge(__VA_ARGS__) << std::endl; EXIT_ABNORMALLY;} while(0)

using Shape = std::vector<int>;

enum { ARRAY_TAKE_MODE_CLIP };

template<typename _Val, typename _Derived> class BasicArray;
template<typename _Val> class Array;
template<typename _Val> class MapArray;
template<typename _Val> class SubArray;

template<typename T>
struct jn_is_array {
private:
    using F = typename std::decay<T>::type;
    template<typename _Val, typename _Derived> static std::true_type check(BasicArray<_Val, _Derived>);
    template<typename _Val> static std::true_type check(Array<_Val>);
    template<typename _Val> static std::true_type check(MapArray<_Val>);
    template<typename _Val> static std::true_type check(SubArray<_Val>);
    static std::false_type check(...);
public:
    enum { value = std::is_same<decltype(check(JN_SELF(F))), std::true_type>::value };
};

#define JN_IS_ARRAY(_type) jn_is_array<_type>::value

inline double square(double n) { return n*n; }

template<typename _T, typename... _Types>
struct JnNumTypes { enum { N = 1 + JnNumTypes<_Types...>::N }; };

template<typename _T>
struct JnNumTypes<_T> { enum { N = 1 }; };

template<typename _T, typename... _Types>
struct JnFirstType { using type = _T; };

namespace array_helper {

#define ARRAY_EACH(shape, ind) \
    for (Shape ind = array_helper::ind_begin(shape); array_helper::ind_in(ind, shape); array_helper::ind_next(ind, shape))

#define ARRAY_EACH_R(shape, ind) \
    for (Shape ind = array_helper::ind_end(shape); array_helper::ind_in(ind, shape); array_helper::ind_prev(ind, shape))

#define ARRAY3_EACH(shape, i, a, b, c) \
    for (int i = 0, a = 0; a < shape[0]; a++) \
    for (int b = 0; b < shape[1]; b++) \
    for (int c = 0; c < shape[2]; c++, i++)

#define ARRAY3_EACH_R(shape, i, a, b, c) \
    for (int i = shape[0] * shape[1] * shape[2] - 1, a = shape[0] - 1; a >= 0; a--) \
    for (int b = shape[1]-1; b >= 0; b--) \
    for (int c = shape[2]-1; c >= 0; c--, i--)

    inline Shape ind_begin(const Shape &shape) {
        return Shape(shape.size(), 0);
    }

    inline Shape ind_end(const Shape &shape) {
        Shape ind = shape;
        for (auto && n : ind) n--;
        return std::move(ind);
    }

    inline bool ind_in(const Shape &ind, const Shape &shape) {
        int dim = int(ind.size());
        for (int i = 0; i < dim; i++) if (ind[i] < 0 || ind[i] >= shape[i]) return false;
        return true;
    }

    inline void ind_next(Shape &ind, const Shape &shape) {
        int dim = int(ind.size());
        ind[dim-1]++;
        bool flag;
        do {
            flag = false;
            for (int i = dim-1; i > 0; i--) {
                if (ind[i] >= shape[i]) {
                    ind[i] = 0;
                    ind[i-1]++;
                    flag = true;
                }
            }
        } while (flag);
    }

    inline void ind_prev(Shape &ind, const Shape &shape) {
        int dim = int(ind.size());
        ind[dim-1]--;
        bool flag;
        do {
            flag = false;
            for (int i = dim-1; i > 0; i--) {
                if (ind[i] < 0) {
                    ind[i] = shape[i]-1;
                    ind[i-1]--;
                    flag = true;
                }
            }
        } while (flag);
    }

    template<typename... _Int>
    inline void set_shape_helper(int n_ind, Shape & shape, int a, _Int ...inds) {
        int dim = int(shape.size());
        shape[dim-n_ind] = a;
        set_shape_helper(n_ind-1, shape, inds...);
    }

    template<>
    inline void set_shape_helper<>(int n_ind, Shape & shape, int a) {
        if (n_ind != 1) DIE("Error in ", __LINE__);
        int dim = int(shape.size());
        shape[dim-n_ind] = a;
    }

    template<typename... _Int>
    inline int index_helper(int n_ind, const Shape & shape, const Shape & coeff, int a, _Int ...inds) {
        int dim = int(shape.size());
        int n = coeff[dim-n_ind];
        //        for (int i = int(dim - n_ind + 1); i < dim; i++) n *= shape[i];
        return a * n + index_helper(n_ind - 1, shape, coeff, inds...);
    }

    template<>
    inline int index_helper<>(int n_ind, const Shape & shape, const Shape & coeff, int a) {
        if (n_ind != 1) DIE("Error in ", __LINE__);
        return a;
    }

    template<typename... _Int>
    inline int index(const Shape & shape, const Shape &coeff, _Int ...inds) {
        return index_helper(int(shape.size()), shape, coeff, inds...);
    }

    inline int index(const Shape &shape, const Shape &coeff, const Shape &ind) {
        int dim = int(shape.size());
        int sum = 0;
        //        int k = 1;
        //        for (int i = dim-1; i >= 0; i--) {
        //            sum += ind[i] * k;
        //            k *= shape[i];
        //        }
        for (int i = 0; i < dim; i++) sum += coeff[i] * ind[i];
        return sum;
    }

} // namespace array_helper

template<typename _NumType> class SubArray;

/**
 * The basic array class.
 */
template<typename _NumType, typename _Derived>
class BasicArray {
public:
    using value_type = _NumType;
    using derived_type = _Derived;
    using type = BasicArray<value_type, derived_type>;

    /**
     * Get ith elements of the shape of the array.
     */
    int shape(int i) const {
        return m_shape[i];
    }

    /**
     * Get the shape of the array.
     */
    const Shape &shape() const {
        return m_shape;
    }

    /**
     * Get ith elements of the coefficients of the array.
     */
    int coeff(int i) const {
        return m_coeff[i];
    }

    /**
     * Get the coefficients of the array.
     */
    const Shape &coeff() const {
        return m_coeff;
    }

    /**
     * Get the ith element of the array.
     */
    virtual value_type &at(int ind) = 0;

    /**
     * Get the ith element of a constant array.
     */
    virtual const value_type &at(int ind) const = 0;

    /** Get the element of a constant array through indices.
     */
    value_type &at(const Shape &ind) {
        return at(array_helper::index(m_shape, m_coeff, ind));
    }

    /**
     * Get the element of a constant array through indices.
     */
    const value_type &at(const Shape &ind) const {
        return at(array_helper::index(m_shape, m_coeff, ind));
    }

    /**
     * Get the element of a constant array through indices.
     */
    template<typename... _Int>
    value_type &at(_Int ...inds) {
        if (JnNumTypes<_Int...>::N != dim()) DIE("Array::at error at ", __LINE__);
        return at(array_helper::index_helper(int(m_shape.size()), m_shape, m_coeff, inds...));
    }

    /**
     * Get the element of a constant array through indices.
     */
    template<typename... _Int>
    const value_type &at(_Int ...inds) const {
        if (JnNumTypes<_Int...>::N != dim()) DIE("Array::at error at ", __LINE__);
        return at(array_helper::index_helper(int(m_shape.size()), m_shape, m_coeff, inds...));
    }

    value_type &operator ()(const Shape &ind) { return at(ind); }
    const value_type &operator ()(const Shape &ind) const { return at(ind); }
    template<typename... _Int> value_type &operator ()(_Int ...inds) { return at(inds...); }
    template<typename... _Int> const value_type &operator ()(_Int ...inds) const { return at(inds...); }
    //template<typename... _Int> value_type &operator [](_Int ...inds) { return at(inds...); }
    //template<typename... _Int> const value_type &operator [](_Int ...inds) const { return at(inds...); }

//    int alloc_size() const { return std::accumulate(m_shape.begin(), m_shape.end(), 1, [](int a, int b){return a*b;}); }
//    int size() const { return size(); }
    int size() const { return m_size; }
    int dim() const { return int(m_shape.size()); }

    /**
     * Reshape the Array.
     */
    derived_type &reshape(const Shape & _shape) {
        if (size() != std::accumulate(_shape.begin(), _shape.end(), 1, [](int a, int b){return a*b;})) {
            DIE("Reshape error in ", __LINE__);
        }
        set_shape(_shape);
        return *(derived_type *)this;
    }

    /**
     * Return the flattened array.
     *
     * This method is equal to reshape({size()}).
     */
    derived_type &flatten() {
        set_shape({size()});
        return *(derived_type *)this;
    }

    /**
     * Set all values.
     *
     * Set the values of all the elements of the array to a same value.
     */
    template<typename _N>
    void set_all(const _N & n) {
        for (int i = 0; i < size(); i++) at(i) = n;
    }

    /**
     * Iterator class of array.
     */
    template<typename _Arr, typename _V>
    class iterator : public std::iterator<std::random_access_iterator_tag, _V> {
    public:
        using array_type = _Arr;
        using array_value_type = _V;
        using difference_type = int;
        using self = iterator<array_type, array_value_type>;

        array_type *array;
        int n;

        //iterator() : array(nullptr), n(-1) {}
        //iterator(array_type *_array) : array(_array), n(0) {}
        iterator(array_type *_array, int _n) : array(_array) { setn(_n); }
        iterator(const self &rhs) : array(rhs.array) { setn(rhs.n); }

        // Operators : misc
    public:
        inline self& operator=(const self &rhs) {setn(rhs.n); array = rhs.array; return *this;}
        inline self& operator+=(int rhs) { setn(n+rhs); return *this;}
        inline self& operator-=(int rhs) { setn(n-rhs); return *this;}
        inline array_value_type& operator*() {return array->at(n); }
        inline array_value_type* operator->() {return &(array->at(n)); }

        // Operators : arithmetic
        inline self& operator++() { setn(n+1); return *this;}
        inline self& operator--() { setn(n-1); return *this;}
        inline self operator++(int) {self tmp(*this); setn(n+1); return tmp;}
        inline self operator--(int) {self tmp(*this); setn(n-1); return tmp;}
        inline int operator-(const self& rhs) const {return n-rhs.n;}
        inline self operator+(int rhs) {return self(array, n+rhs);}
        inline self operator-(int rhs) {return self(array, n-rhs);}
        friend inline self operator+(int lhs, const self& rhs) {return self(rhs.array, lhs+rhs.n);}
        friend inline self operator-(int lhs, const self& rhs) {return self(rhs.array, lhs-rhs.n);}

        // Operators : comparison
        inline bool operator==(const self& rhs) {return n == rhs.n && array == rhs.array;}
        inline bool operator!=(const self& rhs) {return n != rhs.n || array != rhs.array;}
        inline bool operator>=(const self& rhs) {return n >= rhs.n && array == rhs.array;}
        inline bool operator<=(const self& rhs) {return n <= rhs.n && array == rhs.array;}
        inline bool operator> (const self& rhs) {return n >  rhs.n && array == rhs.array;}
        inline bool operator< (const self& rhs) {return n <  rhs.n && array == rhs.array;}

        void setn(int _n) { if (_n < 0 || _n >= array->size()) n = array->size(); else n = _n; }
    };

    /**
     * Return the begin iterator of the array.
     */
    iterator<type, value_type> begin() {
        return iterator<type, value_type>(this, 0);
    }

    /**
     * Return the begin iterator of the array.
     */
    iterator<const type, const value_type> begin() const {
        return iterator<const type, const value_type>(this, 0);
    }

    /**
     * Return the end iterator of the array.
     */
    iterator<type, value_type> end() {
        return iterator<type, value_type>(this, -1);
    }

    /**
     * Return the end iterator of the array.
     */
    iterator<const type, const value_type> end() const {
        return iterator<const type, const value_type>(this, -1);
    }

    /**
     * Return the position of the minimum element.
     */
    int argmin() const {
        int ind = 0;
        value_type min = at(0);
        for (int i = 1; i < size(); i++) if (min > at(i)) {
            min = at(i);
            ind = i;
        }
        return ind;
    }

    /**
     * Return the position of the maximum element.
     */
    int argmax() const {
        int ind = 0;
        value_type max = at(0);
        for (int i = 1; i < size(); i++) if (max < at(i)) {
            max = at(i);
            ind=i;
        }
        return ind;
    }

    /**
     * Return the number of true elements.
     */
    bool count() const {
        int n = 0;
        for (int i = 0; i < size(); i++) {
            if ((at(i))) n++;
        }
        return n;
    }

    /**
     * Return the number elements satisfying a specific condition.
     */
    template<typename _F>
    bool count(_F && f) const {
        int n = 0;
        for (int i = 0; i < size(); i++) {
            if (f(at(i))) n++;
        }
        return n;
    }

    /**
     * Return an array representing if there exists a true element along an axis.
     */
    template<typename _Val = bool>
    Array<_Val> axis_any(int axis) const {
        return axis_any(axis, [](const value_type &v)->bool{return v;});
    }

    /**
     * Return an array representing if there exists
     * an element satisfying a specific condition along an axis.
     */
    template<typename _Val = bool, typename _F>
    Array<_Val> axis_any(int axis, _F && f) const {
        int nd = dim();

        // Set shape
        Shape shape_;
        for (auto && n : shape()) if (n != axis) shape_.push_back(n);
       
        // Set coeff
        Shape coeff_(nd);
        int k = 1;
        for (int i = nd-1; i>=0; i--) {
            if (i == axis) {
                coeff_[i] = 0;
            }
            else {
                coeff_[i] = k;
                k *= shape(i);
            }
        }

        // Set value
        Array<_Val> output(shape_, false);
        int n = 0;
        ARRAY_EACH(shape(), ind) {
            if (f(at(n))) {
                int sum = 0;
                for (int i = 0; i < nd; i++) sum += coeff_[i]*ind[i];
                output[sum] = true;
            }
            n++;
        }
        return std::move(output);
    }

    /** Return if there exists a true element.  */
    bool any() const {
        for (int i = 0; i < size(); i++) {
            if ((at(i))) return true;
        }
        return false;
    }

    /** Return if there exists an element satisfying a specific condition.  */
    template<typename _F>
    bool any(_F && f) const {
        for (int i = 0; i < size(); i++) {
            if (f(at(i))) return true;
        }
        return false;
    }

    /** Return if all the elements are true.  */
    bool all() const {
        for (int i = 0; i < size(); i++) {
            if (!(at(i))) return false;
        }
        return true;
    }

    /** Return if all the elements satisfy a specific condition.  */
    template<typename _F>
    bool all(_F && f) const {
        for (int i = 0; i < size(); i++) {
            if (!f(at(i))) return false;
        }
        return true;
    }

    /** Return the minimum element.  */
    template<typename T = value_type>
    T min() const {
        T min  = at(0);
        for (int i = 1; i < size(); i++) if (min > at(i)) min = at(i);
        return min;
    }

    /** Return the maximum element.  */
    template<typename T = value_type>
    T max() const {
        T max  = at(0);
        for (int i = 1; i < size(); i++) if (max < at(i)) max = at(i);
        return max;
    }

    /** Return the summation.  */
    template<typename T = value_type>
    T sum() const {
        T sum  = 0;
        for (int i = 0; i < size(); i++) sum  += at(i);
        return sum;
    }

    /** Return the summation along an axis.  */
    template<typename T = value_type>
    Array<T> sum(int axis) const {
        int nd = dim();
        Shape shape2(nd-1);
        Shape coeff2(nd);
        int k = 1;
        for (int i = nd-1, j = nd-2; i >= 0; i--) {
            if (i != axis) {
                shape2[j] = shape(i);
                coeff2[i] = k;
                k *= shape(i);
                j--;
            }
            else {
                coeff2[i] = 0;
            }
        }
        Array<T> output(shape2, 0);
        k = 0;
        ARRAY_EACH(shape(), ind) {
            int n = 0;
            for (int i = 0; i < nd; i++) n += ind[i] * coeff2[i];
            output[n] += at(k);
            k++;
        }
        //output.selfDivide(shape(axis));
        return std::move(output);
    }

    /** Return the mean.  */
    template<typename T = value_type>
    T mean() const {
        T mean = 0;
        for (int i = 0; i < size(); i++) mean += at(i);
        return T(mean / (double)size());
    }

    /** Return the variance.  */
    template<typename T = value_type>
    T var() const {
        T m = mean();
        T sum = 0;
        for (int i = 0; i < size(); i++) sum += std::pow(at(i)-m, 2);
        return T(sum / (double)size());
    }

    /** Return the median.  */
    template<typename T = value_type>
    T median()  const {
        // TODO
    }

    /** Return the norm.  */
    template<typename T = value_type>
    T norm()  const {
        T norm = 0;
        for (int i = 0; i < size(); i++) norm += square(at(i));
        return std::sqrt(norm);
    }

    /** Return the trace.  */
    template<typename T = value_type>
    T trace() const {
        int sum = 0;
        for (auto && n : coeff()) sum += n;
        T trace = 0;
        for (int i = 0; i < size(); i += sum) trace += at(i);
        return trace;
    }

    /** Return the transposed array.  */
    template<typename _RT = Array<value_type>>
    _RT transpose() {
        int nd = dim();
        int a = shape(0);
        int b = shape(1);
        ERROR_CHECK(nd!=2, "Transpose operation is only supported by 2-dimensional matrix now!");
        _RT n({b,a});
        for (int i = 0; i < b; i++) for (int j = 0; j < a; j++) {
            n(i, j) = at(j, i);
        }
        return std::move(n);
    }

    /** Transpose the array itself.  */
    derived_type &selfTranspose() {
        int nd = dim();
        int a = shape(0);
        int b = shape(1);
        ERROR_CHECK(a!=b || nd!=2, "Transpose operation is only supported by 2-dimensional square matrix now!");
        for (int i = 0; i < a; i++) for (int j = i+1; j < b; j++) {
            std::swap(at(i,j), at(j,i));
        }
        return *(derived_type *)this;
    }

    /** Reverse the array itself.  */
    derived_type &selfReverse() {
        int l = size();
        for (int i = 0; i < l/2; i++) std::swap(at(i), at(l-1-i)); 
        return *(derived_type*)this;
    }

    /** Apply a function to every element of the array itself.  */
    template<typename _F>
    derived_type &selfApply(_F && f) {
        for (int i = 0; i < size(); i++) at(i) = f(at(i));
        return *(derived_type*)this;
    }

    /** Change the array itself to it's conjugation.  */
    derived_type &selfConj() {
        for (int i = 0; i < size(); i++) at(i) = std::conj(at(i));
        return *(derived_type*)this;
    }

    /** Change the array itself to it's square.  */
    derived_type &selfSquare() {
        for (int i = 0; i < size(); i++) at(i) *= at(i);
        return *(derived_type*)this;
    }

    /** Change the array itself to it's exponation.  */
    derived_type &selfExp   () {
        for (int i = 0; i < size(); i++) at(i) = std::exp(at(i));
        return *(derived_type*)this;
    }

    /** Change the array itself to it's logarithm.  */
    derived_type &selfLog   () {
        for (int i = 0; i < size(); i++) at(i) = std::log(at(i));
        return *(derived_type*)this;
    }

    /** Change the array itself to it's squared root.  */
    derived_type &selfSqrt  () {
        for (int i = 0; i < size(); i++) at(i) = std::sqrt(at(i));
        return *(derived_type*)this;
    }

/**
 * @name Arithmetic methods
 * @{
 */
#define ARRAY_DEF_METHOD(name) \
    template<typename _N> auto name(const _N & n) \
        -> typename std::enable_if<jn_is_array<_N>::value, derived_type &>::type \
    { \
        for (int i = 0; i < size(); i++) at(i) = ARRAY_METHOD(at(i), n(i)); \
        return *(derived_type*)this; \
    } \
    template<typename _N> auto name(const _N & n) \
        -> typename std::enable_if<!jn_is_array<_N>::value, derived_type &>::type \
    { \
        for (int i = 0; i < size(); i++) at(i) = ARRAY_METHOD(at(i), n); \
        return *(derived_type*)this; \
    }

#define ARRAY_METHOD(m, n) m + n
        ARRAY_DEF_METHOD(selfPlus)
#undef ARRAY_METHOD
#define ARRAY_METHOD(m, n) m + n
        ARRAY_DEF_METHOD(operator +=)
#undef ARRAY_METHOD
#define ARRAY_METHOD(m, n) m - n
        ARRAY_DEF_METHOD(selfMinus)
#undef ARRAY_METHOD
#define ARRAY_METHOD(m, n) m - n
        ARRAY_DEF_METHOD(operator -=)
#undef ARRAY_METHOD
#define ARRAY_METHOD(m, n) n - m
        ARRAY_DEF_METHOD(selfMinusedBy)
#undef ARRAY_METHOD
#define ARRAY_METHOD(m, n) m * n
        ARRAY_DEF_METHOD(selfTimes)
#undef ARRAY_METHOD
#define ARRAY_METHOD(m, n) m * n
        ARRAY_DEF_METHOD(operator *=)
#undef ARRAY_METHOD
#define ARRAY_METHOD(m, n) m / n
        ARRAY_DEF_METHOD(selfDivide)
#undef ARRAY_METHOD
#define ARRAY_METHOD(m, n) m / n
        ARRAY_DEF_METHOD(operator /=)
#undef ARRAY_METHOD
#define ARRAY_METHOD(m, n) (n / m)
        ARRAY_DEF_METHOD(selfDividedBy)
#undef ARRAY_METHOD
#define ARRAY_METHOD(m, n) (m && n)
        ARRAY_DEF_METHOD(selfLogicalAnd)
#undef ARRAY_METHOD
#define ARRAY_METHOD(m, n) (m || n)
        ARRAY_DEF_METHOD(selfLogicalOr)
#undef ARRAY_METHOD
#define ARRAY_METHOD(m, n) (m & n)
        ARRAY_DEF_METHOD(selfBitAnd)
#undef ARRAY_METHOD
#define ARRAY_METHOD(m, n) (m | n)
        ARRAY_DEF_METHOD(selfBitOr)
#undef ARRAY_METHOD
#define ARRAY_METHOD(m, n) std::pow(m, n)
        ARRAY_DEF_METHOD(selfPow)
#undef ARRAY_METHOD
#define ARRAY_METHOD(m, n) std::pow(n, m)
        ARRAY_DEF_METHOD(selfPowedBy)
#undef ARRAY_METHOD
#define ARRAY_METHOD(m, n) (m > n)
        ARRAY_DEF_METHOD(selfGt)
#undef ARRAY_METHOD
#define ARRAY_METHOD(m, n) (m < n)
        ARRAY_DEF_METHOD(selfLt)
#undef ARRAY_METHOD
#define ARRAY_METHOD(m, n) (m == n)
        ARRAY_DEF_METHOD(selfEq)
#undef ARRAY_METHOD
#define ARRAY_METHOD(m, n) (m >= n)
        ARRAY_DEF_METHOD(selfGe)
#undef ARRAY_METHOD
#define ARRAY_METHOD(m, n) (m <= n)
        ARRAY_DEF_METHOD(selfLe)
#undef ARRAY_METHOD
#define ARRAY_METHOD(m, n) n
        ARRAY_DEF_METHOD(operator =)
#undef ARRAY_METHOD
#define ARRAY_METHOD(m, n) n
        ARRAY_DEF_METHOD(assign)
#undef ARRAY_METHOD

#undef ARRAY_DEF_METHOD
/**
 * @}
 */

    /**
     * @name Methods that return SubArray
     * @{
     */

    /** Return the sub array. */
    SubArray<value_type> sub(const std::initializer_list<std::initializer_list<int>> &range) {
        return SubArray<value_type>(*this, range);
    }

    /** Return the sub array.  */
    SubArray<const value_type> sub(const std::initializer_list<std::initializer_list<int>> &range) const {
        return SubArray<const value_type>(*this, range);
    }

    /** Return the sub array.  */
    SubArray<value_type> sub(const std::vector<std::vector<int>> &range) {
        return SubArray<value_type>(*this, range);
    }

    /** Return the sub array.  */
    SubArray<const value_type> sub(const std::vector<std::vector<int>> &range) const {
        return SubArray<const value_type>(*this, range);
    }

    /** Return the sub array.  */
    SubArray<value_type> sub(const Array<int> &range) {
        return SubArray<value_type>(*this, range);
    }

    /** Return the sub array.  */
    SubArray<const value_type> sub(const Array<int> &range) const {
        return SubArray<const value_type>(*this, range);
    }

    /** Return the sub array.  */
    SubArray<value_type> sub(int val) {
        return SubArray<value_type>(*this, 0, val);
    }

    /** Return the sub array.  */
    SubArray<const value_type> sub(int val) const {
        return SubArray<const value_type>(*this, 0, val);
    }

    /** Return the sub array.  */
    SubArray<value_type> sub(int axis, int val) {
        return SubArray<value_type>(*this, axis, val);
    }

    /** Return the sub array.  */
    SubArray<const value_type> sub(int axis, int val) const {
        return SubArray<const value_type>(*this, axis, val);
    }

    /** Return the sub array.  */
    SubArray<value_type> take(int val, int axis, int mode = ARRAY_TAKE_MODE_CLIP) {
        return SubArray<value_type>::take(*this, val, axis, mode);
    }

    /** Return the sub array.  */
    SubArray<const value_type> take(int val, int axis, int mode = ARRAY_TAKE_MODE_CLIP) const {
        return SubArray<const value_type>::take(*this, val, axis, mode);
    }

    /** Return the sub array.  */
    template<typename _LS>
    SubArray<value_type> takes(_LS && val, int axis, int mode = ARRAY_TAKE_MODE_CLIP) {
        return SubArray<value_type>::takes(*this, val, axis, mode);
    }

    /** Return the sub array.  */
    template<typename _LS>
    SubArray<const value_type> takes(_LS && val, int axis, int mode = ARRAY_TAKE_MODE_CLIP) const {
        return SubArray<const value_type>::takes(*this, val, axis, mode);
    }

    /// }

    /** Print the identification.  */
    void identify(std::string name) const {
        MPI_LOG << name << " (";
        for (auto && n : shape()) MPI_LOG << n << ' ';
        MPI_LOG << ')'
            << " num:" << size()
            << " argmin:" << argmin()
            << " min:" << min()
            << " argmax:" << argmax()
            << " max:" << max()
            << " sum:" << sum<double>()
            << " mean:" << mean<double>()
            << std::endl;
    }

    /** Pretty print the array.  */
    void print(std::ostream &stream = std::cout) const {
#ifdef USEMPI
        if (MPI_IS_ROOT) {
#endif
        stream << std::right;
        int dim = int(shape().size());
        const auto & w = shape();
        std::vector<int> v(dim);
        int n = 1;
        for (int i = dim - 1; i >= 0; i--) {
            n *= shape(i);
            v[i] = n;
        }
        for (int i = 0; i < size(); i++) {
            for (int j = 0; j < dim; j++) {
                if (i % v[j] == 0) {
                    stream << '[';
                }
                else if (i % v[dim-1] == 0) {
                    stream << ' ';
                    int d = (i % v[j]) / v[j+1];
                    if (d >= 3 && d < w[j]-3) {
                        i += v[j+1]*(w[j]-6)-1;
                        stream << "...\n";
                        goto array_print_for;
                    }
                }
            }
            //if (i % v[dim-1] != v[dim-1]-1) stream << std::setw(5);
            stream << std::setw(5);
            if (i % v[dim-1] >= 3 && i % v[dim-1] < w[dim-1]-3) {
                stream << "...";
                i += w[dim-1]-7;
                continue;
            }
            else {
                stream << at(i);
            }
            for (int j = dim - 1; j >= 0; j--) {
                if (i % v[j] == v[j] - 1) stream << ']';
            }
            if (i % v[dim-1] == v[dim-1]-1) stream << '\n';
            else stream << ' ';
array_print_for:;
        }
#ifdef USEMPI
        } // MPI_IS_ROOT
#endif
    }

    /** Print all the elements of the array.  */
    void printFull(std::ostream &stream = std::cout) const {
#ifdef USEMPI
        if (MPI_IS_ROOT) {
#endif
        int k = m_shape[dim()-1];
        for (int i = 0; i < size(); i++) {
            stream << at(i);
            if (i % k == k-1) stream << "\n";
            else stream << "\t";
        }
#ifdef USEMPI
        } // MPI_IS_ROOT
#endif
    }

protected:
    Shape m_shape;
    Shape m_coeff;
    int m_size;

    /**
     * Set the shape of the Array.
     *
     * This method is equal to reshape
     */
    void set_size() {
        m_size = std::accumulate(m_shape.begin(), m_shape.end(), 1, [](int a, int b){return a*b;});
    }

//    /**
//     * Set the shape of the Array.
//     *
//     * This method is equal to reshape
//     */
//    template<typename... _Int>
//    void set_shape(_Int ...inds) {
//        int dim = JnNumTypes<_Int...>::N;
//        m_shape.resize(dim);
//        array_helper::set_shape_helper(dim, m_shape, inds...);
//        set_size();
//        set_coeff();
//    }

    /**
     * Set the shape of the Array.
     *
     * This method is equal to reshape
     */
    void set_shape(const Shape & _shape) {
        m_shape = _shape;
        set_size();
        set_coeff();
    }

    /**
     * Set the coefficients of the Array.
     */
    void set_coeff() {
        int k = 1;
        int d = dim();
        m_coeff.resize(d);
        for (int i = d-1; i>=0; i--) {
            m_coeff[i] = k;
            k *= m_shape[i];
        }
    }

    template<typename _T, typename _V>
    int get_ind(_T && _inds, _V && _shape) {
        int dim = _inds.size();
        int n = 0;
        int m = 1;
        for (int i = dim - 1; i >= 0; i--) {
            n += _inds[i] * m;
            m *= _shape[i];
        }
        return n;
    }

    template<typename _T, typename _V>
    bool next_ind(_T && _inds, _V && _range) {
        _inds.back()++;
        int dim = _inds.size();
        bool flag;
        do {
            flag = false;
            for (int i = dim - 1; i > 0; i--) {
                if (_inds[i] >= _range[i][1]) {
                    _inds[i] = _range[i][0];
                    _inds[i-1]++;
                    flag = true;
                    break;
                }
            }
        } while (flag);
        return _inds[0] < _range[0][1];
    }

};

/** Operator << reload for array.
 */
template<typename _Derived, typename _Val>
std::ostream &operator <<(std::ostream &stream, const BasicArray<_Val, _Derived> &array) {
    array.print(stream);
    return stream;
}

/** Operator == reload for array.
 */
template<typename _V, typename _D1, typename _D2>
bool operator ==(const BasicArray<_V, _D1> &a1, const BasicArray<_V, _D2> &a2) {
    if (a1.shape() != a2.shape()) return false;
    for (int i = 0; i < a1.size(); i++) {
        if (a1(i) != a2(i)) return false;
    }
    return true;
}

/** Class Array.
 */
template<typename _NumType>
struct Array : public BasicArray<_NumType, Array<_NumType>> {
public:
    using value_type = _NumType;
    using type = Array<value_type>;
    using base_type = BasicArray<value_type, type>;
protected:
	value_type *data;
public:
    virtual value_type &at(int a) { return data[a]; }
    virtual const value_type &at(int a) const { return data[a]; }

	value_type operator [](int ind) const { return data[ind]; }
	value_type &operator [](int ind)  { return data[ind]; }

	const value_type* data_rptr() const {return data;}
	value_type* data_wptr() { return data; }
    /**
     * @name Array constructors
     * @{
     */
    /**
     * Default constructor.
     */
    Array() { init(); }

    /**
     * Copy constructor.
     */
    Array(const type &array) { init(); _assign(array); }

    /**
     * Move constructor.
     */
    Array(type &&array) { init(); swap(array); }

    /**
     * Extended copy constructor.
     */
    template<typename _T, JN_ENABLE(JN_IS_ARRAY(_T))>
    Array(const _T &array) { init(); _assign(array); }

    /**
     * Copy assignment.
     */
    type &operator =(const type &array) { _assign(array); return *this; }

    /**
     * Move assignment.
     */
    type &operator =(type &&array) { swap(array); return *this; }

    /**
     * Extended copy assignment.
     */
    template<typename _T> type &operator =(const _T &array) { _assign(array); return *this; }

    /**
     * Constructor from a shape with the type of initializer_list.
     *
     * @code
     * Array array({2,2});
     * @endcode
     */
    Array(const std::initializer_list<int> & _shape) {
        base_type::set_shape(Shape(_shape));
        data = (value_type*)aMalloc(base_type::size()*sizeof(value_type), 64);
    }

    /**
     * Constructor from a shape.
     *
     * @code
     * Shape shape{2,2};
     * Array array(shape);
     * @endcode
     */
    Array(const Shape & _shape) {
        base_type::set_shape(_shape);
        data = (value_type*)aMalloc(base_type::size()*sizeof(value_type), 64);
    }

//    /**
//     * Constructor from multiple integer parameters.
//     *
//     * @code
//     * Array array(2,2);
//     * @endcode
//     */
//    template<typename... _Int, typename = typename std::enable_if<std::is_same<typename JnFirstType<_Int...>::type, int>::value>::type>
//    Array(_Int ...inds) {
//        base_type::set_shape(inds...);
//        data = (value_type*)aMalloc(base_type::size()*sizeof(value_type), 64);
//    }

    /**
     * Array constructor from shape and data.
     *
     * @code
     * Array array({2,2},{1,2,3,4});
     * @endcode
     */
    Array(const Shape & _shape, const std::initializer_list<value_type> & _data) {
        base_type::set_shape(_shape);
        data = (value_type*)aMalloc(base_type::size()*sizeof(value_type), 64);
        int j = 0; for (auto && n : _data) { data[j] = n; j++; }
    }

    /**
     * Array constructor from shape and value.
     *
     * @code
     * Array array({2,2},0);
     * @endcode
     */
    Array(const Shape & _shape, const value_type & _value) {
        base_type::set_shape(_shape);
        int n = base_type::size();
//        data = (value_type*)aMalloc(n*sizeof(value_type), 64);
        alloc_data(n);
        for (int i = 0; i < n; i++) data[i] = _value;
    }

    /**
     * Array constructor from range with the type of std::vector<std::vector<int>>.
     *
     * @code
     * Array a({3,3},{1,2,3,4,5,6,7,8,9});
     * Array b(a, {{1,2},{2,3}});
     * @endcode
     */
    template<typename _V, typename _D>
    Array(const BasicArray<_V, _D> & array, const std::vector<std::vector<int>> &range) {
        int dim = array.dim();
        if (range.size() != dim) DIE("Array constructor error at ", __LINE__);

        Shape shape(dim);
        for (int i = 0; i < dim; i++) shape[i] = range[i][1] - range[i][0];
        base_type::set_shape(shape);

//        base_type::m_shape.resize(dim);
//        for (int i = 0; i < dim; i++) base_type::m_shape[i] = range[i][1] - range[i][0];
//        base_type::set_size();
//        base_type::set_coeff();
//        data.resize(base_type::size());
        alloc_data(base_type::size());

        std::vector<int> inds(dim);
        for (int i = 0; i < dim; i++) inds[i] = range[i][0];
        int n = 0;
        do {
            data[n] = array(ind(inds, array.shape()));
            n++;
        } while (base_type::next_ind(inds, range));
    }

    /**
     * Array constructor from range with the type of std::vector<std::vector<int>>.
     *
     * @code
     * Array a({3,3},{1,2,3,4,5,6,7,8,9});
     * Arrayi range({3}, {1,3,4,6});
     * Array b(a, range);
     * @endcode
     */
    template<typename _V, typename _D>
    Array(const BasicArray<_V, _D> & array, const Array<int> &range) {
        base_type::set_shape(range.shape());
        data.resize(base_type::size());
        for (int i = 0; i < base_type::size(); i++) data[i] = array(range(i));
    }
    /**
     * @}
     */

    ~Array() {
        clear();
    }

    /**
     * @name Don't use theses methods!!!
     * @{
     */
    void init() {
        base_type::set_shape({0});
        data = nullptr;
    }

//    template<typename... _Int>
//    void realloc(_Int ..._inds) {
//        clear();
//        base_type::set_shape(_inds...);
//        data = new value_type[base_type::size()];
//    }

    void realloc(const Shape &_shape) {
        clear();
        base_type::set_shape(_shape);
        data = (value_type*)aMalloc(base_type::size()*sizeof(value_type), 64);
    }

    void clear() {
        if (data != nullptr) aFree(data);
        init();
    }

    template<typename _T>
    void _assign(const _T &array) {
        if (base_type::size() != array.size()) realloc(array.shape());
        base_type::set_shape(array.shape());
        for (int i = 0; i < base_type::size(); i++) data[i] = array(i);
    }

    void swap(type &array) {
        Shape shape = array.shape();
        base_type::set_shape(shape);
        array.set_shape(base_type::m_shape);
        std::swap(data, array.data);
    }
    /**
     * @}
     */

    /**
     * @name FFT (fast fourior transformation) methods
     * @{
     */
    type &selfFFT()       {       ::resumap::fft(data, data, base_type::dim(), base_type::shape(0)); return *this; }
    type &selfIFFT()      {      ::resumap::ifft(data, data, base_type::dim(), base_type::shape(0)); return *this; }
    type &selfFFTshift()  {  ::resumap::fftshift(data, data, base_type::dim(), base_type::shape(0)); return *this; }
    type &selfFFT1()      {      ::resumap::fft1(data, data,                   base_type::shape(0)); return *this; }
    type &selfIFFT1()     {     ::resumap::ifft1(data, data,                   base_type::shape(0)); return *this; }
    type &selfFFTshift1() { ::resumap::fftshift1(data, data,                   base_type::shape(0)); return *this; }
    type &selfFFT2()      {      ::resumap::fft2(data, data,                   base_type::shape(0)); return *this; }
    type &selfIFFT2()     {     ::resumap::ifft2(data, data,                   base_type::shape(0)); return *this; }
    type &selfFFTshift2() { ::resumap::fftshift2(data, data,                   base_type::shape(0)); return *this; }
    type &selfFFT3()      {      ::resumap::fft3(data, data,                   base_type::shape(0)); return *this; }
    type &selfIFFT3()     {     ::resumap::ifft3(data, data,                   base_type::shape(0)); return *this; }
    type &selfFFTshift3() { ::resumap::fftshift3(data, data,                   base_type::shape(0)); return *this; }
    template<typename T>
    static type       fft(T && m) { type n(m.shape());       ::resumap::fft(m.data, n.data, m.dim(), m.shape(0)); return std::move(n); }
    template<typename T>
    static type      ifft(T && m) { type n(m.shape());      ::resumap::ifft(m.data, n.data, m.dim(), m.shape(0)); return std::move(n); }
    template<typename T>
    static type  fftshift(T && m) { type n(m.shape());  ::resumap::fftshift(m.data, n.data, m.dim(), m.shape(0)); return std::move(n); }
    template<typename T>
    static type      fft1(T && m) { type n(m.shape());      ::resumap::fft1(m.data, n.data,          m.shape(0)); return std::move(n); }
    template<typename T>
    static type     ifft1(T && m) { type n(m.shape());     ::resumap::ifft1(m.data, n.data,          m.shape(0)); return std::move(n); }
    template<typename T>
    static type fftshift1(T && m) { type n(m.shape()); ::resumap::fftshift1(m.data, n.data,          m.shape(0)); return std::move(n); }
    template<typename T>
    static type      fft2(T && m) { type n(m.shape());      ::resumap::fft2(m.data, n.data,          m.shape(0)); return std::move(n); }
    template<typename T>
    static type     ifft2(T && m) { type n(m.shape());     ::resumap::ifft2(m.data, n.data,          m.shape(0)); return std::move(n); }
    template<typename T>
    static type fftshift2(T && m) { type n(m.shape()); ::resumap::fftshift2(m.data, n.data,          m.shape(0)); return std::move(n); }

    template<typename _V>
    static type fft3(const Array< ::std::complex<_V> > & m) {
        type n(m.shape());
        ::resumap::fft3(m.data, n.data, m.shape(0));
        return std::move(n);
    }

    template<typename T>
    static type     ifft3(T && m) { type n(m.shape());     ::resumap::ifft3(m.data, n.data,          m.shape(0)); return std::move(n); }
    template<typename T>
    static type fftshift3(T && m) { type n(m.shape()); ::resumap::fftshift3(m.data, n.data,          m.shape(0)); return std::move(n); }
    /**
     * @}
     */

    /**
     * @name Return initialized array
     * @{
     */
    /**
     * Return an array with all the elements being 0
     */
    static type zero(const Shape &shape) { type m(shape); for (int i = 0; i < m.size(); i++) m(i) = 1; return std::move(m); }

    /**
     * Return an array with all the elements being 1
     */ 
    static type ones(const Shape &shape) { type m(shape); for (int i = 0; i < m.size(); i++) m[i] = 1; return std::move(m); }

    /**
     * Return an array with all the elements being value.
     */
    template<typename _V>
    static type constants(const Shape &shape, const _V &value) {
        type m(shape);
        for (int i = 0; i < m.size(); i++) m(i) = 1;
        return std::move(m);
    }
    /**
     * @}
     */

    /**
     * @name Array arithmetic methods
     *
     * @code
     * Array a({3},{1,2,3});
     * Array b({3},{3,4,5});
     * std::cout << Array::plus(a, b) << std::endl;
     * // [4, 6, 8]
     * @endcode
     *
     * @{
     */
	// TODO : use operator [] instead of operator () !!!
#define ARRAY_DEF_METHOD(name) \
    template<typename _M, typename _N> static auto name(const _M & m, const _N & n) \
    -> typename std::enable_if<jn_is_array<_M>::value && jn_is_array<_N>::value, type>::type{ \
        type r(m.shape()); \
        for (int i = 0; i < r.size(); i++) r[i] = ARRAY_METHOD(m(i), n(i)); \
        return std::move(r);\
    }\
    template<typename _M, typename _N> static auto name(const _M & m, const _N & n) \
    -> typename std::enable_if<jn_is_array<_M>::value && !jn_is_array<_N>::value, type>::type{ \
        type r(m.shape()); \
        for (int i = 0; i < r.size(); i++) r[i] = ARRAY_METHOD(m(i), n); \
        return std::move(r);\
    }\
    template<typename _M, typename _N> static auto name(const _M & m, const _N & n) \
    -> typename std::enable_if<!jn_is_array<_M>::value && jn_is_array<_N>::value, type>::type{ \
        type r(n.shape()); \
        for (int i = 0; i < r.size(); i++) r[i] = ARRAY_METHOD(m, n(i)); \
        return std::move(r);\
    }

#define ARRAY_METHOD(m, n) m + n
    ARRAY_DEF_METHOD(plus)
#undef ARRAY_METHOD

#define ARRAY_METHOD(m, n) m - n
    ARRAY_DEF_METHOD(minus)
#undef ARRAY_METHOD

#define ARRAY_METHOD(m, n) m * n
    ARRAY_DEF_METHOD(times)
#undef ARRAY_METHOD

#define ARRAY_METHOD(m, n) m / n
    ARRAY_DEF_METHOD(divide)
#undef ARRAY_METHOD

#define ARRAY_METHOD(m, n) m > n
    ARRAY_DEF_METHOD(gt)
#undef ARRAY_METHOD

#define ARRAY_METHOD(m, n) m >= n
    ARRAY_DEF_METHOD(ge)
#undef ARRAY_METHOD

#define ARRAY_METHOD(m, n) m < n
    ARRAY_DEF_METHOD(lt)
#undef ARRAY_METHOD

#define ARRAY_METHOD(m, n) m <= n
    ARRAY_DEF_METHOD(le)
#undef ARRAY_METHOD

#define ARRAY_METHOD(m, n) m == n
    ARRAY_DEF_METHOD(eq)
#undef ARRAY_METHOD

#define ARRAY_METHOD(m, n) m & n
    ARRAY_DEF_METHOD(bitAnd)
#undef ARRAY_METHOD

#define ARRAY_METHOD(m, n) m | n
    ARRAY_DEF_METHOD(bitOr)
#undef ARRAY_METHOD

#define ARRAY_METHOD(m, n) m && n
    ARRAY_DEF_METHOD(logicalAnd)
#undef ARRAY_METHOD

#define ARRAY_METHOD(m, n) m || n
    ARRAY_DEF_METHOD(logicalOr)
#undef ARRAY_METHOD

#define ARRAY_METHOD(m, n) std::pow(m, n)
    ARRAY_DEF_METHOD(pow)
#undef ARRAY_METHOD

#undef ARRAY_DEF_METHOD

    /**
     * @}
     */

    /**
     * Return an array with all the elements being mapped by the function.
     */
    template<typename _F, typename _M1>
    static type map(_F && f, _M1 && m1) {
        type r(m1.shape()); 
        for (int i = 0; i < r.size(); i++) r[i] = f(m1[i]); 
        return std::move(r);
    }

    /**
     * Return an array with all the elements being mapped by the function.
     */
    template<typename _F, typename _M1, typename _M2>
    static type map(_F && f, _M1 && m1, _M2 && m2) {
        ERROR_CHECK(m1.shape() != m2.shape(), "The shape of the two arrays should be same!");
        type r(m1.shape()); 
        for (int i = 0; i < r.size(); i++) r(i) = f(m1(i), m2(i)); 
        return std::move(r);
    }

    /**
     * Same as map(f, m1)
     */
    template<typename _M, typename _F>
    static type apply(_M && m1, _F && f) {
        type m2(m1.shape()); 
        for (int i = 0; i < m2.size(); i++) m2(i) = f(m1(i)); 
        return std::move(m2);
    }

    /**
     * Return a new array with all the elements being reversed
     */
    template<typename _M, typename _F> static type reverse(_M && m1) {
        type m2(m1.shape()); 
        int l = m1.size();
        for (int i = 0; i < l; i++) m2(i) = m1(l-1-i); 
        return std::move(m2);
    }

    /**
     * Return a conjugate array.
     */
    template<typename _M> static type conj   (_M && m1) {
        type m2(m1.shape()); 
        for (int i = 0; i < m2.size(); i++) m2(i) = std::conj(m1(i)); 
        return std::move(m2);
    }

    /**
     * Return the real part of the array.
     */
    template<typename _M> static type real   (_M && m1) {
        type m2(m1.shape()); 
        for (int i = 0; i < m2.size(); i++) m2[i] = std::real(m1(i)); 
        return std::move(m2);
    }

    /**
     * Return the imaginary part of the array.
     */
    template<typename _M> static type imag   (_M && m1) {
        type m2(m1.shape()); 
        for (int i = 0; i < m2.size(); i++) m2(i) = std::imag(m1(i)); 
        return std::move(m2);
    }

    /**
     * Return the square of the array.
     */
    template<typename _M> static type square (_M && m1) {
        type m2(m1.shape()); 
        for (int i = 0; i < m2.size(); i++) m2(i) = m1(i)*m1(i); 
        return std::move(m2);
    }

    /**
     * Return the exponation of the array.
     */
    template<typename _M> static type exp    (_M && m1) {
        type m2(m1.shape()); 
        for (int i = 0; i < m2.size(); i++) m2(i) = std::exp(m1(i)); 
        return std::move(m2);
    }

    /**
     * Return the logarithm of the array.
     */
    template<typename _M> static type log    (_M && m1) {
        type m2(m1.shape()); 
        for (int i = 0; i < m2.size(); i++) m2(i) = std::log(m1(i)); 
        return std::move(m2);
    }

    /**
     * Return the squared root of the array.
     */
    template<typename _M> static type sqrt   (_M && m1) {
        type m2(m1.shape()); 
        for (int i = 0; i < m2.size(); i++) m2(i) = std::sqrt(m1(i)); 
        return std::move(m2);
    }

    /**
     * Return the not of the array.
     */
    template<typename _M> static type logicalNot(_M && m1) {
        type m2(m1.shape()); 
        for (int i = 0; i < m2.size(); i++) m2(i) = !m1(i); 
        return std::move(m2);
    }

    /**
     * Return the absolution of the array.
     */
    template<typename _M>
    static type abs(const _M & m1) {
        type m2(m1.shape());
        for (int i = 0; i < m2.size(); i++) m2(i) = ::std::abs(m1(i));
        return std::move(m2);
    }

    /**
     * Return the difference of the array.
     *
     * @code
     * Arrayd a({4}, {1,2,3,4});
     * auto && b = Arrayd::diff(a);
     * @endcode
     * b is:
     * 1 2 3
     */
    template<typename _M> static type diff(const _M & m1) {
        type m2({m1.size()-1});
        for (int i = 0; i < m2.size(); i++) m2(i) = m1(i+1)-m1(i);
        return std::move(m2);
    }

    /**
     * Return the diagonized array.
     */
    template<typename _M> static type diag(const _M & m1) {
        if (m1.dim() < 1) ERROR_REPORT("It's an empty array!");
        if (m1.dim() > 1) {
            int l = *std::min_element(m1.shape().begin(), m1.shape().end());
            type m2({l});
            int k = 1; for (int i = 1; i < l; i++) k *= m1.shape(i);
            for (int i = 0, j = 0; i < l; i++, j += k) {
                m2[i] = m1[k];
            }
            return std::move(m2);
        }
        else {
            int n = m1.size();
            type m2({n, n}, 0);
            for (int i = 0; i < n; i++) m2(i, i) = m1[i];
            return std::move(m2);
        }
    }

    /**
     * Extract array satisfying the condition.
     *
     * @code
     * Arrayb a({4}, {true,false,false,true});
     * Arrayd b({4}, {1,2,3,4});
     * auto && c = Arrayd::extract(a, b);
     * @endcode
     * c is:
     * 1 4
     */
    template<typename _M>
    static type extract(const Array<bool> &condition, _M && m) {
        int n = 0;
        for (int i = 0; i < condition.size(); i++) if (condition[i]) n++;
        type r({n});
        for (int i = 0, j = 0; i < condition.size(); i++) if (condition[i]) {
            r[j] = m[i];
            j++;
        }
        return std::move(r);
    }

    /**
     * Create an array satisfying the range.
     *
     * @code
     * auto && a = Arrayd::range(0, 6);
     * auto && b = Arrayd::range(0, 6, 2);
     * auto && c = Arrayd::range(0, 6, 0, 2);
     * @endcode
     * a: 0 1 2 3 4 5
     * b: 0 2 4
     * c: 0 5
     */
    static type range(double a, double b, double step = 1, int n = 0) {
        if (step != 0 && n == 0) n = int(std::ceil((b-a)/step));
        else if (step == 0 && n != 0) step = (b-a)/(n-1);
        else throw "Array::range error!";

        type m({n});
        for (int i = 0; i < n; i++) m[i] = value_type(a + i * step);
        return std::move(m);
    }

    /**
     * Create an array satisfying the range.
     *
     * range({a,b,step,n}) is equal to range(a, b, step, n).
     */
    static type range(const std::vector<double> &t) {
        return range(t[0], t[1], t[2], int(t[3]));
    }

    /**
     * Create an array satisfying the range.
     *
     * range(n) is equal to range(0, n).
     */
    static type range(int a) {
        type m({a});
        for (int i = 0; i < a; i++) m[i] = value_type(i);
        return std::move(m);
    }

    /**
     * Create an array satisfying the space.
     *
     * linspace(start, stop, num) is equal to range(start, stop, 0, num).
     */
    static type linspace(double start, double stop, int num) {
        type m({num});
        double d = (stop - start) / (num - 1);
        for (int i = 0; i < num; i++) m[i] = start + i * d;
        return std::move(m);
    }

    /**
     * Create indices satisfying the condition.
     */
    static type where(const Array<bool> &condition) {
        std::vector<Shape> ls;
        int i = 0;
        ARRAY_EACH(condition.shape(), ind) {
            if (condition[i]) ls.push_back(ind);
            i++;
        }
        type m({condition.dim(), int(ls.size())});
        i = 0;
        ARRAY_EACH(m.shape(), ind) {
            m[i] = ls[ind[1]][ind[0]];
            i++;
        }
        return std::move(m);
    }

    /**
     * Return the convolve of two arrays.
     */
    template<typename _M, typename _N>
    static type convolve(const _M &m, const _N &n) {
        Array<std::complex<double>> cm = m, cn = n;
        cm.selfFFT();
        return type::real(cn.selfFFT().selfConj().selfTimes(cm).selfIFFT());
    }

protected:
    void alloc_data(int n) {
        data = (value_type*)aMalloc(n*sizeof(value_type), 64);
    }
};

/**
 * Map a pointer to an array
 */
template<typename _NumType>
class MapArray : public BasicArray<_NumType, MapArray<_NumType>> {
public:
    using value_type = _NumType;
    using type = MapArray<value_type>;
    using self_type = MapArray<value_type>;
    using base_type = BasicArray<value_type, type>;

    value_type *data;

    /**
     * Creates a map array.
     *
     * @code
     * double a[4] = {1,2,3,4};
     * MapArray b(a, {2, 2});
     * @endcode
     * The array b is :
     * 1 2
     * 3 4
     */
    MapArray(value_type * _data, const Shape & _shape) {
        data = _data;
        base_type::set_shape(_shape);
    }

//    /**
//     * Creates a map array.
//     *
//     * @code
//     * double a[4] = {1,2,3,4};
//     * MapArray b(a, 2, 2);
//     * @endcode
//     * The array b is :
//     * 1 2
//     * 3 4
//     */
//    template<typename... _Int, typename = typename std::enable_if<std::is_same<typename JnFirstType<_Int...>::type, int>::value>::type>
//    MapArray(value_type * _data, _Int ..._inds) {
//        data = _data;
//        base_type::set_shape(_inds...);
//    }

    virtual value_type &at(int a) { return data[a]; }
    virtual const value_type &at(int a) const { return data[a]; }
};

/**
 * Sub array
 */
template<typename _NumType>
class SubArray : public BasicArray<_NumType, SubArray<_NumType>> {
public:
    using value_type = _NumType;
    using type = SubArray<value_type>;
    using self_type = SubArray<value_type>;
    using base_type = BasicArray<value_type, type>;

    /**
     * Default SubArray constructor.
     */
    SubArray() {
        base_type::set_shape({0});
    }

    /**
     * SubArray assign operator reloading.
     *
     * @code
     * Arrayd a({3,3},{1,2,3,4,5,6,7,8,9});
     * auto && b = a.sub({{0,2},{0,2}});
     * Arrayd c({2,2},{-1,-2,-4,-5});
     * b = c;
     * @endcode
     * Now the array a is:
     * -1 -2  3
     * -4 -5  6
     *  7  8  9
     */
    template<typename _V, typename _D>
    type &operator =(const BasicArray<_V, _D> &m) {
        if (base_type::size() != m.size()) throw "error";
        for (int i = 0; i < base_type::size(); i++) at(i) = m(i);
        return *this;
    }

    template<typename _Array>
    SubArray(_Array && array, int axis, int val) {
        int nd = array.dim();
        assert(axis >= 0 && axis < array.dim());
        Shape shape;
        for (int i = 0; i < nd; i++) if (i != axis) shape.push_back(array.shape(i));
        base_type::set_shape(shape);

        m_data.resize(array.size()/array.shape(axis));
        int i = 0;
        int j = 0;
        ARRAY_EACH(array.shape(), ind) {
            if (ind[axis] == val) {
                m_data[j] = &(array(i));
                j++;
            }
            i++;
        }
    }

    template<typename _Array>
    SubArray(_Array && array, const std::initializer_list<std::initializer_list<int>> &_range) {
        std::vector<std::vector<int>> range;
        for (auto && r : _range) range.push_back(r);
        init(array, range);
    }

    template<typename _Array>
    SubArray(_Array && array, const std::vector<std::vector<int>> &range) {
        init(array, range);
    }

    template<typename _Array>
    void init(_Array && array, const std::vector<std::vector<int>> &range) {
        int dim = array.dim();
        if (range.size() != dim) DIE("BasicArray::map error at ", __LINE__);

        Shape shape(dim);
        for (int i = 0; i < dim; i++) shape[i] = range[i][1] - range[i][0];
        base_type::set_shape(shape);
//        base_type::m_shape.resize(dim);
//        for (int i = 0; i < dim; i++) base_type::m_shape[i] = range[i][1] - range[i][0];
//        base_type::set_coeff();
        //m_data = new value_type *[base_type::size()];
        m_data.resize(base_type::size());

        std::vector<int> inds(dim);
        for (int i = 0; i < dim; i++) inds[i] = range[i][0];
        int n = 0;
        do {
            m_data[n] = &array(base_type::get_ind(inds, array.shape()));
            n++;
        } while (base_type::next_ind(inds, range));
    }

    template<typename _Array>
    SubArray(_Array && array, const Array<int> &range) {
        base_type::set_shape(range.shape());
        m_data.resize(base_type::size());
        for (int i = 0; i < base_type::size(); i++) m_data[i] = &array(range(i));
    }

    virtual value_type &at(int a) { return *(m_data[a]); }
    virtual const value_type &at(int a) const { return *(m_data[a]); }

    template<typename _Array>
    static self_type take(_Array && array, int val, int axis, int mode = ARRAY_TAKE_MODE_CLIP) {
        ERROR_CHECK(axis < 0 || axis >= array.dim(), "ILLEGAL axis!!!");

        self_type sub;
        int nd = array.dim();

        // Set shape
        Shape shape;
        for (int i = 0; i < nd; i++) if (i != axis) shape.push_back(array.shape(i));
        sub.set_shape(shape);

        // Set m_data
        sub.m_data.resize(array.size()/array.shape(axis));
        int i = 0;
        int j = 0;
        ARRAY_EACH(array.shape(), ind) {
            if (ind[axis] == val) {
                sub.m_data[j] = &(array[i]);
                j++;
            }
            i++;
        }

        return std::move(sub);
    }

    template<typename _Array, typename _LS>
    static self_type takes(_Array && array, _LS && ls, int axis, int mode = ARRAY_TAKE_MODE_CLIP) {
        ERROR_CHECK(axis < 0 || axis >= array.dim(), "ILLEGAL axis!!!");

        self_type sub;
        int nd = array.dim();

        // Set shape
        Shape shape = array.shape();
        shape[axis] = ls.size();
        sub.set_shape(shape);

        // Set m_data
        sub.m_data.resize(array.size()/array.shape(axis)*ls.size());
        int i = 0;
        ARRAY_EACH(array.shape(), ind) {
            if (ind[axis] == array.shape(axis) - 1) {
                int d = 0;
                for (auto && n : ls) {
                    if (n >= array.shape(axis)-1) {
                        if (n > array.shape(axis)-1 && mode != ARRAY_TAKE_MODE_CLIP) ERROR_REPORT("Indices error!");
                        Shape v = ind;
                        v[axis] = d;
                        int sum = 0; for (int k = 0; k < nd; k++) sum += sub.coeff(k) * v[k];
                        sub.m_data[sum] = &(array(i));
                    }
                    d++;
                }
            }
            else if (ind[axis] == 0) {
                int d = 0;
                for (auto && n : ls) {
                    if (n <= 0) {
                        if (n < 0 && mode != ARRAY_TAKE_MODE_CLIP) ERROR_REPORT("Indices error!");
                        Shape v = ind;
                        v[axis] = d;
                        int sum = 0; for (int k = 0; k < nd; k++) sum += sub.coeff(k) * v[k];
                        sub.m_data[sum] = &(array(i));
                    }
                    d++;
                }
            }
            else {
                auto it = std::find(ls.begin(), ls.end(), ind[axis]);
                if (it != ls.end()) {
                    int d = std::distance(ls.begin(), it);
                    Shape v = ind;
                    v[axis] = d;
                    int sum = 0; for (int k = 0; k < nd; k++) sum += sub.coeff(k) * v[k];
                    sub.m_data[sum] = &(array(i));
                }
            }
            i++;
        }

        return std::move(sub);
    }

protected:
    std::vector<value_type *> m_data;

};


/**
 * @name Deprecated Class MatrixBase3D, Matrix3D, SubMatrix3D and MapMatrix3D
 * @{
 */

#define MAT3_EACH(mat, i, j, k) for (int i = 0; i < (mat).shape(0); i++) for (int j = 0; j < (mat).shape(1); j++) for (int k = 0; k < (mat).shape(2); k++)

template<typename _NumType> class SubMatrix3D;

template<typename _NumType>
class MatrixBase3D {
public:
    using value_type = _NumType;
    using type = MatrixBase3D<value_type>;

    virtual value_type &at(long int a, long int b, long int c) = 0;
    virtual const value_type &at(long int a, long int b, long int c) const = 0;
    virtual value_type &at(long int a) = 0;
    virtual const value_type &at(long int a) const = 0;
    virtual int shape(int n) const = 0;

    template<typename _V> type &set_all(_V && v) {
        MAT3_EACH(*this, i, j, k) at(i, j, k) = v;
        return *this;
    }

    value_type mean() const {
        value_type n = 0;
        MAT3_EACH(*this, i, j, k) n += at(i, j, k);
        return value_type(n/(shape(0)*shape(1)*shape(2)));
    }

    template<typename T = value_type>
    T sum() const {
        T n = 0;
        MAT3_EACH(*this, i, j, k) n += at(i, j, k);
        return value_type(n);
    }

    value_type min() const {
        value_type min = at(0, 0, 0);
        MAT3_EACH(*this, i, j, k) {
            if (min > at(i, j, k)) min = at(i, j, k);
        }
        return min;
    }

    value_type max() const {
        value_type max = at(0, 0, 0);
        MAT3_EACH(*this, i, j, k) {
            if (max < at(i, j, k)) max = at(i, j, k);
        }
        return max;
    }

    SubMatrix3D<value_type> sub(int a1, int a2, int b1, int b2, int c1, int c2) {
        return SubMatrix3D<value_type>(this, a1, a2, b1, b2, c1, c2);
    }

    template<typename T> void assign(T &&m) {
        MAT3_EACH(*this, i, j, k) at(i, j, k) = m(i, j, k);
    }

    void selfSquare (           ) { MAT3_EACH(*this, i, j, k) at(i, j, k) *= at(i, j, k);      }
    template<typename _N> void selfEq     (const _N &n) { MAT3_EACH(*this, i, j, k) at(i, j, k)  = at(i, j, k) == n; }
    template<typename _N> void selfLt     (const _N &n) { MAT3_EACH(*this, i, j, k) at(i, j, k)  = at(i, j, k) <  n; }
    template<typename _N> void selfLe     (const _N &n) { MAT3_EACH(*this, i, j, k) at(i, j, k)  = at(i, j, k) <= n; }
    template<typename _N> void selfGt     (const _N &n) { MAT3_EACH(*this, i, j, k) at(i, j, k)  = at(i, j, k) >  n; }
    template<typename _N> void selfGe     (const _N &n) { MAT3_EACH(*this, i, j, k) at(i, j, k)  = at(i, j, k) >= n; }
    template<typename _N> void selfPlus   (const _N &n) { MAT3_EACH(*this, i, j, k) at(i, j, k) += n;                }
    template<typename _N> void selfMinus  (const _N &n) { MAT3_EACH(*this, i, j, k) at(i, j, k) -= n;                }
    template<typename _N> void selfTimes  (const _N &n) { MAT3_EACH(*this, i, j, k) at(i, j, k) *= n;                }
    template<typename _N> void selfDivide (const _N &n) { MAT3_EACH(*this, i, j, k) at(i, j, k) /= n;                }
    template<typename _N> void selfPlusM  (const _N &n) { MAT3_EACH(*this, i, j, k) at(i, j, k) += n(i, j, k);       }
    template<typename _N> void selfMinusM (const _N &n) { MAT3_EACH(*this, i, j, k) at(i, j, k) -= n(i, j, k);       }
    template<typename _N> void selfTimesM (const _N &n) { MAT3_EACH(*this, i, j, k) at(i, j, k) *= n(i, j, k);       }
    template<typename _N> void selfDivideM(const _N &n) { MAT3_EACH(*this, i, j, k) at(i, j, k) /= n(i, j, k);       }

    value_type &operator()(int a, int b, int c) { return at(a, b, c); }
    const value_type &operator()(int a, int b, int c) const { return at(a, b, c); }
    value_type &operator()(int a) { return at(a); }
    const value_type &operator()(int a) const { return at(a); }

    friend std::ostream &operator <<(std::ostream &stream, const type &m) {
        MAT3_EACH(m, i, j, k) {
            stream << m(i, j, k);
            if (k == m.shape(2) - 1) stream << std::endl;
            else stream << ' ';
        }
        return stream;
    }

};

template<typename _NumType>
class Matrix3D : public MatrixBase3D<_NumType> {
public:
    using value_type = _NumType;
    using base_type = MatrixBase3D<value_type>;
    using type = Matrix3D<value_type>;

    std::array<int, 3> size;
    value_type *data;

    Matrix3D()                                         { init();                      }
    Matrix3D(int a, int b, int c)                      { init();    resize(a, b, c);  }
    template<typename T> Matrix3D(int a, int b, int c, const T & t)         { init();    resize(a, b, c, t);  }
    template<typename _M> Matrix3D(const _M &m)        { init(); resize(m.shape(0), m.shape(1), m.shape(2));  assign(m);        }
    Matrix3D(type &&m)                                 { init();    swap(m);          }
    template<typename _M> type &operator=(const _M &m) { assign(m); return *this;     }
    type &operator=(type &&m)                          { swap(m);   return *this;     }
    ~Matrix3D()                                        { clear();                     }

    virtual       value_type &at(long int a, long int b, long int c)       { return data[a*size[1]*size[2]+b*size[2]+c]; }
    virtual const value_type &at(long int a, long int b, long int c) const { return data[a*size[1]*size[2]+b*size[2]+c]; }
    virtual value_type &at(long int a) { return data[a]; }
    virtual const value_type &at(long int a) const {return data[a]; }
    virtual int shape(int n) const { return size[n]; }

    void resize(int a, int b, int c) {
        clear();
        size[0] = a;
        size[1] = b;
        size[2] = c;
        data = new value_type[a*b*c];
    }

    template<typename T>
        void resize(int a, int b, int c, const T & t) {
            clear();
            size[0] = a;
            size[1] = b;
            size[2] = c;
            data = new value_type[a*b*c];
            set_all(t);
        }

    static type Zero(int a, int b, int c) { type m(a, b, c); m.set_all(0); return std::move(m); }
    static type Ones(int a, int b, int c) { type m(a, b, c); m.set_all(1); return std::move(m); }
    template<typename T> static type Constant(int a, int b, int c, T && v) {
        type m(a, b, c); m.set_all(v); return std::move(m);
    }

    template<typename _M             > static type square (_M && m1              ) {
        type m2(std::forward<_M>(m1)); m2.selfSquare ( ); return std::move(m2);
    }
    template<typename _M, typename _N> static type eq     (_M && m1, const _N & n) {
        type m2(std::forward<_M>(m1)); m2.selfEq   (n); return std::move(m2);
    }
    template<typename _M, typename _N> static type lt     (_M && m1, const _N & n) {
        type m2(std::forward<_M>(m1)); m2.selfLt   (n); return std::move(m2);
    }
    template<typename _M, typename _N> static type le     (_M && m1, const _N & n) {
        type m2(std::forward<_M>(m1)); m2.selfLe   (n); return std::move(m2);
    }
    template<typename _M, typename _N> static type gt     (_M && m1, const _N & n) {
        type m2(std::forward<_M>(m1)); m2.selfGt   (n); return std::move(m2);
    }
    template<typename _M, typename _N> static type ge     (_M && m1, const _N & n) {
        type m2(std::forward<_M>(m1)); m2.selfGe   (n); return std::move(m2);
    }
    template<typename _M, typename _N> static type plus   (_M && m1, const _N & n) {
        type m2(std::forward<_M>(m1)); m2.selfPlus   (n); return std::move(m2);
    }
    template<typename _M, typename _N> static type minus  (_M && m1, const _N & n) {
        type m2(std::forward<_M>(m1)); m2.selfMinus  (n); return std::move(m2);
    }
    template<typename _M, typename _N> static type times  (_M && m1, const _N & n) {
        type m2(std::forward<_M>(m1)); m2.selfTimes  (n); return std::move(m2);
    }
    template<typename _M, typename _N> static type divide (_M && m1, const _N & n) {
        type m2(std::forward<_M>(m1)); m2.selfDivide (n); return std::move(m2);
    }
    template<typename _M, typename _N> static type plusM  (_M && m1, const _N & n) {
        type m2(std::forward<_M>(m1)); m2.selfPlusM  (n); return std::move(m2);
    }
    template<typename _M, typename _N> static type minusM (_M && m1, const _N & n) {
        type m2(std::forward<_M>(m1)); m2.selfMinusM (n); return std::move(m2);
    }
    template<typename _M, typename _N> static type timesM (_M && m1, const _N & n) {
        type m2(std::forward<_M>(m1)); m2.selfTimesM (n); return std::move(m2);
    }
    template<typename _M, typename _N> static type divideM(_M && m1, const _N & n) {
        type m2(std::forward<_M>(m1)); m2.selfDivideM(n); return std::move(m2);
    }

    template<typename _M> static type abs(const _M & m1) {
        type m2(m1.shape(0), m1.shape(1), m1.shape(2));
        MAT3_EACH(m2, i, j, k) m2(i, j, k) = std::abs(m1(i, j, k));
        return std::move(m2);
    }

protected:
    void init() {
        size = { 0, 0, 0 };
        data = NULL;
    }

    template<typename T> void swap(T && m) {
        std::swap(size, m.size);
        std::swap(data, m.data);
    }

    type &clear() {
        if (data != NULL) delete[] data;
        init();
    }

};

template<typename _NumType>
class MapMatrix3D : public MatrixBase3D<_NumType> {
public:
    using value_type = _NumType;
    using type = MapMatrix3D<value_type>;
    using base_type = MatrixBase3D<value_type>;

    value_type *data;
    std::array<int, 3> size;

    MapMatrix3D() = default;
    MapMatrix3D(const type & m) = default;
    // MapMatrix3D(type && m) = default;
    type &operator =(const type & m) = default;
    // type &operator =(type && m) = default;

    MapMatrix3D(value_type *p, int a, int b, int c) {
        data = p;
        size[0] = a;
        size[1] = b;
        size[2] = c;
    }

    virtual value_type &at(long int a, long int b, long int c) { return (data[a*size[1]*size[2]+b*size[2]+c]); }
    virtual const value_type &at(long int a, long int b, long int c) const { return (data[a*size[1]*size[2]+b*size[2]+c]); }
    virtual value_type &at(long int a) { return data[a]; }
    virtual const value_type &at(long int a) const {return data[a]; }
    virtual int shape(int n) const { return size[n]; }

};

template<typename _NumType>
class SubMatrix3D : public MatrixBase3D<_NumType> {
public:
    using value_type = _NumType;
    using type = SubMatrix3D<value_type>;
    using orig_type = Matrix3D<value_type>;
    using base_type = MatrixBase3D<value_type>;

    SubMatrix3D() = default;
    SubMatrix3D(const type & m) = default;
    // SubMatrix3D(type && m) = default;
    type &operator =(const type & m) = default;
    // type &operator =(type && m) = default;

    SubMatrix3D(base_type *orig, int a1, int a2, int b1, int b2, int c1, int c2) {
        mat = orig;
        for (int i = a1; i < a2; i++) inds[0].push_back(i);
        for (int i = b1; i < b2; i++) inds[1].push_back(i);
        for (int i = c1; i < c2; i++) inds[2].push_back(i);
    }

    virtual value_type &at(long int a, long int b, long int c) { MPI_LOG << inds[0][2] << std::endl; return mat->at(inds[0][a], inds[1][b], inds[2][c]); }
    virtual const value_type &at(long int a, long int b, long int c) const { return mat->at(inds[0][a], inds[1][b], inds[2][c]); }
    virtual value_type &at(long int a) { return at(a/(shape(1)*shape(2)), a/shape(2), a%shape(2)); }
    virtual const value_type &at(long int a) const { return at(a/(shape(1)*shape(2)), a/shape(2), a%shape(2)); }
    virtual int shape(int n) const { return inds[n].size(); }

protected:
    base_type *mat;
    std::array<std::vector<int>, 3> inds;

};

#undef MAT3_EACH

/**
 * @}
 */

/**
 * @name Class Redefinitions
 * @{
 */
using Arrayb  = Array<bool>;
using Arrayc  = Array<char>;
using Arrays  = Array<std::string>;
using Arrayi  = Array<int>;
using Arrayf  = Array<float>;
using Arrayd  = Array<double>;
using Arrayci = Array<std::complex<int>>;
using Arraycf = Array<std::complex<float>>;
using Arraycd = Array<std::complex<double>>;

using MapArrayb  = MapArray<bool>;
using MapArrayc  = MapArray<char>;
using MapArrays  = MapArray<std::string>;
using MapArrayi  = MapArray<int>;
using MapArrayf  = MapArray<float>;
using MapArrayd  = MapArray<double>;
using MapArrayci = MapArray<std::complex<int>>;
using MapArraycf = MapArray<std::complex<float>>;
using MapArraycd = MapArray<std::complex<double>>;

using SubArrayb  = SubArray<bool>;
using SubArrayc  = SubArray<char>;
using SubArrays  = SubArray<std::string>;
using SubArrayi  = SubArray<int>;
using SubArrayf  = SubArray<float>;
using SubArrayd  = SubArray<double>;
using SubArrayci = SubArray<std::complex<int>>;
using SubArraycf = SubArray<std::complex<float>>;
using SubArraycd = SubArray<std::complex<double>>;

using Mat3b     =    Matrix3D<bool                 >;
using Mat3i     =    Matrix3D<int                  >;
using Mat3f     =    Matrix3D<float                >;
using Mat3d     =    Matrix3D<double               >;
using Mat3ci    =    Matrix3D<std::complex<int>    >;
using Mat3cf    =    Matrix3D<std::complex<float>  >;
using Mat3cd    =    Matrix3D<std::complex<double> >;

using MapMat3b  = MapMatrix3D<bool                 >;
using MapMat3i  = MapMatrix3D<int                  >;
using MapMat3f  = MapMatrix3D<float                >;
using MapMat3d  = MapMatrix3D<double               >;
using MapMat3ci = MapMatrix3D<std::complex<int>    >;
using MapMat3cf = MapMatrix3D<std::complex<float>  >;
using MapMat3cd = MapMatrix3D<std::complex<double> >;

using SubMat3b  = SubMatrix3D<bool                 >;
using SubMat3i  = SubMatrix3D<int                  >;
using SubMat3f  = SubMatrix3D<float                >;
using SubMat3d  = SubMatrix3D<double               >;
using SubMat3ci = SubMatrix3D<std::complex<int>    >;
using SubMat3cf = SubMatrix3D<std::complex<float>  >;
using SubMat3cd = SubMatrix3D<std::complex<double> >;

/**
 * @}
 */

