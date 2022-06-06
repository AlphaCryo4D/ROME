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

#include "res_array.h"
#include "res_mkl.h"

namespace linalg {

/**
 * Eigen values decomposition.
 */
struct eig {
    int info;
    Arrayd a, wr, wi, vl, vr;

    eig(const Arrayd &_a) 
        : a(_a)
    {
        int n = a.shape(0);
        int lda = n;
        wr = Arrayd({n});
        wi = Arrayd({n});
        vl = Arrayd({n,n});
        vr = Arrayd({n,n});
        int ldvl = n;
        int ldvr = n;
        info = LAPACKE_dgeev(LAPACK_ROW_MAJOR, 'N', 'V', n, a.data_wptr(), lda, wr.data_wptr(), wi.data_wptr(), vl.data_wptr(), ldvl, vr.data_wptr(), ldvr );
    }
};

/**
 * Singular values decomposition.
 */
struct svd {
    int m;
    int n;
    int lda;
    int ldu;
    int ldvt;
    Arrayd s;
    Arrayd u;
    Arrayd vt;
    Arrayd superb;
    int info;

    svd(Arrayd a) :
        s({std::min(a.shape(0), a.shape(1))}), superb({std::min(a.shape(0), a.shape(1))-1}),
        u({a.shape(0), a.shape(0)}), vt({a.shape(1), a.shape(1)}),
        lda(a.shape(1)), ldu(a.shape(0)), ldvt(a.shape(1)), 
        m(a.shape(0)), n(a.shape(1))
    {
        info = LAPACKE_dgesvd( LAPACK_ROW_MAJOR, 'A', 'A', m, n, a.data_wptr(), lda,
                s.data_wptr(), u.data_wptr(), ldu, vt.data_wptr(), ldvt, superb.data_wptr());
    }

    bool converge() {
        return info <= 0;
    }
};

/**
 * Dot multiplication of matrix.
 */
inline Arrayd dot(const Arrayd &a, const Arrayd &b) {
    int m = a.shape(0);
    int k = a.shape(1);
    int n = b.shape(1);
    double alpha = 1;
    double beta = 0;
    ERROR_CHECK(k != b.shape(0), 
            "The number of columns of a must be equal to the number of rows of b.");
    Arrayd c({m, n});
    cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, 
                            m, n, k, alpha, a.data_rptr(), k, b.data_rptr(), n, beta, c.data_wptr(), n);
    return std::move(c);
}

template<typename _V, typename _A, typename _B>
_V vdot(const _A &a, const _B &b) {
    ERROR_CHECK(a.size() != b.size(), "The size of a must be equal to that of b.");
    _V n = 0;
	int a_size = a.size();
	//#pragma omp simd TODO
    for (int i = 0; i < a_size; i++)
		n += a(i) * b(i);
    return n;
}

inline Arrayd pinv(const Arrayd &a, double rcond=1e-15 ) {
//    a, wrap = _makearray(a)
//    if _isEmpty2d(a):
//        res = empty(a.shape[:-2] + (a.shape[-1], a.shape[-2]), dtype=a.dtype)
//        return wrap(res)
    int m = a.shape(0);
    int n = a.shape(1);

    //a.selfConj();
    svd rt(a);
    double cutoff = rcond * rt.s.max();
    for (auto && n : rt.s) {
        if (n < cutoff) n = 0;
        else n = 1.0 / n;
    }
    Arrayd s({n, m}, 0);
    for (int i = 0; i < std::min(m, n); i++) {
        s(i, i) = rt.s[i];
    }
    //return dot(rt.vt.transpose(), dot(s, rt.u.transpose()));
    return dot(rt.vt.transpose(), dot(s, rt.u.transpose()));
}

} // namespace linalg


