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

#include "resmap_macros.h"

#include "res_array.h"
#include "res_mkl.h"

namespace resumap {

template<typename T>
struct lstsq;

template<>
struct lstsq<std::complex<double>> {
    Arraycd c;
    Arrayd s;
    int rank;

    lstsq(Array<std::complex<double>> &x, Array<std::complex<double>> y, double rcond) {
        assert(x.shape(0) >= x.shape(1));

        Arrayd s({x.shape(0)});
        int rank;

        if (LAPACKE_zgelsd(LAPACK_ROW_MAJOR, x.shape(0), x.shape(1), y.shape(1), 
                    (MKL_Complex16*)x.data_wptr(), x.shape(1), (MKL_Complex16*)y.data_wptr(), y.shape(1), s.data_wptr(), rcond, &rank) > 0) {
            std::cerr << "The algorithm computing SVD failed to converge;\n" << std::endl;
            std::cerr << "the least squares solution could not be computed.\n" << std::endl;
            EXIT_ABNORMALLY;
        }

        c = y.sub({{0, x.shape(1)}, {0, y.shape(1)}});
        s = s;
        rank = rank;
    }
};

template<>
struct lstsq<double> {
    Arrayd c;
    Arrayd s;
    int rank;

    lstsq(Array<double> &x, Array<double> y, double rcond) {
//        MPI_LOG << "x: " << x << std::endl;
//        MPI_LOG << "y: " << y << std::endl;

        ERROR_CHECK(x.shape(0) < x.shape(1), "The number of rows of x must be larger than cols.");

        Arrayd s({x.shape(0)});
        int rank;

        if (LAPACKE_dgelsd(LAPACK_ROW_MAJOR, x.shape(0), x.shape(1), y.shape(1), 
                    x.data_wptr(), x.shape(1), y.data_wptr(), y.shape(1), s.data_wptr(), rcond, &rank) > 0) {
            std::cerr << "The algorithm computing SVD failed to converge;\n" << std::endl;
            std::cerr << "the least squares solution could not be computed.\n" << std::endl;
            EXIT_ABNORMALLY;
        }

//        MPI_LOG << "x: " << x << std::endl;
//        MPI_LOG << "y: " << y << std::endl;

        c = y.sub({{0, x.shape(1)}, {0, y.shape(1)}});
        s = s;
        rank = rank;
    }
};

template<typename _V, typename _X, typename _C>
Array<_V> polyval(const _X &x, const _C &c) {
    Array<_V> y = x;
    int c_size = c.size();
    for (int i = 0; i < y.size(); i++) {
        _V sum = 0;
#if 0 // before optimize
        for (int j = 0; j < c.size(); j++) {
            sum += std::pow(x[i], j) * c[j];
        }
#else
        _V pow_x_i = 1;
        for (int j = 0; j < c_size; j++) {
            sum += pow_x_i * c[j];
            pow_x_i *= x[i];
        }
#endif
        y[i] = sum;
    }
    return std::move(y);
}

template<typename _V>
Array<_V> rollaxis(const Array<_V> &x, int axis, int start) {
    int nd = x.dim();
    int size = x.size();
    int m = x.shape(axis);
    int n = m - start;

    Shape shape(nd);
    Shape coeff(nd);
    for (int i = nd-1, j = nd-1; i >= 0; i--) {
        if (i != axis) {
            shape[j] = x.shape(i);
            coeff[j] = x.coeff(i);
            j--;
        }
    }
    shape[0] = x.shape(axis);
    coeff[0] = x.coeff(axis);

    Array<_V> output(shape);
    int l = 0;
    ARRAY_EACH(shape, ind) {
        int nn = (ind[0]+start) * coeff[0];
        for (int i = 1; i < nd; i++) nn += ind[i] * coeff[i];
        output[l] = x[nn];
        l++;
    }
    return std::move(output);
}

template<typename _V>
Array<_V> polyvander(const Array<_V> &x, int deg) {
    int nx = x.size();
    Arrayd o({nx, deg+1});
    for (int i = 0; i < nx; i++) {
        o(i, 0) = 1;
        for (int j = 1; j < deg+1; j++) {
            o(i, j) = o(i, j-1) * x[i];
        }
    }
    return std::move(o);

// The below codes are translated from numpy, and they are deprecated now.
//    assert(deg >= 0);
//    Shape shape = x.shape();
//    shape.push_front(deg+1);
//    Array<_V> v(shape);
//    v.sub(0) = 1;
//    if (deg > 0) {
//        v.sub(1) = x;
//        for (int i = 2; i < deg+1; i++) {
//            v.sub(i) = v.sub(i-1) * x;
//        }
//    }
//    return rollaxis(v, 0, v.dim());
}

struct polyfit {
    Arrayd c;
    Arrayd s;
    int rank;

    polyfit(const Arrayd &x, Arrayd y, Arrayi deg, const Arrayd &w = Arrayd{}, double rcond=-1) {
        init(x,y,deg,w,rcond);
    }

    polyfit(const Arrayd &x, Arrayd y, int deg, const Arrayd &w = Arrayd{}, double rcond=-1) {
        init(x, y, Arrayi::range(deg+1), w, rcond);
    }

    void init(const Arrayd &x, Arrayd y, Arrayi deg, const Arrayd &w = Arrayd{}, double rcond=-1) {
        // check arguments.
        //    assert((deg>=0).all());
        //    assert(x.dim()==1);
        //    assert(x.size() != 0);
//        MPI_LOG << "x: " << x << std::endl;
//        MPI_LOG << "y: " << y << std::endl;
//        MPI_LOG << "deg: " << deg << std::endl;
        ERROR_CHECK(y.dim() != 2, "y must be a 2-dimensional matrix!");
        //    assert(x.shape(0) == y.shape(0));

        std::sort(deg.begin(), deg.end());
        int nx = x.size();
        int ny = y.size()/nx;
        int order = deg.size();
        int lmax = deg[order-1];
        Arrayd van_full = polyvander(x, lmax);
        Arrayd van({nx, order});
        for (int i = 0; i < nx; i++) for (int j = 0; j < order; j++) van(i, j) = van_full(i, deg[j]);

//        MPI_LOG << "nx: " << nx << std::endl;
//        MPI_LOG << "ny: " << ny << std::endl;
//        MPI_LOG << "order: " << order << std::endl;
//        MPI_LOG << "lmax: " << lmax << std::endl;
//        MPI_LOG << "van: " << van << std::endl;

        // set up the least squares matrices in transposed form
        if (w.size() != 0) {
            for (int i = 0; i < nx; i++) for (int j = 0; j < order; j++) van(i, j) *= w[i];
            for (int i = 0; i < nx; i++) for (int j = 0; j < ny; j++) y(i, j) *= w[i];
        }

//        MPI_LOG << "van: " << van << std::endl;
//        MPI_LOG << "y: " << y << std::endl;

        // set rcond
        if (rcond == -1) rcond = x.size() * DBL_EPSILON;

//        MPI_LOG << "rcond: " << rcond << std::endl;

        // Determine the norms of the design matrix columns.
        Arrayd scl({order}, 0);
        for (int j = 0; j < order; j++) {
            for (int i = 0; i < nx; i++) {
                scl[j] += square(van(i, j));
            }
            if (scl[j] == 0) scl[j] = 1;
            else scl[j] = std::sqrt(scl[j]);
            for (int i = 0; i < nx; i++) {
                van(i, j) /= scl[j];
            }
        }

//        MPI_LOG << "scl: " << scl << std::endl;
//        MPI_LOG << "van: " << van << std::endl;

        // Solve the least squares problem.
        lstsq<double> rt(van, y, rcond);
//        MPI_LOG << "rt.c: " << rt.c << std::endl;

        for (int i = 0; i < order; i++) for (int j = 0; j < ny; j++) {
            rt.c(i, j) /= scl[i];
        }

//        MPI_LOG << "rt.c: " << rt.c << std::endl;

        // Expand c to include non-fitted coefficients which are set to zero
        Arrayd cc({lmax+1, ny}, 0);
        for (int i = 0; i < order; i++) for (int j = 0; j < ny; j++) {
            cc(deg[i], j) = rt.c(i, j);
        }
        rt.c = cc;

//        MPI_LOG << "final rt.c: " << rt.c << std::endl;

        c = rt.c;
        s = rt.s;
        rank = rt.rank;
    }

};

}

