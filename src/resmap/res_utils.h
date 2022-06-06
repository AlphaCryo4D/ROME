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
#include "res_stats.h"

namespace resumap {

inline Arrayd mgrid(const std::vector<std::vector<double>> &ls) {
    int l = int(ls.size());
    std::vector<Arrayd> u(l);
//    std::vector<Arrayd> v(l);
    Shape shape(l);
    for (int i = 0; i < l; i++) {
        u[i] = Arrayd::range(ls[i]);
        shape[i] = u[i].size();
    }
    Shape oshape(shape.size()+1);
    oshape[0] = l;
    for (int i = 0; i < l; i++) oshape[i+1] = shape[i];
    Arrayd v(oshape);
    int size = std::accumulate(shape.begin(), shape.end(), 1, [](int a,int b)->int{return a*b;});
    for (int i = 0; i < l; i++) {
//        v[i].realloc(shape);
        int j = 0;
        ARRAY_EACH(shape, ind) {
            v[i*size+j] = u[i][ind[i]];
//            v[i][j] = u[i][ind[i]];
            j++;
        }
    }
    return std::move(v);
}

template<typename T>
Arrayi unravel_index_helper(T && values, const Shape &shape) {
    int m = int(shape.size());
    int n = int(values.size());
    Arrayi ind({m, n});
    int x = 1;
    int y = 1;
    for (int i = m-1; i >= 0; i--) {
        x *= shape[i];
        int j = 0;
        for (auto && v : values) {
            ind(i, j) = (v % x) / y;
//            std::cout << i << ' ' << j << ' ' << x << ' ' << y << ' ' << ind(i, j) << std::endl;
            j++;
        }
        y *= shape[i];
    }
    return std::move(ind);
}

inline Arrayi unravel_index(int n, const Shape &shape) {
    auto && ind = unravel_index_helper(std::vector<int>{n}, shape);
    ind.reshape({int(shape.size())});
    return std::move(ind);
}

template<typename T>
Arrayi unravel_index(const std::initializer_list<T> && values, const Shape &shape) {
    return unravel_index_helper(values, shape);
}

template<typename T>
Arrayi unravel_index(T && values, const Shape &shape) {
    return unravel_index_helper(values, shape);
}

inline Arrayd hamming(int n) {
    Arrayd v({n});
    for (int i = 0; i < n; i++) {
        v[i] = 0.54 - 0.46 * std::cos(2 * PI * i / (n - 1));
    }
    return std::move(v);
}

template<typename _V1, typename _V2>
Arrayd array_outer_product(_V1 && v1, _V2 && v2) {
    int l = v1.shape(0);
    Arrayd v({l, l});
    for (int i = 0; i < l; i++) {
        for (int j = 0; j < l; j++) v(i, j) = v1[i]*v2[j];
    }
    return std::move(v);
}

template<typename _V1, typename _V2>
Arrayd array_outer_product2(_V1 && v1, _V2 && v2) {
    int a = v1.shape(0), b = v1.shape(1), c = v2.shape(0);
    Arrayd m({a, b, c});
    for (int i = 0; i < a; i++) for (int j = 0; j < b; j++) for (int k = 0; k < c; k++) {
        m(i, j, k) = v1(i, j) * v2(i, k);
    }
    return std::move(m);
}

/**
 *
 * Takes (x, y, z) numpy.mgrid inputs and generates 16 unit normals that point to the edges and 
 * faces of the icosahedron (Alp Kucukelbir, 2013)
 * 
 * See the following citation for an explanation of why this is needed for 3D steerable filters
 * 
 * Konstantinos G Derpanis and Jacob M Gryn. Three-dimensional nth derivative of
 * gaussian separable steerable filters. In Image Processing, 2005. ICIP 2005. IEEE In-
 * ternational Conference on, volume 3. IEEE, 2005.
 * 
 * Parameters
 * ----------
 * (x,y,z): outputs of numpy.mgrid
 * 
 * Returns
 * -------
 * The unit normal matrices oriented towards the edges and faces of the icosahedron
 * 
 */

inline Arrayd make3DsteerableDirections(const Arrayd & x, const Arrayd & y, const Arrayd & z) {
    Shape shape{x.shape(0), x.shape(1), x.shape(2), 16};
    Arrayd dirs(shape, 0);

    Arrayd a({3}, {1, 0, (std::sqrt(5.0)+1)/2.0});
    a.selfDivide(a.norm());
    Arrayd b({3}, {a[0], a[1], -a[2]});
    Arrayd c({3}, {1, (std::sqrt(5.0)+1)/2.0, 2.0/(std::sqrt(5.0)+1)});
    c.selfDivide(c.norm());
    Arrayd d({3}, {c[0], -c[1], c[2]});

//    a.printFull();
//    b.printFull();
//    c.printFull();
//    d.printFull();

    for (int i = 0, j = 0; j < x.size(); j++) {
        // 6 rotations for G2
        // Unit normals to the faces of the dodecahedron
        dirs[i] = x[j] * a[0] + y[j] * a[1] + z[j] * a[2]; i++;
        dirs[i] = x[j] * a[1] + y[j] * a[2] + z[j] * a[0]; i++;
        dirs[i] = x[j] * a[2] + y[j] * a[0] + z[j] * a[1]; i++;
        // Flip sign of golden ratio (arbitrary choice, just stay consistent)
        dirs[i] = x[j] * b[0] + y[j] * b[1] + z[j] * b[2]; i++;
        dirs[i] = x[j] * b[1] + y[j] * b[2] + z[j] * b[0]; i++;
        dirs[i] = x[j] * b[2] + y[j] * b[0] + z[j] * b[1]; i++;
        // 10 rotations for H2
        // Unit normals to the faces of the icosahedron
        dirs[i] = x[j] * c[0] + y[j] * c[1] + z[j] * c[2]; i++;
        dirs[i] = x[j] * c[1] + y[j] * c[2] + z[j] * c[0]; i++;
        dirs[i] = x[j] * c[2] + y[j] * c[0] + z[j] * c[1]; i++;
        // Flip sign of golden ratio (arbitrary choice, just stay consistent)
        dirs[i] = x[j] * d[0] + y[j] * d[1] + z[j] * d[2]; i++;
        dirs[i] = x[j] * d[1] + y[j] * d[2] + z[j] * d[0]; i++;
        dirs[i] = x[j] * d[2] + y[j] * d[0] + z[j] * d[1]; i++;
        // Unit normals to the vertices of the cube
        dirs[i] = 1. / std::sqrt(3.0) * ( x[j]+y[j]+z[j]); i++;
        dirs[i] = 1. / std::sqrt(3.0) * (-x[j]+y[j]+z[j]); i++;
        dirs[i] = 1. / std::sqrt(3.0) * ( x[j]-y[j]+z[j]); i++;
        dirs[i] = 1. / std::sqrt(3.0) * (-x[j]-y[j]+z[j]); i++;
    }

    return std::move(dirs);
}

//template<typename _Array>
struct ArrayWindowRoller {
    Shape shape;
    Shape winSize;
    Shape step;
    int dim;
    std::vector<std::vector<int>> range;

    ArrayWindowRoller(const Shape & _shape, const Shape & _winSize, const Shape & _step = {1, 1, 1}) {
        shape = _shape;
        winSize = _winSize;
        step = _step;
        dim = int(shape.size());
        range.resize(dim);

		for (int i = 0; i < dim; i++) {
			std::vector<int> tmp = { 0, winSize[i] };
			range[i] = tmp;
		}
    }

    bool roll() {
        range[dim-1][0] += step[dim-1];
        range[dim-1][1] += step[dim-1];
        for (int i = dim-1; i > 0; i--) {
            if (range[i][1] > shape[i]) {
                range[i][0] = 0;
                range[i][1] = winSize[i];
                range[i-1][0] += step[i-1];
                range[i-1][1] += step[i-1];
            }
        }
        return range[0][1] <= shape[0];
    }

    Shape pos() const {
        Shape _pos(dim);
        for (int i = 0; i < dim; i++) _pos[i] = range[i][0];
        return std::move(_pos);
    }

};

struct ruben_rt {
    double density;
    int ifault;
    double res;
};

ruben_rt ruben(const Arrayd &weights, double c, Arrayd mult=Arrayd{}, Arrayd delta=Arrayd{}, int mode=1, int maxit=100000, double eps=1e-5);

inline double evaluateRuben(double c, double alpha, const Arrayd &weights) {
    auto evaluated = ruben(weights,c);
    double answer    = std::abs(alpha-evaluated.res);
    return answer;
}

}

