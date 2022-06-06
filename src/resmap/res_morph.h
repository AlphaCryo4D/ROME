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

namespace resumap {

/** Helper of distance transform
 */
struct DistanceTransform {
    Shape shape;
    Shape center;
    int dim;
    int min_size;
    int n;

    template<typename T>
    void init(const Array<T> &in, Arrayi &out) {
        shape = in.shape();
        dim = int(shape.size());

        for (int i = 0; i < out.size(); i++) {
            if (in(i) == 0) out(i) = 0;
            else out(i) = 99999;
        }
    }

    void scan3x3(const Shape &center, Arrayi &out) {
        if (out(center) != 0) {
            int min = min_size;
            Shape ind(dim);
            for (int i = 0; i < dim; i++) ind[i] = center[i]-1;
            bool flag;
            do {
                flag = true;
                for (int i = 0; i < dim; i++) {
                    if (ind[i] < 0 || ind[i] >= shape[i]) {
                        flag = false;
                        break;
                    }
                }
                if (flag) {
                    if (min > out(ind)) min = out(ind);
                }
                ind[dim-1]++;
                do {
                    flag = false;
                    for (int i = dim-1; i > 0; i--) {
                        if (ind[i] > center[i]+1) {
                            ind[i] = center[i]-1;
                            ind[i-1]++;
                            flag = true;
                        }
                    }
                } while (flag);
            } while (ind[0]<=center[0]+1);

            if (out(center) != min+1) {
                n++;
                out(center) = min+1;
            }
        }
    }

    void taxicab(int n, Shape &center, Arrayi &out, int direct) {
        int ind, min;
        if (out[n] != 0) {
            min = out[n];
            for (int i = 0; i < dim; i++) {
                center[i] += direct;
                //ind = array_helper::index(out.shape(), out.coeff(), center);
                if (center[i] >= 0 && center[i] < out.shape(i)) {
                    if (min > out(center)) min = out(center);
                    //if (min > out(ind)) min = out(ind);
                }
                center[i] -= direct;
            }
            if (out[n] > min+1) out[n] = min+1;
        }
    }

    template<typename T> Arrayi operator ()(const Array<T> &in) {
        Arrayi out(in.shape());
        init(in, out);
        int i;
        i = 0;
        ARRAY_EACH(shape, center) {
            //MPI_LOG << "center: "; for (auto && n : center) MPI_LOG << n << ' '; MPI_LOG << std::endl;
            taxicab(i, center, out, -1);
            i++;
        }
        i = in.size()-1;
        ARRAY_EACH_R(shape, center) {
            //MPI_LOG << "center: "; for (auto && n : center) MPI_LOG << n << ' '; MPI_LOG << std::endl;
            taxicab(i, center, out, 1);
            i--;
        }
        return std::move(out);
    }
};

/** Same as scipy.ndimage.morphology.distance_transform_cdt.
 */
template<typename T>
Arrayi distance_transform_cdt(const Array<T> &in) {
    return DistanceTransform()(in);
}

inline void taxicab3(Arrayi &out, int a, int b, int c, int direct) {
    const Shape &shape = out.shape();
    int ind = a*shape[1]*shape[2]+b*shape[2]+c;
    int ind2;
    int min;
    if (out[ind] != 0) {
        min = out[ind];
        if (a+direct>=0 && a+direct<shape[0]) {
            ind2 = (a+direct)*shape[1]*shape[2]+b*shape[2]+c;
            if (min > out[ind2]) min = out[ind2];
        }
        if (b+direct>=0 && b+direct<shape[1]) {
            ind2 = a*shape[1]*shape[2]+(b+direct)*shape[2]+c;
			if (min > out[ind2]) min = out[ind2];
        }
        if (c+direct>=0 && c+direct<shape[2]) {
            ind2 = a*shape[1]*shape[2]+b*shape[2]+c+direct;
            if (min > out[ind2]) min = out[ind2];
        }
        if (out[ind] > min+1) out[ind] = min+1;
    }
}

template<typename T>
Arrayi distance_transform_cdt3(const Array<T> &in) {
    Arrayi out(in.shape());
    for (int i = 0; i < out.size(); i++) {
        if (in(i) == 0) out(i) = 0;
        else out(i) = 99999;
    }
    int dim1 = in.shape(0);
    int dim2 = in.shape(1);
    int dim3 = in.shape(2);
    for (int i = 0; i < dim1; i++) {
        for (int j = 0; j < dim2; j++) {
            for (int k = 0; k < dim3; k++) {
                taxicab3(out, i, j, k, -1);
            }
        }
    }
    for (int i = dim1-1; i >= 0; i--) {
        for (int j = dim2-1; j >= 0; j--) {
            for (int k = dim3-1; k >= 0; k--) {
                taxicab3(out, i, j, k, 1);
            }
        }
    }
    return std::move(out);
}

inline Arrayb binary_dilation(const Arrayb &input, int iterations = 1) {
    return Arrayb(input.shape(), true);
}

inline Arrayb binary_dilation(const Arrayb &input, const Arrayd &structure, int iterations = 1) {
    Shape shape = input.shape();
    int nd = shape.size();
    Arrayb p = input;
    Arrayb q = input;

    std::vector<std::vector<int>> range(nd);
    for (auto && r : range) r = {0, 0};

    std::vector<std::array<int, 2>> l(nd);
    for (int i = 0; i < nd; i++) {
        l[i][0] = structure.shape(i)/2;
        l[i][1] = (structure.shape(i)+1)/2;
    }

    for (int it = 0; it < iterations; it++) {
        int i = 0;
        ARRAY_EACH(shape, ind) {
            if (p[i]) {
                for (int j = 0; j < nd; j++) {
                    range[j][0] = ind[j]-l[j][0];
                    range[j][1] = ind[j]+l[j][1];
                }
                q.sub(range) = structure;
            }
            i++;
        }
        p = q;
    }
    return std::move(q);
}

}

