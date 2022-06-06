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

#include "res_ndimage.h"

namespace cppsci { namespace ndimage {

#define TOLERANCE 1e-15

using npy_intp = int;
using Float64 = double;
using UInt32 = unsigned long;

#define MAXDIM 32

struct NI_Iterator {
    int rank_m1;
    npy_intp dimensions[MAXDIM];
    npy_intp coordinates[MAXDIM];
    npy_intp strides[MAXDIM];
    npy_intp backstrides[MAXDIM];
};

typedef enum {
    NI_EXTEND_FIRST = 0,
    NI_EXTEND_NEAREST = 0,
    NI_EXTEND_WRAP = 1,
    NI_EXTEND_REFLECT = 2,
    NI_EXTEND_MIRROR = 3,
    NI_EXTEND_CONSTANT = 4,
    NI_EXTEND_LAST = NI_EXTEND_CONSTANT,
    NI_EXTEND_DEFAULT = NI_EXTEND_MIRROR
} NI_ExtendMode;

static double map_coordinate(double in, npy_intp len, int mode);
static int zoom_shift(Arrayd &input, Arrayd * zoom_ar, Arrayd * shift_ar, 
        Arrayd & output, int order, int mode, double cval);

/**Convert an extension mode to the corresponding integer code.
  */
static int extend_mode_to_code(std::string mode) {
       if (mode == "nearest") return 0;
       else if (mode == "wrap") return 1;
       else if (mode == "reflect") return 2;
       else if (mode == "mirror") return 3;
       else if (mode == "constant") return 4;
       else throw __FUNCTION__;
}

static void
spline_coefficients(double x, int order, double *result)
{
    int hh;
    double y, start;

    if (order & 1) {
        start = (int)floor(x) - order / 2;
    } else {
        start = (int)floor(x + 0.5) - order / 2;
    }

    for(hh = 0; hh <= order; hh++)  {
        y = fabs(start - x + hh);

        switch(order) {
        case 1:
            result[hh] = y > 1.0 ? 0.0 : 1.0 - y;
            break;
        case 2:
            if (y < 0.5) {
                result[hh] = 0.75 - y * y;
            } else if (y < 1.5) {
                y = 1.5 - y;
                result[hh] = 0.5 * y * y;
            } else {
                result[hh] = 0.0;
            }
            break;
        case 3:
            if (y < 1.0) {
                result[hh] =
                    (y * y * (y - 2.0) * 3.0 + 4.0) / 6.0;
            } else if (y < 2.0) {
                y = 2.0 - y;
                result[hh] = y * y * y / 6.0;
            } else {
                result[hh] = 0.0;
            }
            break;
        case 4:
            if (y < 0.5) {
                y *= y;
                result[hh] = y * (y * 0.25 - 0.625) + 115.0 / 192.0;
            } else if (y < 1.5) {
                result[hh] = y * (y * (y * (5.0 / 6.0 - y / 6.0) - 1.25) +
                                                    5.0 / 24.0) + 55.0 / 96.0;
            } else if (y < 2.5) {
                y -= 2.5;
                y *= y;
                result[hh] = y * y / 24.0;
            } else {
                result[hh] = 0.0;
            }
            break;
        case 5:
            if (y < 1.0) {
                double f = y * y;
                result[hh] =
                    f * (f * (0.25 - y / 12.0) - 0.5) + 0.55;
            } else if (y < 2.0) {
                result[hh] = y * (y * (y * (y * (y / 24.0 - 0.375)
                                                                        + 1.25) -  1.75) + 0.625) + 0.425;
            } else if (y < 3.0) {
                double f = 3.0 - y;
                y = f * f;
                result[hh] = f * y * y / 120.0;
            } else {
                result[hh] = 0.0;
            }
            break;
        }
    }
}

static void spline_filter_set_pars(int order, int &npoles, double (&pole)[2], double &weight) {
    switch (order) {
        case 2:
            npoles = 1;
            pole[0] = sqrt(8.0) - 3.0;
            break;
        case 3:
            npoles = 1;
            pole[0] = sqrt(3.0) - 2.0;
            break;
        case 4:
            npoles = 2;
            pole[0] = sqrt(664.0 - sqrt(438976.0)) + sqrt(304.0) - 19.0;
            pole[1] = sqrt(664.0 + sqrt(438976.0)) - sqrt(304.0) - 19.0;
            break;
        case 5:
            npoles = 2;
            pole[0] = sqrt(67.5 - sqrt(4436.25)) + sqrt(26.25) - 6.5;
            pole[1] = sqrt(67.5 + sqrt(4436.25)) - sqrt(26.25) - 6.5;
            break;
        default:
            break;
    }

    weight = 1.0;
    for(int hh = 0; hh < npoles; hh++) weight *= (1.0 - pole[hh]) * (1.0 - 1.0 / pole[hh]);
}

/* one-dimensional spline filter: */
static void spline_filter1d(double *ln, int len, int order) {
    using namespace std;

    if (order < 0 || order > 5) throw "error";
    if (order == 0 || order == 1) return;

    int npoles = 0;
    double pole[2];
    double weight;
    spline_filter_set_pars(order, npoles, pole, weight);
    int ll, hh;

    /* spline filter: */
    for(ll = 0; ll < len; ll++)
        ln[ll] *= weight;
    for(hh = 0; hh < npoles; hh++) {
        double p = pole[hh];
        int max = (int)ceil(log(TOLERANCE) / log(fabs(p)));
        if (max < len) {
            double zn = p;
            double sum = ln[0];
            for(ll = 1; ll < max; ll++) {
                sum += zn * ln[ll];
                zn *= p;
            }
            ln[0] = sum;
        } else {
            double zn = p;
            double iz = 1.0 / p;
            double z2n = pow(p, (double)(len - 1));
            double sum = ln[0] + z2n * ln[len - 1];
            z2n *= z2n * iz;
            for(ll = 1; ll <= len - 2; ll++) {
                sum += (zn + z2n) * ln[ll];
                zn *= p;
                z2n *= iz;
            }
            ln[0] = sum / (1.0 - zn * zn);
        }
        for(ll = 1; ll < len; ll++)
            ln[ll] += p * ln[ll - 1];
        ln[len-1] = (p / (p * p - 1.0)) * (ln[len-1] + p * ln[len-2]);
        for(ll = len - 2; ll >= 0; ll--)
            ln[ll] = p * (ln[ll + 1] - ln[ll]);
    }
}

void spline_filter(Arrayd &input, int order) {
    ERROR_CHECK(order < 2 || order > 5, "order should be in [2, 5] !");
    int nd = input.dim();
    Arrayi c({nd});
    for (int i = 0; i < nd; i++) {
        int li = input.shape(i);
        int k = li;
        for (int j = nd-1; j >= 0; j--) {
            if (j == i) c[j] = 1;
            else {
                c[j] = k;
                k *= input.shape(j);
            }
        }
        Arrayd temp({input.size()/li, li});
        int m = 0;
        ARRAY_EACH(input.shape(), ind) {
            int n = 0;
            for (int j = 0; j < nd; j++) n += c[j]*ind[j];
            temp[n] = input[m];
            m++;
        }
        for (int j = 0; j < temp.shape(0); j++) {
            spline_filter1d(temp.data_wptr() + j*li, li, order);
        }
        m = 0;
        ARRAY_EACH(input.shape(), ind) {
            int n = 0;
            for (int j = 0; j < nd; j++) n += c[j]*ind[j];
            input[m] = temp[n];
            m++;
        }
    }
}

Arrayd zoom(const Arrayd &input, double _zoom, std::string _mode, int order, double cval, bool prefilter) {
    int nd = input.dim();
    if (nd < 1) throw "input and output rank must be > 0";
    int mode = extend_mode_to_code(_mode);
    Arrayd filtered = input;
    if (prefilter && order > 1) spline_filter(filtered, order);
    else filtered = input;
//    MPI_LOG << "filtered:\n" << filtered << std::endl;

    Arrayd zoom({input.dim()}, _zoom);
    Shape output_shape(nd);
    for (int i = 0; i < nd; i++) {
        output_shape[i] = int(round(input.shape(i) * _zoom));
//        MPI_LOG << input.shape(i) << ' ' << output_shape[i] << ' ' << _zoom << ' ' << std::round(input.shape(i)*_zoom) << std::endl;
    }
    for (int i = 0; i < nd; i++) {
        int d = output_shape[i] - 1;
        if (d < 1e-5) zoom[i] = 1;
        else zoom[i] = (input.shape(i) - 1) / (double)d;
    }
//    MPI_LOG << "zoom:\n";
//    zoom.printFull();
    Arrayd output(output_shape);
    zoom_shift(filtered, &zoom, NULL, output, order, mode, cval);
    return std::move(output);

}

#define NI_ITERATOR_NEXT(iterator, pointer)                         \
do {                                                                   \
    int _ii;                                                          \
    for(_ii = (iterator).rank_m1; _ii >= 0; _ii--)                    \
        if ((iterator).coordinates[_ii] < (iterator).dimensions[_ii]) { \
            (iterator).coordinates[_ii]++;                                \
            pointer += (iterator).strides[_ii];                           \
            break;                                                        \
        } else {                                                        \
            (iterator).coordinates[_ii] = 0;                              \
            pointer -= (iterator).backstrides[_ii];                       \
        }                                                               \
} while(0)

/* go to the next point in two arrays of the same size */
#define NI_ITERATOR_NEXT2(iterator1, iterator2,  pointer1, pointer2)  \
{                                                                     \
    int _ii;                                                            \
    for(_ii = (iterator1).rank_m1; _ii >= 0; _ii--)                     \
        if ((iterator1).coordinates[_ii] < (iterator1).dimensions[_ii]) { \
            (iterator1).coordinates[_ii]++;                                 \
            pointer1 += (iterator1).strides[_ii];                           \
            pointer2 += (iterator2).strides[_ii];                           \
            break;                                                          \
        } else {                                                          \
            (iterator1).coordinates[_ii] = 0;                               \
            pointer1 -= (iterator1).backstrides[_ii];                       \
            pointer2 -= (iterator2).backstrides[_ii];                       \
        }                                                                 \
}

int NI_InitPointIterator(const Arrayd &array, NI_Iterator *iterator)
{
    int ii;

    iterator->rank_m1 = array.dim() - 1;
    for(ii = 0; ii < array.dim(); ii++) {
        /* adapt dimensions for use in the macros: */
        iterator->dimensions[ii] = array.shape(ii) - 1;
        /* initialize coordinates: */
        iterator->coordinates[ii] = 0;
        /* initialize strides: */
        iterator->strides[ii] = sizeof(double)*array.coeff(ii);
        /* calculate the strides to move back at the end of an axis: */
        iterator->backstrides[ii] =
            sizeof(double)*array.coeff(ii) * iterator->dimensions[ii];
    }
    return 1;
}

static int zoom_shift(Arrayd &input, Arrayd * zoom_ar,
        Arrayd * shift_ar, Arrayd & output,
        int order, int mode, double cval)
{
    char *po, *pi;
    npy_intp **zeros = NULL, **offsets = NULL, ***edge_offsets = NULL;
    npy_intp ftmp[MAXDIM], *fcoordinates = NULL, *foffsets = NULL;
    npy_intp jj, hh, kk, filter_size, odimensions[MAXDIM];
    npy_intp idimensions[MAXDIM], istrides[MAXDIM];
    npy_intp size;
    double ***splvals = NULL;
    NI_Iterator io;
    double *zooms = (zoom_ar ? zoom_ar->data_wptr() : NULL);
    double *shifts = (shift_ar ? shift_ar->data_wptr() : NULL);
    int rank = 0, qq;

    bool error = false;

    for(kk = 0; kk < input.dim(); kk++) {
        idimensions[kk] = input.shape(kk);
        istrides[kk] = sizeof(double)*input.coeff(kk);
        odimensions[kk] = output.shape(kk);
    }
    rank = input.dim();

    /* if the mode is 'constant' we need some temps later: */
    if (mode == NI_EXTEND_CONSTANT) {
        zeros = (npy_intp**)aMalloc(rank * sizeof(npy_intp*), 64);
        if (!zeros) {
            error = true;
            goto exit;
        }
        for(jj = 0; jj < rank; jj++)
            zeros[jj] = NULL;
        for(jj = 0; jj < rank; jj++) {
            zeros[jj] = (npy_intp*)aMalloc(odimensions[jj] * sizeof(npy_intp), 64);
            if(!zeros[jj]) {
                error = true;
                goto exit;
            }
        }
    }

    /* store offsets, along each axis: */
    offsets = (npy_intp**)aMalloc(rank * sizeof(npy_intp*), 64);
    /* store spline coefficients, along each axis: */
    splvals = (double***)aMalloc(rank * sizeof(double**), 64);
    /* store offsets at all edges: */
    edge_offsets = (npy_intp***)aMalloc(rank * sizeof(npy_intp**), 64);
    if (!offsets || !splvals || !edge_offsets) {
        error = true;
        goto exit;
    }
    for(jj = 0; jj < rank; jj++) {
        offsets[jj] = NULL;
        splvals[jj] = NULL;
        edge_offsets[jj] = NULL;
    }
    for(jj = 0; jj < rank; jj++) {
        offsets[jj] = (npy_intp*)aMalloc(odimensions[jj] * sizeof(npy_intp), 64);
        splvals[jj] = (double**)aMalloc(odimensions[jj] * sizeof(double*), 64);
        edge_offsets[jj] = (npy_intp**)aMalloc(odimensions[jj] * sizeof(npy_intp*), 64);
        if (!offsets[jj] || !splvals[jj] || !edge_offsets[jj]) {
            error = true;
            goto exit;
        }
        for(hh = 0; hh < odimensions[jj]; hh++) {
            splvals[jj][hh] = NULL;
            edge_offsets[jj][hh] = NULL;
        }
    }

    /* precalculate offsets, and offsets at the edge: */
    for(jj = 0; jj < rank; jj++) {
        double shift = 0.0, zoom = 0.0;
        if (shifts)
            shift = shifts[jj];
        if (zooms)
            zoom = zooms[jj];
        for(kk = 0; kk < odimensions[jj]; kk++) {
            double cc = (double)kk;
            if (shifts)
                cc += shift;
            if (zooms)
                cc *= zoom;
            cc = map_coordinate(cc, idimensions[jj], mode);
            if (cc > -1.0) {
                int start;
                if (zeros && zeros[jj])
                    zeros[jj][kk] = 0;
                if (order & 1) {
                    start = (int)floor(cc) - order / 2;
                } else {
                    start = (int)floor(cc + 0.5) - order / 2;
                }
                offsets[jj][kk] = istrides[jj] * start;
                if (start < 0 || start + order >= idimensions[jj]) {
                    edge_offsets[jj][kk] = (npy_intp*)aMalloc((order + 1) * sizeof(npy_intp), 64);
                    if (!edge_offsets[jj][kk]) {
                        error = true;
                        goto exit;
                    }
                    for(hh = 0; hh <= order; hh++) {
                        int idx = start + hh;
                         int len = idimensions[jj];
                        if (len <= 1) {
                            idx = 0;
                        } else {
                            int s2 = 2 * len - 2;
                            if (idx < 0) {
                                idx = s2 * (int)(-idx / s2) + idx;
                                idx = idx <= 1 - len ? idx + s2 : -idx;
                            } else if (idx >= len) {
                                idx -= s2 * (int)(idx / s2);
                                if (idx >= len)
                                    idx = s2 - idx;
                            }
                        }
                        edge_offsets[jj][kk][hh] = istrides[jj] * (idx - start);
                    }
                }
                if (order > 0) {
                    splvals[jj][kk] = (double*)aMalloc((order + 1) * sizeof(double), 64);
                    if (!splvals[jj][kk]) {
                        error = true;
                        goto exit;
                    }
                    spline_coefficients(cc, order, splvals[jj][kk]);
                }
            } else {
                zeros[jj][kk] = 1;
            }
        }
    }

    filter_size = 1;
    for(jj = 0; jj < rank; jj++)
        filter_size *= order + 1;

    if (!NI_InitPointIterator(output, &io))
        goto exit;

    //pi = (void *)PyArray_DATA(input);
    //po = (void *)PyArray_DATA(output);
    pi = reinterpret_cast<char *>(input.data_wptr());
    po = reinterpret_cast<char *>(output.data_wptr());

    /* store all coordinates and offsets with filter: */
    fcoordinates = (npy_intp*)aMalloc(rank * filter_size * sizeof(npy_intp), 64);
    foffsets = (npy_intp*)aMalloc(filter_size * sizeof(npy_intp), 64);
    if (!fcoordinates || !foffsets) {
        error = true;
        goto exit;
    }

    for(jj = 0; jj < rank; jj++)
        ftmp[jj] = 0;
    kk = 0;
    for(hh = 0; hh < filter_size; hh++) {
        for(jj = 0; jj < rank; jj++)
            fcoordinates[jj + hh * rank] = ftmp[jj];
        foffsets[hh] = kk;
        for(jj = rank - 1; jj >= 0; jj--) {
            if (ftmp[jj] < order) {
                ftmp[jj]++;
                kk += istrides[jj];
                break;
            } else {
                ftmp[jj] = 0;
                kk -= istrides[jj] * order;
            }
        }
    }
    size = 1;
    for(qq = 0; qq < output.dim(); qq++)
        size *= output.shape(qq);
    for(kk = 0; kk < size; kk++) {
        double t = 0.0;
        int edge = 0, oo = 0, zero = 0;

        for(hh = 0; hh < rank; hh++) {
            if (zeros && zeros[hh][io.coordinates[hh]]) {
                /* we use constant border condition */
                zero = 1;
                break;
            }
            oo += offsets[hh][io.coordinates[hh]];
            if (edge_offsets[hh][io.coordinates[hh]])
                edge = 1;
        }

        if (!zero) {
            npy_intp *ff = fcoordinates;
            t = 0.0;
            for(hh = 0; hh < filter_size; hh++) {
                int idx = 0;
                double coeff = 0.0;

                if (edge) {
                        /* use precalculated edge offsets: */
                    for(jj = 0; jj < rank; jj++) {
                        if (edge_offsets[jj][io.coordinates[jj]])
                            idx += edge_offsets[jj][io.coordinates[jj]][ff[jj]];
                        else
                            idx += ff[jj] * istrides[jj];
                    }
                    idx += oo;
                } else {
                    /* use normal offsets: */
                    idx += oo + foffsets[hh];
                }
                coeff = *(double*)(pi + idx);
                /* calculate interpolated value: */
                for(jj = 0; jj < rank; jj++)
                    if (order > 0)
                        coeff *= splvals[jj][io.coordinates[jj]][ff[jj]];
                t += coeff;
                ff += rank;
            }
        } else {
            t = cval;
        }
        /* store output: */
        *(double*)po = (double)t;
        NI_ITERATOR_NEXT(io, po);
    }

 exit:
    if (zeros) {
        for(jj = 0; jj < rank; jj++)
            if (zeros[jj])
                aFree(zeros[jj]);
        aFree(zeros);
    }
    if (offsets) {
        for(jj = 0; jj < rank; jj++)
            if (offsets[jj])
                aFree(offsets[jj]);
        aFree(offsets);
    }
    if (splvals) {
        for(jj = 0; jj < rank; jj++) {
            if (splvals[jj]) {
                for(hh = 0; hh < odimensions[jj]; hh++)
                    if (splvals[jj][hh])
                        aFree(splvals[jj][hh]);
                aFree(splvals[jj]);
            }
        }
        aFree(splvals);
    }
    if (edge_offsets) {
        for(jj = 0; jj < rank; jj++) {
            if (edge_offsets[jj]) {
                for(hh = 0; hh < odimensions[jj]; hh++)
                    if (edge_offsets[jj][hh])
                        aFree(edge_offsets[jj][hh]);
                aFree(edge_offsets[jj]);
            }
        }
        aFree(edge_offsets);
    }
    if (foffsets)
        aFree(foffsets);
    if (fcoordinates)
        aFree(fcoordinates);
    return !error;
}

/* map a coordinate outside the borders, according to the requested
     boundary condition: */
static double map_coordinate(double in, npy_intp len, int mode)
{
    if (in < 0) {
        switch (mode) {
        case NI_EXTEND_MIRROR:
            if (len <= 1) {
                in = 0;
            } else {
                npy_intp sz2 = 2 * len - 2;
                in = sz2 * (npy_intp)(-in / sz2) + in;
                in = in <= 1 - len ? in + sz2 : -in;
            }
            break;
        case NI_EXTEND_REFLECT:
            if (len <= 1) {
                in = 0;
            } else {
                npy_intp sz2 = 2 * len;
                if (in < -sz2)
                    in = sz2 * (npy_intp)(-in / sz2) + in;
                in = in < -len ? in + sz2 : -in - 1;
            }
            break;
        case NI_EXTEND_WRAP:
            if (len <= 1) {
                in = 0;
            } else {
                npy_intp sz = len - 1;
                // Integer division of -in/sz gives (-in mod sz)
                // Note that 'in' is negative
                in += sz * ((npy_intp)(-in / sz) + 1);
            }
            break;
        case NI_EXTEND_NEAREST:
            in = 0;
            break;
        case NI_EXTEND_CONSTANT:
            in = -1;
            break;
        }
    } else if (in > len-1) {
        switch (mode) {
        case NI_EXTEND_MIRROR:
            if (len <= 1) {
                in = 0;
            } else {
                npy_intp sz2 = 2 * len - 2;
                in -= sz2 * (npy_intp)(in / sz2);
                if (in >= len)
                    in = sz2 - in;
            }
            break;
        case NI_EXTEND_REFLECT:
            if (len <= 1) {
                in = 0;
            } else {
                npy_intp sz2 = 2 * len;
                in -= sz2 * (npy_intp)(in / sz2);
                if (in >= len)
                    in = sz2 - in - 1;
            }
            break;
        case NI_EXTEND_WRAP:
            if (len <= 1) {
                in = 0;
            } else {
                npy_intp sz = len - 1;
                in -= sz * (npy_intp)(in / sz);
            }
            break;
        case NI_EXTEND_NEAREST:
            in = len - 1;
            break;
        case NI_EXTEND_CONSTANT:
            in = -1;
            break;
        }
    }

    return in;
}

/* initialize iteration over a lower sub-space: */
static int NI_SubspaceIterator(NI_Iterator *iterator, UInt32 axes)
{
    int ii, last = 0;

    for(ii = 0; ii <= iterator->rank_m1; ii++) {
        if (axes & (((UInt32)1) << ii)) {
            if (last != ii) {
                iterator->dimensions[last] = iterator->dimensions[ii];
                iterator->strides[last] = iterator->strides[ii];
                iterator->backstrides[last] = iterator->backstrides[ii];
            }
            ++last;
        }
    }
    iterator->rank_m1 = last - 1;
    return 1;
}

/* initialize iteration over array lines: */
static int NI_LineIterator(NI_Iterator *iterator, int axis)
{
    UInt32 axes = ((UInt32)1) << axis;
    return NI_SubspaceIterator(iterator, ~axes);
}

static void basic_geometric_transform(/*const*/ Arrayd &input, 
        int (*map)(npy_intp*, double*, int, int, void*),
        void* map_data, /*const*/ Arrayd &coordinates,
        /*const*/ Arrayd &matrix_ar, /*const*/ Arrayd &shift_ar, 
        Arrayd &output, int order, int mode, double cval)
{
    using namespace std;

    std::string error;

    char *po, *pi, *pc = NULL;
    npy_intp **edge_offsets = NULL, **data_offsets = NULL, filter_size;
    npy_intp ftmp[MAXDIM], *fcoordinates = NULL, *foffsets = NULL;
    npy_intp cstride = 0, kk, hh, ll, jj;
    npy_intp size;
    double **splvals = NULL, icoor[MAXDIM];
    npy_intp idimensions[MAXDIM], istrides[MAXDIM];
    NI_Iterator io, ic;
    double *matrix = matrix_ar.data_wptr();
    double *shift = shift_ar.data_wptr();
//    Float64 *matrix = matrix_ar ? (Float64*)PyArray_DATA(matrix_ar) : NULL;
//    Float64 *shift = shift_ar ? (Float64*)PyArray_DATA(shift_ar) : NULL;
    int irank = 0, orank, qq;

//    MPI_LOG << "input" << std::endl;
//    input.printFull();
//    MPI_LOG << "coordinates" << std::endl;
//    coordinates.printFull();

    for(kk = 0; kk < input.dim(); kk++) {
        idimensions[kk] = input.shape(kk);
        istrides[kk] = input.coeff(kk)*sizeof(double);
    }
    irank = input.dim();
    orank = output.dim();

    /* if the mapping is from array coordinates: */
    if (coordinates.size() != 0) {
        /* initialze a line iterator along the first axis: */
        if (!NI_InitPointIterator(coordinates, &ic))
            goto exit;
        cstride = ic.strides[0];
        if (!NI_LineIterator(&ic, 0))
            goto exit;
//        pc = (void *)(PyArray_DATA(coordinates));
        pc = reinterpret_cast<char *>(coordinates.data_wptr());
    }

    /* offsets used at the borders: */
    edge_offsets = (npy_intp**)aMalloc(irank * sizeof(npy_intp*), 64);
    data_offsets = (npy_intp**)aMalloc(irank * sizeof(npy_intp*), 64);
    if (!edge_offsets || !data_offsets) {
        error = "no memory";
        goto exit;
    }
    for(jj = 0; jj < irank; jj++)
        data_offsets[jj] = NULL;
    for(jj = 0; jj < irank; jj++) {
        data_offsets[jj] = (npy_intp*)aMalloc((order + 1) * sizeof(npy_intp), 64);
        if (!data_offsets[jj]) {
            error = "no memory";
            goto exit;
        }
    }
    /* will hold the spline coefficients: */
    splvals = (double**)aMalloc(irank * sizeof(double*), 64);
    if (!splvals) {
        error = "no memory";
        goto exit;
    }
    for(jj = 0; jj < irank; jj++)
        splvals[jj] = NULL;
    for(jj = 0; jj < irank; jj++) {
        splvals[jj] = (double*)aMalloc((order + 1) * sizeof(double), 64);
        if (!splvals[jj]) {
            error = "no memory";
            goto exit;
        }
    }

    filter_size = 1;
    for(jj = 0; jj < irank; jj++)
        filter_size *= order + 1;

    /* initialize output iterator: */
    if (!NI_InitPointIterator(output, &io))
        goto exit;

    /* get data pointers: */
//    pi = (void *)PyArray_DATA(input);
//    po = (void *)PyArray_DATA(output);
    pi = reinterpret_cast<char *>(input.data_wptr());
    po = reinterpret_cast<char *>(output.data_wptr());

    /* make a table of all possible coordinates within the spline filter: */
    fcoordinates = (npy_intp*)aMalloc(irank * filter_size * sizeof(npy_intp), 64);
    /* make a table of all offsets within the spline filter: */
    foffsets = (npy_intp*)aMalloc(filter_size * sizeof(npy_intp), 64);
    if (!fcoordinates || !foffsets) {
        error = "no memory";
        goto exit;
    }
    for(jj = 0; jj < irank; jj++)
        ftmp[jj] = 0;
    kk = 0;
    for(hh = 0; hh < filter_size; hh++) {
        for(jj = 0; jj < irank; jj++)
            fcoordinates[jj + hh * irank] = ftmp[jj];
        foffsets[hh] = kk;
        for(jj = irank - 1; jj >= 0; jj--) {
            if (ftmp[jj] < order) {
                ftmp[jj]++;
                kk += istrides[jj];
                break;
            } else {
                ftmp[jj] = 0;
                kk -= istrides[jj] * order;
            }
        }
    }

    size = 1;
    for(qq = 0; qq < output.dim(); qq++)
        size *= output.shape(qq);
    for(kk = 0; kk < size; kk++) {
        double t = 0.0;
        int constant = 0, edge = 0, offset = 0;
        if (map) {
            /* call mappint functions: */
            if (!map(io.coordinates, icoor, orank, irank, map_data)) {
                error = "unknown error in mapping function";
                goto exit;
            }
        } else if (matrix) {
            /* do an affine transformation: */
            Float64 *p = matrix;
            for(hh = 0; hh < irank; hh++) {
                icoor[hh] = 0.0;
                for(ll = 0; ll < orank; ll++)
                    icoor[hh] += io.coordinates[ll] * (*p++);
                icoor[hh] += shift[hh];
            }
        } else if (coordinates.size() != 0) {
            /* mapping is from an coordinates array: */
            char *p = pc;
            for(int _hh = 0; _hh < irank; _hh++) {
                icoor[_hh] = *(double *)p;
                p += cstride;
            }                                                            
        }
        /* iterate over axes: */
        for(hh = 0; hh < irank; hh++) {
            /* if the input coordinate is outside the borders, map it: */
            double cc = map_coordinate(icoor[hh], idimensions[hh], mode);
            if (cc > -1.0) {
                /* find the filter location along this axis: */
                int start;
                if (order & 1) {
                    start = (int)floor(cc) - order / 2;
                } else {
                    start = (int)floor(cc + 0.5) - order / 2;
                }
                /* get the offset to the start of the filter: */
                offset += istrides[hh] * start;
                if (start < 0 || start + order >= idimensions[hh]) {
                    /* implement border mapping, if outside border: */
                    edge = 1;
                    edge_offsets[hh] = data_offsets[hh];
                    for(ll = 0; ll <= order; ll++) {
                        int idx = start + ll;
                        int len = idimensions[hh];
                        if (len <= 1) {
                            idx = 0;
                        } else {
                            int s2 = 2 * len - 2;
                            if (idx < 0) {
                                idx = s2 * (int)(-idx / s2) + idx;
                                idx = idx <= 1 - len ? idx + s2 : -idx;
                            } else if (idx >= len) {
                                idx -= s2 * (int)(idx / s2);
                                if (idx >= len)
                                    idx = s2 - idx;
                            }
                        }
                        /* calculate and store the offests at this edge: */
                        edge_offsets[hh][ll] = istrides[hh] * (idx - start);
                    }
                } else {
                    /* we are not at the border, use precalculated offsets: */
                    edge_offsets[hh] = NULL;
                }
                spline_coefficients(cc, order, splvals[hh]);
            } else {
                /* we use the constant border condition: */
                constant = 1;
                break;
            }
        }

        if (!constant) {
            npy_intp *ff = fcoordinates;
            t = 0.0;
            for(hh = 0; hh < filter_size; hh++) {
                double coeff = 0.0;
                int idx = 0;

                if (edge) {
                    for(ll = 0; ll < irank; ll++) {
                        if (edge_offsets[ll])
                            idx += edge_offsets[ll][ff[ll]];
                        else
                            idx += ff[ll] * istrides[ll];
                    }
                } else {
                    idx = foffsets[hh];
                }
                idx += offset;
                coeff = *(double *)(pi+idx);
                /* calculate the interpolated value: */
                for(ll = 0; ll < irank; ll++)
                    if (order > 0)
                        coeff *= splvals[ll][ff[ll]];
                t += coeff;
                ff += irank;
            }
        } else {
            t = cval;
        }
        /* store output value: */
        *(double *)po = t;
//        CASE_INTERP_OUT(po, t, Float64);
        if (coordinates.size() != 0) {
            NI_ITERATOR_NEXT2(io, ic, po, pc);
        } else {
            NI_ITERATOR_NEXT(io, po);
        }
    }

 exit:
    if (edge_offsets)
        aFree(edge_offsets);
    if (data_offsets) {
        for(jj = 0; jj < irank; jj++)
            aFree(data_offsets[jj]);
        aFree(data_offsets);
    }
    if (splvals) {
        for(jj = 0; jj < irank; jj++)
            aFree(splvals[jj]);
        aFree(splvals);
    }
    if (foffsets)
        aFree(foffsets);
    if (fcoordinates)
        aFree(fcoordinates);

    if (error.empty()) return;
    else throw error;
}

Arrayd map_coordinates(const Arrayd &input, /*const*/ Arrayd &coordinates,
        int order, std::string mode, double cval, bool prefilter)
{
    ERROR_CHECK(order < 0 || order > 5, "spline order not supported");

    int nd = coordinates.dim()-1;
    Shape output_shape(nd);
    for (int i = 0; i < nd; i++) output_shape[i] = coordinates.shape(i+1);
    ERROR_CHECK(input.dim() < 1 || output_shape.size() < 1, "input and output rank must be > 0");

    ERROR_CHECK(coordinates.shape(0) != input.dim(), "invalid shape for coordinate array");

    int mode_code = extend_mode_to_code(mode);

    Arrayd filtered = input;
    if (prefilter && order > 1) spline_filter(filtered, order);

    Arrayd output(output_shape);
    auto matrix_ar = Arrayd{};
    auto shift_ar = Arrayd{};
    basic_geometric_transform(filtered, NULL, NULL, coordinates, matrix_ar, shift_ar, output, order, mode_code, cval);
    return std::move(output);

}

}} // namespace cppsci::ndimage

