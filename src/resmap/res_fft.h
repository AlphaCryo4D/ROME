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

#include "fftw/fftw3.h"
#include "omp.h"

namespace resumap {

/**
 * N-dimensional fast fourier transform.
 */
template<typename T, typename V>
void fft(T in, V out, int dim, int l, int direct = -1) {
    int *v = (int*)aMalloc(dim*sizeof(int), 64);
    for (int i = 0; i < dim; i++) v[i] = l;
    fftw_plan plan = fftw_plan_dft(dim, v, reinterpret_cast<fftw_complex*>(in), reinterpret_cast<fftw_complex*>(out), direct, FFTW_ESTIMATE);
    aFree(v);
    fftw_execute(plan);
    fftw_destroy_plan(plan);
}

/**
 * 1-dimensional fast fourier transform.
 */
template<typename T, typename V>
void fft1(T in, V out, int l, int direct = -1) {
    fftw_plan plan = fftw_plan_dft_1d(l, reinterpret_cast<fftw_complex*>(in), reinterpret_cast<fftw_complex*>(out), direct, FFTW_ESTIMATE);
    fftw_execute(plan);
    fftw_destroy_plan(plan);
}

/**
 * 2-dimensional fast fourier transform.
 */
template<typename T, typename V>
void fft2(T in, V out, int l, int direct = -1) {
    fftw_plan plan = fftw_plan_dft_2d(l, l, reinterpret_cast<fftw_complex*>(in), reinterpret_cast<fftw_complex*>(out), direct, FFTW_ESTIMATE);
    fftw_execute(plan);
    fftw_destroy_plan(plan);
}

/**
 * 3-dimensional fast fourier transform.
 */
template<typename T, typename V>
void fft3(T in, V out, int l, int direct = -1) {
    fftw_plan plan = fftw_plan_dft_3d(l, l, l, reinterpret_cast<fftw_complex*>(in), reinterpret_cast<fftw_complex*>(out), direct, FFTW_ESTIMATE);
    fftw_execute(plan);
    fftw_destroy_plan(plan);
}

/**
 * N-dimensional inverse fast fourier transform.
 */
template<typename T, typename V>
void ifft(T in, V out, int dim, int l) {
    fft(in, out, dim, l, 1);
    int k = 1;
    for (int i = 0; i < dim; i++) k *= l;
    for (int i = 0; i < k; i++) out[i] /= k;
}

/**
 * 1-dimensional inverse fast fourier transform.
 */
template<typename T, typename V>
void ifft1(T in, V out, int l) {
    fft(in, out, l, 1);
    for (int i = 0; i < l; i++) out[i] /= l;
}

/**
 * 2-dimensional inverse fast fourier transform.
 */
template<typename T, typename V>
void ifft2(T in, V out, int l) {
    fft2(in, out, l, 1);
    for (int i = 0; i < l*l; i++) out[i] /= l*l;
}

/**
 * 3-dimensional inverse fast fourier transform.
 */
template<typename T, typename V>
void ifft3(T in, V out, int l) {
    fft3(in, out, l, 1);
    for (int i = 0; i < l*l*l; i++) out[i] /= l*l*l;
}

/**
 * N-dimensional shift of the result of fast fourier transform.
 */
template<typename T, typename V>
void fftshift(T in, V out, int dim, int l) {
    int b = (int)(l / 2);
    std::vector<int> v(dim, 0);
    std::vector<int> k(dim, 1);
    for (int i = dim-2; i>=0; i--) k[i]=k[i+1]*l;
    auto ind = [&dim, &k](const std::vector<int> &v) -> int {
        int sum = 0;
        for (int i = 0; i < dim; i++) sum += v[i]*k[i];
        return sum;
    };
    auto symm = [&b](std::vector<int> v) -> std::vector<int> {
        for (auto && n : v) {
            if (n < b) n += b;
            else n -=b;
        }
        return std::move(v);
    };
    auto next = [&dim, &l](std::vector<int> &v, int lim) -> bool {
        v[dim-1]++;
        for (int i = dim-1; i>=1; i--) {
            if (v[i] > l) {
                v[i] = 0;
                v[i-1]++;
            }
        }
        return v[0] < lim && v[0] >= 0;
    };

    if (in != out) do { out[ind(v)] = in[ind(symm(v))];          } while (next(v, l)); 
    else           do { std::swap(in[ind(v)], in[ind(symm(v))]); } while (next(v, b));
}


/**
 * 1-dimensional shift of the result of fast fourier transform.
 */
template<typename T, typename V>
void fftshift1(T in, V out, int l) {
    int b = (int)(l / 2);

    if (in != out) {
        for (int i = 0; i < b; i++) out[i] = in[i+b];
        for (int i = b; i < l; i++) out[i] = in[i-b];
    }
    else {
        for (int i = 0; i < b; i++) std::swap(in[i], in[i+b]);
    }
}

/**
 * 2-dimensional shift of the result of fast fourier transform.
 */
template<typename T, typename V>
void fftshift2(T in, V out, int l) {
    int b = (int)(l / 2);

    if (in != out) {
        for (int i = 0; i < b; i++) {
            for (int j = 0; j < b; j++) out[i*l+j] = in[(i+b)*l+j+b];
            for (int j = b; j < l; j++) out[i*l+j] = in[(i+b)*l+j-b];
        }
        for (int i = b; i < l; i++) {
            for (int j = 0; j < b; j++) out[i*l+j] = in[(i-b)*l+j+b];
            for (int j = b; j < l; j++) out[i*l+j] = in[(i-b)*l+j-b];
        }
    }
    else {
        for (int i = 0; i < b; i++) {
            for (int j = 0; j < b; j++) std::swap(in[i*l+j], in[(i+b)*l+j+b]);
            for (int j = b; j < l; j++) std::swap(in[i*l+j], in[(i+b)*l+j-b]);
        }
    }
}

/**
 * 3-dimensional shift of the result of fast fourier transform.
 */
template<typename T, typename V>
void fftshift3(T in, V out, int l) {
    int b = (int)(l / 2);
    int l2 = l*l;

    if (in != out) {
        for (int k = 0; k < l; k++) {
            int kp = k + b;
            if (kp < 0) kp += l;
            else if (kp >= l) kp -= l;

            for (int i = 0; i < l; i++) {
                int ip = i + b;
                if (ip < 0) ip += l;
                else if (ip >= l) ip -= l;

                for (int j = 0; j < l; j++) {
                    int jp = j + b;
                    if (jp < 0) jp += l;
                    else if (jp >= l) jp -= l;

                    out[kp*l2+ip*l+jp] = in[k*l2+i*l+j];

                }
            }
        }
    }
    else {
        for (int i = 0; i < l/2; i++) {
            for (int j = 0; j < l/2; j++) {
                for (int k = 0; k < l/2; k++) {
                    std::swap(in[i*l2+j*l+k], out[(i+b)*l2+(j+b)*l+(k+b)]);
                }
                for (int k = l/2; k < l; k++) {
                    std::swap(in[i*l2+j*l+k], out[(i+b)*l2+(j+b)*l+(k-b)]);
                }
            }
            for (int j = l/2; j < l; j++) {
                for (int k = 0; k < l/2; k++) {
                    std::swap(in[i*l2+j*l+k], out[(i+b)*l2+(j-b)*l+(k+b)]);
                }
                for (int k = l/2; k < l; k++) {
                    std::swap(in[i*l2+j*l+k], out[(i+b)*l2+(j-b)*l+(k-b)]);
                }
            }
        }
    }
}

}

