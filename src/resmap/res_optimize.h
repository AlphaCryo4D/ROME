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

struct bracket_rt {
    double xa, xb, xc, fa, fb, fc;
    int funcalls;
};

struct minimize_scalar_rt {
    double xmin;
    double fval;
    int iter = 0;
    int funcalls = 0;
};

template<typename _F>
bracket_rt bracket(_F && func, double xa=0.0, double xb=1.0, double grow_limit=110.0, int maxiter=1000) {
    double _gold = 1.618034;
    double _verysmall_num = 1e-21;
    double fa = func(xa);
    double fb = func(xb);
    if (fa < fb) {                      // Switch so fa > fb
        std::swap(xa, xb);
        std::swap(fa, fb);
    }
    double xc = xb + _gold * (xb - xa);
    double fc = func(xc);
    int funcalls = 3;
    int iter = 0;
    double denom;
    double fw;
    double tmp1, tmp2;
    double val;
    double w, wlim;
    while (fc < fb) {
        tmp1 = (xb - xa) * (fb - fc);
        tmp2 = (xb - xc) * (fb - fa);
        val = tmp2 - tmp1;
        if (std::abs(val) < _verysmall_num) denom = 2.0 * _verysmall_num;
        else denom = 2.0 * val;
        w = xb - ((xb - xc) * tmp2 - (xb - xa) * tmp1) / denom;
        wlim = xb + grow_limit * (xc - xb);
        if (iter > maxiter) {
            throw "Too many iterations.";
        }
        iter += 1;
        if ((w - xc) * (xb - w) > 0.0) {
            fw = func(w);
            funcalls += 1;
            if (fw < fc) {
                xa = xb;
                xb = w;
                fa = fb;
                fb = fw;
                return {xa, xb, xc, fa, fb, fc, funcalls};
            }
            else if (fw > fb) {
                xc = w;
                fc = fw;
                return {xa, xb, xc, fa, fb, fc, funcalls};
            }
            w = xc + _gold * (xc - xb);
            fw = func(w);
            funcalls += 1;
        }
        else if ((w - wlim)*(wlim - xc) >= 0.0) {
            w = wlim;
            fw = func(w);
            funcalls += 1;
        }
        else if ((w - wlim)*(xc - w) > 0.0) {
            fw = func(w);
            funcalls += 1;
            if (fw < fc) {
                xb = xc;
                xc = w;
                w = xc + _gold * (xc - xb);
                fb = fc;
                fc = fw;
                fw = func(w);
                funcalls += 1;
            }
        }
        else {
            w = xc + _gold * (xc - xb);
            fw = func(w);
            funcalls += 1;
        }
        xa = xb;
        xb = xc;
        xc = w;
        fa = fb;
        fb = fc;
        fc = fw;
    }
    return {xa, xb, xc, fa, fb, fc, funcalls};
}

template<typename _RT = double, typename _V = double>
struct Brent : public minimize_scalar_rt {
    using func_t = std::function<_RT(_V)>;

    func_t func;
    double tol = 1.48e-8;
    int maxiter = 500;
    double _mintol = 1.0e-11;
    double _cg = 0.3819660;

    std::vector<double> brack;

    Brent(func_t _func, double _tol = 1.48e-18, int _maxiter = 500) 
        : func(_func), tol(_tol), maxiter(_maxiter)
    {}

    void set_bracket(std::vector<double> _brack) {
        brack = _brack;
    }

    bracket_rt get_bracket_info() {
        // BEGIN core bracket_info code 
        // carefully DOCUMENT any CHANGES in core 
        if (brack.empty()) {
            return bracket(func);
        }
        else if (brack.size() == 2) {
            return bracket(func, brack[0], brack[1]);
        }
        else if (brack.size() == 3) {
            double xa = brack[0], xb = brack[1], xc = brack[2];
            if (xa > xc) {
                std::swap(xc, xa);
            }
            if (xa >= xb || xb >= xc) {
                throw "Not a bracketing interval.";
            }
            double fa = func(xa);
            double fb = func(xb);
            double fc = func(xc);
            if (fb >= fa || fb >= fc) {
                throw "Not a bracketing interval.";
            }
            funcalls = 3;
            return {xa, xb, xc, fa, fb, fc, funcalls};
        }
        else {
            throw "Bracketing interval must be length 2 or 3 sequence.";
        }
        // END core bracket_info code 

    }

    void optimize() {
        // set up for optimization
        auto rt = get_bracket_info();
        double xa = rt.xa, xb = rt.xb, xc = rt.xc, fa = rt.fa, fb = rt.fb, fc = rt.fc;
        //////////////////////////////////
        //BEGIN CORE ALGORITHM
        //we are making NO CHANGES in this
        //////////////////////////////////
        double x = xb, w = xb, v = xb;
        double fw = func(x);
        double fv = fw;
        double fx = fw;
        double a = std::min(xa, xc), b = std::max(xa, xc);
        double deltax = 0.0;
        int funcalls = 1;
        int iter = 0;
        double rat;
        double tol1, tol2, xmid, tmp1, tmp2, p, dx_temp, u, fu;
        while (iter < maxiter) {
            tol1 = tol * std::abs(x) + _mintol;
            tol2 = 2.0 * tol1;
            xmid = 0.5 * (a + b);
            // check for convergence
            if (std::abs(x - xmid) < (tol2 - 0.5 * (b - a))) break;
            if (std::abs(deltax) <= tol1) {
                if (x >= xmid) deltax = a - x;       // do a golden section step
                else deltax = b - x;
                rat = _cg * deltax;
            }
            else {                              // do a parabolic step
                tmp1 = (x - w) * (fx - fv);
                tmp2 = (x - v) * (fx - fw);
                p = (x - v) * tmp2 - (x - w) * tmp1;
                tmp2 = 2.0 * (tmp2 - tmp1);
                if (tmp2 > 0.0) p = -p;
                tmp2 = std::abs(tmp2);
                dx_temp = deltax;
                deltax = rat;
                // check parabolic fit
                if ((p > tmp2 * (a - x)) && (p < tmp2 * (b - x)) &&
                        (std::abs(p) < std::abs(0.5 * tmp2 * dx_temp))) {
                    rat = p * 1.0 / tmp2;        // if parabolic step is useful.
                    u = x + rat;
                    if ((u - a) < tol2 || (b - u) < tol2) {
                        if (xmid - x >= 0) rat = tol1;
                        else rat = -tol1;
                    }
                }
                else {
                    if (x >= xmid) deltax = a - x;  // if it's not do a golden section step
                    else deltax = b - x;
                    rat = _cg * deltax;
                }
            }

            if (std::abs(rat) < tol1) {            // update by at least tol1
                if (rat >= 0) u = x + tol1;
                else u = x - tol1;
            }
            else {
                u = x + rat;
            }
            fu = func(u);      // calculate new output value
            funcalls += 1;

            if (fu > fx) {                 // if it's bigger than current
                if (u < x) a = u;
                else b = u;
                if ((fu <= fw) || (w == x)) {
                    v = w;
                    w = u;
                    fv = fw;
                    fw = fu;
                }
                else if ((fu <= fv) || (v == x) || (v == w)) {
                    v = u;
                    fv = fu;
                }
            }
            else {
                if (u >= x) a = x;
                else b = x;
                v = w;
                w = x;
                x = u;
                fv = fw;
                fw = fx;
                fx = fu;
            }

            iter += 1;
        }
        //////////////////////////////////////////////////////////////////
        //END CORE ALGORITHM
        //////////////////////////////////////////////////////////////////

        this->xmin = x;
        this->fval = fx;
        this->iter = iter;
        this->funcalls = funcalls;
    }

};

template<typename _RT, typename _V>
minimize_scalar_rt minimize_scalar_brent(std::function<_RT(_V)> func, 
        std::vector<double> brack, double xtol = 1.48e-8, int maxiter = 500)
{
    double tol = xtol;

    Brent<_RT, _V> brent(func, tol, maxiter);
    brent.set_bracket(brack);
    brent.optimize();
    return brent;
    
}

//template<typename _RT, typename _V>
inline minimize_scalar_rt minimize_scalar(std::function<double(double)> func, std::string method = "brent", 
        std::vector<double> brack = {}, double xtol = 1.48e-8, int maxiter = 500)
{
    if (method != "brent") {
        ERROR_REPORT("Only brent method is supported now!");
    }
    return minimize_scalar_brent(func, brack, xtol, maxiter);
}


