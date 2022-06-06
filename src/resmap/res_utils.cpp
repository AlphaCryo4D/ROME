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

#include "res_utils.h"

namespace resumap {

ruben_rt ruben(const Arrayd &weights, double c, Arrayd mult, Arrayd delta, int mode, int maxit, double eps) {
    int n;
    double dnsty;
    int ifault;
    double res;
    double tol;
    double bbeta;
    double summ;
    double sum1;
    int k;
    double hold;
    double ao;
    double z;
    int i;
    double lans, dans, pans;
    double prbty;
    double eps2;
    double aoinv;

    // Initialize parameters
    n     = weights.size();
    Arrayd gamma(weights.shape(), 0);
    Arrayd theta(weights.shape(), 1);
    std::list<double> alist, blist;

    // If no multiplicities are given, assume 1 for all
    if (mult.size() == 0) {
        mult = Arrayd(weights.shape(), 1);
    }

    // If no non-centralities are given, assume 0 for all
    if (delta.size() == 0) {
        delta = Arrayd(weights.shape(), 0);
    }

    // Basic error checking
    if ((n<1) || (c<=0) || (maxit<1) || (eps<=0.0)) {
        dnsty = 0.0;
        ifault  = 2;
        res     = -2.0;
        return {dnsty, ifault, res};
    }
    else {
        tol = -200.0;

        bbeta = weights.min();
        summ = weights.max();

        //JN_INFOV(bbeta);
        //JN_INFOV(summ);

        // Some more error checking
        if (bbeta <= 0 || summ <= 0 || (mult<1).any() || (delta<0).any()) {
            dnsty  = 0.0;
            ifault = -1;
            res    = -7.0;
            return {dnsty, ifault, res};
        }

        // Calculate BetaB
        if (mode > 0.0) bbeta = mode*bbeta;
        else bbeta = 2.0/(1.0/bbeta + 1.0/summ);

        k = 0;
        summ = 1.0;
        sum1 = 0.0;
        for (int ii = 0; ii < n; ii++) {
            hold       = bbeta/weights[ii];
            gamma[ii]   = 1.0 - hold;
            summ       = summ*(std::pow(hold,mult[ii])); //this is ok -- A.K.
            sum1       = sum1 + delta[ii];
            k          = k + mult[ii];
            // theta[ii]   = 1.0
        }

        ao = std::exp(0.5*(std::log(summ)-sum1));

        //JN_INFOV(ao);

        if (ao <= 0.0) {
            dnsty  = 0.0;
            ifault = 1;
            res    = 0.0;
            return {dnsty, ifault, res};
        }
        else { // evaluate probability and density of chi-squared on k degrees of freedom. 
            z = c/bbeta;

            //JN_INFOV(z);

            if (k % 2 == 0) {
                i    = 2;
                lans = -0.5*z;
                dans = std::exp(lans);
                pans = 1.0 - dans;
            }
            else {
                i    = 1;
                lans = -0.5*(z+std::log(z)) - std::log(std::sqrt(PI/2.0));
                dans = std::exp(lans);
                // pans = normcdf(sqrt(z),0,1) - normcdf(-1*sqrt(z),0,1)
                pans = norm_cdf(std::sqrt(z)) - norm_cdf(-1*std::sqrt(z));
            }

            //JN_INFOV(lans);
            //JN_INFOV(dans);
            //JN_INFOV(pans);

            k = k-2;
            for (int j = i; j < int(k+2); j += 2) {
                if (lans < tol) {
                    lans = lans + std::log(z/j);
                    dans = std::exp(lans);
                }
                else {
                    dans = dans*z/j;
                }
                pans = pans - dans;
            }

            //JN_INFOV(lans);
            //JN_INFOV(dans);
            //JN_INFOV(pans);

            // Evaluate successive terms of expansion
            prbty = pans;
            dnsty = dans;
            eps2  = eps/ao;
            aoinv = 1.0/ao;
            summ  = aoinv - 1.0;

            //JN_INFOV(prbty);
            //JN_INFOV(dnsty);
            //JN_INFOV(eps);
            //JN_INFOV(eps2);
            //JN_INFOV(aoinv);
            //JN_INFOV(summ);

            ifault = 4;
            for (int m = 1; m < maxit; m++) {

                //JN_INFOV(m);

                sum1 = 0.5*(theta*gamma*mult + m*delta*(theta-(theta*gamma))).sum();

                //JN_INFOA(theta);
                //JN_INFOA(gamma);
                //JN_INFOA(mult);

                //JN_INFOV(sum1);

                theta = theta*gamma;

                // b[m] = sum1
                blist.push_back(sum1);
                if (m>1) {
                    auto rt1 = blist.begin();
                    auto rt2 = alist.rbegin();
                    for (; std::next(rt1) != blist.end() && rt2 != alist.rend(); rt1++, rt2++) {
                        sum1 += *rt1 * *rt2;
                    }
                }
                //MPI_LOG << "alist "; for (auto && n : alist) MPI_LOG << n << ' '; MPI_LOG << std::endl;
                //MPI_LOG << "blist "; for (auto && n : blist) MPI_LOG << n << ' '; MPI_LOG << std::endl;
                //JN_INFOV(sum1);
                //if (m>1) sum1 = sum1 + linalg::dot(blist[:-1],alist[::-1]);

                sum1 = sum1/m;
                // a[m] = sum1
                alist.push_back(sum1);
                k    = k + 2;
                if (lans < tol) {
                    lans = lans + std::log(z/k);
                    dans = std::exp(lans);
                }
                else {
                    dans = dans*z/k;
                }

                pans  = pans - dans;
                summ  = summ - sum1;
                dnsty = dnsty + dans*sum1;
                sum1  = pans*sum1;
                prbty = prbty + sum1;
                if (prbty < -aoinv) {
                    dnsty  = 0.0;
                    ifault = 3;
                    res    = -3.0;
                    return {dnsty, ifault, res};
                }

                //JN_INFOV(pans);
                //JN_INFOV(summ);
                //JN_INFOV(dnsty);
                //JN_INFOV(sum1);
                //JN_INFOV(prbty);

                if (std::abs(pans*summ) < eps && std::abs(sum1) < eps2) {
                    ifault = 0;
                    break;
                }
            }

            dnsty  = ao*dnsty/(bbeta+bbeta);
            prbty  = ao*prbty;
            if (prbty<0.0 || prbty>1.0) {
                ifault = ifault + 5;
                res = 1e10;
            }
            else {
                if (dnsty < 0.0) ifault = ifault + 6;
                res = prbty;
            }
        }
            
        return {dnsty, ifault, res};
    }
}

}

