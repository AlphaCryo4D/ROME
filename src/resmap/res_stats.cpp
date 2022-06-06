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

#include "res_stats.h"

namespace resumap {

//const double MAXNUM = 1.79769313486231570815E308; /* 2**1024*(1-MACHEP) */
//const double PI = 3.14159265358979323846; /* pi */
//const double PIO2 = 1.57079632679489661923;   /* pi/2 */
//const double SQRT2 = 1.41421356237309504880;  /* sqrt(2) */
const double SQRTH = 7.07106781186547524401E-1;   /* sqrt(2)/2 */
//const double LOG2E = 1.4426950408889634073599;    /* 1/log(2) */
//const double SQ2OPI = 7.9788456080286535587989E-1;    /* sqrt( 2/pi ) */
//const double LOGE2 = 6.93147180559945309417E-1;   /* log(2) */
//const double LOGSQ2 = 3.46573590279972654709E-1;  /* log(2)/2 */
//const double THPIO4 = 2.35619449019234492885; /* 3*pi/4 */
//const double TWOOPI = 6.36619772367581343075535E-1;   /* 2/pi */

// !!! not used
//#define DOMAIN      1   /* argument domain error */
//#define SING        2   /* argument singularity */
//#define OVERFLOW    3   /* overflow range error */
//#define UNDERFLOW   4   /* underflow range error */
//#define TLOSS       5   /* total loss of precision */
//#define PLOSS       6   /* partial loss of precision */
//#define TOOMANY         7   /* too many iterations */
//#define MAXITER        500

#define EDOM        33
#define ERANGE      34

double MAXLOG = 7.09782712893383996732E2;   /* log(MAXNUM) */
                                                      /* double MINLOG = -7.44440071921381262314E2; *//* log(2**-1074) */
double MINLOG = -7.451332191019412076235E2; /* log(2**-1075) */

static double P[] = {
    2.46196981473530512524E-10,
    5.64189564831068821977E-1,
    7.46321056442269912687E0,
    4.86371970985681366614E1,
    1.96520832956077098242E2,
    5.26445194995477358631E2,
    9.34528527171957607540E2,
    1.02755188689515710272E3,
    5.57535335369399327526E2
};

static double Q[] = {
    /* 1.00000000000000000000E0, */
    1.32281951154744992508E1,
    8.67072140885989742329E1,
    3.54937778887819891062E2,
    9.75708501743205489753E2,
    1.82390916687909736289E3,
    2.24633760818710981792E3,
    1.65666309194161350182E3,
    5.57535340817727675546E2
};

static double R[] = {
    5.64189583547755073984E-1,
    1.27536670759978104416E0,
    5.01905042251180477414E0,
    6.16021097993053585195E0,
    7.40974269950448939160E0,
    2.97886665372100240670E0
};

static double S[] = {
    /* 1.00000000000000000000E0, */
    2.26052863220117276590E0,
    9.39603524938001434673E0,
    1.20489539808096656605E1,
    1.70814450747565897222E1,
    9.60896809063285878198E0,
    3.36907645100081516050E0
};

static double T[] = {
    9.60497373987051638749E0,
    9.00260197203842689217E1,
    2.23200534594684319226E3,
    7.00332514112805075473E3,
    5.55923013010394962768E4
};

static double U[] = {
    /* 1.00000000000000000000E0, */
    3.35617141647503099647E1,
    5.21357949780152679795E2,
    4.59432382970980127987E3,
    2.26290000613890934246E4,
    4.92673942608635921086E4
};

static inline double polevl(double x, double coef[], int N)
{
    double ans;
    int i;
    double *p;

    p = coef;
    ans = *p++;
    i = N;

    do
        ans = ans * x + *p++;
    while (--i);

    return (ans);
}

static inline double p1evl(double x, double coef[], int N)
{
    double ans;
    double *p;
    int i;

    p = coef;
    ans = x + *p++;
    i = N - 1;

    do
        ans = ans * x + *p++;
    while (--i);

    return (ans);
}

static inline double erfc(double a)
{
    double p, q, x, y, z;

    if (a < 0.0)
        x = -a;
    else
        x = a;

    if (x < 1.0)
        return (1.0 - erf(a));

    z = -a * a;

    if (z < -MAXLOG) {
under:
        //mtherr("erfc", UNDERFLOW);
        if (a < 0)
            return (2.0);
        else
            return (0.0);
    }

    z = std::exp(z);

    if (x < 8.0) {
        p = polevl(x, P, 8);
        q = p1evl(x, Q, 8);
    }
    else {
        p = polevl(x, R, 5);
        q = p1evl(x, S, 6);
    }
    y = (z * p) / q;

    if (a < 0)
        y = 2.0 - y;

    if (y == 0.0)
        goto under;

    return y;
}

static inline double erf(double x)
{
    double y, z;

    if (std::fabs(x) > 1.0)
        return (1.0 - erfc(x));
    z = x * x;

    y = x * polevl(z, T, 4) / p1evl(z, U, 5);
    return y;

}

inline double ndtr(double a)
{
    double x, y, z;

    x = a * SQRTH;
    z = std::fabs(x);

    if (z < SQRTH)
	y = 0.5 + 0.5 * erf(x);

    else {
	y = 0.5 * erfc(z);

	if (x > 0)
	    y = 1.0 - y;
    }

    return y;
}

double norm_cdf(double a) {
    return ndtr(a);
}

}

