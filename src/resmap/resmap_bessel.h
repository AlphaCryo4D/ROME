/***************************************************************************
 *
 * Authors: "Yongbei(Glow) Ma,Jiayi (Timmy) Wu, Youdong (Jack) Mao"
 * Dana-Farber Cancer Institute, Harvard Medical School and Peking University
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

#ifndef BESSEL_H_
#define BESSEL_H_

#include <cmath>
#include <iostream>
#include <algorithm>

#include "resmap_error.h"
#include "resmap_macros.h"

double kaiser_Fourier_value(double w, double a, double alpha, int m);

// Bessel functions --------------------------------------------------------
double bessj0(double x);
double bessj3_5(double x);
double bessj1_5(double x);

double bessi0(double x);
double bessi1(double x);
double bessi0_5(double x);
double bessi1_5(double x);
double bessi2(double x);
double bessi3(double x);
double bessi2_5(double x);
double bessi3_5(double x);
double bessi4(double x);

#endif
