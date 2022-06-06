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

#pragma once

/* Clang/LLVM. ---------------------------------------------- */
#if defined(__clang__)
#  define JN_CLANG

/* Intel ICC/ICPC. ------------------------------------------ */
#elif defined(__ICC) || defined(__INTEL_COMPILER)
#  define JN_ICC

/* GNU GCC/G++. --------------------------------------------- */
#elif defined(__GNUC__) || defined(__GNUG__)
#  define JN_GCC

/* Hewlett-Packard C/aC++. ---------------------------------- */
#elif defined(__HP_cc) || defined(__HP_aCC)
#  define JN_HPCC

/* IBM XL C/C++. -------------------------------------------- */
#elif defined(__IBMC__) || defined(__IBMCPP__)
#  define JN_IBMC

/* Microsoft Visual Studio. --------------------------------- */
#elif defined(_MSC_VER)
#  define JN_MSC

/* Portland Group PGCC/PGCPP. ------------------------------- */
#elif defined(__PGI)

/* Oracle Solaris Studio. ----------------------------------- */
#elif defined(__SUNPRO_C) || defined(__SUNPRO_CC)

#endif


