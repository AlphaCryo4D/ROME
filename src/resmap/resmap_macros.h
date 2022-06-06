/***************************************************************************
 *
 * Intel® Parallel Computing Center for Structural Biology
 * Principal Investigator : Youdong (Jack) Mao (Youdong_Mao@dfci.harvard.edu)
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
/***************************************************************************
 *
 * Authors:     Carlos Oscar S. Sorzano (coss@cnb.csic.es)
 *
 * Unidad de  Bioinformatica of Centro Nacional de Biotecnologia , CSIC
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
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA
 * 02111-1307  USA
 *
 *  All comments concerning this program package may be sent to the
 *  e-mail address 'xmipp@cnb.csic.es'
 ***************************************************************************/
#ifndef MACROS_H_
#define MACROS_H_

#include "./resmap_util.h"

//#define FLOAT_PRECISION
// precision
#ifdef FLOAT_PRECISION
#define FDOUBLE float
#else
#define FDOUBLE double
#endif

#ifndef MINFLOAT
#define MINFLOAT -1e30
#endif
#ifndef MAXFLOAT
#define MAXFLOAT  1e30
#endif

//#ifndef FLT_EPSILON
//#define FLT_EPSILON 1.19209e-07
//#endif

#ifndef PI
#define PI 3.14159265358979323846
#endif

/** Equal accuracy In a comparison if two values are closer than this epsilon they are said to
    be the same. Actually set to 1e-6 **/
#ifdef FLOAT_PRECISION
#define XMIPP_EQUAL_ACCURACY 1e-4
#else
#define XMIPP_EQUAL_ACCURACY 1e-6
#endif

#define ROME_EQUAL_ACCURACY XMIPP_EQUAL_ACCURACY

/** Sign of
 *
 * Valid for any kind of number (int, short, float, etc). It returns +1 or -1
 *
 * @code
 * if (SGN(x) == -1)
 *     std::cout << "x is negative" << std::endl;
 * @endcode
 */
#ifndef SGN
#define SGN(x) (((x) >= 0) ? 1 : -1)
#endif

/** Sign of, considering 0 as 0
 *
 * Valid for any kind of number (int, short, float, etc). It returns +1 if the
 * number is positive, -1 if the number is negative, and 0 if the number is 0.
 *
 * @code
 * if (SGN0(x) == -1)
 *     std::cout << "x is negative" << std::endl;
 * @endcode
 */
#ifndef SGN0
#define SGN0(x) (((x) >= 0) ? (((x) == 0) ? 0:1) : -1)
#endif


#if defined(_WIN32)

//c++ in vs2012 has no round() function
#ifndef round
#define round(x) (((x) > 0) ? (int)((x) + 0.5) : (int)((x) - 0.5))
#endif

#endif


/** Return the fractional part of a value
 *
 * The fractional part of 3.7 is 0.7 and of -3.7 is -0.7.
 */
//#define FRACTION(x) ((x) - (int)(x))

/** Clip in a saturation fashion
 *
 * CLIP is a macro which acts like a saturation curve, a value x is "clipped" to
 * a range defined by x0 and xF, for example the output values for the following
 * x and CLIP(x,-2,2) would be
 *
 * @code
 * x = ... -8 -7 -6 -5 -4 -3 -2 -1 0 1 2 3 4 5 6 7 8 ...
 * output = ... -2 -2 -2 -2 -2 -2 -2 -1 0 1 2 2 2 2 2 2 2 ...
 * @endcode
 */
#define CLIP(x, x0, xF) (((x) < (x0)) ? (x0) : (((x) > (xF)) ? (xF) : (x)))

/** Wrapping for integers
 *
 * intWRAP performs a wrapping in the integer set, when the cycle is finsihed it
 * begins again. For example, for intWRAP(x,-2,2) would be
 *
 * @code
 * x = ... -8 -7 -6 -5 -4 -3 -2 -1  0  1  2  3  4  5  6  7  8 ...
 * output = ...  2 -2 -1  0  1  2 -2 -1  0  1  2 -2 -1  0  1  2 -2 ...
 * @endcode
 */
#define intWRAP(x, x0, xF) (((x) >= (x0) && (x) <= (xF)) ? (x) : ((x) < (x0)) \
                            ? ((x) - (int)(((x) - (x0) + 1) / ((xF) - (x0) + 1) - 1) * \
                               ((xF) - (x0) + 1)) : ((x) - (int)(((x) - (xF) - 1) / ((xF) - (x0) + 1) \
                                                                 + 1) * ((xF) - (x0) + 1)))

/** Wrapping for real numbers
 *
 * realWRAP is used to keep a floating number between a range with a wrapping
 * fashion. For instance, it is used in trigonometry to say that an angle of
 * 5*PI is the same as PI, ie, to keep an angle in the range 0...2*PI
 *
 * @code
 * Corrected_angle = realWRAP(angle, 0, 2*PI);
 * @endcode
 */
#define realWRAP(x, x0, xF) (((x) >= (x0) && (x) <= (xF)) ? (x) : ((x) < (x0)) \
                             ? ((x) - (int)(((x) - (x0)) / ((xF) - (x0)) - 1) * ((xF) - (x0))) : \
                             ((x) - (int)(((x) - (xF)) / ((xF) - (x0)) + 1) * ((xF) - (x0))))

/** Degrees to radians **/
#define DEG2RAD(d) ((d) * PI / 180)

/** Radians to degrees **/
#define RAD2DEG(r) ((r) * 180 / PI)

/** ArcCosine in degrees **/
#define ACOSD(x) acos((x)) * 180. / PI

/** big-endian and little-endian conveter**/
#define SWAP32(x) ( (((x)&0x000000FF)<<24) | (((x)&0x0000FF00)<<8) | (((x)&0x00FF0000)>>8) | (((x)&0xFF000000)>>24) )


#endif
