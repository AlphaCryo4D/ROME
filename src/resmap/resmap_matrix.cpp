/***************************************************************************
 *
 * Authors: "Yongbei(Glow) Ma,Jiayi (Timmy) Wu, Youdong (Jack) Mao"
 * Dana-Farber Cancer Institute, Harvard Medical School and Peking University
 * Authors: "Brett, Bevin"
 * Intel Corporation
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

#include "resmap_util.h"		// used for building precompiled headers on Windows

#include "resmap_matrix.h"

/***************************************************************************
 *
 * Author: "Sjors H.W. Scheres"
 * MRC Laboratory of Molecular Biology
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

/* from transformations.cpp ---------------------------------------------------- */
/* Rotation 3D around the system axes -------------------------------------- */
void rotation3DMatrix(FDOUBLE ang, char axis, Matrix2D< FDOUBLE > &result,
                      bool homogeneous/* = true*/)
{
    if (homogeneous)
    {
        result.initZeros(4,4);
        result(3,3) = 1;
    }
    else
        result.initZeros(3,3);
    
    FDOUBLE cosine, sine;
    ang = DEG2RAD(ang);
    cosine = cos(ang);
    sine = sin(ang);
    
    switch (axis)
    {
        case 'Z':
            result(0, 0) = cosine;
            result(0, 1) = -sine;
            result(1, 0) = sine;
            result(1, 1) = cosine;
            result(2, 2) = 1;
            break;
        case 'Y':
            result(0, 0) = cosine;
            result(0, 2) = -sine;
            result(2, 0) = sine;
            result(2, 2) = cosine;
            result(1, 1) = 1;
            break;
        case 'X':
            result(1, 1) = cosine;
            result(1, 2) = -sine;
            result(2, 1) = sine;
            result(2, 2) = cosine;
            result(0, 0) = 1;
            break;
        default:
            ERROR_REPORT("rotation3DMatrix: Unknown axis");
    }
}
/* Align a vector with Z axis */
void alignWithZ(const Matrix1D< FDOUBLE > &axis, Matrix2D< FDOUBLE >& result,
                bool homogeneous/* = true*/)
{
    if (axis.size() != 3)
        ERROR_REPORT("alignWithZ: Axis is not in R3");
    if (homogeneous)
    {
        result.initZeros(4,4);
        result(3, 3) = 1;
    }
    else
        result.initZeros(3,3);
    Matrix1D<FDOUBLE>  Axis(axis);
    Axis.selfNormalize();
    
    // Compute length of the projection on YZ plane
    FDOUBLE proj_mod = sqrt(YY(Axis) * YY(Axis) + ZZ(Axis) * ZZ(Axis));
    if (proj_mod > XMIPP_EQUAL_ACCURACY)
    {   // proj_mod!=0
        // Build Matrix result, which makes the turning axis coincident with Z
        result(0, 0) = proj_mod;
        result(0, 1) = -XX(Axis) * YY(Axis) / proj_mod;
        result(0, 2) = -XX(Axis) * ZZ(Axis) / proj_mod;
        result(1, 0) = 0;
        result(1, 1) = ZZ(Axis) / proj_mod;
        result(1, 2) = -YY(Axis) / proj_mod;
        result(2, 0) = XX(Axis);
        result(2, 1) = YY(Axis);
        result(2, 2) = ZZ(Axis);
    }
    else
    {
        // I know that the Axis is the X axis
        result(0, 0) = 0;
        result(0, 1) = 0;
        result(0, 2) = -1;
        result(1, 0) = 0;
        result(1, 1) = 1;
        result(1, 2) = 0;
        result(2, 0) = 1;
        result(2, 1) = 0;
        result(2, 2) = 0;
    }
}
/* Rotation 3D around any axis -------------------------------------------- */
void rotation3DMatrix(FDOUBLE ang, const Matrix1D<FDOUBLE> &axis,
                      Matrix2D<FDOUBLE> &result, bool homogeneous/* = true*/)
{
    // Compute a matrix which makes the turning axis coincident with Z
    // And turn around this axis
    Matrix2D<FDOUBLE> A,R;
    alignWithZ(axis,A,homogeneous);
    rotation3DMatrix(ang, 'Z', R, homogeneous);
    result=A.transpose() * R * A;
}