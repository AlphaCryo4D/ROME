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

#include "resmap_util.h"		// used for building precompiled headers on Windows

#include "resmap_ctf.h"

/* Read -------------------------------------------------------------------- */
//void CTF::read(MetaDataTable &MD1, MetaDataTable &MD2, long int objectID)
//{
//    
//    if (!MD1.getValue(EMDL_CTF_VOLTAGE, kV, objectID))
//        if (!MD2.getValue(EMDL_CTF_VOLTAGE, kV, objectID))
//            kV=200;
//    
//    if (!MD1.getValue(EMDL_CTF_DEFOCUSU, DeltafU, objectID))
//        if (!MD2.getValue(EMDL_CTF_DEFOCUSU, DeltafU, objectID))
//            DeltafU=0;
//    
//    if (!MD1.getValue(EMDL_CTF_DEFOCUSV, DeltafV, objectID))
//        if (!MD2.getValue(EMDL_CTF_DEFOCUSV, DeltafV, objectID))
//            DeltafV=DeltafU;
//    
//    if (!MD1.getValue(EMDL_CTF_DEFOCUS_ANGLE, azimuthal_angle, objectID))
//        if (!MD2.getValue(EMDL_CTF_DEFOCUS_ANGLE, azimuthal_angle, objectID))
//            azimuthal_angle=0;
//    
//    if (!MD1.getValue(EMDL_CTF_CS, Cs, objectID))
//        if (!MD2.getValue(EMDL_CTF_CS, Cs, objectID))
//            Cs=0;
//    
//    if (!MD1.getValue(EMDL_CTF_BFACTOR, Bfac, objectID))
//        if (!MD2.getValue(EMDL_CTF_BFACTOR, Bfac, objectID))
//            Bfac=0;
//    
//    if (!MD1.getValue(EMDL_CTF_SCALEFACTOR, scale, objectID))
//        if (!MD2.getValue(EMDL_CTF_SCALEFACTOR, scale, objectID))
//            scale=1;
//    
//    if (!MD1.getValue(EMDL_CTF_Q0, Q0, objectID))
//        if (!MD2.getValue(EMDL_CTF_Q0, Q0, objectID))
//            Q0=0;
//    
//    initialise();
//}
void CTF::setValues(double _defU, double _defV, double _defAng, double _voltage,
                    double _Cs, double _Q0, double _Bfac, double _scale)
{
    kV              = _voltage;
    DeltafU         = _defU;
    DeltafV         = _defV;
    azimuthal_angle = _defAng;
    Cs              = _Cs;
    Bfac            = _Bfac;
    scale           = _scale;
    Q0              = _Q0;
    
    initialise();
}
/* Read from 1 MetaDataTable ----------------------------------------------- */
//void CTF::read(MetaDataTable &MD)
//{
//    MetaDataTable MDempty;
//    MDempty.addObject(); // add one empty object
//    read(MD, MDempty);
//}
//
///* Write ------------------------------------------------------------------- */
//void CTF::write(std::ostream &out)
//{
//    MetaDataTable MD;
//    MD.addObject();
//    MD.setValue(EMDL_CTF_VOLTAGE, kV);
//    MD.setValue(EMDL_CTF_DEFOCUSU, DeltafU);
//    MD.setValue(EMDL_CTF_DEFOCUSV, DeltafV);
//    MD.setValue(EMDL_CTF_DEFOCUS_ANGLE, azimuthal_angle);
//    MD.setValue(EMDL_CTF_CS, Cs);
//    MD.setValue(EMDL_CTF_BFACTOR, Bfac);
//    MD.setValue(EMDL_CTF_SCALEFACTOR, scale);
//    MD.setValue(EMDL_CTF_Q0, Q0);
//    MD.write(out);
//}

/* Default values ---------------------------------------------------------- */
void CTF::clear()
{
    kV = 200;
    DeltafU = DeltafV = azimuthal_angle = 0;
    Cs = Bfac = 0;
    Q0 = 0;
    scale = 1;
}

/* Initialise the CTF ------------------------------------------------------ */
void CTF::initialise()
{
    
    // Change units
    double local_Cs = Cs * 1e7;
    double local_kV = kV * 1e3;
    rad_azimuth = DEG2RAD(azimuthal_angle);
    
    // Average focus and deviation
    defocus_average   = -(DeltafU + DeltafV) * 0.5;
    defocus_deviation = -(DeltafU - DeltafV) * 0.5;
    
    // lambda=h/sqrt(2*m*e*kV)
    //    h: Planck constant
    //    m: electron mass
    //    e: electron charge
    // lambda=0.387832/sqrt(kV*(1.+0.000978466*kV)); // Hewz: Angstroms
    // lambda=h/sqrt(2*m*e*kV)
    lambda=12.2643247 / sqrt(local_kV * (1. + local_kV * 0.978466e-6)); // See http://en.wikipedia.org/wiki/Electron_diffraction
    
    // Helpful constants
    // ICE: X(u)=-PI/2*deltaf(u)*lambda*u^2+PI/2*Cs*lambda^3*u^4
    //          = K1*deltaf(u)*u^2         +K2*u^4
    K1 = PI / 2 * 2 * lambda;
    K2 = PI / 2 * local_Cs * lambda * lambda * lambda;
    K3 = sqrt(1-Q0*Q0);
    K4 = -Bfac / 4.;
    
    if (Q0 < 0. || Q0 > 1.)
        ERROR_REPORT("CTF::initialise ERROR: AmplitudeContrast Q0 cannot be smaller than zero or larger than one!");
    
    if (fabs(DeltafU) < 1e-6 && fabs(DeltafV) < 1e-6 && fabs(Q0) < 1e-6 && fabs(Cs) < 1e-6)
        ERROR_REPORT("CTF::initialise: ERROR: CTF initialises to all-zero values. Was a correct STAR file provided?");
    
}

/* Generate a complete CTF Image ------------------------------------------------------ */
void CTF::getFftwImage(double* result, int dim,int orixdim, int oriydim, double angpix,
                       bool do_abs, bool do_only_flip_phases, bool do_intact_until_first_peak, bool do_damping)
{
    
    double xs = (double)orixdim * angpix;
    double ys = (double)oriydim * angpix;
    int hdim = dim/2+1;

    for (int i = 0; i < dim; i++){
        for (int j = 0; j < hdim; j++)
        {
            long int ip = (i < hdim) ? i : i - dim;
            long int jp = j;
            double x = (double)jp / xs;
            double y = (double)ip / ys;
            result[i*hdim+j] = getCTF(x, y, do_abs, do_only_flip_phases, do_intact_until_first_peak, do_damping);
        }
    }
}

void CTF::getFftwImage(float* result, int dim,int orixdim, int oriydim, double angpix,
                       bool do_abs, bool do_only_flip_phases, bool do_intact_until_first_peak, bool do_damping)
{
    
    double xs = (double)orixdim * angpix;
    double ys = (double)oriydim * angpix;
    int hdim = dim/2+1;
    
    for (int i = 0; i < dim; i++){
        for (int j = 0; j < hdim; j++)
        {
            int ip = (i < hdim) ? i : i - dim;
            int jp = j;
            double x = (double)jp / xs;
            double y = (double)ip / ys;
            result[i*hdim+j] = getCTF(x, y, do_abs, do_only_flip_phases, do_intact_until_first_peak, do_damping);
        }
    }
}