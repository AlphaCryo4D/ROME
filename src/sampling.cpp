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
#include "util.h"		// used for building precompiled headers on Windows

#include "sampling.h"

HealpixSampler::HealpixSampler()
{
    is_3D = false;
    healpix_order = -1;
    fn_sym = "C1";
    pgGroup = pgOrder = 0;
    random_perturbation = 0;
    perturbation_factor = 0.5;
}

HealpixSampler::~HealpixSampler()
{
    translations_x  .resize(0);
    translations_y  .resize(0);
    psi_angles      .resize(0);
    rot_angles      .resize(0);
    tilt_angles     .resize(0);
    healpix_order   = -1;
    single_translations		.resize(0);
    single_translations_x_index	.resize(0);
    single_translations_y_index	.resize(0);
}


void HealpixSampler::initialize(FDOUBLE _offset_step /*= -1.*/,FDOUBLE _offset_range /*= -1.*/,FDOUBLE _psi_step /*= -1*/,
                                int _order /*= -1*/,std::string _fn_sym /*= "C1"*/,int prior_mode /*= NOPRIOR*/)
{
    orientational_prior_mode = prior_mode;
    
    // --------- initialize offset,rotate,healpix step -----------
    if (_offset_step > 0. && _offset_range >= 0.)
    {
        offset_step = _offset_step;
        offset_range = _offset_range;
    }
    else{
        offset_step = 2;
        offset_range = 10;
    }
    
    if(_order > 0){
        healpix_order = _order;
        assert(healpix_order < 5);
        is_3D = true;
    }
    
    if(_psi_step == -1){
        if(is_3D)
            psi_step = 360. / (6 * round(std::pow(2., healpix_order)));
        else
            psi_step = 10;
    }
    else
        psi_step = _psi_step;
    
    // ----------- initialize rot_angles,tilt_angles  --------------
    if (is_3D) {
#ifdef SAMPLING3D
        healpix_base.Set(healpix_order, NEST);
        fn_sym = _fn_sym;
        
        // symmetry
        // Set up symmetry
        SymList SL;
        SL.isSymmetryGroup(fn_sym, pgGroup, pgOrder);
        SL.read_sym_file(fn_sym);
        
        // Precalculate (3x3) symmetry matrices
        Matrix2D<FDOUBLE>  L(4, 4), R(4, 4);
        Matrix2D<FDOUBLE>  Identity(3,3);
        Identity.initIdentity();
        R_repository.clear();
        L_repository.clear();
        R_repository.push_back(Identity);
        L_repository.push_back(Identity);
        for (int isym = 0; isym < SL.SymsNo(); isym++)
        {
            SL.get_matrices(isym, L, R);
            R.resize(3, 3);
            L.resize(3, 3);
            R_repository.push_back(R);
            L_repository.push_back(L);
        }
#endif
    }
    else{
        assert(_fn_sym=="C1");
        fn_sym = "C1"; // This may not be set yet if restarting a 2D run....
    }
    
    setTranslations();
    
    setOrientations();
    
    resetRandomlyPerturbedSampling();
}

void HealpixSampler::setTranslations(FDOUBLE _offset_step, FDOUBLE _offset_range)
{
    if (_offset_step > 0. && _offset_range >= 0.)
    {
        offset_step = _offset_step;
        offset_range = _offset_range;
    }
    else
    {
        assert(offset_step>0);
        assert(offset_range>=0);
    }
    
    int maxr = ceil(offset_range / offset_step);
    
    translations_x  .resize((2*maxr+1)*(2*maxr+1));
    translations_y  .resize((2*maxr+1)*(2*maxr+1));
    int trans_ind_x = 0;
    nrTrans = 0;
    for (int ix = -maxr; ix <= maxr; ix++,trans_ind_x++)
    {
        FDOUBLE xoff = ix * offset_step;
        single_translations.push_back(xoff);
        int trans_ind_y = 0;
        for (int iy = -maxr; iy <= maxr; iy++,trans_ind_y++)
        {
            FDOUBLE yoff = iy * offset_step;
            if (xoff*xoff + yoff*yoff <= offset_range * offset_range){
                translations_x[nrTrans] = xoff;
                translations_y[nrTrans] = yoff;
                nrTrans++;
                //
                single_translations_x_index.push_back(trans_ind_x);
                single_translations_y_index.push_back(trans_ind_y);
            }
        }
    }
    //
    // #define DEBUG_SAMPLING
#ifdef DEBUG_SAMPLING
    std::cout<<std::setw(10)<<"itrans"<<std::setw(15)<<"translations_x"<<std::setw(15)<<"translations_y"<<std::endl;
    for(int itrans = 0;itrans < nrTrans;itrans++){
        std::cout<<std::setw(10)<<itrans<<std::setw(15)<<translations_x[itrans]<<std::setw(15)<<translations_y[itrans]<<std::endl;
        double trans_x = single_translations[single_trans_x_index[itrans]];
        double trans_y = single_translations[single_trans_y_index[itrans]];
        if(trans_x!=translations_x[itrans] || trans_y!=translations_y[itrans]){
            std::cout<<"!!!!"<<std::setw(10)<<itrans<<std::setw(15)<<trans_x<<std::setw(15)<<trans_y<<" "<<std::endl;
            ERROR_REPORT("diff..sampling");
        }
    }
#endif
}

void HealpixSampler::setOrientations(int _order, FDOUBLE _psi_step)
{
    // 2D in-plane angles
    // By default in 3D case: use more-or-less same psi-sampling as the 3D healpix object
    // By default in 2D case: use 5 degree
    if (_psi_step > 0.)
        psi_step = _psi_step;
    
    // ---------- initialize psi_angles   -----------------
    nrPsi = 0;
    int nr_psi = ceil(360./psi_step);
    psi_step = 360./(double)nr_psi;
    
    psi_angles  .resize(nr_psi);
    for (int ipsi = 0; ipsi < nr_psi; ipsi++)
    {
        psi_angles[nrPsi] = ipsi * psi_step;
        nrPsi++;
    }
    
    // Setup the HealPix object
    // For adaptive oversampling only precalculate the COARSE sampling!
    if (_order >= 0)
    {
        assert(_order<=13);// Npix = 805306368,up_bound of integer is 2147483648...angle is 25".8
        healpix_base.Set(_order, NEST);
        healpix_order = _order;
    }
    
    // ----------- initialize rot_angles,tilt_angles  --------------
    if (is_3D) {
#ifdef SAMPLING3D
        //
        nrPix = healpix_base.Npix();
        double zz, phi;
        
        rot_angles  	.resize(nrPix);
        tilt_angles 	.resize(nrPix);
        directions_ipix	.resize(nrPix);
        for (int ipix = 0; ipix < nrPix; ipix++)
        {
            healpix_base.pix2ang_z_phi(ipix, zz, phi);
            rot_angles[ipix] = rad2deg(phi);
            tilt_angles[ipix] = acosd(zz);
            checkDirection(rot_angles[ipix], tilt_angles[ipix]);
            directions_ipix[ipix] = ipix;
        }
        
//#define DEBUG_SAMPLING
#ifdef  DEBUG_SAMPLING
        writeAllOrientationsToBild("./orients_all.bild", "1 0 0 ", 0.020);
#endif
        // Now remove symmetry-related pixels
        // TODO check size of healpix_base.max_pixrad
        removeSymmetryEquivalentPoints(0.5 * RAD2DEG(healpix_base.max_pixrad()));
        nrPix = directions_ipix.size();
        
#ifdef  DEBUG_SAMPLING
        writeAllOrientationsToBild("./orients_sym.bild", "0 1 0 ", 0.021);
#endif
        
#endif
    }
    else{
        nrPix = 1;
    }
}

void HealpixSampler::resetRandomlyPerturbedSampling()
{
    // Actual instance of random perturbation
    // Add to the random perturbation from the last iteration, so it keeps changing strongly...
    random_perturbation += dontShare_Random_generator.rnd_unif(0.5*perturbation_factor, perturbation_factor);
    random_perturbation = realWrap(random_perturbation, -perturbation_factor, perturbation_factor);
}

size_t HealpixSampler::NrDir(int oversampling_order /*= 0*/) const
{
    if (oversampling_order == 0)
        return nrPix;
    else
        return round(std::pow(2., oversampling_order * 2))*nrPix;
}

size_t HealpixSampler::NrPsi(int oversampling_order /*= 0*/) const
{
    if (oversampling_order == 0)
        return nrPsi;
    else
        return round(pow(2., oversampling_order)) * nrPsi;
}

size_t HealpixSampler::NrTrans(int oversampling_order /*= 0*/) const
{
    if (oversampling_order == 0)
        return nrTrans;// translations_x.size();
    else
        return round(pow(2., oversampling_order * 2)) * nrTrans;//translations_x.size();   
}

size_t HealpixSampler::NrPoints(int oversampling_order) const
{
    if (is_3D)
        return NrPsi(oversampling_order) * NrTrans(oversampling_order) * NrDir(oversampling_order);
    else
        return NrPsi(oversampling_order) * NrTrans(oversampling_order);
}

double HealpixSampler::getTranslationalSampling(int adaptive_oversampling) const
{
    return offset_step / std::pow(2., adaptive_oversampling);
}

double HealpixSampler::getAngularSampling(int adaptive_oversampling /*= 0*/) const
{
    if (is_3D) {
        int order =  healpix_order + adaptive_oversampling;
        return 360. / (6 * round(std::pow(2., order)));
    }
    else
        return psi_step / pow(2., adaptive_oversampling);
}

size_t HealpixSampler::oversamplingFactorOrientations(int oversampling_order) const
{
	assert(oversampling_order < 4);
    if (is_3D)
        return round(std::pow(2., oversampling_order * 3));
    else
        return round(std::pow(2., oversampling_order));
}

int HealpixSampler::oversamplingFactorTranslations(int oversampling_order) const
{
	assert(oversampling_order < 4);
	return round(pow(2., oversampling_order * 2));
}

void HealpixSampler::getTranslations(int itrans,int oversampling_order,FDOUBLE* over_trans_x,FDOUBLE* over_trans_y) const
{
    assert(oversampling_order < 4);
    
    if (oversampling_order == 0)
    {
        over_trans_x[0] = translations_x[itrans];
        over_trans_y[0] = translations_y[itrans];
    }
    else
    {
        int nr_oversamples = round(std::pow(2., oversampling_order));
        for (int iover_trans_y = 0; iover_trans_y < nr_oversamples; iover_trans_y++)
        {
            for (int iover_trans_x = 0; iover_trans_x < nr_oversamples; iover_trans_x++)
            {
                FDOUBLE over_yoff = translations_y[itrans] - 0.5 * offset_step + (0.5 + iover_trans_y) * offset_step / nr_oversamples;
                FDOUBLE over_xoff = translations_x[itrans] - 0.5 * offset_step + (0.5 + iover_trans_x) * offset_step / nr_oversamples;
                
                int iover_trans = iover_trans_y*nr_oversamples+iover_trans_x;
                
                over_trans_x[iover_trans] = over_xoff;
                over_trans_y[iover_trans] = over_yoff;
            }
        }
    }
    
    if (fabs(random_perturbation) > 0.)
    {
        double myperturb = random_perturbation * offset_step;
        int nr_over_trans = oversamplingFactorTranslations(oversampling_order);
        for (int iover_trans = 0; iover_trans < nr_over_trans; iover_trans++)
        {
            over_trans_x[iover_trans] += myperturb;
            over_trans_y[iover_trans] += myperturb;
        }
    }
}

void HealpixSampler::getAllTranslationsAndOverTrans(int oversampling_order,Vector1d& trans_x,Vector1d& trans_y,
                                                    Vector1d& trans_x_over,Vector1d& trans_y_over) const
{
    assert(oversampling_order < 4);
    
    for (int itrans = 0; itrans < nrTrans; itrans++)
    {
        trans_x[itrans] = translations_x[itrans];
        trans_y[itrans] = translations_y[itrans];
    }
    //
    if (oversampling_order > 0)
    {
        int nr_oversamples = round(std::pow(2., oversampling_order));
        for (int iover_trans_y = 0; iover_trans_y < nr_oversamples; iover_trans_y++)
        {
            for (int iover_trans_x = 0; iover_trans_x < nr_oversamples; iover_trans_x++)
            {
                double over_yoff = - 0.5 * offset_step + (0.5 + iover_trans_y) * offset_step / nr_oversamples;
                double over_xoff = - 0.5 * offset_step + (0.5 + iover_trans_x) * offset_step / nr_oversamples;
                
                int iover_trans = iover_trans_y*nr_oversamples+iover_trans_x;
                
                trans_x_over[iover_trans] = over_xoff;
                trans_y_over[iover_trans] = over_yoff;
            }
        }
    }
    else
    {
        trans_x_over[0] = 0;
        trans_y_over[0] = 0;
    }
    //
    if (fabs(random_perturbation) > 0.)
    {
        double myperturb = random_perturbation * offset_step;
#ifdef TTTT // where to add perturbation
        for (int itrans = 0; itrans < nrTrans; itrans++) {
            trans_x[itrans] += myperturb;
            trans_y[itrans] += myperturb;
        }
#else
        if (oversampling_order > 0){
            int nr_oversamples = round(std::pow(2., oversampling_order));
            int exp_nr_over_trans = nr_oversamples*nr_oversamples;
            for (int iover_trans = 0; iover_trans < exp_nr_over_trans; iover_trans++) {
                trans_x_over[iover_trans] += myperturb;
                trans_y_over[iover_trans] += myperturb;
            }
        }
        else{
            trans_x_over[0] += myperturb;
            trans_y_over[0] += myperturb;
        }
#endif
    }
}

void HealpixSampler::getAllSingleTranslationsAndOverTrans(int oversampling_order,std::map<double,int>& exp_positive_shift_index,
                                                          std::vector<XYShift>& exp_trans_xyshift,Vector1d& trans_x_over,Vector1d& trans_y_over) const
{
    assert(oversampling_order < 4);
    exp_positive_shift_index.clear();
    int index = 0;
    for (int iSingleTrans = 0; iSingleTrans < single_translations.size(); iSingleTrans++)
    {
        double shift = single_translations[iSingleTrans];
        if (shift > 0) { exp_positive_shift_index.insert(std::make_pair(shift, index));index++;}
    }
    exp_trans_xyshift.resize(0);
    for (int itrans = 0; itrans < nrTrans; itrans++) {
        int index_x = single_translations_x_index[itrans];
        int index_y = single_translations_y_index[itrans];
        double shiftx = single_translations[index_x];
        double shifty = single_translations[index_y];
        exp_trans_xyshift.push_back(
			std::make_pair(
				Shift(fabs(shiftx),shiftx>0?1:(shiftx!=0?-1:0)), 
				Shift(fabs(shifty),shifty>0?1:(shifty!=0?-1:0)) ));
    }
    //
    if (oversampling_order > 0)
    {
        int nr_oversamples = round(std::pow(2., oversampling_order));
        for (int iover_trans_y = 0; iover_trans_y < nr_oversamples; iover_trans_y++)
        {
            for (int iover_trans_x = 0; iover_trans_x < nr_oversamples; iover_trans_x++)
            {
                double over_yoff = - 0.5 * offset_step + (0.5 + iover_trans_y) * offset_step / nr_oversamples;
                double over_xoff = - 0.5 * offset_step + (0.5 + iover_trans_x) * offset_step / nr_oversamples;
                
                int iover_trans = iover_trans_y*nr_oversamples+iover_trans_x;
                
                trans_x_over[iover_trans] = over_xoff;
                trans_y_over[iover_trans] = over_yoff;
            }
        }
    }
    else
    {
        trans_x_over[0] = 0;
        trans_y_over[0] = 0;
    }
    if (fabs(random_perturbation) > 0.)
    {
        double myperturb = random_perturbation * offset_step;
#ifdef TTTT // where to add perturbation
        for (int itrans = 0; itrans < nrTrans; itrans++) {
            trans_x[itrans] += myperturb;
            trans_y[itrans] += myperturb;
        }
#else
        if (oversampling_order > 0){
            int nr_oversamples = round(std::pow(2., oversampling_order));
            int exp_nr_over_trans = nr_oversamples*nr_oversamples;
            for (int iover_trans = 0; iover_trans < exp_nr_over_trans; iover_trans++) {
                trans_x_over[iover_trans] += myperturb;
                trans_y_over[iover_trans] += myperturb;
            }
        }
        else{
            trans_x_over[0] += myperturb;
            trans_y_over[0] += myperturb;
        }
#endif
    }
}

void HealpixSampler::getOrientations2D(int ipsi,int oversampling_order,FDOUBLE* over_psi) const
{
    assert(oversampling_order < 4);

    if (oversampling_order == 0)
    {
        over_psi[0] = psi_angles[ipsi];
    }
    else
    {
        // for 2D sampling, only push back oversampled psi rotations
        int nr_ipsi_over = round(pow(2., oversampling_order));
        for (int ipsi_over = 0; ipsi_over < nr_ipsi_over; ipsi_over++)
        {
            double overpsi = psi_angles[ipsi] - 0.5 * psi_step + (0.5 + ipsi_over) * psi_step / nr_ipsi_over;
            over_psi[ipsi_over] = overpsi;
        }
    }
    
    // Random perturbation
    if (fabs(random_perturbation) > 0.)
    {
        double myperturb = random_perturbation * getAngularSampling();
        
        int nr_ipsi_over = round(pow(2., oversampling_order));
        for (int ipsi_over = 0; ipsi_over < nr_ipsi_over; ipsi_over++)
            over_psi[ipsi_over] += myperturb;
    }
}

void HealpixSampler::addPerturbation(int oversampling_order,FDOUBLE* over_rot,FDOUBLE* over_tilt,FDOUBLE* over_psi) const
{
    // Random perturbation
    if (fabs(random_perturbation) > 0.)
    {
        FDOUBLE A[3][3],R[3][3],AR[3][3];
        double myperturb = random_perturbation * getAngularSampling();
        int nr_over = oversamplingFactorOrientations(oversampling_order);
        for (int iover = 0; iover < nr_over; iover++) {
            Euler_angles2matrix(over_rot[iover],over_tilt[iover],over_psi[iover],A);
            Euler_angles2matrix(myperturb, myperturb, myperturb, R);
#define A_MULTIPLY_R(i,j) AR[i][j] = A[i][0]*R[0][j]+A[i][1]*R[1][j]+A[i][2]*R[2][j];
            A_MULTIPLY_R(0,0) A_MULTIPLY_R(0,1) A_MULTIPLY_R(0,2)
            A_MULTIPLY_R(1,0) A_MULTIPLY_R(1,1) A_MULTIPLY_R(1,2)
            A_MULTIPLY_R(2,0) A_MULTIPLY_R(2,1) A_MULTIPLY_R(2,2)
            Euler_matrix2angles(AR,over_rot[iover],over_tilt[iover],over_psi[iover]);
        }// end loop iover
    }
}

void HealpixSampler::getOrientations3D(int idir,int ipsi,int oversampling_order,FDOUBLE* over_psi,
                                       FDOUBLE* over_rot,FDOUBLE* over_tilt) const
{
    assert(oversampling_order < 4);
	//
#ifdef SAMPLING3D
    if (oversampling_order == 0)
    {
        over_psi[0] = psi_angles[ipsi];
        over_tilt[0] = tilt_angles[idir];
        over_rot[0] = rot_angles[idir];
    }
    else{
        Healpix_Base HealPixOver(oversampling_order + healpix_order, NEST);
        int fact = HealPixOver.Nside()/healpix_base.Nside();
        int nr_ipsi_over = round(std::pow(2., oversampling_order));
        int nr_idir_over = round(std::pow(2., oversampling_order*2));
        int x, y, face;
        FDOUBLE rot, tilt;
        int ipix = directions_ipix[idir];
        int idir_over = 0;
        // Get x, y and face for the original, coarse grid
        healpix_base.nest2xyf(ipix, x, y, face);
        // Loop over the oversampled Healpix pixels on the fine grid
        for (int j = fact * y; j < fact * (y+1); ++j)
        {
            for (int i = fact * x; i < fact * (x+1); ++i)
            {
                int overpix = HealPixOver.xyf2nest(i, j, face);
                // this one always has to be double (also for SINGLE_PRECISION CALCULATIONS) for call to external library
                double zz, phi;
                HealPixOver.pix2ang_z_phi(overpix, zz, phi);
                rot = rad2deg(phi);
                tilt = acosd(zz);
                
                // The geometrical considerations about the symmetry below require that rot = [-180,180] and tilt [0,180]
                checkDirection(rot, tilt);
                
                for (int ipsi_over = 0; ipsi_over < nr_ipsi_over; ipsi_over++)
                {
                    double overpsi = psi_angles[ipsi] - 0.5 * psi_step + (0.5 + ipsi_over) * psi_step / nr_ipsi_over;
                    int iover = idir_over*nr_ipsi_over+ipsi_over;
                    assert(iover >= 0);
                    over_rot[iover] = rot;
                    over_tilt[iover] = tilt;
                    over_psi[iover] = overpsi;
                }// end loop ipsi_over
                idir_over++;
            }// end loop i
        } // end loop j
        assert(idir_over == nr_idir_over);
    }// oversampling_order != 0
#endif
    //
    addPerturbation(oversampling_order, over_rot, over_tilt, over_psi);
}

// 2D
void HealpixSampler::getAllOrientations(int oversampling_order,Vector2d& my_psi) const
{
    assert(oversampling_order < 4);
    // TODO
    int nr_over = oversamplingFactorOrientations(oversampling_order);
    for (int ipsi = 0; ipsi < nrPsi; ipsi++) {
        getOrientations2D(ipsi,oversampling_order,my_psi[ipsi].data());
    }
}

// 3D
void HealpixSampler::getAllOrientations(int oversampling_order,Vector2d& my_psi,Vector2d& my_rot,Vector2d& my_tilt) const
{
    assert(oversampling_order < 4);
    int nr_over = oversamplingFactorOrientations(oversampling_order);
#ifdef SAMPLING3D
    if (oversampling_order == 0)
    {
        for (int idir = 0; idir < nrPix; idir++) {
            my_tilt[idir][0] = tilt_angles[idir];
            my_rot[idir][0] = rot_angles[idir];
        }
        for (int ipsi = 0; ipsi < nrPsi; ipsi++) {
            my_psi[ipsi][0] = psi_angles[ipsi];
        }
    }
    else
    {
        Healpix_Base HealPixOver(oversampling_order + healpix_order, NEST);
        int fact = HealPixOver.Nside()/healpix_base.Nside();
        int nr_idir_over = round(std::pow(2., oversampling_order*2));
        int x, y, face;
        FDOUBLE rot, tilt;
        for (int idir = 0; idir < nrPix; idir++)
        {
            auto over_rot = my_rot[idir].data();
            auto over_tilt = my_tilt[idir].data();
            int ipix = directions_ipix[idir];
            int idir_over = 0;
            // Get x, y and face for the original, coarse grid
            healpix_base.nest2xyf(ipix, x, y, face);
            // Loop over the oversampled Healpix pixels on the fine grid
            for (int j = fact * y; j < fact * (y+1); ++j)
            {
                for (int i = fact * x; i < fact * (x+1); ++i)
                {
                    int overpix = HealPixOver.xyf2nest(i, j, face);
                    // this one always has to be double (also for SINGLE_PRECISION CALCULATIONS) for call to external library
                    double zz, phi;
                    HealPixOver.pix2ang_z_phi(overpix, zz, phi);
                    rot = rad2deg(phi);
                    tilt = acosd(zz);
                    
                    // The geometrical considerations about the symmetry below require that rot = [-180,180] and tilt [0,180]
                    checkDirection(rot, tilt);
                    //
                    over_rot[idir_over] = rot;
                    over_tilt[idir_over] = tilt;
                    //
                    idir_over++;
                }// end loop i
            } // end loop j
            assert(idir_over == nr_idir_over);
        }
        //
        int nr_ipsi_over = round(std::pow(2., oversampling_order));
        assert(nr_over==nr_idir_over*nr_ipsi_over);
        for (int ipsi = 0; ipsi < nrPsi; ipsi++)
        {
            auto over_psi = my_psi[ipsi].data();
            for (int ipsi_over = 0; ipsi_over < nr_ipsi_over; ipsi_over++)
            {
                double overpsi = psi_angles[ipsi] - 0.5 * psi_step + (0.5 + ipsi_over) * psi_step / nr_ipsi_over;
                over_psi[ipsi_over] = overpsi;
            }// end loop ipsi_over
        }
    }// oversampling_order != 0
#endif
}

void HealpixSampler::getAllTranslations(int oversampling_order, FDOUBLE* my_translations_x, FDOUBLE* my_translations_y) const
{
    assert(oversampling_order < 4);
    //
    int nr_over_trans = oversamplingFactorTranslations(oversampling_order);
    for (int itrans = 0; itrans < nrTrans; itrans++)
        getTranslations(itrans,oversampling_order,
                        my_translations_x+itrans*nr_over_trans,
                        my_translations_y+itrans*nr_over_trans);
}

void HealpixSampler::getAllTranslations(int oversampling_order,Vector2d& my_translations_x,Vector2d& my_translations_y) const
{
    assert(oversampling_order < 4);
    //
    int nr_over_trans = oversamplingFactorTranslations(oversampling_order);
    for (int itrans = 0; itrans < nrTrans; itrans++)
        getTranslations(itrans,oversampling_order,
                        my_translations_x[itrans].data(),
                        my_translations_y[itrans].data());
}

//
void HealpixSampler::writeOutSampling(std::string fn_sampling)
{
    std::ofstream  samplingFile;
    samplingFile.open((fn_sampling+".star").c_str(), std::ios::out);
    {
        samplingFile << std::endl;
        samplingFile << "data_sampling_general" <<std::endl;
        samplingFile << std::endl;
#define COUTMETADATA(v1,v2) samplingFile << "_rln" << std::left<<std::setw(30) << v1 << std::right<<std::setw(18) << v2 <<std::endl;
        //
        COUTMETADATA("Is3DSampling"				, is_3D					)
        COUTMETADATA("Is3DTranslationalSampling", 0						)
        COUTMETADATA("HealpixOrder"				, healpix_order			)
        COUTMETADATA("SymmetryGroup"			, fn_sym				)
        COUTMETADATA("TiltAngleLimit"			, -91					)
        COUTMETADATA("PsiStep"					, psi_step				)
        COUTMETADATA("OffsetRange"				, offset_range			)
        COUTMETADATA("OffsetStep"				, offset_step			)
        COUTMETADATA("SamplingPerturbInstance"	, random_perturbation	)
        COUTMETADATA("SamplingPerturbFactor"	, 0.5					)
        //
#undef COUTMETADATA
        samplingFile << std::endl << std::endl;
    }
    // data_sampling_directions
    {
        SmallMetataDataTable data_sampling_directions("data_sampling_directions");
        data_sampling_directions.appendName({"AngleRot","AngleTilt"});
        data_sampling_directions.appendType({ElemTypeDouble,ElemTypeDouble});
        data_sampling_directions.appendElem({rot_angles.data(),tilt_angles.data()}, nrPix);
        data_sampling_directions.print(samplingFile);
    }
    samplingFile.close();
}

void HealpixSampler::readFromSampling(std::string fn_sampling,bool debug_flag /*= false*/)
{
    ifstreamCheckingExistence samplingFile(fn_sampling.c_str());
    {
        bool startingRead = false;
        std::string line;
        while (true) {
            if (startingRead) {
                double doubleTmp;std::string stringTmp;
#define CINMETADATADOUBLE(V) samplingFile >> line;samplingFile >> doubleTmp;if(debug_flag) {MASTERNODE std::cout<<std::setw(30)<<line<<" "<<doubleTmp<<std::endl;}
#define CINMETADATASTR(V) samplingFile >> line;samplingFile >> stringTmp;if(debug_flag) {MASTERNODE std::cout<<std::setw(30)<<line<<" "<<stringTmp<<std::endl;}
                //
                CINMETADATADOUBLE(	"Is3DSampling"				)
                CINMETADATADOUBLE(	"Is3DTranslationalSampling"	)
                CINMETADATADOUBLE(	"HealpixOrder"				)
                CINMETADATASTR(		"SymmetryGroup"				)
                CINMETADATADOUBLE(	"TiltAngleLimit"			)
                CINMETADATADOUBLE(	"PsiStep"					)
                CINMETADATADOUBLE(	"OffsetRange"				)
                CINMETADATADOUBLE(	"OffsetStep"				)
                CINMETADATADOUBLE(	"SamplingPerturbInstance"	);//random_perturbation = doubleTmp;
                CINMETADATADOUBLE(	"SamplingPerturbFactor"		);ERROR_CHECK(0.5!=doubleTmp, "Set sampling SamplingPerturbFactor.");
                //
#undef CINMETADATADOUBLE
#undef CINMETADATASTR
                break;
            }
            else{
                getline(samplingFile,line);
                if(debug_flag) MASTERNODE std::cout<<line<<std::endl;
                if ((line.find("data_sampling_general") !=std::string::npos) ){
                    startingRead = true;
                    getline(samplingFile,line);assert(line=="");// escape a empty line
                }
                ERROR_CHECK(samplingFile.eof(), "end of sampling file,can not find data_sampling_general.");
            }
        }
    }
    samplingFile.close();
}

// --------------------------------------------------------------------------------------------------- //

void SamplingGrid::initialize(HealpixSampler& sampler3d,int _adaptive_oversampling)
{
    nr_trans = sampler3d.NrTrans();
    nr_dir = sampler3d.NrDir();
    nr_psi = sampler3d.NrPsi();
    // size_t nr_orientation = exp_nr_dir*exp_nr_psi;
    // int exp_nr_over_rot_max   = sampler3d.oversamplingFactorOrientations(adaptive_oversampling);
    adaptive_oversampling = _adaptive_oversampling;
    int exp_nr_over_trans_max = sampler3d.oversamplingFactorTranslations(adaptive_oversampling);
    
    int nr_ipsi_over = round(std::pow(2., adaptive_oversampling));
    int nr_idir_over = round(std::pow(2., adaptive_oversampling*2));
    exp_over_psi	.resize(nr_psi, std::vector<FDOUBLE>(nr_ipsi_over,0));
    exp_over_rot	.resize(nr_dir, std::vector<FDOUBLE>(nr_idir_over,0));
    exp_over_tilt	.resize(nr_dir, std::vector<FDOUBLE>(nr_idir_over,0));
    thread_over_psi	.init(omp_get_max_threads(), nr_ipsi_over*nr_idir_over);thread_over_psi.fill_with_first_touch(0.);
    thread_over_rot	.init(omp_get_max_threads(), nr_ipsi_over*nr_idir_over);thread_over_rot.fill_with_first_touch(0.);
    thread_over_tilt.init(omp_get_max_threads(), nr_ipsi_over*nr_idir_over);thread_over_tilt.fill_with_first_touch(0.);
    //
    exp_over_trans_x.resize(nr_trans, std::vector<FDOUBLE>(exp_nr_over_trans_max,0));
    exp_over_trans_y.resize(nr_trans, std::vector<FDOUBLE>(exp_nr_over_trans_max,0));
    exp_trans_x		.resize(nr_trans);
    exp_trans_y		.resize(nr_trans);
    exp_trans_x_over.resize(exp_nr_over_trans_max);
    exp_trans_y_over.resize(exp_nr_over_trans_max);
    exp_trans_xyshift.resize(nr_trans);
    nr_over_trans = 0;
    nr_over_rot = 0;
}

void SamplingGrid::finalize()
{
    exp_over_psi	.resize(0);
    exp_over_tilt	.resize(0);
    exp_over_rot	.resize(0);
    thread_over_psi	.fini();
    thread_over_rot	.fini();
    thread_over_tilt.fini();
    //
    exp_over_trans_x.resize(0);
    exp_over_trans_y.resize(0);
    exp_trans_x		.resize(0);
    exp_trans_y		.resize(0);
    exp_trans_x_over.resize(0);
    exp_trans_y_over.resize(0);
    exp_positive_shift_index.clear();
    exp_trans_xyshift.resize(0);
    nr_trans = 0;nr_dir = 0;nr_psi = 0;
    nr_over_trans = 0;nr_over_rot = 0;
    adaptive_oversampling = 0;current_oversampling = 0;
}

void SamplingGrid::computeGrid2D(HealpixSampler& sampler2d,int _current_oversampling)
{
    assert(current_oversampling <= adaptive_oversampling);
    current_oversampling = _current_oversampling;
    sampler2d.getAllOrientations(current_oversampling, exp_over_psi);
    sampler2d.getAllTranslations(current_oversampling, exp_over_trans_x, exp_over_trans_y);
    sampler2d.getAllTranslationsAndOverTrans(current_oversampling, exp_trans_x, exp_trans_y, exp_trans_x_over, exp_trans_y_over);
    sampler2d.getAllSingleTranslationsAndOverTrans(current_oversampling, exp_positive_shift_index,
                                                   exp_trans_xyshift, exp_trans_x_over, exp_trans_y_over);
    nr_over_rot = sampler2d.oversamplingFactorOrientations(current_oversampling);
    nr_over_trans = sampler2d.oversamplingFactorTranslations(current_oversampling);
}

void SamplingGrid::computeGrid3D(HealpixSampler& sampler3d,int _current_oversampling)
{
    assert(current_oversampling <= adaptive_oversampling);
    current_oversampling = _current_oversampling;
    sampler3d.getAllOrientations(current_oversampling, exp_over_psi, exp_over_rot, exp_over_tilt);
    sampler3d.getAllTranslations(current_oversampling, exp_over_trans_x, exp_over_trans_y);
    sampler3d.getAllTranslationsAndOverTrans(current_oversampling, exp_trans_x, exp_trans_y, exp_trans_x_over, exp_trans_y_over);
    sampler3d.getAllSingleTranslationsAndOverTrans(current_oversampling, exp_positive_shift_index,
                                                   exp_trans_xyshift, exp_trans_x_over, exp_trans_y_over);
    nr_over_rot = sampler3d.oversamplingFactorOrientations(current_oversampling);
    nr_over_trans = sampler3d.oversamplingFactorTranslations(current_oversampling);
}

void SamplingGrid::selectOrientationsForLocalSearch(HealpixSampler& sampler3d,FDOUBLE sigma_rot,FDOUBLE sigma_tilt,FDOUBLE sigma_psi,
                                                    const MetaDataElem* exp_metadata,int nr_images)
{
    //assert(nr_images==1);
    assert(sigma_rot==sigma_tilt);
    assert(sigma_tilt==sigma_psi);
    pdf_orientation_for_nonzero_orientations.clear();
   
    pointer_dir_nonzeroprior.resize(nr_images);
    directions_prior		.resize(nr_images);
    pointer_psi_nonzeroprior.resize(nr_images);
    psi_prior				.resize(nr_images);
    std::vector<FDOUBLE> significant_images(nr_images,0);
    for (int iimage = 0; iimage < nr_images; iimage++)
    {
        // First try if there are some fixed prior angles
        FDOUBLE prior_rot = exp_metadata[iimage].ROT_PRIOR;
        FDOUBLE prior_tilt = exp_metadata[iimage].TILT_PRIOR;
        FDOUBLE prior_psi = exp_metadata[iimage].PSI_PRIOR;
        
        // If there were no defined priors (i.e. their values were 999.), then use the "normal" angles
        if (prior_rot > 998.99 && prior_rot < 999.01)
            prior_rot = exp_metadata[iimage].ROT;
        if (prior_tilt > 998.99 && prior_tilt < 999.01)
            prior_tilt = exp_metadata[iimage].TILT;
        if (prior_psi > 998.99 && prior_psi < 999.01)
            prior_psi = exp_metadata[iimage].PSI;
        
        pointer_dir_nonzeroprior[iimage].resize(0);directions_prior[iimage].resize(0);
        pointer_psi_nonzeroprior[iimage].resize(0);psi_prior[iimage].resize(0);
        sampler3d.selectOrientationsWithNonZeroPriorProbability(prior_rot, prior_tilt, prior_psi,
                                                                sqrt(sigma_rot), sqrt(sigma_tilt), sqrt(sigma_psi),
                                                                pointer_dir_nonzeroprior[iimage], directions_prior[iimage],
                                                                pointer_psi_nonzeroprior[iimage], psi_prior[iimage]);
        
        //
        for (const auto& directions_prior_p : directions_prior[iimage])
            ERROR_CHECK(directions_prior_p <= 0, "directions_prior[iimage] smaller than zero???");
        for (const auto& psi_prior_p : psi_prior[iimage])
            ERROR_CHECK(psi_prior_p <= 0, "psi_prior[iimage] smaller than zero???");
            
//#define DEBUG_LOCAL_SEARCH
#ifdef DEBUG_LOCAL_SEARCH
#define PRINT(V) std::cout<<#V<<" : ";for (const auto& v : V) std::cout<<v<<" ";std::cout<<std::endl;
        PRINT(pointer_dir_nonzeroprior[iimage])
        PRINT(directions_prior[iimage])
        PRINT(pointer_psi_nonzeroprior[iimage])
        PRINT(psi_prior[iimage])
#undef PRINT
#endif
        //
        for (int i_idir = 0; i_idir < pointer_dir_nonzeroprior[iimage].size(); i_idir++) {
            for (int i_ipsi = 0; i_ipsi < pointer_psi_nonzeroprior[iimage].size(); i_ipsi++) {
                auto orientation = std::make_pair(pointer_dir_nonzeroprior[iimage][i_idir], pointer_psi_nonzeroprior[iimage][i_ipsi]);
                auto pdf_orientation = directions_prior[iimage][i_idir]*psi_prior[iimage][i_ipsi];
                auto it = pdf_orientation_for_nonzero_orientations.find(orientation);
                if(it == pdf_orientation_for_nonzero_orientations.end()){
                    // insert newer
                    auto newit = pdf_orientation_for_nonzero_orientations.insert(it, std::make_pair(orientation,significant_images));
                    newit->second[iimage] = pdf_orientation;
                }
                else{
                    // modify
                    assert(it->second[iimage]==0);
                    it->second[iimage] = pdf_orientation;
                }
            }
        }
    }
#ifdef DEBUG_LOCAL_SEARCH
    for (const auto& orientationAndImage : pdf_orientation_for_nonzero_orientations) {
        std::cout<<"{ " <<orientationAndImage.first.first<<","<<orientationAndImage.first.second<<" } : ";
        for (const auto& p : orientationAndImage.second) {
            std::cout<<p<<" , ";
        }
        std::cout<<std::endl;
    }
#endif
}

void SamplingGrid::setOrientation(HealpixSampler& sampler3d,int _current_oversampling,int idir,int ipsi,
                                  FDOUBLE* &over_rot,FDOUBLE* &over_tilt,FDOUBLE* &over_psi)
{
    assert(current_oversampling==_current_oversampling);
    auto tid = omp_get_thread_num();
    over_psi = thread_over_psi[tid].wptrAll();
    over_rot = thread_over_rot[tid].wptrAll();
    over_tilt = thread_over_tilt[tid].wptrAll();
    int nr_over = sampler3d.oversamplingFactorOrientations(_current_oversampling);
    if (nr_over==1)
    {
        over_rot[0] = exp_over_rot[idir][0];
        over_tilt[0] = exp_over_tilt[idir][0];
        over_psi[0] = exp_over_psi[ipsi][0];
        sampler3d.addPerturbation(current_oversampling, over_rot, over_tilt, over_psi);
    }
    else
    {
        int nr_idir_over = round(std::pow(2., current_oversampling*2));
        int nr_ipsi_over = round(std::pow(2., current_oversampling));
        int iover = 0;
        for (int idir_over = 0; idir_over < nr_idir_over; idir_over++) {
            for (int ipsi_over = 0; ipsi_over < nr_ipsi_over; ipsi_over++) {
                over_psi[iover] = exp_over_psi[ipsi][ipsi_over];
                over_rot[iover] = exp_over_rot[idir][idir_over];
                over_tilt[iover] = exp_over_tilt[idir][idir_over];
                iover++;
            }
        }
        sampler3d.addPerturbation(current_oversampling, over_rot, over_tilt, over_psi);
    }
}

void SamplingGrid::testGetShift()
{
    //std::cout<<"Starting compare get shift.."<<std::endl;
    for (int itrans = 0; itrans < nr_trans; itrans++) {
        for (int iover_trans = 0; iover_trans < nr_over_trans; iover_trans++) {
            FDOUBLE shiftx1,shifty1;
            FDOUBLE shiftx2,shifty2,shiftx2Over,shifty2Over;
            getShiftxy(shiftx1, shifty1, itrans, iover_trans);
            getShiftxy(shiftx2, shifty2, shiftx2Over, shifty2Over, itrans, iover_trans);
            //std::cout<<"itrans = "<<itrans<<",iover_trans = "<<iover_trans<<std::endl;
            if (fabs(shiftx1-shiftx2-shiftx2Over) > 1e-12 || fabs(shifty1-shifty2-shifty2Over) > 1e-12) {
                std::cout<<"!!! "<<shiftx1<<" != "<<shiftx2<<" + "<<shiftx2Over<<","<<shifty1<<" != "<<shifty2<<" + "<<shifty2Over<<std::endl;
                ERROR_REPORT("Wrong getShiftxy..");
            }
            else{
                //std::cout<<shiftx1<<" = "<<shiftx2<<" + "<<shiftx2Over<<","<<shifty1<<" = "<<shifty2<<" + "<<shifty2Over<<std::endl;
            }
        }
    }
}

///

namespace GTMSampling{
    // info is var_num*3 matrix each row store basic info(start,end,order(indexBits))
    // notice ,if we declare an pointer out of this function and then allocate the memmory for X in this function
    // the allocated memory will be free...
    double* uniformSampling(const double *infos,int var_num,int &K,int &L){
        
        L = var_num;
        
        K = 1;
        int bit_sum = 0;
        for(int var = 0;var < L;var++){
            bit_sum += infos[var*3+2];
            K *= (1 << (int)infos[var*3+2]);
        }
        // K = pow(2,bit_sum)
        // std::cout<<"K = "<<K<<",L = "<<L<<std::endl;
        // std::cout<<"need space:"<<(double)K*L*8./1024./1024.<<" MB"<<std::endl;
        
        ERROR_CHECK(bit_sum > 32, "bit sum out of range(32).");
        
        double *X = vNew(double,K*L);
        int subIndex = 0;
        int tempIndex;
        
        for(int index = 0;index < K;index++){
            tempIndex = index;
            
            for(int var = 0;var < L;var++){
                subIndex = tempIndex & ((1 << (int)infos[var*3+2]) - 1);
                
                X[index*L+var] = infos[var*3+0] + (infos[var*3+1] - infos[var*3+0])/((1 << (int)infos[var*3+2]) - 1)*subIndex;
                tempIndex = tempIndex >> (int)infos[var*3+2];
            }
        }
        
        std::cout<<"end sampling."<<std::endl;
        return X;
    }
    
    
    double* sample1L(const double *infos,int &K,int &L){
        
        L = 1;
        K = (int)infos[2];
        
        double *X = vNew(double, K*L);
        
        for(int i = 0;i < K;i++)
            X[i] = infos[0] + (infos[1] - infos[0])/(K-1)*i;
        
        return X;
    }
    
    
    //   example:infos[6] = {1,3,3,1,2,2}
    //   result:
    //   1 1
    //   2 1
    //   3 1
    //   1 2
    //   2 2
    //   3 2
    //
    double* sample2L(const double *infos,int &K,int &L,bool inv){
        
        L = 2;
        K = (int)infos[2]*(int)infos[5];
        
        double *X = vNew(double,K*L);
        int i,j;
        
        if(inv){
            for(int k = 0;k < K;k++){
                i = k / (int)infos[5];
                j = k % (int)infos[5];
                X[k*L+0] = infos[0] + (infos[1] - infos[0])/((int)infos[2]-1)*i;
                X[k*L+1] = infos[3] + (infos[4] - infos[3])/((int)infos[5]-1)*j;
            }
        }
        else
        {
            for(int k = 0;k < K;k++){
                i = k / (int)infos[2];
                j = k % (int)infos[2];
                X[k*L+0] = infos[0] + (infos[1] - infos[0])/((int)infos[2]-1)*j;
                X[k*L+1] = infos[3] + (infos[4] - infos[3])/((int)infos[5]-1)*i;
            }
        }
        
        return X;
    }
    
    
    //   example:infos[9] = {1,3,3,1,2,2,1,2,2}
    //   result:
    //   1               1               1
    //   1               1               2
    //   1               2               1
    //   1               2               2
    //   2               1               1
    //   2               1               2
    //   2               2               1
    //   2               2               2
    //   3               1               1
    //   3               1               2
    //   3               2               1
    //   3               2               2
    //
    double* sample3L(const double *infos,int &K,int &L,bool inv){
        
        L = 3;
        K = (int)infos[2]*(int)infos[5]*(int)infos[8];
        
        double *X = vNew(double,K*L);
        int i1,i2,i3,i;
        if(inv){
            for(int k = 0;k < K;k++){
                i1 = k / ((int)infos[5]*(int)infos[8]);
                
                i = k % ((int)infos[5]*(int)infos[8]);
                i2 = i / (int)infos[8];
                i3 = i % (int)infos[8];
                
                X[k*L+0] = infos[0] + (infos[1] - infos[0])/((int)infos[2]-1)*i1;
                X[k*L+1] = infos[3] + (infos[4] - infos[3])/((int)infos[5]-1)*i2;
                X[k*L+2] = infos[6] + (infos[7] - infos[6])/((int)infos[8]-1)*i3;
            }
        }
        else
        {
            for(int k = 0;k < K;k++){
                i1 = k / ((int)infos[5]*(int)infos[2]);
                
                i = k % ((int)infos[5]*(int)infos[2]);
                i2 = i / (int)infos[2];
                i3 = i % (int)infos[2];
                
                X[k*L+0] = infos[0] + (infos[1] - infos[0])/((int)infos[2]-1)*i3;
                X[k*L+1] = infos[3] + (infos[4] - infos[3])/((int)infos[5]-1)*i2;
                X[k*L+2] = infos[6] + (infos[7] - infos[6])/((int)infos[8]-1)*i1;
            }
        }
        
        return X;
    }
    
    double* sample4L(const double *infos,int &K,int &L,bool inv){
        
        L = 4;
        K = (int)infos[2]*(int)infos[5]*(int)infos[8]*(int)infos[11];
        
        double *X = vNew(double,K*L);
        int i1,i2,i3,i4,i;
        if(inv){
            for(int k = 0;k < K;k++){
                i1 = k / ((int)infos[5]*(int)infos[8]*(int)infos[11]);
                
                i = k % ((int)infos[5]*(int)infos[8]*(int)infos[11]);
                i2 = i / ((int)infos[8]*(int)infos[11]);
                
                i = i % ((int)infos[8]*(int)infos[11]);
                i3 = i / (int)infos[11];
                i4 = i % (int)infos[11];
                
                X[k*L+0] = infos[0] + (infos[1] - infos[0])/((int)infos[2]-1)*i1;
                X[k*L+1] = infos[3] + (infos[4] - infos[3])/((int)infos[5]-1)*i2;
                X[k*L+2] = infos[6] + (infos[7] - infos[6])/((int)infos[8]-1)*i3;
                X[k*L+3] = infos[9] + (infos[10] - infos[9])/((int)infos[11]-1)*i4;
            }
        }
        else
        {
            for(int k = 0;k < K;k++){
                i4 = k / ((int)infos[8]*(int)infos[5]*(int)infos[2]);
                
                i = k % ((int)infos[8]*(int)infos[5]*(int)infos[2]);
                i3 = i / ((int)infos[5]*(int)infos[2]);
                
                i = i % ((int)infos[5]*(int)infos[2]);
                i2 = i / (int)infos[2];
                i1 = i % (int)infos[2];
                
                X[k*L+0] = infos[0] + (infos[1] - infos[0])/((int)infos[2]-1)*i1;
                X[k*L+1] = infos[3] + (infos[4] - infos[3])/((int)infos[5]-1)*i2;
                X[k*L+2] = infos[6] + (infos[7] - infos[6])/((int)infos[8]-1)*i3;
                X[k*L+3] = infos[9] + (infos[10] - infos[9])/((int)infos[11]-1)*i4;
            }
        }
        
        return X;
    }
    
    double* sample5L(const double *infos,int &K,int &L,bool inv){
        
        L = 5;
        K = (int)infos[2]*(int)infos[5]*(int)infos[8]*(int)infos[11]*(int)infos[14];
        
        double *X = vNew(double,K*L);
        int i1,i2,i3,i4,i5,i;
        if(inv){
            for(int k = 0;k < K;k++){
                i1 = k / ((int)infos[5]*(int)infos[8]*(int)infos[11]*(int)infos[14]);
                
                i = k % ((int)infos[5]*(int)infos[8]*(int)infos[11]*(int)infos[14]);
                i2 = i / ((int)infos[8]*(int)infos[11]*(int)infos[14]);
                
                i = i % ((int)infos[8]*(int)infos[11]*(int)infos[14]);
                i3 = i / ((int)infos[11]*(int)infos[14]);
                
                i = i % ((int)infos[11]*(int)infos[14]);
                i4 = i / (int)infos[14];
                i5 = i % (int)infos[14];
                
                X[k*L+0] = infos[0] + (infos[1] - infos[0])/((int)infos[2]-1)*i1;
                X[k*L+1] = infos[3] + (infos[4] - infos[3])/((int)infos[5]-1)*i2;
                X[k*L+2] = infos[6] + (infos[7] - infos[6])/((int)infos[8]-1)*i3;
                X[k*L+3] = infos[9] + (infos[10] - infos[9])/((int)infos[11]-1)*i4;
                X[k*L+4] = infos[12] + (infos[13] - infos[12])/((int)infos[14]-1)*i5;
            }
        }
        else
        {
            for(int k = 0;k < K;k++){
                i5 = k / ((int)infos[11]*(int)infos[8]*(int)infos[5]*(int)infos[2]);
                
                i = k % ((int)infos[11]*(int)infos[8]*(int)infos[5]*(int)infos[2]);
                i4 = i / ((int)infos[8]*(int)infos[5]*(int)infos[2]);
                
                i = i % ((int)infos[8]*(int)infos[5]*(int)infos[2]);
                i3 = i / ((int)infos[5]*(int)infos[2]);
                
                i = i % ((int)infos[5]*(int)infos[2]);
                i2 = i / (int)infos[2];
                i1 = i % (int)infos[2];
                
                X[k*L+0] = infos[0] + (infos[1] - infos[0])/((int)infos[2]-1)*i1;
                X[k*L+1] = infos[3] + (infos[4] - infos[3])/((int)infos[5]-1)*i2;
                X[k*L+2] = infos[6] + (infos[7] - infos[6])/((int)infos[8]-1)*i3;
                X[k*L+3] = infos[9] + (infos[10] - infos[9])/((int)infos[11]-1)*i4;
                X[k*L+4] = infos[12] + (infos[13] - infos[12])/((int)infos[14]-1)*i5;
            }
        }
        
        return X;
    }
    
    // #include "./Healpix_2.15a/healpix_base.h"
    // double *sampleHealpix(int order,int &K,int &L){
    
    // 	std::cout<<"staring s sampling."<<std::endl;
    
    //     const double Pi = 3.1415926535897;
    
    //     Healpix_Base healpix_base;
    //     int healpix_order = order;   //N_side = 2^order
    //     healpix_base.Set(healpix_order, RING);   //initialize,RING//NEST
    //     int Npix = healpix_base.Npix();
    
    //     K = Npix;
    //     L = 2;
    
    //     ///////////////////some unnecessary test codes/////////////////////////////////
    //     std::cout<<"Order = "<<healpix_order<<std::endl;
    //     std::cout<<"Nside = "<<healpix_base.Nside()<<std::endl;
    //     std::cout<<"Npix = "<<Npix<<std::endl;
    
    //     double *X = (double*)aMalloc(sizeof(double)*2*Npix,64);
    
    //     double zz,phi;
    //     for(int ipix = 0;ipix < Npix;ipix++){
    //       healpix_base.pix2ang_z_phi(ipix, zz, phi);
    //       X[2*ipix] = (acos(zz)/*-Pi/2*/);
    //       X[2*ipix+1] = phi;
    //     }
    
    //     std::cout<<"end sampling."<<std::endl;
    
    //     return X;
    // }
}

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
// The way symmetry is handled was copied from Xmipp.
// The original disclaimer is copied below
/***************************************************************************
 *
 * Authors:     Roberto Marabini
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
/* Calculate an angular distance between two sets of Euler angles */
FDOUBLE HealpixSampler::calculateAngularDistance(FDOUBLE rot1, FDOUBLE tilt1, FDOUBLE psi1,
                                                 FDOUBLE rot2, FDOUBLE tilt2, FDOUBLE psi2)
{
    if (is_3D)
    {
        Matrix1D<FDOUBLE>  direction1(3), direction1p(3), direction2(3);
        Euler_angles2direction(rot1, tilt1, direction1.data().data());
        Euler_angles2direction(rot2, tilt2, direction2.data().data());
        
        // Find the symmetry operation where the Distance based on Euler axes is minimal
        FDOUBLE min_axes_dist = 3600.;
        FDOUBLE rot2p, tilt2p, psi2p;
        Matrix2D<FDOUBLE> E1(3,3), E2(3,3);
        Matrix1D<FDOUBLE> v1(3), v2(3);
        for (int j = 0; j < R_repository.size(); j++)
        {
            Euler_apply_transf(L_repository[j], R_repository[j], rot2, tilt2, psi2, rot2p, tilt2p, psi2p);
            // Distance based on Euler axes
            Euler_angles2matrix(rot1, tilt1, psi1, E1);
            Euler_angles2matrix(rot2p, tilt2p, psi2p, E2);
            FDOUBLE axes_dist = 0;
            for (int i = 0; i < 3; i++)
            {
                E1.getRow(i, v1);
                E2.getRow(i, v2);
                axes_dist += ACOSD(CLIP(dotProduct(v1, v2), -1., 1.));
            }
            axes_dist /= 3.;
            
            if (axes_dist < min_axes_dist)
                min_axes_dist = axes_dist;
            
        }// for all symmetry operations j
        
        return min_axes_dist;
    }
    else
    {
        FDOUBLE diff = fabs(psi2 - psi1);
        return realWRAP(diff, 0., 360.);
    }
}
void HealpixSampler::selectOrientationsWithNonZeroPriorProbability(FDOUBLE prior_rot, FDOUBLE prior_tilt, FDOUBLE prior_psi,
                                                                   FDOUBLE sigma_rot, FDOUBLE sigma_tilt, FDOUBLE sigma_psi,
                                                                   std::vector<int> &pointer_dir_nonzeroprior, std::vector<FDOUBLE> &directions_prior,
                                                                   std::vector<int> &pointer_psi_nonzeroprior, std::vector<FDOUBLE> &psi_prior,
                                                                   FDOUBLE sigma_cutoff)
{
    pointer_dir_nonzeroprior.clear();
    directions_prior.clear();
    
    if (is_3D)
    {
        // Loop over all directions
        FDOUBLE sumprior = 0.;
        // Keep track of the closest distance to prevent 0 orientations
        FDOUBLE best_ang = 9999.;
        int best_idir = -999;
        for (int idir = 0; idir < nrPix; idir++)
        {
            // Any prior involving rot and/or tilt.
            if (sigma_rot > 0. || sigma_tilt > 0. )
            {
                Matrix1D<FDOUBLE> prior_direction(3), my_direction(3), sym_direction(3), best_direction(3);
                // Get the direction of the prior
                Euler_angles2direction(prior_rot, prior_tilt, prior_direction.data().data());
                
                // Get the current direction in the loop
                Euler_angles2direction(rot_angles[idir], tilt_angles[idir], my_direction.data().data());
                
                // Loop over all symmetry operators to find the operator that brings this direction nearest to the prior
                FDOUBLE best_dotProduct = dotProduct(prior_direction, my_direction);
                best_direction = my_direction;
                for (int j = 0; j < R_repository.size(); j++)
                {
                    sym_direction =  L_repository[j] * (my_direction.transpose() * R_repository[j]).transpose();
                    FDOUBLE my_dotProduct = dotProduct(prior_direction, sym_direction);
                    if (my_dotProduct > best_dotProduct)
                    {
                        best_direction = sym_direction;
                        best_dotProduct = my_dotProduct;
                    }
                }
                
                if (sigma_rot > 0. && sigma_tilt > 0.)
                {
                    
                    FDOUBLE diffang = acosd( dotProduct(best_direction, prior_direction) );
                    if (diffang > 180.) diffang = fabs(diffang - 360.);
                    
                    // Only consider differences within sigma_cutoff * sigma_rot
                    if (diffang < sigma_cutoff * sigma_rot)
                    {
                        // TODO!!! If tilt is zero then any rot will be OK!!!!!
                        FDOUBLE prior = gaussian1D(diffang, sigma_rot, 0.);
                        pointer_dir_nonzeroprior.push_back(idir);
                        directions_prior.push_back(prior);
                        sumprior += prior;
                    }
                    
                    // Keep track of the nearest direction
                    if (diffang < best_ang)
                    {
                        best_idir = idir;
                        best_ang = diffang;
                    }
                }
                else if (sigma_rot > 0.)
                {
                    FDOUBLE best_rot, best_tilt;
                    
                    Euler_direction2angles(best_direction.data().data(), best_rot, best_tilt);
                    FDOUBLE diffrot = fabs(best_rot - prior_rot);
                    if (diffrot > 180.) diffrot = fabs(diffrot - 360.);
                    
                    // Only consider differences within sigma_cutoff * sigma_rot
                    if (diffrot < sigma_cutoff * sigma_rot)
                    {
                        FDOUBLE prior = gaussian1D(diffrot, sigma_rot, 0.);
                        pointer_dir_nonzeroprior.push_back(idir);
                        directions_prior.push_back(prior);
                        sumprior += prior;
                    }
                    
                    // Keep track of the nearest direction
                    if (diffrot < best_ang)
                    {
                        best_idir = idir;
                        best_ang = diffrot;
                    }
                    
                }
                else if (sigma_tilt > 0.)
                {
                    
                    FDOUBLE best_rot, best_tilt;
                    
                    Euler_direction2angles(best_direction.data().data(), best_rot, best_tilt);
                    FDOUBLE difftilt = fabs(best_tilt - prior_tilt);
                    if (difftilt > 180.) difftilt = fabs(difftilt - 360.);
                    
                    // Only consider differences within sigma_cutoff * sigma_tilt
                    if (difftilt < sigma_cutoff * sigma_tilt)
                    {
                        FDOUBLE prior = gaussian1D(difftilt, sigma_tilt, 0.);
                        pointer_dir_nonzeroprior.push_back(idir);
                        directions_prior.push_back(prior);
                        sumprior += prior;
                    }
                    
                    // Keep track of the nearest direction
                    if (difftilt < best_ang)
                    {
                        best_idir = idir;
                        best_ang = difftilt;
                    }
                    
                }
                
            } // end if any prior involving rot and/or tilt
            else
            {
                // If no prior on the directions: just add all of them
                pointer_dir_nonzeroprior.push_back(idir);
                directions_prior.push_back(1.);
                sumprior += 1.;
            }
            
        } // end for idir
        
        
        //Normalise the prior probability distribution to have sum 1 over all psi-angles
        for (int idir = 0; idir < directions_prior.size(); idir++)
            directions_prior[idir] /= sumprior;
        
        // If there were no directions at all, just select the single nearest one:
        if (directions_prior.size() == 0)
        {
            if (best_idir < 0)
                ERROR_REPORT("HealpixSampler::selectOrientationsWithNonZeroPriorProbability BUG: best_idir < 0");
            pointer_dir_nonzeroprior.push_back(best_idir);
            directions_prior.push_back(1.);
        }
        
#ifdef  DEBUG_SAMPLING
        writeNonZeroPriorOrientationsToBild("orients_local.bild", prior_rot, prior_tilt, pointer_dir_nonzeroprior, "0 0 1", 0.023);
        std::cerr << " directions_prior.size()= " << directions_prior.size() << " pointer_dir_nonzeroprior.size()= " << pointer_dir_nonzeroprior.size() << std::endl;
        std::cerr << " sumprior= " << sumprior << std::endl;
        char c;
        std::cerr << "Written orients_local.bild for prior on angles ("<<prior_rot<<","<<prior_tilt<<") Press any key to continue.." << std::endl;
        std::cin >> c;
#endif
        
        
    }
    else
    {
        pointer_dir_nonzeroprior.push_back(0);
        directions_prior.push_back(1.);
    }
    
    
    // Psi-angles
    pointer_psi_nonzeroprior.clear();
    psi_prior.clear();
    
    FDOUBLE sumprior = 0.;
    FDOUBLE best_diff = 9999.;
    int best_ipsi = -999;
    for (int ipsi = 0; ipsi < nrPsi; ipsi++)
    {
        if (sigma_psi > 0.)
        {
            FDOUBLE diffpsi = fabs(psi_angles[ipsi] - prior_psi);
            if (diffpsi > 180.) diffpsi = fabs(diffpsi - 360.);
            
            // Only consider differences within sigma_cutoff * sigma_psi
            if (diffpsi < sigma_cutoff * sigma_psi)
            {
                FDOUBLE prior = gaussian1D(diffpsi, sigma_psi, 0.);
                pointer_psi_nonzeroprior.push_back(ipsi);
                psi_prior.push_back(prior);
                sumprior += prior;
                
                // TMP DEBUGGING
                if (prior == 0.)
                {
                    std::cerr << " psi_angles[ipsi]= " << psi_angles[ipsi] << " prior_psi= " << prior_psi << " orientational_prior_mode= " << orientational_prior_mode << std::endl;
                    std::cerr << " diffpsi= " << diffpsi << " sigma_cutoff= " << sigma_cutoff << " sigma_psi= " << sigma_psi << std::endl;
                    ERROR_REPORT("prior on psi is zero!");
                }
                
            }
            // Keep track of the nearest sampling point
            if (diffpsi < best_diff)
            {
                best_ipsi = ipsi;
                best_diff = diffpsi;
            }
        }
        else
        {
            pointer_psi_nonzeroprior.push_back(ipsi);
            psi_prior.push_back(1.);
            sumprior += 1.;
        }
    }
    // Normalise the prior probability distribution to have sum 1 over all psi-angles
    for (int ipsi = 0; ipsi < psi_prior.size(); ipsi++)
        psi_prior[ipsi] /= sumprior;
    
    // If there were no directions at all, just select the single nearest one:
    if (psi_prior.size() == 0)
    {
        if (best_ipsi < 0)
            ERROR_REPORT("HealpixSampling::selectOrientationsWithNonZeroPriorProbability BUG: best_ipsi < 0");
        pointer_psi_nonzeroprior.push_back(best_ipsi);
        psi_prior.push_back(1.);
    }
    
    
#ifdef  DEBUG_SAMPLING
    std::cerr << " psi_angles.size()= " << psi_angles.size() << " psi_step= " << psi_step << std::endl;
    std::cerr << " psi_prior.size()= " << psi_prior.size() << " pointer_psi_nonzeroprior.size()= " << pointer_psi_nonzeroprior.size() << " sumprior= " << sumprior << std::endl;
#endif
    
    
}

void HealpixSampler::writeAllOrientationsToBild(std::string fn_bild, std::string rgb, FDOUBLE size)
{
    std::ofstream out;
    out.open (fn_bild.c_str());
    if (!out)
        ERROR_REPORT( (std::string)"HealpixSampling::writeAllOrientationsToBild: Cannot write file: " + fn_bild);
    
    
    out << ".color 1 0 0 \n";
    out << ".arrow 0 0 0 1 0 0 0.01 \n";
    out << ".color 0 1 0 \n";
    out << ".arrow 0 0 0  0 1 0 0.01 \n";
    out << ".color 0 0 1 \n";
    out << ".arrow 0 0 0 0 0 1 0.01 \n";
    
    
    Matrix1D<FDOUBLE> v(3);
    out << ".color " << rgb << std::endl;
    
    for (unsigned long int ipix = 0; ipix < rot_angles.size(); ipix++)
    {
        Euler_angles2direction(rot_angles[ipix], tilt_angles[ipix], v.data().data());
        out <<  ".sphere " << XX(v) << " " << YY(v) << " " << ZZ(v) <<  " " << std::to_string((long long)size) << std::endl;
    }
    
    out.close();
}

void HealpixSampler::writeBildFileOrientationalDistribution(VectorOfFDOUBLE &pdf_direction,
                                                             std::string &fn_bild, FDOUBLE R, FDOUBLE offset, FDOUBLE Rmax_frac, FDOUBLE width_frac)
{
    if (!is_3D)
        return;

    if (pdf_direction.size() != rot_angles.size())
        ERROR_REPORT("HealpixSampling::writeBildFileOrientationalDistribution XSIZE(pdf_direction) != rot_angles.size()!");
    
    
    FDOUBLE pdfmax, pdfmin, pdfmean, pdfsigma;
    computeStats(pdf_direction.rptrAll(), pdf_direction.size(), pdfmean, pdfsigma, pdfmin, pdfmax);
    
    std::ofstream fh_bild;
    fh_bild.open(fn_bild.c_str(), std::ios::out);
    if (!fh_bild)
        ERROR_REPORT("HealpixSampling::writeBildFileOrientationalDistribution: cannot open " + fn_bild);
    
    // 2 * PI * R = 360 degrees, 2*radius should cover angular sampling at width_frac=1
    FDOUBLE width = width_frac * PI*R*(getAngularSampling()/360.);
    Matrix1D<FDOUBLE> v(3);
    
    for (int iang = 0; iang < rot_angles.size(); iang++)
    {
        FDOUBLE pdf = pdf_direction[iang];
        
        // Don't make a cylinder for pdf==0
        if (pdf > 0.)
        {
            // Colour from blue to red according to deviations from sigma_pdf
            FDOUBLE colscale = (pdf - pdfmean) / pdfsigma;
            colscale = std::min(colscale, (FDOUBLE)5.);
            colscale = std::max(colscale, (FDOUBLE)-1.);
            colscale /= 6.;
            colscale += 1./6.; // colscale ranges from 0 (-5 sigma) to 1 (+5 sigma)
            
            // The length of the cylinder will depend on the pdf_direction
            FDOUBLE Rp = R + Rmax_frac * R * pdf / pdfmax;
            
            Euler_angles2direction(rot_angles[iang], tilt_angles[iang], v.data().data());
            
            // Don't include cylinders with zero length, as chimera will complain about that....
            if (fabs((R - Rp) * XX(v)) > 0.01 ||
                fabs((R - Rp) * YY(v)) > 0.01 ||
                fabs((R - Rp) * ZZ(v)) > 0.01)
            {
                // The width of the cylinders will be determined by the sampling:
                fh_bild << ".color " << colscale << " 0 " << 1. - colscale << std::endl;
                fh_bild << ".cylinder "
                << R  * XX(v) + offset << " "
                << R  * YY(v) + offset << " "
                << R  * ZZ(v) + offset << " "
                << Rp * XX(v) + offset << " "
                << Rp * YY(v) + offset << " "
                << Rp * ZZ(v) + offset << " "
                << width
                <<"\n";
            }
        }
        
    }
    
    // Close and write file to disc
    fh_bild.close();
}

void HealpixSampler::removeSymmetryEquivalentPoints(FDOUBLE max_ang)
{
    // Maximum distance
    FDOUBLE cos_max_ang = cos(DEG2RAD(max_ang));
    FDOUBLE my_dotProduct;
    Matrix1D<FDOUBLE>  direction(3), direction1(3);
    std::vector< Matrix1D<FDOUBLE> > directions_vector;
    
    // Calculate all vectors and fill directions_vector
    for (long int i = 0; i < rot_angles.size(); i++)
    {
        Euler_angles2direction(rot_angles[i], tilt_angles[i], direction.data().data());
        directions_vector.push_back(direction);
    }
    
    // First call to conventional remove_redundant_points
    removeSymmetryEquivalentPointsGeometric(pgGroup, pgOrder, directions_vector);
    
#ifdef  DEBUG_SAMPLING
    writeAllOrientationsToBild("orients_sym0.bild", "0 1 0", 0.021);
#endif
    
    // Only correct the seams (i.e. the borders of the asymmetric units) for small numbers of directions
    // For large numbers, the sampling is very fine and the probability distributions are probably delta functions anyway
    // Large numbers take long times to calculate...
    // Only a small fraction of the points at the border of the AU is thrown away anyway...
    if (rot_angles.size() < 4000)
    {
        // Create no_redundant vectors
        std::vector <Matrix1D<FDOUBLE> > no_redundant_directions_vector;
        std::vector <FDOUBLE> no_redundant_rot_angles;
        std::vector <FDOUBLE> no_redundant_tilt_angles;
        std::vector <int> no_redundant_directions_ipix;
        
        // Then check all points versus each other
        for (long int i = 0; i < rot_angles.size(); i++)
        {
            
            direction1=directions_vector[i];
            bool uniq = true;
            
            //for (long int k = 0; k < no_redundant_directions_vector.size(); k++)
            // i is probably closer to latest additions: loop backwards over k....
            for (long int k = no_redundant_directions_vector.size() -1; k >= 0; k--)
            {
                for (int j = 0; j < R_repository.size(); j++)
                {
                    auto dir_mul_R = no_redundant_directions_vector[k].transpose() * R_repository[j];
                    direction =  L_repository[j] * ( (dir_mul_R).transpose() );
                    //Calculate distance
                    my_dotProduct = dotProduct(direction,direction1);
                    if (my_dotProduct > cos_max_ang)
                    {
                        uniq = false;
                        break;
                    }
                }// for j
                if (!uniq) break;
            } // for k
            
            if (uniq)
            {
                no_redundant_directions_vector.push_back(directions_vector[i]);
                no_redundant_rot_angles.push_back(rot_angles[i]);
                no_redundant_tilt_angles.push_back(tilt_angles[i]);
                no_redundant_directions_ipix.push_back(directions_ipix[i]);
            }
        } // for i
        
        // Now overwrite the rot/tilt_angles and directions_vectors with their no_redundant counterparts
        rot_angles = no_redundant_rot_angles;
        tilt_angles = no_redundant_tilt_angles;
        directions_ipix = no_redundant_directions_ipix;
    }
}

void HealpixSampler::removeSymmetryEquivalentPointsGeometric(const int symmetry,int sym_order,
                                                             std::vector <Matrix1D<FDOUBLE> >  &directions_vector)
{
    Matrix2D<FDOUBLE>  L(4, 4), R(4, 4);
    Matrix2D<FDOUBLE>  aux(3, 3);
    Matrix1D<FDOUBLE>  row1(3), row2(3), row(3);
    
    std::vector <Matrix1D<FDOUBLE> > no_redundant_directions_vector;
    std::vector <FDOUBLE> no_redundant_rot_angles;
    std::vector <FDOUBLE> no_redundant_tilt_angles;
    std::vector <int> no_redundant_directions_ipix;
    
    FDOUBLE my_dotProduct;
    if (symmetry == pg_CN)
    {//OK
        for (long int i = 0; i < rot_angles.size(); i++)
        {
            if (rot_angles[i] >= (-180. / sym_order) &&
                rot_angles[i] <= (180. / sym_order))
            {
                no_redundant_rot_angles.push_back(rot_angles[i]);
                no_redundant_tilt_angles.push_back(tilt_angles[i]);
                no_redundant_directions_vector.push_back(directions_vector[i]);
                no_redundant_directions_ipix.push_back(directions_ipix[i]);
            }
        }// for i
    }
    else if (symmetry == pg_CI  ||
             symmetry == pg_CS )
    {//OK
        for (long int i = 0; i < rot_angles.size(); i++)
        {
            if (tilt_angles[i] <= 90)
            {
                no_redundant_rot_angles.push_back(rot_angles[i]);
                no_redundant_tilt_angles.push_back(tilt_angles[i]);
                no_redundant_directions_vector.push_back(directions_vector[i]);
                no_redundant_directions_ipix.push_back(directions_ipix[i]);
            }
        }// for i
    }
    else if (symmetry  == pg_CNV )
    {//OK
        for (long int i = 0; i < rot_angles.size(); i++)
        {
            if (rot_angles[i] >=    0. / sym_order &&
                rot_angles[i] <=  180. / sym_order)
            {
                no_redundant_rot_angles.push_back(rot_angles[i]);
                no_redundant_tilt_angles.push_back(tilt_angles[i]);
                no_redundant_directions_vector.push_back(directions_vector[i]);
                no_redundant_directions_ipix.push_back(directions_ipix[i]);
            }
        }// for i
    }
    else if (symmetry  == pg_CNH )
    {//OK
        for (long int i = 0; i < rot_angles.size(); i++)
        {
            if (rot_angles[i] >= -180. / sym_order &&
                rot_angles[i] <=  180. / sym_order &&
                tilt_angles[i] <=    90.
                )
            {
                no_redundant_rot_angles.push_back(rot_angles[i]);
                no_redundant_tilt_angles.push_back(tilt_angles[i]);
                no_redundant_directions_vector.push_back(directions_vector[i]);
                no_redundant_directions_ipix.push_back(directions_ipix[i]);
            }
        }// for i
    }
    else if (symmetry  == pg_SN )
    {//OK
        for (long int i = 0; i < rot_angles.size(); i++)
        {
            if (rot_angles[i] >= -180.*2. / sym_order &&
                rot_angles[i] <=  180.*2. / sym_order &&
                tilt_angles[i] <=    90.
                )
            {
                no_redundant_rot_angles.push_back(rot_angles[i]);
                no_redundant_tilt_angles.push_back(tilt_angles[i]);
                no_redundant_directions_vector.push_back(directions_vector[i]);
                no_redundant_directions_ipix.push_back(directions_ipix[i]);
            }
        }// for i
    }
    else if (symmetry  == pg_DN )
    {
        for (long int i = 0; i < rot_angles.size(); i++)
        {
            if (rot_angles[i] >= -180. / (sym_order) + 90. &&
                rot_angles[i] <=  180. / (sym_order) + 90. &&
                tilt_angles[i] <=    90.
                )
            {
                no_redundant_rot_angles.push_back(rot_angles[i]);
                no_redundant_tilt_angles.push_back(tilt_angles[i]);
                no_redundant_directions_vector.push_back(directions_vector[i]);
                no_redundant_directions_ipix.push_back(directions_ipix[i]);
            }
        }// for i
    }
    else if (symmetry  == pg_DNV )
    {
        for (long int i = 0; i < rot_angles.size(); i++)
        {
            if (rot_angles[i] >=   90.  &&
                rot_angles[i] <=  180. / (sym_order) + 90. &&
                tilt_angles[i] <=    90.
                )
            {
                no_redundant_rot_angles.push_back(rot_angles[i]);
                no_redundant_tilt_angles.push_back(tilt_angles[i]);
                no_redundant_directions_vector.push_back(directions_vector[i]);
                no_redundant_directions_ipix.push_back(directions_ipix[i]);
            }
        }// for i
    }
    else if (symmetry  == pg_DNH )
    {
        for (long int i = 0; i < rot_angles.size(); i++)
        {
            if (rot_angles[i] >=   90. &&
                rot_angles[i] <=  180. / (sym_order) + 90. &&
                tilt_angles[i] <=   90.
                )
            {
                no_redundant_rot_angles.push_back(rot_angles[i]);
                no_redundant_tilt_angles.push_back(tilt_angles[i]);
                no_redundant_directions_vector.push_back(directions_vector[i]);
                no_redundant_directions_ipix.push_back(directions_ipix[i]);
            }
        }// for i
    }
    else if (symmetry  == pg_T )
    {//OK
        Matrix1D<FDOUBLE>  _3_fold_axis_1_by_3_fold_axis_2(3);
        _3_fold_axis_1_by_3_fold_axis_2 = vectorR3((FDOUBLE)-0.942809, (FDOUBLE)0., (FDOUBLE)0.);
        _3_fold_axis_1_by_3_fold_axis_2.selfNormalize();
        Matrix1D<FDOUBLE>  _3_fold_axis_2_by_3_fold_axis_3(3);
        _3_fold_axis_2_by_3_fold_axis_3 = vectorR3((FDOUBLE)0.471405, (FDOUBLE)0.272165, (FDOUBLE)0.7698);
        _3_fold_axis_2_by_3_fold_axis_3.selfNormalize();
        Matrix1D<FDOUBLE>  _3_fold_axis_3_by_3_fold_axis_1(3);
        _3_fold_axis_3_by_3_fold_axis_1 = vectorR3((FDOUBLE)0.471404, (FDOUBLE)0.816497, (FDOUBLE)0.);
        _3_fold_axis_3_by_3_fold_axis_1.selfNormalize();
        for (long int i = 0; i < rot_angles.size(); i++)
        {
            if (rot_angles[i] >=     90. &&
                rot_angles[i] <=   150. ||
                rot_angles[i] ==     0
                )
                if (
                    dotProduct(directions_vector[i], _3_fold_axis_1_by_3_fold_axis_2) >= 0 &&
                    dotProduct(directions_vector[i], _3_fold_axis_2_by_3_fold_axis_3) >= 0 &&
                    dotProduct(directions_vector[i], _3_fold_axis_3_by_3_fold_axis_1) >= 0
                    )
                {
                    no_redundant_rot_angles.push_back(rot_angles[i]);
                    no_redundant_tilt_angles.push_back(tilt_angles[i]);
                    no_redundant_directions_vector.push_back(directions_vector[i]);
                    no_redundant_directions_ipix.push_back(directions_ipix[i]);
                }
        }// for i
    }
    else if (symmetry  == pg_TD )
    {//OK
        Matrix1D<FDOUBLE>  _2_fold_axis_1_by_3_fold_axis_2(3);
        _2_fold_axis_1_by_3_fold_axis_2 = vectorR3((FDOUBLE)-0.942809, (FDOUBLE)0., (FDOUBLE)0.);
        _2_fold_axis_1_by_3_fold_axis_2.selfNormalize();
        Matrix1D<FDOUBLE>  _3_fold_axis_2_by_3_fold_axis_5(3);
        _3_fold_axis_2_by_3_fold_axis_5 = vectorR3((FDOUBLE)0.471405, (FDOUBLE)0.272165, (FDOUBLE)0.7698);
        _3_fold_axis_2_by_3_fold_axis_5.selfNormalize();
        Matrix1D<FDOUBLE>  _3_fold_axis_5_by_2_fold_axis_1(3);
        _3_fold_axis_5_by_2_fold_axis_1 = vectorR3((FDOUBLE)0., (FDOUBLE)0.471405, (FDOUBLE)-0.666667);
        _3_fold_axis_5_by_2_fold_axis_1.selfNormalize();
        for (long int i = 0; i < rot_angles.size(); i++)
        {
            //           if ( rot_angles[i]>=     120. &&
            //                 rot_angles[i]<=   150. ||
            //                 rot_angles[i]==     0
            //              )
            if (
                dotProduct(directions_vector[i], _2_fold_axis_1_by_3_fold_axis_2) >= 0 &&
                dotProduct(directions_vector[i], _3_fold_axis_2_by_3_fold_axis_5) >= 0 &&
                dotProduct(directions_vector[i], _3_fold_axis_5_by_2_fold_axis_1) >= 0
                )
            {
                no_redundant_rot_angles.push_back(rot_angles[i]);
                no_redundant_tilt_angles.push_back(tilt_angles[i]);
                no_redundant_directions_vector.push_back(directions_vector[i]);
                no_redundant_directions_ipix.push_back(directions_ipix[i]);
            }
        }// for i
    }
    else if (symmetry  == pg_TH )
    {//OK
        Matrix1D<FDOUBLE>  _3_fold_axis_1_by_2_fold_axis_1(3);
        _3_fold_axis_1_by_2_fold_axis_1 = vectorR3((FDOUBLE)-0.816496, (FDOUBLE)0., (FDOUBLE)0.);
        _3_fold_axis_1_by_2_fold_axis_1.selfNormalize();
        Matrix1D<FDOUBLE>  _2_fold_axis_1_by_2_fold_axis_2(3);
        _2_fold_axis_1_by_2_fold_axis_2 = vectorR3((FDOUBLE)0.707107, (FDOUBLE)0.408248, (FDOUBLE)-0.57735);
        _2_fold_axis_1_by_2_fold_axis_2.selfNormalize();
        Matrix1D<FDOUBLE>  _2_fold_axis_2_by_3_fold_axis_1(3);
        _2_fold_axis_2_by_3_fold_axis_1 = vectorR3((FDOUBLE)-0.408248, (FDOUBLE)-0.707107, (FDOUBLE)0.);
        _2_fold_axis_2_by_3_fold_axis_1.selfNormalize();
        for (long int i = 0; i < rot_angles.size(); i++)
        {
            //           if ( rot_angles[i]>=     120. &&
            //                 rot_angles[i]<=   150. ||
            //                 rot_angles[i]==     0
            //              )
            if (
                dotProduct(directions_vector[i], _3_fold_axis_1_by_2_fold_axis_1) >= 0 &&
                dotProduct(directions_vector[i], _2_fold_axis_1_by_2_fold_axis_2) >= 0 &&
                dotProduct(directions_vector[i], _2_fold_axis_2_by_3_fold_axis_1) >= 0
                )
            {
                no_redundant_rot_angles.push_back(rot_angles[i]);
                no_redundant_tilt_angles.push_back(tilt_angles[i]);
                no_redundant_directions_vector.push_back(directions_vector[i]);
                no_redundant_directions_ipix.push_back(directions_ipix[i]);
            }
        }// for i
    }
    else if (symmetry  == pg_O )
    {//OK
        Matrix1D<FDOUBLE>  _3_fold_axis_1_by_3_fold_axis_2(3);
        _3_fold_axis_1_by_3_fold_axis_2 = vectorR3((FDOUBLE)0., (FDOUBLE)-1., (FDOUBLE)1.);
        _3_fold_axis_1_by_3_fold_axis_2.selfNormalize();
        Matrix1D<FDOUBLE>  _3_fold_axis_2_by_4_fold_axis(3);
        _3_fold_axis_2_by_4_fold_axis = vectorR3((FDOUBLE)1., (FDOUBLE)1., (FDOUBLE)0.);
        _3_fold_axis_2_by_4_fold_axis.selfNormalize();
        Matrix1D<FDOUBLE>  _4_fold_axis_by_3_fold_axis_1(3);
        _4_fold_axis_by_3_fold_axis_1 = vectorR3((FDOUBLE)-1., (FDOUBLE)1., (FDOUBLE)0.);
        _4_fold_axis_by_3_fold_axis_1.selfNormalize();
        for (long int i = 0; i < rot_angles.size(); i++)
        {
            if ((rot_angles[i] >=   45. &&
                 rot_angles[i] <=  135. &&
                 tilt_angles[i] <=  90.) ||
                rot_angles[i] ==  0.
                )
                if (
                    dotProduct(directions_vector[i], _3_fold_axis_1_by_3_fold_axis_2) >= 0 &&
                    dotProduct(directions_vector[i], _3_fold_axis_2_by_4_fold_axis) >= 0 &&
                    dotProduct(directions_vector[i], _4_fold_axis_by_3_fold_axis_1) >= 0
                    )
                {
                    no_redundant_rot_angles.push_back(rot_angles[i]);
                    no_redundant_tilt_angles.push_back(tilt_angles[i]);
                    no_redundant_directions_vector.push_back(directions_vector[i]);
                    no_redundant_directions_ipix.push_back(directions_ipix[i]);
                }
        }// for i
    }
    else if (symmetry  == pg_OH )
    {//OK
        Matrix1D<FDOUBLE>  _3_fold_axis_1_by_3_fold_axis_2(3);
        _3_fold_axis_1_by_3_fold_axis_2 = vectorR3((FDOUBLE)0., (FDOUBLE)-1., (FDOUBLE)1.);
        _3_fold_axis_1_by_3_fold_axis_2.selfNormalize();
        Matrix1D<FDOUBLE>  _3_fold_axis_2_by_4_fold_axis(3);
        _3_fold_axis_2_by_4_fold_axis = vectorR3((FDOUBLE)1., (FDOUBLE)1., (FDOUBLE)0.);
        _3_fold_axis_2_by_4_fold_axis.selfNormalize();
        Matrix1D<FDOUBLE>  _4_fold_axis_by_3_fold_axis_1(3);
        _4_fold_axis_by_3_fold_axis_1 = vectorR3((FDOUBLE)-1., (FDOUBLE)1., (FDOUBLE)0.);
        _4_fold_axis_by_3_fold_axis_1.selfNormalize();
        for (long int i = 0; i < rot_angles.size(); i++)
        {
            if (rot_angles[i] >=   90. &&
                rot_angles[i] <=  135. &&
                tilt_angles[i] <=  90.)
                if (
                    dotProduct(directions_vector[i], _3_fold_axis_1_by_3_fold_axis_2) >= 0 &&
                    dotProduct(directions_vector[i], _3_fold_axis_2_by_4_fold_axis) >= 0 &&
                    dotProduct(directions_vector[i], _4_fold_axis_by_3_fold_axis_1) >= 0
                    )
                {
                    no_redundant_rot_angles.push_back(rot_angles[i]);
                    no_redundant_tilt_angles.push_back(tilt_angles[i]);
                    no_redundant_directions_vector.push_back(directions_vector[i]);
                    no_redundant_directions_ipix.push_back(directions_ipix[i]);
                }
        }// for i
    }
    else if (symmetry  == pg_I || symmetry  == pg_I2)
    {//OK
        Matrix1D<FDOUBLE>  _5_fold_axis_1_by_5_fold_axis_2(3);
        _5_fold_axis_1_by_5_fold_axis_2 = vectorR3((FDOUBLE)0., (FDOUBLE)1., (FDOUBLE)0.);
        _5_fold_axis_1_by_5_fold_axis_2.selfNormalize();
        Matrix1D<FDOUBLE>  _5_fold_axis_2_by_3_fold_axis(3);
        _5_fold_axis_2_by_3_fold_axis = vectorR3((FDOUBLE)-0.4999999839058737,
                                                 (FDOUBLE)-0.8090170074556163,
                                                 (FDOUBLE)0.3090169861701543);
        _5_fold_axis_2_by_3_fold_axis.selfNormalize();
        Matrix1D<FDOUBLE>  _3_fold_axis_by_5_fold_axis_1(3);
        _3_fold_axis_by_5_fold_axis_1 = vectorR3((FDOUBLE)0.4999999839058737,
                                                 (FDOUBLE)-0.8090170074556163,
                                                 (FDOUBLE)0.3090169861701543);
        _3_fold_axis_by_5_fold_axis_1.selfNormalize();
        for (long int i = 0; i < rot_angles.size(); i++)
        {
            if (
                dotProduct(directions_vector[i], _5_fold_axis_1_by_5_fold_axis_2) >= 0 &&
                dotProduct(directions_vector[i], _5_fold_axis_2_by_3_fold_axis) >= 0 &&
                dotProduct(directions_vector[i], _3_fold_axis_by_5_fold_axis_1) >= 0
                )
            {
                no_redundant_rot_angles.push_back(rot_angles[i]);
                no_redundant_tilt_angles.push_back(tilt_angles[i]);
                no_redundant_directions_vector.push_back(directions_vector[i]);
                no_redundant_directions_ipix.push_back(directions_ipix[i]);
            }
        }// for i
    }
    else if (symmetry  == pg_I1)
    {//OK
        FDOUBLE _A[3][3];
        Euler_angles2matrix(0, 90, 0, _A);
        Matrix2D<FDOUBLE>  A(_A);
        Matrix1D<FDOUBLE>  _5_fold_axis_1_by_5_fold_axis_2(3);
        _5_fold_axis_1_by_5_fold_axis_2 = A * vectorR3((FDOUBLE)0., (FDOUBLE)1., (FDOUBLE)0.);
        _5_fold_axis_1_by_5_fold_axis_2.selfNormalize();
        Matrix1D<FDOUBLE>  _5_fold_axis_2_by_3_fold_axis(3);
        _5_fold_axis_2_by_3_fold_axis = A * vectorR3((FDOUBLE)-0.4999999839058737,
                                                     (FDOUBLE)-0.8090170074556163,
                                                     (FDOUBLE)0.3090169861701543);
        _5_fold_axis_2_by_3_fold_axis.selfNormalize();
        Matrix1D<FDOUBLE>  _3_fold_axis_by_5_fold_axis_1(3);// TODO : check this??double to float????
        _3_fold_axis_by_5_fold_axis_1 = A * vectorR3((FDOUBLE)0.4999999839058737,
                                                     (FDOUBLE)-0.8090170074556163,
                                                     (FDOUBLE)0.3090169861701543);
        _3_fold_axis_by_5_fold_axis_1.selfNormalize();
        
        for (long int i = 0; i < rot_angles.size(); i++)
        {
            if (
                dotProduct(directions_vector[i], _5_fold_axis_1_by_5_fold_axis_2) >= 0 &&
                dotProduct(directions_vector[i], _5_fold_axis_2_by_3_fold_axis) >= 0 &&
                dotProduct(directions_vector[i], _3_fold_axis_by_5_fold_axis_1) >= 0
                )
            {
                no_redundant_rot_angles.push_back(rot_angles[i]);
                no_redundant_tilt_angles.push_back(tilt_angles[i]);
                no_redundant_directions_vector.push_back(directions_vector[i]);
                no_redundant_directions_ipix.push_back(directions_ipix[i]);
            }
        }// for i
    }
    else if (symmetry  == pg_I3)
    {//OK
        FDOUBLE _A[3][3];
        Euler_angles2matrix(0,31.7174745559,0, _A);
        Matrix2D<FDOUBLE>  A(_A);
        Matrix1D<FDOUBLE>  _5_fold_axis_1_by_5_fold_axis_2(3);
        _5_fold_axis_1_by_5_fold_axis_2 = A * vectorR3((FDOUBLE)0., (FDOUBLE)1., (FDOUBLE)0.);
        _5_fold_axis_1_by_5_fold_axis_2.selfNormalize();
        Matrix1D<FDOUBLE>  _5_fold_axis_2_by_3_fold_axis(3);
        _5_fold_axis_2_by_3_fold_axis = A * vectorR3((FDOUBLE)-0.4999999839058737,
                                                     (FDOUBLE)-0.8090170074556163,
                                                     (FDOUBLE)0.3090169861701543);
        _5_fold_axis_2_by_3_fold_axis.selfNormalize();
        Matrix1D<FDOUBLE>  _3_fold_axis_by_5_fold_axis_1(3);
        _3_fold_axis_by_5_fold_axis_1 = A * vectorR3((FDOUBLE)0.4999999839058737,
                                                     (FDOUBLE)-0.8090170074556163,
                                                     (FDOUBLE)0.3090169861701543);
        _3_fold_axis_by_5_fold_axis_1.selfNormalize();
        
        for (long int i = 0; i < rot_angles.size(); i++)
        {
            if (
                dotProduct(directions_vector[i], _5_fold_axis_1_by_5_fold_axis_2) >= 0 &&
                dotProduct(directions_vector[i], _5_fold_axis_2_by_3_fold_axis) >= 0 &&
                dotProduct(directions_vector[i], _3_fold_axis_by_5_fold_axis_1) >= 0
                )
            {
                no_redundant_rot_angles.push_back(rot_angles[i]);
                no_redundant_tilt_angles.push_back(tilt_angles[i]);
                no_redundant_directions_vector.push_back(directions_vector[i]);
                no_redundant_directions_ipix.push_back(directions_ipix[i]);
            }
        }// for i
    }
    else if (symmetry  == pg_I4)
    {//OK
        FDOUBLE _A[3][3];
        Euler_angles2matrix(0,-31.7174745559,0, _A);
        Matrix2D<FDOUBLE>  A(_A);
        Matrix1D<FDOUBLE>  _5_fold_axis_1_by_5_fold_axis_2(3);
        _5_fold_axis_1_by_5_fold_axis_2 = A * vectorR3((FDOUBLE)0., (FDOUBLE)0., (FDOUBLE)1.);
        _5_fold_axis_1_by_5_fold_axis_2.selfNormalize();
        Matrix1D<FDOUBLE>  _5_fold_axis_2_by_3_fold_axis(3);
        _5_fold_axis_2_by_3_fold_axis = A * vectorR3((FDOUBLE)0.187592467856686,
                                                     (FDOUBLE)-0.303530987314591,
                                                     (FDOUBLE)-0.491123477863004);
        _5_fold_axis_2_by_3_fold_axis.selfNormalize();
        Matrix1D<FDOUBLE>  _3_fold_axis_by_5_fold_axis_1(3);
        _3_fold_axis_by_5_fold_axis_1 = A * vectorR3((FDOUBLE)0.187592467856686,
                                                     (FDOUBLE)0.303530987314591,
                                                     (FDOUBLE)-0.491123477863004);
        _3_fold_axis_by_5_fold_axis_1.selfNormalize();
        
        for (long int i = 0; i < rot_angles.size(); i++)
        {
            if (
                dotProduct(directions_vector[i],
                           _5_fold_axis_2_by_3_fold_axis)   <= 0 &&
                dotProduct(directions_vector[i],
                           _3_fold_axis_by_5_fold_axis_1)   <= 0 &&
                dotProduct(directions_vector[i],
                           _5_fold_axis_1_by_5_fold_axis_2) <= 0
                )
            {
                no_redundant_rot_angles.push_back(rot_angles[i]);
                no_redundant_tilt_angles.push_back(tilt_angles[i]);
                no_redundant_directions_vector.push_back(directions_vector[i]);
                no_redundant_directions_ipix.push_back(directions_ipix[i]);
            }
        }// for i
    }
    else if (symmetry  == pg_I5)
    {//OK
        std::cerr << "ERROR: Symmetry pg_I5 not implemented" << std::endl;
		EXIT_ABNORMALLY;
    }
    else if (symmetry  == pg_IH || symmetry  == pg_I2H)
    {//OK
        Matrix1D<FDOUBLE>  _5_fold_axis_1_by_5_fold_axis_2(3);
        _5_fold_axis_1_by_5_fold_axis_2 = vectorR3((FDOUBLE)0., (FDOUBLE)1., (FDOUBLE)0.);
        _5_fold_axis_1_by_5_fold_axis_2.selfNormalize();
        Matrix1D<FDOUBLE>  _5_fold_axis_2_by_3_fold_axis(3);
        _5_fold_axis_2_by_3_fold_axis = vectorR3((FDOUBLE)-0.4999999839058737,
                                                 (FDOUBLE)-0.8090170074556163,
                                                 (FDOUBLE)0.3090169861701543);
        _5_fold_axis_2_by_3_fold_axis.selfNormalize();
        Matrix1D<FDOUBLE>  _3_fold_axis_by_5_fold_axis_1(3);
        _3_fold_axis_by_5_fold_axis_1 = vectorR3((FDOUBLE)0.4999999839058737,
                                                 (FDOUBLE)-0.8090170074556163,
                                                 (FDOUBLE)0.3090169861701543);
        _3_fold_axis_by_5_fold_axis_1.selfNormalize();
        Matrix1D<FDOUBLE>  _3_fold_axis_by_2_fold_axis(3);
        _3_fold_axis_by_2_fold_axis =  vectorR3((FDOUBLE)1.,(FDOUBLE)0.,(FDOUBLE)0.);
        _3_fold_axis_by_2_fold_axis.selfNormalize();
        for (long int i = 0; i < rot_angles.size(); i++)
        {
            if (
                dotProduct(directions_vector[i], _5_fold_axis_1_by_5_fold_axis_2) >= 0 &&
                dotProduct(directions_vector[i], _5_fold_axis_2_by_3_fold_axis) >= 0 &&
                dotProduct(directions_vector[i], _3_fold_axis_by_2_fold_axis) >= 0
                )
            {
                no_redundant_rot_angles.push_back(rot_angles[i]);
                no_redundant_tilt_angles.push_back(tilt_angles[i]);
                no_redundant_directions_vector.push_back(directions_vector[i]);
                no_redundant_directions_ipix.push_back(directions_ipix[i]);
            }
        }// for i
    }
    else if (symmetry  == pg_I1H)
    {//OK
        FDOUBLE _A[3][3];
        Euler_angles2matrix(0, 90, 0, _A);
        Matrix2D<FDOUBLE>  A(_A);
        Matrix1D<FDOUBLE>  _5_fold_axis_1_by_5_fold_axis_2(3);
        _5_fold_axis_1_by_5_fold_axis_2 = A * vectorR3((FDOUBLE)0., (FDOUBLE)1., (FDOUBLE)0.);
        _5_fold_axis_1_by_5_fold_axis_2.selfNormalize();
        Matrix1D<FDOUBLE>  _5_fold_axis_2_by_3_fold_axis(3);
        _5_fold_axis_2_by_3_fold_axis = A * vectorR3((FDOUBLE)-0.4999999839058737,
                                                     (FDOUBLE)-0.8090170074556163,
                                                     (FDOUBLE)0.3090169861701543);
        _5_fold_axis_2_by_3_fold_axis.selfNormalize();
        Matrix1D<FDOUBLE>  _3_fold_axis_by_5_fold_axis_1(3);
        _3_fold_axis_by_5_fold_axis_1 = A * vectorR3((FDOUBLE)0.4999999839058737,
                                                     (FDOUBLE)-0.8090170074556163,
                                                     (FDOUBLE)0.3090169861701543);
        _3_fold_axis_by_5_fold_axis_1.selfNormalize();
        Matrix1D<FDOUBLE>  _3_fold_axis_by_2_fold_axis(3);
        _3_fold_axis_by_2_fold_axis =  A * vectorR3((FDOUBLE)1.,(FDOUBLE)0.,(FDOUBLE)0.);
        _3_fold_axis_by_2_fold_axis.selfNormalize();
        for (long int i = 0; i < rot_angles.size(); i++)
        {
            if (
                dotProduct(directions_vector[i], _5_fold_axis_1_by_5_fold_axis_2) >= 0 &&
                dotProduct(directions_vector[i], _5_fold_axis_2_by_3_fold_axis) >= 0 &&
                dotProduct(directions_vector[i], _3_fold_axis_by_2_fold_axis) >= 0
                )
            {
                no_redundant_rot_angles.push_back(rot_angles[i]);
                no_redundant_tilt_angles.push_back(tilt_angles[i]);
                no_redundant_directions_vector.push_back(directions_vector[i]);
                no_redundant_directions_ipix.push_back(directions_ipix[i]);
            }
        }// for i
    }
    else if (symmetry  == pg_I3H)
    {//OK
        FDOUBLE _A[3][3];
        Euler_angles2matrix(0,31.7174745559,0, _A);
        Matrix2D<FDOUBLE>  A(_A);
        Matrix1D<FDOUBLE>  _5_fold_axis_1_by_5_fold_axis_2(3);
        _5_fold_axis_1_by_5_fold_axis_2 = A * vectorR3((FDOUBLE)0., (FDOUBLE)0., (FDOUBLE)1.);
        _5_fold_axis_1_by_5_fold_axis_2.selfNormalize();
        Matrix1D<FDOUBLE>  _5_fold_axis_2_by_3_fold_axis(3);
        _5_fold_axis_2_by_3_fold_axis = A * vectorR3((FDOUBLE)0.187592467856686,
                                                     (FDOUBLE)-0.303530987314591,
                                                     (FDOUBLE)-0.491123477863004);
        _5_fold_axis_2_by_3_fold_axis.selfNormalize();
        Matrix1D<FDOUBLE>  _3_fold_axis_by_5_fold_axis_1(3);
        _3_fold_axis_by_5_fold_axis_1 = A * vectorR3((FDOUBLE)0.187592467856686,
                                                     (FDOUBLE)0.303530987314591,
                                                     (FDOUBLE)-0.491123477863004);
        _3_fold_axis_by_5_fold_axis_1.selfNormalize();
        Matrix1D<FDOUBLE>  _3_fold_axis_by_2_fold_axis(3);
        _3_fold_axis_by_2_fold_axis = vectorR3((FDOUBLE)0.,(FDOUBLE)1.,(FDOUBLE)0.);
        _3_fold_axis_by_2_fold_axis.selfNormalize();
        for (long int i = 0; i < rot_angles.size(); i++)
        {
            if (
                dotProduct(directions_vector[i],
                           _5_fold_axis_2_by_3_fold_axis)   >= 0 &&
                dotProduct(directions_vector[i],
                           _3_fold_axis_by_5_fold_axis_1)   >= 0 &&
                dotProduct(directions_vector[i],
                           _5_fold_axis_1_by_5_fold_axis_2) >= 0 &&
                dotProduct(directions_vector[i],
                           _3_fold_axis_by_2_fold_axis)     >= 0
                )
            {
                no_redundant_rot_angles.push_back(rot_angles[i]);
                no_redundant_tilt_angles.push_back(tilt_angles[i]);
                no_redundant_directions_vector.push_back(directions_vector[i]);
                no_redundant_directions_ipix.push_back(directions_ipix[i]);
            }
        }// for i
    }
    else if (symmetry  == pg_I4H)
    {//OK
        FDOUBLE _A[3][3];
        Euler_angles2matrix(0,-31.7174745559,0, _A);
        Matrix2D<FDOUBLE>  A(_A);
        Matrix1D<FDOUBLE>  _5_fold_axis_1_by_5_fold_axis_2(3);
        _5_fold_axis_1_by_5_fold_axis_2 = A * vectorR3((FDOUBLE)0., (FDOUBLE)0., (FDOUBLE)1.);
        _5_fold_axis_1_by_5_fold_axis_2.selfNormalize();
        Matrix1D<FDOUBLE>  _5_fold_axis_2_by_3_fold_axis(3);
        _5_fold_axis_2_by_3_fold_axis = A * vectorR3((FDOUBLE)0.187592467856686,
                                                     (FDOUBLE)-0.303530987314591,
                                                     (FDOUBLE)-0.491123477863004);
        _5_fold_axis_2_by_3_fold_axis.selfNormalize();
        Matrix1D<FDOUBLE>  _3_fold_axis_by_5_fold_axis_1(3);
        _3_fold_axis_by_5_fold_axis_1 = A * vectorR3((FDOUBLE)0.187592467856686,
                                                     (FDOUBLE)0.303530987314591,
                                                     (FDOUBLE)-0.491123477863004);
        _3_fold_axis_by_5_fold_axis_1.selfNormalize();
        Matrix1D<FDOUBLE>  _3_fold_axis_by_2_fold_axis(3);
        _3_fold_axis_by_2_fold_axis = vectorR3((FDOUBLE)0.,(FDOUBLE)1.,(FDOUBLE)0.);
        _3_fold_axis_by_2_fold_axis.selfNormalize();
        for (long int i = 0; i < rot_angles.size(); i++)
        {
            if (
                dotProduct(directions_vector[i],
                           _5_fold_axis_2_by_3_fold_axis)   <= 0 &&
                dotProduct(directions_vector[i],
                           _3_fold_axis_by_5_fold_axis_1)   <= 0 &&
                dotProduct(directions_vector[i],
                           _5_fold_axis_1_by_5_fold_axis_2) <= 0 &&
                dotProduct(directions_vector[i],
                           _3_fold_axis_by_2_fold_axis)     >= 0
                )
            {
                no_redundant_rot_angles.push_back(rot_angles[i]);
                no_redundant_tilt_angles.push_back(tilt_angles[i]);
                no_redundant_directions_vector.push_back(directions_vector[i]);
                no_redundant_directions_ipix.push_back(directions_ipix[i]);
            }
        }// for i
    }
    else if (symmetry  == pg_I5H)
    {//OK
        std::cerr << "ERROR: pg_I5H Symmetry not implemented" << std::endl;
		EXIT_ABNORMALLY;
    }
    else
    {
        std::cerr << "ERROR: Symmetry " << symmetry  << "is not known" << std::endl;
		EXIT_ABNORMALLY;
    }
    
    
    // Now overwrite the rot/tilt_angles and directions_vectors with their no_redundant counterparts
    rot_angles = no_redundant_rot_angles;
    tilt_angles = no_redundant_tilt_angles;
    directions_vector = no_redundant_directions_vector;
    directions_ipix = no_redundant_directions_ipix;
    
    
}

