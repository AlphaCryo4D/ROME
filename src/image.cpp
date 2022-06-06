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
 /***************************************************************************
 *
 * Authors: "Yongbei(Glow) Ma,Jiayi (Timmy) Wu, Youdong (Jack) Mao"
 * Dana-Farber Cancer Institute, Harvard Medical School and Peking University
 * Bevin Brett
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

#include "image.h"

bool ShiftImageInFourierTransformNew_useSinCosTable = true;

namespace TabulatedSinCos {
	Table::Table() {
		for (int i = 0; i < nr_element; i++) {
			double xx = double(i) * radiansPerEntry;
			tabulated[i].sin = sin(xx);
			tabulated[i].cos = cos(xx);
		}
#ifdef L2_CACHE_MODELING
		epoch = 1;
		for (int i = 0; i < nr_element; i++)
			mostRecentAccess[i] = 0;
#endif
	}
	int Table::initNotingSeqAcc() const {
		int result = 0;
#ifdef L2_CACHE_MODELING
		result = ++epoch;
#endif
		return result;
	}
	void Table::finiNotingSeqAcc(int initResult) const {
#ifdef L2_CACHE_MODELING
		// It might be worth compressing these to the cache line boundary
		// and also recognizing strided (or at least unit stride) accesses
		//
		for (int i = 0; i < nr_element; i++) {
			int mra = mostRecentAccess[i].load();
			if (mra >= initResult)
				L2CacheModel::seqAcc("TabulatedSinCos table", 0, &mostRecentAccess[i], 1);
		}
#endif
	}
	const Table table;
}

IntPerformanceCounter ShiftImageInFourierTransformNew_performanceCounter("ShiftImageInFourierTransformNew_performanceCounter");
IntPerformanceCounter ShiftImageInFourierTransformNew_init_performanceCounter("ShiftImageInFourierTransformNew_init_performanceCounter");

IntPerformanceCounter ShiftImageInFourierTransform1_performanceCounter("ShiftImageInFourierTransform1_performanceCounter");
IntPerformanceCounter ShiftImageInFourierTransform2_performanceCounter("ShiftImageInFourierTransform2_performanceCounter");
IntPerformanceCounter ShiftImageInFourierTransform4_performanceCounter("ShiftImageInFourierTransform4_performanceCounter");

IntPerformanceCounter ShiftImageInFourierTransformSimple_performanceCounter("ShiftImageInFourierTransformSimple_performanceCounter");


void windowFourierTransform(SOAComplexReadonly& in,int in_size,SOAComplexDouble& out,int out_size)
{
    int in_Fsize = in_size/2 + 1;
    int out_Fsize = out_size/2 + 1;
    int out_Fsize2 = out_Fsize*out_size;
    // Check size of the input array
    if (in_size == out_size)
    {
        for (int i = 0; i < out_Fsize2; i++) {
            out.real[i] = in.real[i];
            out.imag[i] = in.imag[i];
        }
        return;
    }
    
    
    for (int i = 0; i < out_Fsize2; i++) {
        out.real[i] = 0;
        out.imag[i] = 0;
    }
    
    if (out_Fsize > in_Fsize)
    {
        int max_r2 = (in_Fsize -1) * (in_Fsize - 1);
        
        for (int i = 0; i<in_size; i++){
            
            int ip = (i < in_Fsize) ? i : i - in_size;
            int in_i = ((ip < 0) ? (ip + in_size) : (ip));
            int out_i = ((ip < 0) ? (ip + out_size) : (ip));
            // max_r2 should be larger than ip*ip
            int jp2 = max_r2-ip*ip;
            if (jp2 < 0) continue;
            int j_max = sqrt(jp2)+1;
            j_max = j_max < in_Fsize?j_max:in_Fsize;
            
#pragma ivdep
            for (int j = 0; j<j_max; j++)
            {
                // Make sure windowed FT has nothing in the corners, otherwise we end up with an asymmetric FT!
                int jp = j;
                int in_j = jp;
                int out_j = jp;
                
                // if (ip*ip + jp*jp <= max_r2)
                out.real[out_i*out_Fsize+out_j] = in.real[in_i*in_Fsize+in_j];
                out.imag[out_i*out_Fsize+out_j] = in.imag[in_i*in_Fsize+in_j];
            }
        }
    }
    else
    {
        for ( int i = 0; i<out_size; i++){
            
            int ip = (i < out_Fsize) ? i : i - out_size;
            int in_i = ((ip < 0) ? (ip + in_size) : (ip));
            int out_i = ((ip < 0) ? (ip + out_size) : (ip));

#pragma ivdep
            for ( int j = 0; j<out_Fsize; j++)
            {
                
                int jp = j;
                int in_j = jp;
                int out_j = jp;
                
                out.real[out_i*out_Fsize+out_j] = in.real[in_i*in_Fsize+in_j];
                out.imag[out_i*out_Fsize+out_j] = in.imag[in_i*in_Fsize+in_j];
            }
        }
    }
}


void getSpectrum(const double* Min,int Min_size,double* spectrum,int spectrum_type)
{
    int Min_Fsize = Min_size/2 + 1;
    int Min_Fsize2 = Min_size*(Min_size/2+1);
    
    SOAComplexDouble Faux;
    Faux.real = (double*)aMalloc(sizeof(double)*Min_Fsize2,64);
    Faux.imag = (double*)aMalloc(sizeof(double)*Min_Fsize2,64);
    double* count = (double*)aMalloc(sizeof(double)*Min_size,64);
    
    for (int i = 0; i < Min_size; i++) {
        count[i] = spectrum[i] = 0.;
    }
    
    FFTWTransformer transformer(Min_size,Min_size);
    
    transformer.FourierTransform(Min, Faux);
    
    for (int i = 0; i<Min_size; i++){
        
        int ip = (i < Min_Fsize) ? i : i - Min_size;
        
        for (int j = 0, jp = 0; j< Min_Fsize; j++, jp = j)
        {
            int idx = round(sqrt(ip*ip + jp*jp));
            
            if (spectrum_type == AMPLITUDE_SPECTRUM)
                spectrum[idx] += sqrt(Faux.real[i*Min_Fsize+j]*Faux.real[i*Min_Fsize+j]+Faux.imag[i*Min_Fsize+j]*Faux.imag[i*Min_Fsize+j]);
            else
                spectrum[idx] += Faux.real[i*Min_Fsize+j]*Faux.real[i*Min_Fsize+j]+Faux.imag[i*Min_Fsize+j]*Faux.imag[i*Min_Fsize+j];
            
            count[idx] += 1.;
        }
    }
    
    for (long int i = 0; i < Min_size; i++)
        if (count[i] > 0.)
            spectrum[i] /= count[i];
    
    aFree(count);
    aFree(Faux.real);
    aFree(Faux.imag);
}


void getSpectrum(const float* Min,int Min_size,float* spectrum,int spectrum_type)
{
    int Min_Fsize = Min_size/2 + 1;
    int Min_Fsize2 = Min_size*(Min_size/2+1);
    
    SOAComplexFloat Faux;
    Faux.real = (float*)aMalloc(sizeof(float)*Min_Fsize2,64);
    Faux.imag = (float*)aMalloc(sizeof(float)*Min_Fsize2,64);
    float* count = (float*)aMalloc(sizeof(float)*Min_size,64);
    
    for (int i = 0; i < Min_size; i++) {
        count[i] = spectrum[i] = 0.;
    }
    
    FFTWFTransformer transformer(Min_size,Min_size);
    
    transformer.FourierTransform(Min, Faux);
    
    for (int i = 0; i < Min_size; i++){
        
        int ip = (i < Min_Fsize) ? i : i - Min_size;
        
        for (int j = 0, jp = 0; j< Min_Fsize; j++, jp = j)
        {
            int idx = round(sqrt(ip*ip + jp*jp));
            
            if (spectrum_type == AMPLITUDE_SPECTRUM)
                spectrum[idx] += sqrt(Faux.real[i*Min_Fsize+j]*Faux.real[i*Min_Fsize+j]+Faux.imag[i*Min_Fsize+j]*Faux.imag[i*Min_Fsize+j]);
            else
                spectrum[idx] += Faux.real[i*Min_Fsize+j]*Faux.real[i*Min_Fsize+j]+Faux.imag[i*Min_Fsize+j]*Faux.imag[i*Min_Fsize+j];
            
            count[idx] += 1.;
        }
    }
    
    for (long int i = 0; i < Min_size; i++)
        if (count[i] > 0.)
            spectrum[i] /= count[i];
    
    aFree(count);
    aFree(Faux.real);
    aFree(Faux.imag);
}


void inverse(double A[][3])
{
    double Aref[3][3];
    
    Aref[0][0] =   A[2][2]*A[1][1]-A[2][1]*A[1][2];
    Aref[0][1] = -(A[2][2]*A[0][1]-A[2][1]*A[0][2]);
    Aref[0][2] =   A[1][2]*A[0][1]-A[1][1]*A[0][2];
    Aref[1][0] = -(A[2][2]*A[1][0]-A[2][0]*A[1][2]);
    Aref[1][1] =   A[2][2]*A[0][0]-A[2][0]*A[0][2];
    Aref[1][2] = -(A[1][2]*A[0][0]-A[1][0]*A[0][2]);
    Aref[2][0] =   A[2][1]*A[1][0]-A[2][0]*A[1][1];
    Aref[2][1] = -(A[2][1]*A[0][0]-A[2][0]*A[0][1]);
    Aref[2][2] =   A[1][1]*A[0][0]-A[1][0]*A[0][1];
    double tmp = A[0][0] * Aref[0][0]+A[1][0] * Aref[0][1] + A[2][0] * Aref[0][2];
    for (int i = 0; i < 3; i++)
        for (int j = 0; j < 3; j++)
            Aref[i][j] = Aref[i][j]/tmp;
    
    for (int i = 0; i < 3; i++)
        for (int j = 0; j < 3; j++)
            A[i][j] = Aref[i][j];
}

// template <typename T>
inline void decenter(const SOAComplexDouble &Min,int Min_size, SOAComplexDouble &Mout, int Mout_size,int my_rmax2)
{
    
    // Mout should already have the right size
    // Initialize to zero,because we not sure the data type,so initialize zero outside this function!!!!!!!
    
    int Min_size_shift = -(int)((float) (Min_size) / 2.0);
    int Min_Fsize = Min_size/2+1;
    int Mout_Fsize = Mout_size/2+1;

    for (int i = 0; i<Mout_size; i++){
        
        int ip = (i < Mout_Fsize) ? i : i - Mout_size;
        //my_rmax2 should be larger than ip*ip
        int jp2 = my_rmax2-ip*ip;
        if (jp2 < 0) continue;
        int j_max = sqrt(jp2)+1;
        j_max = j_max<Mout_Fsize?j_max:Mout_Fsize;
        ip = ip - Min_size_shift;
        
#pragma ivdep
        for (int j = 0; j<j_max; j++)
        {
            int jp = j;
            //if (ip*ip + jp*jp <= my_rmax2)
            Mout.real[i*Mout_Fsize+j] = Min.real[ip*Min_Fsize+jp];//only y dim shift
            Mout.imag[i*Mout_Fsize+j] = Min.imag[ip*Min_Fsize+jp];
        }
    }
}

void shiftImageInFourierTransformUnitTestCorrectness() {
	std::cout << "shiftImageInFourierTransformUnitTestCorrectness" << std::endl;
	static const int inout_size = 20;
    static const int inout_Fsize2 = inout_size*(inout_size/2+1);

	ShiftImageInFourierTransform<float,float> t(inout_size, 0.5, 0.5, 1.0);

	class V { public: std::vector<float> v; V() : v(inout_Fsize2) {} };
	auto compare = [&](V&lhs, V&rhs) {
		for (int i = 0; i < inout_Fsize2; i++) {
			auto l = lhs.v[i];
			auto r = rhs.v[i];
			auto m = std::max(std::abs(l), std::abs(r));
			if (m == 0) continue;
			if (std::abs(l-r) * 10 < m) continue;
			std::cerr << "shiftImageInFourierTransform got wrong answer" << std::endl;
		}
	};

	V vfin_real    [4];
	V vfin_imag    [4];
	V vfout_realOld[4];
	V vfout_imagOld[4];
	V vfout_realNew[4];
	V vfout_imagNew[4];
	for (int j = 0; j < 4; j++) {
		V & fin_real     = vfin_real    [j];
		V & fin_imag     = vfin_imag    [j];
		V & fout_realOld = vfout_realOld[j];
		V & fout_imagOld = vfout_imagOld[j];
		V & fout_realNew = vfout_realNew[j];
		V & fout_imagNew = vfout_imagNew[j];
		for (int i = 0; i < inout_Fsize2; i++) fin_real.v[i] = fin_imag.v[i] = float(i*j);

		shiftImageInFourierTransform(fin_real.v, fin_imag.v, fout_realOld.v, fout_imagOld.v, inout_size, 0.5, 0.5, 1.0);
		t.transform(fin_real.v.data(), fin_imag.v.data(), fout_realNew.v.data(), fout_imagNew.v.data());
		compare(fout_realOld, fout_realNew);
		compare(fout_imagOld, fout_imagNew);
	}

	t.transform4(
		vfin_real[0].v.data(), vfin_imag[0].v.data(), vfout_realNew[0].v.data(), vfout_imagNew[0].v.data(),
		vfin_real[1].v.data(), vfin_imag[1].v.data(), vfout_realNew[1].v.data(), vfout_imagNew[1].v.data(),
		vfin_real[2].v.data(), vfin_imag[2].v.data(), vfout_realNew[2].v.data(), vfout_imagNew[2].v.data(),
		vfin_real[3].v.data(), vfin_imag[3].v.data(), vfout_realNew[3].v.data(), vfout_imagNew[3].v.data());

	for (int j = 0; j < 4; j++) {
		V & fin_real     = vfin_real    [j];
		V & fin_imag     = vfin_imag    [j];
		V & fout_realOld = vfout_realOld[j];
		V & fout_imagOld = vfout_imagOld[j];
		V & fout_realNew = vfout_realNew[j];
		V & fout_imagNew = vfout_imagNew[j];
		compare(fout_realOld, fout_realNew);
		compare(fout_imagOld, fout_imagNew);
	}

	std::cout << "shiftImageInFourierTransformUnitTestCorrectness end" << std::endl;
}

void __declspec(noinline) shiftImageInFourierTransformUnitTestPerformance() {
	std::cout << "shiftImageInFourierTransformUnitTestPerformance" << std::endl;

	static const int inout_size = 30;
    static const int inout_Fsize2 = inout_size*(inout_size/2+1);

	ShiftImageInFourierTransform<double,double> t(inout_size, 0.5, 0.5, 1.0);

	class V { 
		public: double* v; 
		V() : v((double*)aMalloc(sizeof(double)*inout_Fsize2,64)) {} 
		~V() { aFree(v); }
	};
	auto compare = [&](V&lhs, V&rhs) {
		for (int i = 0; i < inout_Fsize2; i++) {
			auto l = lhs.v[i];
			auto r = rhs.v[i];
			auto m = std::max(std::abs(l), std::abs(r));
			if (m == 0) continue;
			if (std::abs(l-r) * 10 < m) continue;
			std::cerr << "shiftImageInFourierTransform got wrong answer" << std::endl;
		}
	};

	static const int numberOfV_guess = (128000 / (sizeof(double)*inout_Fsize2*6)) & ~3;	// must be a multiple of 4 for this test
	static const int numberOfV = numberOfV_guess>0 ? numberOfV_guess : 4;
	std::cout << "    Doing numberOfV:" << numberOfV << std::endl;

	struct HeapAlloc {
		V vfin_real    [numberOfV];
		V vfin_imag    [numberOfV];
		V vfout_realOld[numberOfV];
		V vfout_imagOld[numberOfV];
		V vfout_realNew[numberOfV];
		V vfout_imagNew[numberOfV];
	}* heapAlloc = sNew(HeapAlloc);
	auto & vfin_real     = heapAlloc->vfin_real    ;
	auto & vfin_imag     = heapAlloc->vfin_imag    ;
	auto & vfout_realOld = heapAlloc->vfout_realOld;
	auto & vfout_imagOld = heapAlloc->vfout_imagOld;
	auto & vfout_realNew = heapAlloc->vfout_realNew;
	auto & vfout_imagNew = heapAlloc->vfout_imagNew;

	for (int j = 0; j < numberOfV; j++) {
		V & fin_real     = vfin_real    [j];
		V & fin_imag     = vfin_imag    [j];
		V & fout_realOld = vfout_realOld[j];
		V & fout_imagOld = vfout_imagOld[j];
		V & fout_realNew = vfout_realNew[j];
		V & fout_imagNew = vfout_imagNew[j];
		for (int i = 0; i < inout_Fsize2; i++) fin_real.v[i] = fin_imag.v[i] = float(i*j);
	}

	static const int repeats = 100000;
	for (int trial = 0; trial < 4; trial++) {
		auto const t0_start = dtime();
		for (int r = 0; r < repeats; r++)
		for (int j = 0; j < numberOfV; j+=1) {
			t.transform(
				&vfin_real[j+0].v[0], &vfin_imag[j+0].v[0], &vfout_realNew[j+0].v[0], &vfout_imagNew[j+0].v[0]);
		}
		auto const t0_finish = dtime();

		auto const t1_start = dtime();
		for (int r = 0; r < repeats; r++)
		for (int j = 0; j+1 < numberOfV; j+=2) {
			t.transform2(
				&vfin_real[j+0].v[0], &vfin_imag[j+0].v[0], &vfout_realNew[j+0].v[0], &vfout_imagNew[j+0].v[0],
				&vfin_real[j+1].v[0], &vfin_imag[j+1].v[0], &vfout_realNew[j+1].v[0], &vfout_imagNew[j+1].v[0]);
		}
		auto const t1_finish = dtime();

		auto const t2_start = dtime();
		for (int r = 0; r < repeats; r++)
		for (int j = 0; j+3 < numberOfV; j+=4) {
			t.transform4(
				&vfin_real[j+0].v[0], &vfin_imag[j+0].v[0], &vfout_realNew[j+0].v[0], &vfout_imagNew[j+0].v[0],
				&vfin_real[j+1].v[0], &vfin_imag[j+1].v[0], &vfout_realNew[j+1].v[0], &vfout_imagNew[j+1].v[0],
				&vfin_real[j+2].v[0], &vfin_imag[j+2].v[0], &vfout_realNew[j+2].v[0], &vfout_imagNew[j+2].v[0],
				&vfin_real[j+3].v[0], &vfin_imag[j+3].v[0], &vfout_realNew[j+3].v[0], &vfout_imagNew[j+3].v[0]);
		}
		auto const t2_finish = dtime();
		double sum(0);
		for (int j = 0; j < numberOfV; j++) {
			V & fin_imag     = vfin_imag    [j];
			V & fout_realOld = vfout_realOld[j];
			V & fout_imagOld = vfout_imagOld[j];
			V & fin_real     = vfin_real    [j];
			V & fout_realNew = vfout_realNew[j];
			V & fout_imagNew = vfout_imagNew[j];
			for (int i = 0; i < inout_Fsize2; i++) sum += fin_real.v[i] + fin_imag.v[i];
		}
		std::cout << "trial:" << trial << std::endl
			<< "     transform :" << t0_finish-t0_start << std::endl
			<< "     transform2:" << t1_finish-t1_start << std::endl
			<< "     transform4:" << t2_finish-t2_start << std::endl
			<< " produced << sum:" << sum << std::endl;
	}

	sDelete(heapAlloc);

	std::cout << "shiftImageInFourierTransformUnitTestPerformance end" << std::endl;
}

void testImageModule()
{
    // debug setFourierTransforms
    // TODO testing all image operation in setFourierTransforms...
    //    #define DEBUG_IMAGESTRANS
    // test windowFourierTransform
    // test softMaskOutsideMap
    // test ....
#ifdef DEBUG_IMAGESTRANS
    // write out mask image
    Mrcs::MrcsImages *mask_image = new Mrcs::MrcsImages(exp_image,ori_size,1);
    mask_image->write(fn+"_mask.mrcs");
    delete mask_image;
    std::cout<<"writeout mask particle to "<<fn+"_mask.mrcs"<<std::endl;
    //
    Mrcs::MrcsImages *mask_center_image = new Mrcs::MrcsImages(tempImage,ori_size,1);
    mask_center_image->write(fn+"_mask_center.mrcs");
    delete mask_center_image;
    std::cout<<"writeout mask center particle to "<<fn+"_mask_center.mrcs"<<std::endl;
    //
    Mrcs::MrcsFImages *mask_specturm_image = new Mrcs::MrcsFImages(tempFimage_real,tempFimage_imag,ori_size,1);
    mask_specturm_image->write(fn+"_mask_specturm.mrcs");
    delete mask_specturm_image;
    std::cout<<"writeout nomask fft specturm of particle to "<<fn+"_mask_specturm.mrcs"<<std::endl;
    // write out windows image
    Mrcs::MrcsFImages *mask_win_image = new Mrcs::MrcsFImages(exp_Fimage_mask_fine_real,exp_Fimage_mask_fine_imag,current_size,1);
    mask_win_image->write(fn+"_nomask_win.mrcs");
    delete mask_win_image;
    std::cout<<"writeout mask fft windows of particle to "<<fn+"_mask_win.mrcs,press any key to continue : ";
    std::cin>>debug_temp;
#endif
    
    //
#ifdef DEBUG_IMAGESTRANS
    // write out translate image
    Mrcs::MrcsImages *nomask_image = new Mrcs::MrcsImages(exp_image,ori_size,1);
    nomask_image->write(fn+"_nomask.mrcs");
    delete nomask_image;
    std::cout<<"translate the particle,x : "<<exp_old_offsetx[iimage]<<",y : "<<exp_old_offsety[iimage]<<std::endl;
    std::cout<<"writeout the translated particle to "<<fn+"_nomask.mrcs"<<std::endl;
    // write out centered image
    Mrcs::MrcsImages *nomask_center_image = new Mrcs::MrcsImages(tempImage,ori_size,1);
    nomask_center_image->write(fn+"_nomask_center.mrcs");
    delete nomask_center_image;
    std::cout<<"writeout the translated particle to "<<fn+"_nomask_center.mrcs"<<std::endl;
    // write out specturm of fft image
    Mrcs::MrcsFImages *nomask_specturm_image = new Mrcs::MrcsFImages(tempFimage_real,tempFimage_imag,ori_size,1);
    nomask_specturm_image->write(fn+"_nomask_specturm.mrcs");
    delete nomask_specturm_image;
    std::cout<<"writeout nomask fft specturm of particle to "<<fn+"_nomask_specturm.mrcs"<<std::endl;
    // write out windows image
    Mrcs::MrcsFImages *nomask_win_image = new Mrcs::MrcsFImages(exp_Fimage_nomask_real,exp_Fimage_nomask_imag,current_size,1);
    nomask_win_image->write(fn+"_nomask_win.mrcs");
    delete nomask_win_image;
    std::cout<<"writeout nomask fft windows of particle to "<<fn+"_nomask_win.mrcs,press any key to continue : ";
    std::cin>>debug_temp;
#endif
    
    //
#ifdef DEBUG_IMAGESTRANS
    std::string output_path;
    int debug_temp;
    std::cout<<"starting debugging setFourierTransforms(),where do you want to output the result : ";
    std::cin>>output_path;
    std::string fn = output_path+"/"+std::to_string((long long)(exp_first_image+iimage));
    // write out original particle
    Mrcs::MrcsImages *ori_image = new Mrcs::MrcsImages(exp_image,ori_size,1);
    ori_image->write(fn+"_original.mrcs");
    delete ori_image;
    // write out normalize image
    Mrcs::MrcsImages *normalize_image = new Mrcs::MrcsImages(tempImage,ori_size,1);
    normalize_image->write(fn+"_normalize.mrcs");
    delete normalize_image;
    std::cout<<"writeout original import particle,"<<fn+"_original.mrcs"<<std::endl;
    std::cout<<"write out normalized particle,"<<fn+"_normalize.mrcs,press any key to continue : ";
    std::cin>>debug_temp;
#endif
    
    //
    
}
