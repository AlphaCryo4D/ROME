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

#include "spider.h"

namespace Spider{
    
    void readDatHead(std::string fn_dat,DatHead& dathead,std::string ordering)
    {
        
        FILE* filehandle = fopen(fn_dat.c_str(),"rb");
        
        if (filehandle == NULL)
            ERROR_REPORT("Fail to open "+fn_dat);

        // pointer to the dat head structure
        float *buffer = &dathead.nz;
        
        if(fread((char*)buffer,sizeof(DatHead),1,filehandle) == NULL)
            ERROR_REPORT("Fail to read "+fn_dat+" Head.");
        
        if (ordering == "ieee-be")
            for (int i = 0; i < 256; i++)
                ((unsigned int *)buffer)[i] = ntohl(((unsigned int *)buffer)[i]);
        
        fclose(filehandle);
    }
    
    void readDatData(std::string fn_dat,DatHead& dathead,FloatImages& datImages,std::string ordering){

        // check the dathead
        int N  = int((dathead.istack > 0)?dathead.maxim:dathead.nz);
        int ny = int(dathead.ny);
        int nx = int(dathead.nx);
        
        int imageSize = datImages.imageSide();
        int numberOfImages = datImages.nr_images();
        
        ERROR_CHECK(imageSize != ny || imageSize != nx || N != numberOfImages,
                    "ListOfImage data is not much dathead");
        
        int iform = dathead.iform;
        if(iform != 1)
            ERROR_REPORT(fn_dat+" file' iform is incorrect.");
        
        long offset = dathead.labbyt;

        FILE* filehandle = fopen(fn_dat.c_str(),"rb");
        if (filehandle == NULL)
            ERROR_REPORT("Fail to open "+fn_dat);
        
        int imageSize2 = imageSize*imageSize;
        std::vector<float> buffer(imageSize2);
        for (size_t iimage = 0; iimage < numberOfImages; iimage++) {
            // If it is a stack,it must skip another image's head or its data
            if (dathead.istack > 0){
                // TODO,not know why this one get wrong result
                // offset = dathead.labbyt + (iimage+1)*dathead.labbyt + iimage*ny*nx*sizeof(float);
                offset = dathead.labbyt + (iimage+1)*dathead.labbyt;
                offset += iimage*ny*nx*sizeof(float);
            }
            // warning :
            if(fseek(filehandle,offset,SEEK_SET) == -1)
                ERROR_REPORT("offset out of range for "+fn_dat);
            
            if(fread((char*)(&buffer[0]),sizeof(float)*imageSize2,1,filehandle) == NULL)
                ERROR_REPORT("Fail to read "+fn_dat);
            
            if(ordering == "ieee-be")
                for (size_t i = 0; i < imageSize2; i++)
                    ((unsigned int *)&buffer[0])[i] = ntohl(((unsigned int *)&buffer[0])[i]);

            auto* image_data = datImages.image_ptr(iimage);
            for (int i = 0; i < imageSize2; i++) image_data[i] = buffer[i];
        }
        
        fclose(filehandle);
        
    }
    
    void writeDatData(std::string filename,float* data,int nz,int ny,int nx,bool is_stack,std::string ordering)
    {
        DatImages listOfDatImages(data,nx, nz);
        
        writeDatData(filename,listOfDatImages,ordering);
    }
    
    void writeDatData(std::string filename,FloatImages& datImages, std::string ordering){
        
        int imageSize = datImages.imageSide();
        int numberOfImages = datImages.nr_images();
        bool is_stack = numberOfImages>1?true:false;
        
        DatHead datHead;
        // Notice the head nz always 1!
        datHead.nz = 1;
        datHead.ny = imageSize;
        datHead.nx = imageSize;
        
        datHead.iform = 1;
        
        int lenbyt = datHead.nx * 4;
        int labrec = 1024 / lenbyt;
        if(1024 % lenbyt != 0) labrec = labrec + 1;
        int labbyt = labrec * lenbyt;
        
        // Multipy by nx will be head+data
        datHead.irec = labrec + datHead.ny;
        datHead.labrec = labrec;
        datHead.labbyt = labbyt;
        datHead.lenbyt = lenbyt;
        
        // I donot know what's this means
        datHead.unused3 = -1;
        
        // Set time,ignore this
        char date[12] = "02-OCT-2015";
        char time[9] = "18:25:55";
        strcpy(datHead.cdat,date);
        strcpy(datHead.ctim, time);
        
        if (is_stack) {
            datHead.istack = 2;
            // This is the real nz!
            datHead.maxim = numberOfImages;
        }
        
        FILE* filehandle = fopen((filename+".dat").c_str(),"wb");
        if (filehandle == NULL)
            ERROR_REPORT("Fail to open "+filename);
        
        // Write head first
        // float* headBuffer = (float*)aMalloc(datHead.labbyt,64);
        // memset(headBuffer, 0, datHead.labbyt);
        std::vector<float> headBuffer(datHead.labbyt);
        for(auto& v : headBuffer) v = 0.;
        
        if(ordering == "ieee-be"){
            float* datHead_aux = &datHead.nz;
            for (int i = 0; i < 256; i++) {
                ((unsigned int *)&headBuffer[0])[i] = ntohl(((unsigned int *)datHead_aux)[i]);
            }
        }
        
        fwrite(&headBuffer[0], datHead.labbyt, 1, filehandle);
        
        // write image one by one
        int imageSize2 = imageSize*imageSize;
        std::vector<float> dataBuffer(imageSize2);
        
        for (size_t iimage = 0; iimage < numberOfImages; iimage++) {
            
            // copy the data to buffer
            auto* image_data = datImages.image_ptr(iimage);
            for (int i = 0; i < imageSize2; i++) dataBuffer[i] = image_data[i];
            
            if(ordering == "ieee-be"){
                for (int i = 0; i < imageSize2; i++)
                    ((unsigned int *)&dataBuffer[0])[i] = ntohl(((unsigned int *)&dataBuffer[0])[i]);
            }
            
            if (is_stack) {
                
                datHead.istack = 0;
                datHead.unused3 = 0;
                datHead.maxim = 0;
                datHead.imgnum = iimage+1;
                if(ordering == "ieee-be"){
                    float* datHead_aux = &datHead.nz;
                    for (int i = 0; i < 256; i++)
                        ((unsigned int *)&headBuffer[0])[i] = ntohl(((unsigned int *)datHead_aux)[i]);
                }
                
                if(fwrite(&headBuffer[0], datHead.labbyt, 1, filehandle) == NULL)
                    ERROR_REPORT("Fail to write "+filename);
                
            }
            
            if(fwrite(&dataBuffer[0], sizeof(float)*imageSize2, 1, filehandle) == NULL)
                ERROR_REPORT("Fail to write "+filename);
            
        }
        
        fclose(filehandle);
        
    }
    
    template<typename T>
    void FQ_Q(T* data,float* B,int nx,int ny,fftwf_complex* Fdata,double filter)
    {
        bool spider_sign = true;
        
        int n2x = nx*2;
        int n2y = ny*2;
        // int LSD = n2x + 2;
        for (int i = 0; i < n2x*n2y; i++) B[i] = 0;
        
        // read image in B(nx*ny)
        for (int j = 0; j < ny; j++) {
            for (int i = 0; i < nx; i++) {
                B[j*n2x+i] = data[j*nx+i];
            }
        }
        double ave = 0;
        for (int i = 0; i < nx; i++) {
            ave = ave + B[i*n2y+0] + B[i*n2y+ny-1];
        }
        for (int j = 1; j < ny-1; j++) {
            ave = ave + B[0*n2y+j] + B[(nx-1)*n2y+j];
        }
        ave /= (2*(nx+ny)-4);
        
        for (int j = 0; j < ny; j++)
            for (int i = nx; i < n2x; i++)
                B[j*n2x+i] = ave;
        
        for (int j = nx; j < n2x; j++)
            for (int i = 0; i < n2y; i++)
                B[j*n2x+i] = ave;
        
        fftwf_plan forward_plan,backward_plan;
        forward_plan = fftwf_plan_dft_r2c_2d(n2y,n2x,B,Fdata,FFTW_ESTIMATE);
        fftwf_execute(forward_plan);
        fftwf_destroy_plan(forward_plan);
        
        if (spider_sign) {
            for (int i = 0; i < n2y*(n2x/2+1); i++) {
                Fdata[i][1] = -1*Fdata[i][1];
            }
        }
        
        double parm1 = filter;
        double parm2 = parm1;
        if (filter < 0.0 || filter > 1.0) parm1 = 0.5*parm1/(nx/2);
        
        int nr2    = n2y / 2;
        double x1    = (n2x/2)*(n2x/2);
        double y1     = nr2*nr2;
        double parm   = parm1*parm1;
        double parm22 = parm2*parm2;
        
        for (int j = 1; j <= n2y; j++) {
            
            int iy = j-1;
            if(iy > nr2) iy = iy - n2y;
            
            for (int i = 1; i <= (n2x/2+1); i++) {
                
                int ix = i-1;
                if(0.25*((ix*ix)/x1/parm + (iy*iy)/y1/parm22) > 1.0){
                    Fdata[(j-1)*(n2x/2+1)+i-1][0] = 0.0;
                    Fdata[(j-1)*(n2x/2+1)+i-1][1] = 0.0;
                }
                
            }
        }
        
        if (spider_sign) {
            for (int i = 0; i < n2y*(n2x/2+1); i++) {
                Fdata[i][1] = -1*Fdata[i][1];
            }
        }
        
        backward_plan = fftwf_plan_dft_c2r_2d(n2y,n2x,Fdata,B,FFTW_ESTIMATE);
        fftwf_execute(backward_plan);
        fftwf_destroy_plan(backward_plan);
        
        for (int j = 0; j < ny; j++) {
            for (int i = 0; i < nx; i++) {
                data[j*nx+i] = B[j*n2x+i];
            }
        }
    }
    
    void applyFilter(FloatImages& fLoatImages,double filter){
        int imageSize = fLoatImages.imageSide();
        int numberOfImage = fLoatImages.nr_images();
        int image2size = 2*imageSize;
        // some temporary data
        float*         image_data_tmp  = (float*        )aMalloc(sizeof(float        )*image2size*image2size      ,64);
        fftwf_complex* image_Fdata_tmp = (fftwf_complex*)aMalloc(sizeof(fftwf_complex)*image2size*(image2size/2+1),64);
        
        for (int iimage = 0; iimage < numberOfImage; iimage++) {
            FQ_Q(fLoatImages.image_ptr(iimage),image_data_tmp,imageSize,imageSize,image_Fdata_tmp,filter);
        }
        
        // free temporary data
        aFree(image_data_tmp);
        aFree(image_Fdata_tmp);
    }
    
    float FBS2(float x,float y,int nx,int ny,const float* data,const float* X1,const float* Y1,const float* XY2)
    {
        int i,j,i2,j2,i3,j3;
        float A0,A1,A2,A3,ADX,BDX,DADX,DBDX;
        
        i = floor(x);
        j = floor(y);
        
        float dx    = x - i;
        float dy    = y - j;
        float dxsqr = dx*dx;
        float dxcub = dx*dx*dx;
        
        bool chkbound = true;
        if (chkbound) {
            i2 = modulo(i-1,nx);
            j2 = modulo(j-1,ny);
            i3 = modulo(i,nx);
            j3 = modulo(j,ny);
        }
        else{
            i2 = i - 1;
            j2 = j - 1;
            i3 = i;
            j3 = j;
        }
        
        A0   = data[j2*nx+i2];
        A1 = X1[j2*nx+i2];
        A2   = 3*( data[j2*nx+i3]-A0) -2*A1 - X1[j2*nx+i3];
        A3   = 2*(-data[j2*nx+i3]+A0)  + A1 + X1[j2*nx+i3];
        ADX  = A0 + A1*dx + A2*dxsqr + A3*dxcub;
        
        A0   = data[j3*nx+i2];
        A1   = X1[j3*nx+i2];
        A2   = 3*( data[j3*nx+i3]-A0) - 2*A1 - X1[j3*nx+i3];
        A3   = 2*(-data[j3*nx+i3]+A0) +   A1 + X1[j3*nx+i3];
        BDX  = A0 + A1*dx + A2*dxsqr + A3*dxcub;
        
        A0   = Y1[j2*nx+i2];
        A1   = XY2[j2*nx+i2];
        A2   = 3*( Y1[j2*nx+i3]-A0) -2*A1 - XY2[j2*nx+i3];
        A3   = 2*(-Y1[j2*nx+i3]+A0) +  A1 + XY2[j2*nx+i3];
        DADX = A0 + A1*dx + A2*dxsqr + A3*dxcub;
        
        A0   = Y1[j3*nx+i2];
        A1   = XY2[j3*nx+i2];
        A2   = 3*( Y1[j3*nx+i3]-A0) - 2*A1 - XY2[j3*nx+i3];
        A3   = 2*(-Y1[j3*nx+i3]+A0) +   A1 + XY2[j3*nx+i3];
        DBDX = A0 + A1*dx + A2*dxsqr + A3*dxcub;
        
        A2   = 3*(BDX - ADX) - 2*DADX - DBDX;
        A3   = 2*(ADX - BDX) +   DADX + DBDX;
        
        float fbs2 = ADX + DADX * dy + A2 * dy*dy + A3 * dy*dy*dy;
        
        // std::cout<<ADX<<" "<<BDX<<" "<<DADX<<" "<<DBDX<<" "<<FBS2<<" ";
        // exit(1);
        return fbs2;
    }
    
    // TODO : setup a RT ST plan,then do
    void RT_SF(const float* data,float *out_data,int nx,int ny,float angle,float shiftx,float shifty){
        
        bool spider_sign = true;
        bool spider_scale = true;
        float scale = 1.;
        float spider_PI2 = 6.28318530717958647692;
        float spider_PI = 3.14159265358979323846;
        
        /******  FBS2_PREP  ******/
        
        fftwf_plan forward_plan,backward_plan;
        
        float *data2 = (float*)aMalloc(sizeof(float)*ny*nx,64);
        memcpy(data2, data, sizeof(float)*ny*nx);
        fftwf_complex* Fdata = (fftwf_complex*) aMalloc(sizeof(fftwf_complex)*ny*(nx/2+1), 64);
        fftwf_complex* Ftemp = (fftwf_complex*) aMalloc(sizeof(fftwf_complex)*ny*(nx/2+1), 64);
        float *WX  = (float*)aMalloc(sizeof(float)*(nx/2+1),64);
        float *WY  = (float*)aMalloc(sizeof(float)*ny,      64);
        float *X1  = (float*)aMalloc(sizeof(float)*ny*nx,   64);
        float *Y1  = (float*)aMalloc(sizeof(float)*ny*nx,   64);
        float *XY2 = (float*)aMalloc(sizeof(float)*ny*nx,   64);
        
        forward_plan = fftwf_plan_dft_r2c_2d(ny,nx,data2,Fdata,FFTW_ESTIMATE);
        fftwf_execute(forward_plan);
        fftwf_destroy_plan(forward_plan);
        
        if (spider_sign) {
            for (int i = 0; i < ny*(nx/2+1); i++) {
                Fdata[i][1] = -1*Fdata[i][1];
            }
        }
        
        float A4 = spider_PI2/nx;
        for (int i = 0; i < (nx/2+1); i++) {
            WX[i] = i*A4;
        }
        
        A4 = spider_PI2/ny;
        for (int i = 0; i < ny/2+1; i++) {
            WY[i] = i*A4;
        }
        for (int i = ny/2+1; i < ny; i++) {
            WY[i] = (i-ny)*A4;
        }
        
        
        for (int j = 0; j < ny; j++) {
            for (int i = 0; i < nx/2+1; i++) {
                Ftemp[j*(nx/2+1)+i][0] = Fdata[j*(nx/2+1)+i][1]*WX[i];
                Ftemp[j*(nx/2+1)+i][1] = -Fdata[j*(nx/2+1)+i][0]*WX[i];
            }
        }
        if (spider_sign) {
            for (int i = 0; i < ny*(nx/2+1); i++) {
                Ftemp[i][1] = -1*Ftemp[i][1];
            }
        }
        
        backward_plan = fftwf_plan_dft_c2r_2d(ny,nx,Ftemp,X1,FFTW_ESTIMATE);
        fftwf_execute(backward_plan);
        fftwf_destroy_plan(backward_plan);
        
        if (spider_scale) {
            float pix = 1./(ny*nx);
            for (int i = 0; i < ny*nx; i++) {
                X1[i] *= pix;
            }
        }
        
        
        for (int j = 0; j < ny; j++) {
            for (int i = 0; i < nx/2+1; i++) {
                Ftemp[j*(nx/2+1)+i][0] = Fdata[j*(nx/2+1)+i][1]*WY[j];
                Ftemp[j*(nx/2+1)+i][1] = -Fdata[j*(nx/2+1)+i][0]*WY[j];
            }
        }
        
        if (spider_sign) {
            for (int i = 0; i < ny*(nx/2+1); i++) {
                Ftemp[i][1] = -1*Ftemp[i][1];
            }
        }
        
        backward_plan = fftwf_plan_dft_c2r_2d(ny,nx,Ftemp,Y1,FFTW_ESTIMATE);
        fftwf_execute(backward_plan);
        fftwf_destroy_plan(backward_plan);
        
        if (spider_scale) {
            float pix = 1./(ny*nx);
            for (int i = 0; i < ny*nx; i++) {
                Y1[i] *= pix;
            }
        }
        
        
        for (int j = 0; j < ny; j++) {
            for (int i = 0; i < nx/2+1; i++) {
                Ftemp[j*(nx/2+1)+i][0] = Fdata[j*(nx/2+1)+i][1]*WY[j];
                Ftemp[j*(nx/2+1)+i][1] = -Fdata[j*(nx/2+1)+i][0]*WY[j];
            }
        }
        
        if (spider_sign) {
            for (int i = 0; i < ny*(nx/2+1); i++) {
                Ftemp[i][1] = -1*Ftemp[i][1];
            }
        }
        
        backward_plan = fftwf_plan_dft_c2r_2d(ny,nx,Ftemp,XY2,FFTW_ESTIMATE);
        fftwf_execute(backward_plan);
        fftwf_destroy_plan(backward_plan);
        
        if (spider_scale) {
            float pix = 1./(ny*nx);
            for (int i = 0; i < ny*nx; i++) {
                XY2[i] *= pix;
            }
        }
        
        /******   RTSF ********/
        
        float theta = angle * spider_PI / 180;
        float costh = cos(theta);
        float sinth = sin(theta);
        
        int cx    = nx / 2 + 1;
        int cy    = ny / 2 + 1;
        
        float shx = modulo(shiftx,(float)nx);
        float shy = modulo(shifty,(float)ny);
        // std::cout<<shx<<" "<<shiftx<<" "<<shy<<" "<<shifty<<std::endl;
        
        if (scale == 1.) {
            
            float fy0    = - shy - cy;
            float fy1    = - shy + ny - cy;
            float fy2    = - shy - ny - cy;
            
            float fx0    = - shx - cx;
            float fx1    = - shx + nx - cx;
            float fx2    = - shx - nx - cx;
            
            float shypny = shy + ny;
            float shxpnx = shx + nx;
            
            for (int iy = 1; iy <= ny; iy++) {
                float fy = iy + fy0;
                if((iy-1) < shy) fy = iy + fy1;
                if((iy-1) >= shypny) fy = iy + fy2;
                float ycod =  costh*fy+cy;
                float ysid = -sinth*fy+cx;
                
                for (int ix = 1; ix <= nx; ix++) {
                    float fx = ix+fx0;
                    if ((ix-1) < shx) fx = ix + fx1;
                    if ((ix-1) >= shxpnx) fx = ix + fx2;
                    
                    float xold = costh*fx + ysid;
                    float yold = sinth*fx + ycod;
                    
                    out_data[(iy-1)*nx+(ix-1)] = FBS2(xold,yold,nx,ny,data,X1,Y1,XY2);
                }
            }
        }
        
        aFree(Fdata);
        aFree(Ftemp);
        aFree(data2);
        aFree(WX);
        aFree(WY);
        aFree(X1);
        aFree(Y1);
        aFree(XY2);
    }
    
    void UnitTest(std::string fn_mrcs){
        std::cout<<"testing Spider---------"<<std::endl;
        DatHead datHead;
        std::cout<<"sizeof(DatHead) = "<<sizeof(DatHead)<<",shoud be 1024=256*4."<<std::endl;
        datHead.print();
        
        std::cout<<"reading the image data from mrcs file.---------------"<<std::endl;
        Mrcs::MrcsHead mrcsHead;
        std::cout<<"read mrcs head.--------"<<std::endl;
        Mrcs::readMrcsHead(fn_mrcs, mrcsHead);
        mrcsHead.print();
        
        std::cout<<"reading mrcs image data------------------"<<std::endl;
        int N = mrcsHead.NS;
        int size = mrcsHead.NC;
        DatImages listOfImages(size,N);
        Mrcs::readMrcsData(fn_mrcs, mrcsHead, listOfImages);
        
        std::cout<<"write the image data as *.dat,then use the spider to read this data.-------------------"<<std::endl;
        writeDatData(fn_mrcs, listOfImages);

        std::cout<<"read the new dat head.---------------------------------"<<fn_mrcs<<std::endl;
        readDatHead(fn_mrcs+".dat", datHead);
        datHead.print();
        
        std::cout<<"testing the DatHead copy function -----------------------"<<fn_mrcs<<std::endl;
        DatHead datHead2;
        datHead2 = datHead;
        datHead2.print();

        std::cout<<"reading the dat data-------------------------------------"<<fn_mrcs<<std::endl;
        readDatData(fn_mrcs+".dat", datHead, listOfImages);
        
        std::cout<<"write back the mrcs file---------------------------------"<<std::endl;
        Mrcs::writeMrcsData(fn_mrcs, listOfImages);
    }
    
    void UnitTest2(std::string fn_mrcs){
        std::cout<<"testing image operation-----------------------------------"<<std::endl;
        std::cout<<"reading the image data from mrcs file.--------------------"<<std::endl;
        Mrcs::MrcsHead mrcsHead;
        Mrcs::readMrcsHead(fn_mrcs, mrcsHead);
        int N = mrcsHead.NS;
        int size = mrcsHead.NC;
        Mrcs::MrcsImages listOfImages_original(size, N);
        Mrcs::readMrcsData(fn_mrcs, mrcsHead, listOfImages_original);
        std::cout<<"create two copy of image data,the new and old one,for running new and old version of code.-----"<<std::endl;
        DatImages listOfImages_old(size, N);
        DatImages listOfImages_new(size, N);
        std::cout<<"testing apply fliter-------------------------"<<std::endl;
        for (int i = 1; i < 10; i++) {
            double fliter = static_cast <double> (rand()) / (static_cast <double> (RAND_MAX/5));
            std::cout<<fliter<<std::endl;
            // the new version of code
            memcpy(listOfImages_new.image_ptr(0), listOfImages_original.image_ptr(0), sizeof(float)*N*size*size);
            applyFilter(listOfImages_new, fliter);
            // the old version
            memcpy(listOfImages_old.image_ptr(0), listOfImages_original.image_ptr(0), sizeof(float)*N*size*size);
            for (int iimage = 0; iimage < N; iimage++) {
                // need to change the data type from float to double then test
                spider_fq_q(listOfImages_old.image_ptr(iimage), size, size, fliter);
            }
            // compare the result
            auto* image_data_new = listOfImages_new.image_ptr(0);
            auto* image_data_old = listOfImages_old.image_ptr(0);
            for (int i = 0; i < N*size*size; i++) {
                if (image_data_new[i] != image_data_old[i]) {
                    std::cout<<"wrong.apply fliter...."<<std::endl;
					EXIT_ABNORMALLY;
                }
            }
        }
        std::cout<<"apply flite 0.05 to see what happen.....-----------"<<std::endl;
        memcpy(listOfImages_new.image_ptr(0), listOfImages_original.image_ptr(0), sizeof(float)*N*size*size);
        applyFilter(listOfImages_new, 0.05);
        Mrcs::writeMrcsData(fn_mrcs+"_fliter", listOfImages_new);
        std::cout<<"testing done..... for apply fliter----------------"<<std::endl;
        
        std::cout<<"testing rotation and shift the image----------------"<<std::endl;
        // temp data
        float* new_image_in = (float*)aMalloc(sizeof(float)*size*size,64);
        float* new_image_out = (float*)aMalloc(sizeof(float)*size*size,64);
        float* old_image_in = (float*)aMalloc(sizeof(float)*size*size,64);
        float* old_image_out = (float*)aMalloc(sizeof(float)*size*size,64);
        for (int i = 0; i < 10; i++) {
            memcpy(listOfImages_new.image_ptr(0), listOfImages_original.image_ptr(0), sizeof(float)*N*size*size);
            memcpy(listOfImages_old.image_ptr(0), listOfImages_original.image_ptr(0), sizeof(float)*N*size*size);
            for (int iimage = 0; iimage < N; iimage++) {
                double rotation = static_cast <double> (rand()) / (static_cast <double> (RAND_MAX/1000)) - 500.0;
                double shiftx = static_cast <double> (rand()) / (static_cast <double> (RAND_MAX/1000)) - 500.0;
                double shifty = static_cast <double> (rand()) / (static_cast <double> (RAND_MAX/1000)) - 500.0;
                std::cout<<rotation<<" "<<shiftx<<" "<<shifty<<std::endl;
                // new
                for (int i = 0; i < size*size; i++) new_image_in[i] = (listOfImages_new.image_ptr(iimage))[i];
                RT_SF(new_image_in, new_image_out, size, size, rotation, shiftx, shifty);
                // old
                for (int i = 0; i < size*size; i++) old_image_in[i] = (listOfImages_old.image_ptr(iimage))[i];
                spider_rtsf(old_image_in, old_image_out, size, size, rotation, shiftx, shifty);
                // compare
                for (int i = 0; i < size*size; i++) {
                    if (new_image_out[i] != old_image_out[i]) {
                        std::cout<<"wrong...rotate and shit image...."<<std::endl;
						EXIT_ABNORMALLY;
                    }
                }
            }
        }
        std::cout<<"random shift and rotate image to see what happen-------------"<<std::endl;
        memcpy(listOfImages_new.image_ptr(0), listOfImages_original.image_ptr(0), sizeof(float)*N*size*size);
        for (int iimage = 0; iimage < N; iimage++) {
            double rotation = static_cast <double> (rand()) / (static_cast <double> (RAND_MAX/1000)) - 500.0;
            double shiftx = static_cast <double> (rand()) / (static_cast <double> (RAND_MAX/1000)) - 500.0;
            double shifty = static_cast <double> (rand()) / (static_cast <double> (RAND_MAX/1000)) - 500.0;
            for (int i = 0; i < size*size; i++) new_image_in[i] = (listOfImages_new.image_ptr(iimage))[i];
            RT_SF(new_image_in, new_image_out, size, size, rotation, shiftx, shifty);
            for (int i = 0; i < size*size; i++) (listOfImages_new.image_ptr(iimage))[i] = new_image_out[i];
        }
        Mrcs::writeMrcsData(fn_mrcs+"_rt_sf", listOfImages_new);
        aFree(new_image_in);
        aFree(new_image_out);
        aFree(old_image_in);
        aFree(old_image_out);
        
    }
}


