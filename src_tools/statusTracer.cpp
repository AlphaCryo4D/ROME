/***************************************************************************
 *
 * Dana-Farber Cancer Institute, Harvard Medical School and Peking University
 * Intel® Parallel Computing Center for Structural Biology
 *
 * Authors: "Yong Bei Ma(galowma@gmail.com)"
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
#include "../src/resmap/resmap_util.h"		// used for building precompiled headers on Windows

#include "statusTracer.h"

#define IS_MASTER_NODE (0 == MPI::COMM_WORLD.Get_rank())

void StatusTracer::backupStatus(std::string status_fn)
{
	std::string backup_fn = status_fn + "_backup";
	std::ofstream outBinaryFile;
	outBinaryFile.open((backup_fn + ".back").c_str(), std::ios::out | std::ios::binary | std::ios::trunc);
    if(IS_MASTER_NODE) std::cout << "--- write backup data to " << backup_fn + ".back" << std::endl;

	ERROR_CHECK(outBinaryFile.fail(), "Failed opening " + backup_fn + ".back to write backup");

	std::ofstream outTxtFile;
	outTxtFile.open((backup_fn + ".txt").c_str(), std::ios::out);
	std::cout.precision(20);
	auto newline = [&](std::ostream& os){
		static size_t newline = 0;
		if (newline++ % 5 == 0) os << std::endl;
	};
	{
		outTxtFile << "(what, length, type)" << std::endl;

		for (auto& data : doublePtr)
		{
			outBinaryFile.write((char*)data.data, (data.len)*sizeof(double));
			outTxtFile << " (" << data.what << ", " << data.len << ", " << "double) ";
			newline(outTxtFile);
		}
        for (auto& data : floatPtr)
        {
            outBinaryFile.write((char*)data.data, (data.len)*sizeof(float));
            outTxtFile << " (" << data.what << ", " << data.len << ", " << "int) ";
            newline(outTxtFile);
        }
		for (auto& data : intPtr)
		{
			outBinaryFile.write((char*)data.data, (data.len)*sizeof(int));
			outTxtFile << " (" << data.what << ", " << data.len << ", " << "int) ";
			newline(outTxtFile);
		}
        for (auto& data : boolPtr)
        {
            outBinaryFile.write((char*)data.data, (data.len)*sizeof(bool));
            outTxtFile << " (" << data.what << ", " << data.len << ", " << "bool) ";
            newline(outTxtFile);
        }
	}

    if(IS_MASTER_NODE) std::cout << "--- complete write" << std::endl;

	outBinaryFile.close();
	outTxtFile.close();
	ERROR_CHECK(outBinaryFile.fail(), "Failed writing or closing " + backup_fn + ".back to write backup");
}
// read the tracing data to _backup.model file
void StatusTracer::recoveryStatus(std::string status_fn)
{
	std::string backup_fn = status_fn + "_backup";
	std::ifstream inBinaryFile;
	inBinaryFile.open((backup_fn + ".back").c_str(), std::ios::in | std::ios::binary);
    if(IS_MASTER_NODE) std::cout << "--- read backup data from " << backup_fn + ".back , (name,sum)" << std::endl;

	ERROR_CHECK(inBinaryFile.fail(), "Failed opening " + backup_fn + ".back to read backup");

	std::cout.precision(20);

	{
		for (auto& data : doublePtr) {
			inBinaryFile.read((char*)data.data, (data.len)*sizeof(double));
		}
        for (auto& data : floatPtr) {
            inBinaryFile.read((char*)data.data, (data.len)*sizeof(float));
        }
		for (auto& data : intPtr) {
			inBinaryFile.read((char*)data.data, (data.len)*sizeof(int));
		}
        for (auto& data : boolPtr) {
            inBinaryFile.read((char*)data.data, (data.len)*sizeof(bool));
        }
	}

    if(IS_MASTER_NODE) std::cout << "--- complete read" << std::endl;

	inBinaryFile.close();
	ERROR_CHECK(inBinaryFile.fail(), "Failed writing or closing " + backup_fn + ".back to write backup");
}
//
void StatusTracer::checkStatus(std::string guidance_fn,bool do_recovery /*= true*/)
{
	std::string backup_fn = guidance_fn + "_backup";
	std::ifstream inBinaryFile;
	inBinaryFile.open((backup_fn + ".back").c_str(), std::ios::in | std::ios::binary);
    std::cout << "--- check backup data from " << backup_fn + ".back,starting comparing" << std::endl;
	
	ERROR_CHECK(inBinaryFile.fail(), "Failed opening " + backup_fn + ".back to read backup");

	std::cout.precision(20);
    std::ofstream compare_result(backup_fn+"_compare_result.txt");
    compare_result<<"do compare..."<<std::endl;
	{
		//
		int maxLength = 0;
		for (auto& data : doublePtr) maxLength = std::max(maxLength, data.len);
		double* dataTmp = mallocDoubleAligned(maxLength);
		for (auto& data : doublePtr) {
			inBinaryFile.read((char*)dataTmp, (data.len)*sizeof(double));
			if (!checkvec(compare_result, dataTmp, data.data, data.len)){
				compare_result << " " << data.what << " difference ###" << std::endl;
                if(do_recovery) {
					for (int i = 0; i < data.len; i++) data.data[i] = dataTmp[i];
                }
			}
		}
		freeAligned(dataTmp);
        //
        maxLength = 0;
        for (auto& data : floatPtr) maxLength = std::max(maxLength, data.len);
        float* dataTmp2Float = mallocFloatAligned(maxLength);
        double* dataTmp2Double = mallocDoubleAligned(maxLength);
        for (auto& data : floatPtr) {
            if(data.compareFloatWithDouble)
            {
                inBinaryFile.read((char*)dataTmp2Double, (data.len)*sizeof(double));
                if (!checkDiffTypeVec(compare_result, dataTmp2Double, data.data, data.len)){
                    compare_result << " " << data.what << " difference ###" << std::endl;
                    if(do_recovery) {
                    	for (int i = 0; i < data.len; i++) data.data[i] = dataTmp2Double[i];
                    }
                }
            }
            else
            {
                inBinaryFile.read((char*)dataTmp2Float, (data.len)*sizeof(float));
                if (!checkvec(compare_result, dataTmp2Float, data.data, data.len)){
                    compare_result << " " << data.what << " difference ###" << std::endl;
                    if(do_recovery) {
                    	for (int i = 0; i < data.len; i++) data.data[i] = dataTmp2Float[i];
                    }
                }
            }
        }
        freeAligned(dataTmp2Float);
        freeAligned(dataTmp2Double);
        //
        maxLength = 0;
        for (auto& data : intPtr) maxLength = std::max(maxLength, data.len);
        int* dataTmp3 = mallocIntAligned(maxLength);
        for (auto& data : intPtr) {
            inBinaryFile.read((char*)dataTmp3, (data.len)*sizeof(int));
            if (!checkvec(compare_result, dataTmp3, data.data, data.len)){
                compare_result << " " << data.what << " difference ###" << std::endl;
                if(do_recovery) {
                	for (int i = 0; i < data.len; i++) data.data[i] = dataTmp3[i];
                }
            }
        }
        freeAligned(dataTmp3);
        //
        maxLength = 0;
        for (auto& data : boolPtr) maxLength = std::max(maxLength, data.len);
        bool* dataTmp4 = mallocBoolAligned(maxLength);
        for (auto& data : boolPtr) {
            inBinaryFile.read((char*)dataTmp4, (data.len)*sizeof(bool));
            if (!checkvec(compare_result, dataTmp4, data.data, data.len)){
                compare_result << " " << data.what << " difference ###" << std::endl;
                if(do_recovery) {
                    for (int i = 0; i < data.len; i++) data.data[i] = dataTmp4[i];
                }
            }
        }
        freeAligned(dataTmp4);
	}

    std::cout << "--- complete check" << std::endl;

	inBinaryFile.close();
	ERROR_CHECK(inBinaryFile.fail(), "Failed writing or closing " + backup_fn + ".back to write backup");

}
