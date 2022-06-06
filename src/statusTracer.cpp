#include "util.h"		// used for building precompiled headers on Windows

#include "statusTracer.h"

#include <iostream>
#include "mpi.h"
#include "error.h"

void StatusTracer::backupStatus(std::string status_fn)
{
	std::string backup_fn = status_fn + "_backup";
	std::ofstream outBinaryFile;
	outBinaryFile.open((backup_fn + ".back").c_str(), std::ios::out | std::ios::binary | std::ios::trunc);
	MASTERNODE std::cout << "--- write backup data to " << backup_fn + ".back" << std::endl;

	ERROR_CHECK(outBinaryFile.fail(), "Failed opening " + backup_fn + ".back to write backup");

	std::ofstream outTxtFile;
	outTxtFile.open((backup_fn + ".txt").c_str(), std::ios::out);
	std::cout.precision(20);
	auto newline = [&](std::ostream& os){
		static size_t newline = 0;
		if (newline++ % 5 == 0) os << std::endl;
	};
	MASTERNODE
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
	}

	MASTERNODE std::cout << "--- complete write" << std::endl;

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
	MASTERNODE std::cout << "--- read backup data from " << backup_fn + ".back , (name,sum)" << std::endl;

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
	}

	MASTERNODE std::cout << "--- complete read" << std::endl;

	inBinaryFile.close();
	ERROR_CHECK(inBinaryFile.fail(), "Failed writing or closing " + backup_fn + ".back to write backup");
}
//
void StatusTracer::checkStatus(std::string guidance_fn)
{
	std::string backup_fn = guidance_fn + "_backup";
	std::ifstream inBinaryFile;
	inBinaryFile.open((backup_fn + ".back").c_str(), std::ios::in | std::ios::binary);
	MASTERNODE std::cout << "--- check backup data from " << backup_fn + ".back,starting comparing" << std::endl;

	ERROR_CHECK(inBinaryFile.fail(), "Failed opening " + backup_fn + ".back to read backup");

	std::cout.precision(20);

	MASTERNODE
	{
		//
		int maxLength = 0;
		for (auto& data : doublePtr) maxLength = std::max(maxLength, data.len);
		double* dataTmp = mallocDoubleAligned(maxLength);
		for (auto& data : doublePtr) {
			inBinaryFile.read((char*)dataTmp, (data.len)*sizeof(double));
			if (!checkvec(std::cout, dataTmp, data.data, data.len)){
				std::cout << " " << data.what << " difference ###" << std::endl;
				for (int i = 0; i < data.len; i++)
					data.data[i] = dataTmp[i];
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
                if (!checkDiffTypeVec(std::cout, dataTmp2Double, data.data, data.len)){
                    std::cout << " " << data.what << " difference ###" << std::endl;
                    for (int i = 0; i < data.len; i++)
                        data.data[i] = dataTmp2Double[i];
                }
            }
            else
            {
                inBinaryFile.read((char*)dataTmp2Float, (data.len)*sizeof(float));
                if (!checkvec(std::cout, dataTmp2Float, data.data, data.len)){
                    std::cout << " " << data.what << " difference ###" << std::endl;
                    for (int i = 0; i < data.len; i++)
                        data.data[i] = dataTmp2Float[i];
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
            if (!checkvec(std::cout, dataTmp3, data.data, data.len)){
                std::cout << " " << data.what << " difference ###" << std::endl;
                for (int i = 0; i < data.len; i++)
                    data.data[i] = dataTmp3[i];
            }
        }
        freeAligned(dataTmp3);
	}

	MASTERNODE std::cout << "--- complete check" << std::endl;

	inBinaryFile.close();
	ERROR_CHECK(inBinaryFile.fail(), "Failed writing or closing " + backup_fn + ".back to write backup");

}
