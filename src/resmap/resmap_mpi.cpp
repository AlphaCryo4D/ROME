/***************************************************************************
 *
 * Intel® Parallel Computing Center for Structural Biology
 * Principal Investigator : Youdong (Jack) Mao (Youdong_Mao@dfci.harvard.edu)
 * Dana-Farber Cancer Institute, Harvard Medical School and Peking University
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

#include "resmap_util.h"		// used for building precompiled headers on Windows

#include "resmap_mpi.h"

#ifdef FAKEMPI
namespace MPI {
    Datatype DOUBLE,FLOAT,INT,BOOL;
    Op SUM;
	Intracomm COMM_WORLD;
};
#endif

// Divides a number into most equally groups
int divide_equally(int N, int size, int rank,int &first,int &last)
{
    int jobs_per_worker = N / size;
    int jobs_resting = N % size;
    if (rank < jobs_resting)
    {
        first = rank * (jobs_per_worker + 1);
        last = first + jobs_per_worker;
    }
    else
    {
        first = rank * jobs_per_worker + jobs_resting;
        last = first + jobs_per_worker - 1;
    }
    return last - first + 1;
}

//
int divide_equally_pool(int N, int nr_pool, int size, int rank,int &first,int &last){
    int nr_pool_block_start,nr_pool_bolck_end;
    int nr_pool_block_num = ceil((double)N/nr_pool);
    //divded the images for MPI
    divide_equally(nr_pool_block_num,size,rank, nr_pool_block_start, nr_pool_bolck_end);
    first = nr_pool_block_start*nr_pool;
    last = std::min((nr_pool_bolck_end+1)*nr_pool-1,N-1);
    return last-first+1;
}

// In which group from divide_equally is myself?
int divide_equally_which_group(int N, int size,int myself)
{
    int first, last;
    for (int rank = 0; rank < size; rank++)
    {
        divide_equally(N, size, rank, first, last);
        if (myself >= first && myself <= last)
            return rank;
    }
    return -1;
}


int divideIndex(int N,int nodes,int node,int &n_start,int &n_end){
	int sub_N = ceil((double)N/nodes);
	n_start = 0 + node*sub_N;
	n_end = (n_start+sub_N) < N?(n_start+sub_N):N;
	return (n_end-n_start);
}

//#if 0
///**** it will tranpose matrix(N*D) between cluster,Mat is (N/nodes)xD cross each node,new Mat is N*(D/nodes) ****/
//void transposeMatCluster(double* Mat,double *newMat,int N,int D,int node,int nodes){
//
//	//create buffer for each node,it shoud be same size (larger enough size to keep total data).
//	int sub_N = ceil((double)N/nodes);
//	int sub_D = ceil((double)D/nodes);
//	double *buffer = (double*)_mm_malloc(sizeof(double)*sub_N*sub_D,64);
//	
//	//create final transpose date,it may has different size in final node
//	int my_d_interval,d_start,d_end;
//	my_d_interval = divideIndex(D,nodes,node,d_start,d_end);
//	//double *newMatbuffer = (double*)_mm_malloc(sizeof(double)*my_d_interval*N,64);
//	//intilize some other
//	int my_n_interval,n_start,n_end;
//	my_n_interval = divideIndex(N,nodes,node,n_start,n_end);
//	int *rcount = (int*)_mm_malloc(sizeof(int)*nodes,64);
//	int *displs = (int*)_mm_malloc(sizeof(int)*nodes,64);
//
//	//each node will get different piece of D
//	for(int dist_node = 0;dist_node < nodes;dist_node++){
//		int dist_d_interval = divideIndex(D,nodes,dist_node,d_start,d_end);
//		//copy the data to buffer
//		for(int n = 0;n < my_n_interval;n++)
//			memcpy(buffer+n*dist_d_interval,Mat+n*D+d_start,sizeof(double)*1*dist_d_interval);
//
//		//gatehr the data
//		int count = my_n_interval*dist_d_interval;
//		//MPI::COMM_WORLD.Barrier();
//		MPI::COMM_WORLD.Gather(&count,1,MPI::INT,rcount,1,MPI::INT,dist_node);
//		if(node == dist_node){
//			//std::cout<<rcount[0]<<" "<<rcount[1]<<std::endl;
//			displs[0] = 0;
//			for(int i = 1;i < nodes;i++)
//				displs[i] = displs[i-1]+rcount[i-1];
//		}
//		//MPI::COMM_WORLD.Barrier();
//		MPI::COMM_WORLD.Gatherv(buffer,my_n_interval*dist_d_interval,MPI::DOUBLE,newMat,rcount,displs,MPI::DOUBLE,dist_node);
//		//tranpose the data(N*my_d_interval) to data(my_d_interval*N)
//		if(node == dist_node)
//			mkl_dimatcopy('R','T',N,my_d_interval,1,newMat,my_d_interval,N);
//		MPI::COMM_WORLD.Barrier();
//	}
//
//	_mm_free(buffer);
//	//_mm_free(newMatbuffer);
//	_mm_free(rcount);
//	_mm_free(displs);
//}
//
//
///**** it will tranpose matrix(N*D) between cluster,Mat is (N/nodes)xD cross each node,new Mat is N*(D/nodes) ****/
//void transposeMatCluster(float* Mat,float *newMat,int N,int D,int node,int nodes){
//
//	//create buffer for each node,it shoud be same size (larger enough size to keep total data).
//	int sub_N = ceil((double)N/nodes);
//	int sub_D = ceil((double)D/nodes);
//	float *buffer = (float*)_mm_malloc(sizeof(float)*sub_N*sub_D,64);
//	
//	//create final transpose date,it may has different size in final node
//	int my_d_interval,d_start,d_end;
//	my_d_interval = divideIndex(D,nodes,node,d_start,d_end);
//	//double *newMatbuffer = (double*)_mm_malloc(sizeof(double)*my_d_interval*N,64);
//	//intilize some other
//	int my_n_interval,n_start,n_end;
//	my_n_interval = divideIndex(N,nodes,node,n_start,n_end);
//	int *rcount = (int*)_mm_malloc(sizeof(int)*nodes,64);
//	int *displs = (int*)_mm_malloc(sizeof(int)*nodes,64);
//
//	//each node will get different piece of D
//	for(int dist_node = 0;dist_node < nodes;dist_node++){
//	//int dist_node = 1;
//		int dist_d_interval = divideIndex(D,nodes,dist_node,d_start,d_end);
//		//copy the data to buffer
//		for(int n = 0;n < my_n_interval;n++)
//			memcpy(buffer+n*dist_d_interval,Mat+n*D+d_start,sizeof(float)*1*dist_d_interval);
//
//		//gatehr the data
//		int count = my_n_interval*dist_d_interval;
//		//MPI::COMM_WORLD.Barrier();
//		MPI::COMM_WORLD.Gather(&count,1,MPI::INT,rcount,1,MPI::INT,dist_node);
//		if(node == dist_node){
//			//std::cout<<rcount[0]<<" "<<rcount[1]<<std::endl;
//			displs[0] = 0;
//			for(int i = 1;i < nodes;i++){
//				displs[i] = displs[i-1]+rcount[i-1];
//				//std::cout<<displs[i]<<" ";
//			}
//			//std::cout<<"size = "<<dist_d_interval*N<<std::endl;
//			//double temp = 0;
//			//for(int i = 0;i < count;i++)
//			//	temp += buffer[i];
//			//std::cout<<"buffer = "<<temp<<std::endl;
//		}
//		//MPI::COMM_WORLD.Barrier();
//
//		//double temp1 = 0;
//		//for(int i = 0;i < my_n_interval*dist_d_interval;i++){
//		//	temp1 += buffer[i];
//		//}
//		//for(int i = 0;i < nodes;i++){
//		//	if(node == i) std::cout<<node<<" "<<temp1<<" "<<std::endl<<std::flush;
//		//	MPI::COMM_WORLD.Barrier();
//		//	usleep(100);
//		//}
//
//		MPI::COMM_WORLD.Gatherv(buffer,my_n_interval*dist_d_interval,MPI::FLOAT,newMat,rcount,displs,MPI::FLOAT,dist_node);
//		
//		//if(node == dist_node){
//		//	double temp = 0;
//		//	for(int i = 0;i < N*my_d_interval;i++)
//		//		temp += newMat[i];
//		//	std::cout<<"in tranpsoe "<<temp<<std::endl;
//		//}
//		//tranpose the data(N*my_d_interval) to data(my_d_interval*N)
//		if(node == dist_node)
//			mkl_simatcopy('R','T',N,my_d_interval,1,newMat,my_d_interval,N);
//		MPI::COMM_WORLD.Barrier();
//	}
//
//	_mm_free(buffer);
//	//_mm_free(newMatbuffer);
//	_mm_free(rcount);
//	_mm_free(displs);
//}
//
///**** it will tranpose matrix(N*D) between cluster,Mat is (N/nodes)xD cross each node,new Mat is N*(D/nodes) ****/
//void transposeMatClusterIO(float* Mat,float *newMat,int N,int D,int node,int nodes){
//
//	int n_interval,n_start,n_end;
//	n_interval = divideIndex(N,nodes,node,n_start,n_end);
//	int d_interval,d_start,d_end;
//	d_interval = divideIndex(D,nodes,node,d_start,d_end);
//
//	std::string tempFileName = "transposeMatClusterIO_temp_Mat.data";
//	int amode = MPI::MODE_RDWR| MPI::MODE_CREATE | MPI::MODE_DELETE_ON_CLOSE;
//	MPI::File datafile = MPI::File::Open(MPI::COMM_WORLD,tempFileName.c_str(),amode,MPI::INFO_NULL);
//	//write the Mat data to disk first
//	datafile.Write_at(sizeof(float)*n_start*D,Mat,n_interval*D,MPI::FLOAT);
//
//	//read the data
//	MPI::Datatype myFileType = MPI::FLOAT.Create_vector(N,d_interval,D);
//	myFileType.Commit();
//	datafile.Set_view(d_start*sizeof(float),MPI::FLOAT,myFileType,"native",MPI::INFO_NULL);
//
//	datafile.Read_all(newMat,N*d_interval,MPI::FLOAT);
//
//	mkl_simatcopy('R','T',N,d_interval,1,newMat,d_interval,N);
//
//	datafile.Close();
//}
//
//
///**** gather specific rows of Mat(NxD) from each node,the rowIndex store the specific row index and 
//Mat size is (N/nodes)xD in each node ,the newMat size is number_of_rowxD ****/
//void gaterOneMatCluster(float* Mat,float *newMat,int N,int D,int *rowIndex,int rowIndexLength,int dist_node,int node,int nodes){
//	int *rcount = (int*)_mm_malloc(sizeof(int)*nodes,64);
//	int *displs = (int*)_mm_malloc(sizeof(int)*nodes,64);
//
//	int count,bufferSize;
//	int my_n_interval,n_start,n_end;
//	my_n_interval = divideIndex(N,nodes,node,n_start,n_end);
//
//	//calculate hwo many rows in this node
//	bufferSize = 0;
//	for(int i = 0;i < rowIndexLength;i++)
//		if(n_start <= rowIndex[i] && rowIndex[i] < n_end) bufferSize++;
//
//	float *buffer = (float*)_mm_malloc(sizeof(float)*bufferSize*D,64);
//	//std::cout<<"node = "<<node<<",bufferSize = "<<bufferSize<<" "<<n_start<<" "<<n_end<<std::endl;
//	//copy the data to buffer
//	count = 0;
//	for(int i = 0;i < rowIndexLength;i++)
//		if(n_start <= rowIndex[i] && rowIndex[i] < n_end){
//			memcpy(buffer+count*D,Mat+(rowIndex[i]-n_start)*D,sizeof(float)*1*D);
//			count++;
//		}
//
//	//gather the data to dist_node
//	count = bufferSize*D;
//	MPI::COMM_WORLD.Gather(&count,1,MPI::INT,rcount,1,MPI::INT,dist_node);
//	if(node == dist_node){
//		displs[0] = 0;
//		for(int i = 1;i < nodes;i++)
//			displs[i] = displs[i-1]+rcount[i-1];
//	}
//	MPI::COMM_WORLD.Gatherv(buffer,bufferSize*D,MPI::FLOAT,newMat,rcount,displs,MPI::FLOAT,dist_node);
//	MPI::COMM_WORLD.Barrier();
//
//	_mm_free(rcount);
//	_mm_free(displs);
//	_mm_free(buffer);
//}
//
//
///**** using the file IO to gather specific row of Matrix for each node,
//be careful !!!!!!,make sure your rowIndex for each node is form small to large   ****/
//void gaterMatClusterIO(float* Mat,float *newMat,int N,int D,int *rowIndex,int *rowIndexLength,int node,int nodes,bool updateTemp)
//{
//	int n_interval,n_start,n_end;
//	n_interval = divideIndex(N,nodes,node,n_start,n_end);
//
//	std::string tempFileName = "gaterMatClusterIO_temp_Mat.data";
//	int amode;
//	if(updateTemp) 
//		amode = MPI::MODE_RDWR| MPI::MODE_CREATE;
//	else
//		amode = MPI::MODE_RDWR;
//	MPI::File datafile = MPI::File::Open(MPI::COMM_WORLD,tempFileName.c_str(),amode,MPI::INFO_NULL);
//	//write the Mat data to disk first
//	if(updateTemp) datafile.Write_at(sizeof(float)*n_start*D,Mat,n_interval*D,MPI::FLOAT);
//
//	
//	//get displacements
//	int count = rowIndexLength[node];
//	int displacements[count];
//	int rowIndex_offest = 0;
//	for(int i = 0;i < node;i++)
//		rowIndex_offest += rowIndexLength[i];
//	for(int i = 0;i < count;i++)
//		displacements[i] = rowIndex[rowIndex_offest+i]*D;
//
//	//read the data
//
//	MPI::Datatype myFileType = MPI::FLOAT.Create_indexed_block(count,D,displacements);
//	myFileType.Commit();
//	datafile.Set_view(0,MPI::FLOAT,myFileType,"native",MPI::INFO_NULL);
//
//	datafile.Read_all(newMat,count*D,MPI::FLOAT);
//	
//	datafile.Close();
//
//}
//
///**** divided the Matrix(N*D) in cluster,each Mat in node will have different number of Row for Matrix
//the newMat is equally divided by nodes****/
//void divideMatClusterIO(float* Mat,float *newMat,int *subN,int D,int node,int nodes){
//
//	std::string tempFileName = "divideMatClusterIO_temp_Mat.data";
//	int amode = MPI::MODE_RDWR| MPI::MODE_CREATE | MPI::MODE_DELETE_ON_CLOSE;
//	MPI::File datafile = MPI::File::Open(MPI::COMM_WORLD,tempFileName.c_str(),amode,MPI::INFO_NULL);
//	//write the Mat data to disk first
//	int subN_start = 0;
//	for(int i = 0;i < node;i++)
//		subN_start += subN[i];
//	datafile.Write_at(sizeof(float)*subN_start*D,Mat,subN[node]*D,MPI::FLOAT);
//
//	//divide the Mat by nodes
//	int N = 0;
//	for(int i = 0;i < nodes;i++)
//		N += subN[i];
//	int n_interval,n_start,n_end;
//	n_interval = divideIndex(N,nodes,node,n_start,n_end);
//
//
//	//read the data
//	MPI::Datatype myFileType = MPI::FLOAT.Create_contiguous(n_interval);
//	myFileType.Commit();
//	datafile.Set_view(sizeof(float)*n_start*D,MPI::FLOAT,myFileType,"native",MPI::INFO_NULL);
//
//	datafile.Read_all(newMat,n_interval*D,MPI::FLOAT);
//
//	datafile.Close();
//}
//#endif