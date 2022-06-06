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

#ifndef MPI_H_
#define MPI_H_

#include <cmath>
#include <iostream>
#include <cstring> /*memcpy*/
#include "mkl.h"

#include "./macros.h"

// image division strategy for multi-nodes will be same as one-node.
// for Map2d_optimizer_old,Map3d_optimizer_old
//#define SAME_IMAGE_DIVISION

// build with MPI
//#define USEMPI

#ifdef USEMPI
#include "./util_heap_undefs.h"
#include <mpi.h>
#include "./util_heap_defs.h"
#endif


#ifdef USEMPI
# define NODEONLY(__node,__flag) if(__flag == __node)
# define NODE0ONLY NODEONLY(node,0)
# define MASTERNODE if (0 == MPI::COMM_WORLD.Get_rank()) // when not define node
#else
# define NODEONLY(__node,__flag) if(true)
# define NODE0ONLY NODEONLY(node,0)
# define MASTERNODE
#endif

// Support compile-time checking of the MPI code
//
#ifndef USEMPI
#define FAKEMPI
#endif

#ifdef FAKEMPI
//#define USEMPI UsingFakeMPI
static const int MPI_MAX_PROCESSOR_NAME = 256;
namespace MPI {
    //
    static void Init() {}
    static void Finalize() {}
    static void Get_processor_name(char*,int) {}
    //
    class Datatype {
    public:
        static Datatype Create_contiguous(int) {}
        void Commit() {}
        void Free() {}
        int Get_size() const {}
        bool operator==(Datatype&) const {}
    };
    extern Datatype DOUBLE,FLOAT,INT,BOOL;
    //
    class User_function;
    class Op {
    public:
        void Init(User_function*, bool) {}
        void Free() {}
    };
    extern Op SUM;
    class Intracomm {
    public:
        static void Barrier() {}
        static void Bcast(void*, int, Datatype, int) {}
        static void Send(const void*, int, Datatype, int, int) {}
        static void Recv(void*, int, Datatype, int, int) {}
        static void Reduce(void*, void*, int, Datatype, Op, int) {}
        static void Allreduce(void*, void*, int, Datatype, Op) {}
        static void Gatherv(void*, int, Datatype, void*, int*, int*, Datatype, int) {}
        static void Reduce_scatter(void*, void*, int*, Datatype, Op) {}
        static Intracomm Split(int,int) {}
        int Get_size() {return 1;}
        int Get_rank() {return 0;}
    };
    extern Intracomm COMM_WORLD;
};
#endif

#if defined(FLOAT_PRECISION)
#define MPI_FDOUBLE MPI::FLOAT
#else
#define MPI_FDOUBLE MPI::DOUBLE 
#endif


// Divides a number into most equally groups
int divide_equally(int N, int size, int rank,int &first,int &last);

// Devides a number into most equally groups
// consider the nr_pool in map,this will give same Devided Strategy for different nodes
int divide_equally_pool(int N, int nr_pool, int size, int rank,int &first,int &last);

// In which group from divide_equally is myself?
int divide_equally_which_group(int N, int size,int myself);

// devide N by nodes
int divideIndex(int N,int nodes,int node,int &n_start,int &n_end);


//#if 0
//#include <mpi.h>
//
///**** out-place tranposition of Mat(NxD) in cluster,the Mat size is (N/node)xD in each node,
//and newMat size is Dx(D/nodes) ****/
//void transposeMatCluster(double* Mat,double *newMat,int N,int D,int node,int nodes);
//void transposeMatCluster(float* Mat,float *newMat,int N,int D,int node,int nodes);
//
///**** using the file IO to tranpose matrix in cluster ****/
//void transposeMatClusterIO(float* Mat,float *newMat,int N,int D,int node,int nodes);
//
///**** gather specific rows of Mat(NxD) from each node,the rowIndex store the specific row index and 
//Mat size is (N/nodes)xD in each node ,the newMat size is number_of_rowxD ****/
//void gaterOneMatCluster(float* Mat,float *newMat,int N,int D,int *rowIndex,int rowIndexLength,int dist_node,int node,int nodes);
//
///**** using the file IO to gather specific row of Matrix for each node,
//be careful !!!!!!,make sure your rowIndex for each node is form small to large   ****/
//void gaterMatClusterIO(float* Mat,float *newMat,int N,int D,int *rowIndex,int *rowIndexLength,int node,int nodes,bool updateTemp = true);
//
///**** divided the Matrix(N*D) in cluster,each Mat in node will have different number of Row for Matrix
//the newMat is equally divided by nodes****/
//void divideMatClusterIO(float* Mat,float *newMat,int *subN,int D,int node,int nodes);
//#endif

#endif