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

#ifndef MPI_DEBUG_H_
#define MPI_DEBUG_H_

#include <fstream>

class MPIDebug {
    std::ofstream writeStream;
public:
    MPIDebug(){}
    ~MPIDebug() {
        writeStream.close();
    }
    void init() {
        int i = 0;//mpiRank();
        // writeStream.open("./debug/debug_node_"+std::to_string((long long)i)+".txt");
    }
    std::ofstream& getStream() {
        return writeStream;
    }
};

static MPIDebug mpidebug;

#endif /* defined(MPI_DEBUG_H_) */
