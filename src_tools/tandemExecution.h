/***************************************************************************
 *
 * Dana-Farber Cancer Institute, Harvard Medical School and Peking University
 * Intel® Parallel Computing Center for Structural Biology
 *
 * Authors: "Bevin R Brett(bevin_brett@hotmail.com) 2012-09-21"
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
#pragma once
//
//  Tandem Execution is used to help debug two implementations of an algorithm.
//
//  Consider two implementations, using differing data structures, but going through the same basic steps.
//
//  Run them on one thread each.
//  Have the code execute waypoints at places where the intermediate results should be the same.
//  Run each thread to the first waypoint, have it wait for the other to get there.
//  Compare their intermediate results, to make sure they are the same.
//  Continue until done...

#include "../src/resmap/resmap_util.h"

namespace TandemExecution {

    typedef class Algorithms* AlgorithmsP;
    typedef class Algorithm * AlgorithmP;

    class Algorithms {
        friend class Algorithm;
        class Internals;
        Internals* _internals;
    public:
        Algorithms();
        ~Algorithms();

        size_t size();
        AlgorithmP algorithm(unsigned int i);

        virtual bool compare(unsigned int tag) = 0;

        void run();
    };

    class Algorithm {
        AlgorithmsP _algorithms;
    public:
        Algorithm(AlgorithmsP algorithms);
        virtual void run() = 0;
        void waypoint(const char* message,unsigned int tag = 0,bool finish = false);        // true is for internal use only
    };

    void UnitTest();
};

