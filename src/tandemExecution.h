#pragma once

//  Author: Bevin R Brett  2012-09-21
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

#include "./util.h"

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

