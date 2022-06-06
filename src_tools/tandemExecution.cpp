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

#include "../src/resmap/resmap_util.h"		// used for building precompiled headers on Windows

#include "tandemExecution.h"


namespace TandemExecution {

	// Thread support
	//
	typedef class Thread* ThreadP;

	class Synchronize {
		Lock& lock;
	public:
		 Synchronize(Lock& lock) : lock(lock) { lock.acquire(__FILE__,__LINE__); }
		~Synchronize()						  { lock.release(); }
	};

#if defined(_WIN32)
	class Thread {
		DWORD   threadId;
		static const int numberOfThreads = 1;
		HANDLE hThreadArray[numberOfThreads];
	public:
		Thread() {}
		virtual ~Thread() {}
		void start() {
			hThreadArray[0] =
				CreateThread(
					NULL,                   // default security attributes
					0,                      // use default stack size
					runner,					// thread function name
					this,					// argument to thread function
					0,                      // use default creation flags
					&threadId);				// returns the thread identifier
		}
		void join() {
			WaitForMultipleObjects(numberOfThreads, hThreadArray, TRUE, INFINITE);
		}
	private:
		virtual void run() = 0;
		static DWORD WINAPI runner( LPVOID lpParam ) {
			Thread* thread = (Thread*)(lpParam);
			thread->run();
			return 0;
		}
	};
#else
    
#if 0
	 template <typename TYPE, void (TYPE::*_RunThread)() >
	 void* _thread_t(void* param) {
	 	TYPE* This = (TYPE*)param;
	 	This->_RunThread();
	 	return NULL;
	 }
	
	 class Thread {
	 	pthread_t _tid;
	 public:
	 	Thread() {}
	 	virtual ~Thread() {}
	 	void start() {
	 		pthread_create(&_tid, NULL, _thread_t<Thread, &Thread::_RunThread>, this);
	 	}
	 	void join() {
	 		pthread_join(_tid,NULL);
	 	}
	 	void _RunThread() {
	 		this->run();
	 	}
	 private:
	 	virtual void run() = 0;
	 };
#else
	class Thread {
		static const int numberOfThreads = 1;
		pthread_t pThreadArray[numberOfThreads];
	public:
		Thread() {}
		virtual ~Thread() {}
		void start() {
			pthread_create(
				&pThreadArray[0],			//thread_t
				NULL,
				runner,				//thread function
				this);				//thread function argument
		}
		void join() {
			pthread_join(pThreadArray[0],NULL);
		}
	private:
		virtual void run() = 0;
		static void *runner(void *Param)
		{
			Thread* thread = (Thread*)(Param);
		    thread->run();
		    return NULL;
		}
	};
#endif
    
#endif

	// Some glue...
	//
	typedef class AlgorithmThread : public Thread {
		AlgorithmP _algorithm;
	public:
		AlgorithmThread(AlgorithmP algorithm) : _algorithm(algorithm) {}
		~AlgorithmThread() { _algorithm = 0L; }
		virtual void run() { 
			_algorithm->run();
			_algorithm->waypoint("",0,true);
		}
	}* AlgorithmThreadP;


	// Algorithms
	//
	typedef std::vector<AlgorithmP>		  VectorOfAlgorithmP;
	typedef std::vector<AlgorithmThreadP> VectorOfAlgorithmThreadP;

	class Algorithms::Internals {
	public:
		VectorOfAlgorithmP		 _algorithmVector;
		VectorOfAlgorithmThreadP _threadVector;

		Lock					 _barrierEnter;
		Lock					 _barrierExit;
		unsigned int			 _alive;
		unsigned int			 _wild;
		unsigned int			 _captured;
        unsigned int             _tag;
		Internals() : _alive(0), _wild(0), _captured(0) ,_tag(0){}
	};

	Algorithms::Algorithms() : _internals(sNew(Algorithms::Internals)) {
	}

	Algorithms::~Algorithms() {
		sDelete(_internals);
	}
		
	size_t Algorithms::size() {
		return _internals->_algorithmVector.size();
	}

	Algorithm* Algorithms::algorithm(unsigned int i) {
		return _internals->_algorithmVector[i];
	}
		
	void Algorithms::run() {
		// Don't let them change state
		_internals->_barrierEnter.acquire(__FILE__,__LINE__);
		_internals->_barrierExit.acquire(__FILE__, __LINE__);

		// start all the threads, they will either finish or wait at the barrier
		for (unsigned int i = 0; i < size(); i++) {
			AlgorithmThreadP t = 
#include "../src/resmap/resmap_util_heap_undefs.h"
				sNewA(AlgorithmThread,(_internals->_algorithmVector[i]));
#include "../src/resmap/resmap_util_heap_defs.h"
			_internals->_threadVector.push_back(t);
			_internals->_alive++;
			_internals->_wild++;
			t->start();
		}

		while (_internals->_alive> 0) {
			// Allow all through the barrierEnter
			while (_internals->_wild > 0) {
				_internals->_barrierEnter.release();
				yield(true);
				_internals->_barrierEnter.acquire(__FILE__, __LINE__);
			}
			// Compare all that need to be compared
			if (_internals->_captured > 1) {
				compare(_internals->_tag);
			}
			// Allow all through the barrierExit
			while (_internals->_captured > 0) {
				_internals->_barrierExit.release();
				yield(true);
				_internals->_barrierExit.acquire(__FILE__, __LINE__);
			}
		}

		// Allow all through the barrierExit
		_internals->_barrierExit.release();
		yield(true);

		// join all the threads
		for (unsigned int i = 0; i < size(); i++) {
			AlgorithmThreadP & t = _internals->_threadVector[i];
			t->join();
			sDelete(t); 
		}
	}


	// Algorithm
	//
	Algorithm::Algorithm(Algorithms* algorithms) : _algorithms(algorithms) {
		_algorithms->_internals->_algorithmVector.push_back(this);
	}

	void Algorithm::waypoint(const char* message,unsigned int tag,bool finish) {
		auto i = _algorithms->_internals;

		i->_barrierEnter.acquire(__FILE__,__LINE__);
            std::cerr<<"starting : "<<message<<",tag : "<<tag<<std::endl<<std::flush;
            i->_tag = tag;
			i->_wild--;
			if (finish) i->_alive--; else i->_captured++;
		i->_barrierEnter.release();

		if (finish) return;

		i->_barrierExit.acquire(__FILE__, __LINE__);
            std::cerr<<"finishing : "<<message<<",tag : "<<tag<<std::endl<<std::flush;
			i->_captured--;
			i->_wild++;
		i->_barrierExit.release();
	}

	// UnitTest
	//
	class Algorithms01 : public Algorithms {
	public:

		class Algorithm0 : public Algorithm {
		public:
			Algorithm0(AlgorithmsP algorithms) : Algorithm(algorithms) {}
			void run() {
                waypoint("0",0);
                waypoint("1",1);
                waypoint("2",2);
                waypoint("3",3);
			}
		}* _algorithm0;

		class Algorithm1 : public Algorithm {
		public:
			Algorithm1(AlgorithmsP algorithms) : Algorithm(algorithms) {}
			void run() {
                waypoint("0",0);
                waypoint("1",1);
                waypoint("2",2);
                waypoint("3",3);
			}
		}* _algorithm1;

		Algorithms01() {
#include "../src/resmap/resmap_util_heap_undefs.h"
			_algorithm0 = sNewA(Algorithm0,(this));
			_algorithm1 = sNewA(Algorithm1,(this));
#include "../src/resmap/resmap_util_heap_defs.h"
		};

		~Algorithms01() {
			sDelete(_algorithm0);
			sDelete(_algorithm1);
		}
        
		virtual bool compare(unsigned int tag) {
			return true;
		}
	};

	void UnitTest() {
		std::cout << "TandemExecution::UnitTest" << std::endl;
		Algorithms01 algorithms;
		algorithms.run();
		std::cout << "TandemExecution::UnitTest finished" << std::endl;
	}
};

