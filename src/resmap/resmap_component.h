/***************************************************************************
 *
 * Intel® Parallel Computing Center for Structural Biology
 * Principal Investigator : Youdong (Jack) Mao (Youdong_Mao@dfci.harvard.edu)
 * Dana-Farber Cancer Institute, Harvard Medical School and Peking University
 *
 * Authors: "Jian Wang(wj_hust08@hust.edu.cn) Yong Bei Ma(galowma@gmail.com)"
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

#ifndef COMPONENT_H_
#define COMPONENT_H_

#include <iostream>
#include <map>
#include <functional>
#include <string>

#include "resmap_mpi.h"
//#define JN_PP_CAT(a, b) a##b
//#define JN_PP_STRING(a) #a

//#define REGISTER_ROME_COMPONENT(i)\
//int JN_PP_CAT(rome_, i)(int argc, char **argv);\
//ROMEComponent JN_PP_CAT(component_, i)(JN_PP_STRING(i), JN_PP_CAT(rome_, i));\
//int JN_PP_CAT(rome_, i)(int argc, char **argv)


class Components {
protected:
    ::std::map< ::std::string, ::std::function< int(int, char **)> > methods;
    virtual void show_help() = 0;
public:
    void run(int argc, char **argv) {
        auto &m = methods;
        if (argc > 1 && m.find(argv[1]) != m.end()) {
            m[argv[1]](argc, argv);
        }
        else {
            show_help();
        }
    }
};

class RomeResComponents : public Components {
private:
    RomeResComponents() {}
    void show_help(){
        if (MPI_IS_ROOT)
        {
            std::cout << "#### single volume mode" << std::endl;
            std::cout << " ./bin/rome_res -res -i <MRC_FILENAME> [-minRes <MINIMUM_RESOLUTION>] ";
            std::cout << "[-maxRes <MAX_RESOLUTION>] [-stepRes <STEP_OF_RESOLUTION>] " <<std::endl;
            std::cout << "example1 : ./bin/rome_res -res -i dataset/single-100.mrc " <<std::endl;
            std::cout << "example2 : ./bin/rome_res -res -i dataset/single-100.mrc -minRes 5 -maxRes 15 -stepRes 1.0 "<<std::endl;
            std::cout << "#### split volumes mode" << std::endl;
            std::cout << "./bin/rome_res -res -i1 <MRC_FILENAME1> -i2 <MRC_FILENAME2> [-minRes <MINIMUM_RESOLUTION>] ";
            std::cout << "[-maxRes <MAX_RESOLUTION>] [-stepRes <STEP_OF_RESOLUTION>]" <<std::endl;
            std::cout << "example1 : ./bin/rome_res -res -i1 dataset/split1-100.mrc -i2 dataset/split2-100.mrc "<<std::endl;
            std::cout << "example2 : ./bin/rome_res -res -i1 dataset/split1-100.mrc -i2 dataset/split2-100.mrc -minRes 10 -maxRes 20 -stepRes 1.0"<<std::endl;
        }

    }
    static RomeResComponents &instance() {
        static RomeResComponents components;
        return components;
    }
public:
    template<typename F>
    static void addComponent(const ::std::string &name, F &&f) {
        instance().methods[name] = f;
    }
    static void runComponent(int argc, char **argv) {
        instance().run(argc, argv);
    }
};

// TODO
class RomeToolComponents : public Components {
    
};
#endif
