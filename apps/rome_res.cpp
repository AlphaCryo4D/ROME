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

#include "../src/resmap/resmap_mpi.h"
#include "../src/resmap/resmap_mrcs.h"
#include "../src/resmap/resmap_option.h"
#include "../src/resmap/resmap_component.h"

#include "../src/resmap/res_calculator.h"
#include "../src/resmap/res_array.h"

int rome_desc(int argc, char **argv)
{
    Option option;
    
    // general option
    option.addOption("-desc", "run desc",	"0");
    option.addOption("-i"   , "MRC file.", 	"" );
    
    option.readCommandLine(argc, argv);
    
    bool run_desc = option.getBoolOption("-desc");
    ERROR_CHECK(!run_desc, "Wrong -desc option.");
    std::string inFileName  = option.getStrOption("-i");
    
    // Read volume
    Mrcs::MrcVolume volume;
    volume.read(inFileName);
    
    return 0;
}

int rome_res(int argc, char **argv)
{
    Option option;

    // general option
    option.addOption("-res"				, "run res"				, "0"	 );
    option.addOption("-i"               , "MRC file."           , ""     );
    option.addOption("-i1"              , "MRC file."           , ""     );
    option.addOption("-i2"              , "MRC file."           , ""     );
    option.addOption("-vxSize"          , "voxcel size"         , "0.0"  );
    option.addOption("-pValue"          , "pValue"              , "0.05" );
    option.addOption("-minRes"          , "min resolution"      , "0.0"  );
    option.addOption("-maxRes"          , "max resolution"      , "5.0"  );
    option.addOption("-stepRes"         , "step of resolution"  , "1.0"  );
    option.addOption("-variance"        , "variance"            , "0.0"  );
    // option.addOption("-noiseDiagnostics", "noise"               , "false");
    
    if (argc < 2) {
        option.printHelp();
        std::cerr << "example : " << std::endl;
        std::cerr << std::endl;
        EXIT_ABNORMALLY;
    }

    option.readCommandLine(argc, argv);

    bool run_res = option.getBoolOption("-res");
    ERROR_CHECK(!run_res, "Wrong -res option.");
    
    option.printValue();

    double t1 = dtime();

    resumap::inputFileName    = option.getStrOption  ("-i"               );
    resumap::inputFileName1   = option.getStrOption  ("-i1"              );
    resumap::inputFileName2   = option.getStrOption  ("-i2"              );
    resumap::vxSize           = option.getFloatOption("-vxSize"          );
    resumap::pValue           = option.getFloatOption("-pValue"          );
    resumap::minRes           = option.getFloatOption("-minRes"          );
    resumap::maxRes           = option.getFloatOption("-maxRes"          );
    resumap::stepRes          = option.getFloatOption("-stepRes"         );
    resumap::variance         = option.getFloatOption("-variance"        );
    // resumap::noiseDiagnostics = option.getBoolOption ("-noiseDiagnostics");

    resumap::calculate();

    double t2 = dtime();

    MPI_LOG << "Local resolution calculation costs : " << t2 - t1 << " seconds." << std::endl;

    return 0;
}

int rome_zoom(int argc, char **argv)
{
    Option option;

    // general option
    option.addOption("-zoom"			, "run zoom"				, "0"	);
    option.addOption("-i"               , "MRC file."               , ""   );
    option.addOption("-n"               , "New size after zooming." , ""  );
    option.addOption("-o"               , "Output file."            , ""   );
    
    option.readCommandLine(argc, argv);

    bool run_zoom = option.getBoolOption("-zoom");
    ERROR_CHECK(!run_zoom, "Wrong -zoom option.");
    
    std::string inFileName  = option.getStrOption("-i");
    int size = option.getIntOption("-n");
    std::string outFileName  = option.getStrOption("-o");

    // Read volume
    Mrcs::MrcVolume volume;
    volume.read(inFileName);

    // Zoom
    Arrayf data = cppsci::ndimage::zoom(MapArrayf(volume.data, {volume.size, volume.size, volume.size}), double(size)/volume.size, "reflect");

    // Set new mrc head
    auto head = *(volume.head);
    head.NC = size;
    head.NR = size;
    head.NS = size;
    head.NX = size;
    head.NY = size;
    head.NZ = size;

    // Write new volume
    writeMrcData(outFileName, data.data_wptr(), data.shape(0), head);

    return 0;
}

// compare twn mrc file
int rome_compare(int argc, char **argv)
{
    Option option;
    
    // general option
    option.addOption("-compare"			, "run compare" , "0" );
    option.addOption("-i1"               , "MRC file."  , ""  );
    option.addOption("-i2"               , "MRC file2." , ""  );
    
    option.readCommandLine(argc, argv);
    
    bool run_compare = option.getBoolOption("-compare");
    ERROR_CHECK(!run_compare, "Wrong -compare option.");
    
    std::string FileName1  = option.getStrOption("-i1");
    std::string FileName2  = option.getStrOption("-i2");
    
    // Read volume
    Mrcs::MrcVolume volume1;
    volume1.read(FileName1);
    Mrcs::MrcVolume volume2;
    volume2.read(FileName2);
    
    // compare
    auto head1 = *(volume1.head);
    auto head2 = *(volume2.head);
    ERROR_CHECK(head1.NC != head2.NC, "NC is diff.");
    ERROR_CHECK(head1.NR != head2.NR, "NR is diff.");
    ERROR_CHECK(head1.NS != head2.NS, "NS is diff.");
    ERROR_CHECK(head1.NX != head2.NX, "NX is diff.");
    ERROR_CHECK(head1.NY != head2.NY, "NY is diff.");
    ERROR_CHECK(head1.NZ != head2.NZ, "NZ is diff.");
    
    for (int i = 0; i < head1.NX*head1.NY*head1.NZ; i++) {
        if (volume1.data[i] != volume2.data[i]) {
            std::cout<<volume1.data[i]<<" "<<volume2.data[i]<<std::endl;
        }
    }
    
    return 0;
}

int main(int argc, char **argv)
{
    MPI_INITIALIZE(argc, argv);
    //
    RomeResComponents::addComponent("-desc", rome_desc);
	RomeResComponents::addComponent("-res", rome_res);
    RomeResComponents::addComponent("-zoom", rome_zoom);
    RomeResComponents::addComponent("-compare", rome_compare);
    
    RomeResComponents::runComponent(argc, argv);
    
    MPI_FINALIZE;
}

