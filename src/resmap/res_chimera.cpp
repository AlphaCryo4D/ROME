/***************************************************************************
 *
 * Intel® Parallel Computing Center for Structural Biology
 * Principal Investigator : Youdong (Jack) Mao (Youdong_Mao@dfci.harvard.edu)
 * Dana-Farber Cancer Institute, Harvard Medical School and Peking University
 *
 * Authors: "Jian Wang(wj_hust08@hust.edu.cn)"
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

#include "res_chimera.h"

namespace resumap {

    std::string createChimeraScript(std::string inputFileName, double Mbegin, double Mmax, int N, bool animated) {

        std::string fname, ext, basename, basenameRAW;

        std::tie(fname, ext) = os::path::splitext(inputFileName);
        basename = os::path::basename(inputFileName);
        std::tie(basenameRAW, ext) = os::path::splitext(basename);

        std::string newFileName = fname + "_resmap_chimera.cmd";

        std::ofstream f(newFileName.c_str());

        f << "# Color a map by local resolution computed by ResMap\n";
        f << "#   Initial version of script courtesy of Tom Goddard, UCSF.\n\n";

        f << "# Open both volumes and hide the ResMap volume.\n";
        f << "set bg_color white\n";
        f << "open #0 " << basename << "\n";
        f << "open #1 " << basenameRAW << "_resmap" << ext << "\n";
        f << "volume #1 hide\n\n";

        // Define color mapping
        std::vector<std::string> colorsChimera {"blue","cyan","green","yellow","orange","red"};
        auto && valuesChimera = Arrayd::range(Mbegin, Mmax, 0, int(colorsChimera.size()));

        // Create a string to input to Chimera's scolor command
        std::stringstream stream;
        stream << std::fixed << std::setprecision(2);
        for (int i = 0; i < colorsChimera.size(); i++) {
            stream << valuesChimera[i] << ',' << colorsChimera[i] << ':';
        }
        stream << valuesChimera[valuesChimera.size()-1]+0.01 << ',' << "gray";
        std::string colorStr = stream.str();

        f << "# Color the original map with the values from the ResMap output.\n";
        f << "scolor #0 volume #1 cmap " << colorStr << "\n\n\n\n";

        f << "# OPTIONAL: ResMap Slice Animation.\n\n";

        f << "# Show midway slice of the original map with contour level below the minimum map value.\n";
        if (!animated) f << "# ";
        f << "volume #0 planes z," << N/2 << " step 1 level -1 style surface\n\n";

        f << "# Show a smooth transparent contour surface indicating the structure boundaries.\n";
        if (!animated) f << "# ";
        f << "vop gaussian #0 sDev 5 model #2\n";
        if (!animated) f << "# ";
        f << "volume #2 level 0.02 step 1 color .9,.7,.7,.5\n\n";

        f << "# Zoom out a bit and tilt to a nice viewing angle.\n";
        if (!animated) f << "# ";
        f << "turn x -45\n";
        if (!animated) f << "# ";
        f << "turn y -30\n";
        if (!animated) f << "# ";
        f << "turn z -30\n";
        if (!animated) f << "# ";
        f << "scale 0.5\n\n";

        f << "# Cycle through planes from N/2 to 4N/5 up to N/5 and back to N/2.\n";
        if (!animated) f << "# ";
        f << "volume #0 planes z," << N/2 << ',' << 4*N/5 << ",0.25\n";
        if (!animated) f << "# ";
        f << "wait " << 4*int(4*N/5 - N/2) << '\n';
        if (!animated) f << "# ";
        f << "volume #0 planes z," << 4*N/5 << ',' << N/5 << ",0.25\n";
        if (!animated) f << "# ";
        f << "wait " << 4*int(4*N/5 - N/5) << "\n";
        if (!animated) f << "# ";
        f << "volume #0 planes z," << N/5 << ',' << N/2 << ",0.25\n";

        f.close();

        return newFileName;

    }


}

