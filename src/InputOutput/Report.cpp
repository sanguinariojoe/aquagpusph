/*
 *  This file is part of AQUAgpusph, a free CFD program based on SPH.
 *  Copyright (C) 2012  Jose Luis Cercos Pita <jl.cercos@upm.es>
 *
 *  AQUAgpusph is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  (at your option) any later version.
 *
 *  AQUAgpusph is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with AQUAgpusph.  If not, see <http://www.gnu.org/licenses/>.
 */

/** @file
 * @brief Base class for all the report file managers.
 * (See Aqua::InputOutput::Report for details)
 */

#include <sstream>
#include <iostream>
#include <fstream>
#include <InputOutput/Report.h>
#include <AuxiliarMethods.h>

namespace Aqua{ namespace InputOutput{

Report::Report()
{
}

Report::~Report()
{
}

void Report::file(const std::string& basename,
                  const unsigned int& startindex)
{
    if(basename.find("%d") == std::string::npos){
        // We cannot replace nothing in the basename, just test if the file
        // does not exist
        std::ifstream f(basename.c_str());
        if(f.good()){
            f.close();
            throw std::runtime_error("Bad file name");
        }

        file(basename);
        return;
    }

    unsigned int i = startindex;
    while(true){
        std::ostringstream number_str;
        number_str.str("");
        number_str << i;
        std::string newname = replaceAllCopy(basename, "%d", number_str.str());

        std::ifstream f(newname.c_str());
        if(!f.good()){
            file(newname);
            break;
        }
        f.close();

        i++;
    }
}

}}  // namespace
