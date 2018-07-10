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

#include <sys/time.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>
#include <sstream>
#include <iostream>

#include <InputOutput/Report.h>
#include <AuxiliarMethods.h>

namespace Aqua{ namespace InputOutput{

Report::Report()
{
}

Report::~Report()
{
}

void Report::file(std::string filename)
{
    _output_file = filename;
}

void Report::file(std::string basename, unsigned int startindex)
{
    FILE *f;
    unsigned int i = startindex;
    std::string newname;
    std::ostringstream number_str;

    if(basename.find("%d") == std::string::npos){
        // We cannot replace nothing in the basename, just test if the file
        // does not exist
        f = fopen(basename.c_str(), "r");
        if(f){
            // The file already exist, so we cannot operate
            fclose(f);
            throw std::invalid_argument("Invalid file name pattern");
        }

        file(basename);
        return;
    }

    while(true){
        number_str.str("");
        number_str << i;
        newname = replaceAllCopy(basename, "%d", number_str.str());

        f = fopen(newname.c_str(), "r");
        if(!f)
            break;
        fclose(f);
        i++;
    }

    file(newname);
}

}}  // namespace
