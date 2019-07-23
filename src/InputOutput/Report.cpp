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
#include <InputOutput/Logger.h>
#include <InputOutput/Report.h>
#include <AuxiliarMethods.h>

#ifdef HAVE_MPI
#include <mpi.h>
#endif

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

    // Start replacing all the old-school formatting string instances by the new
    // one, based on a more intelligible variable name
    newname = replaceAllCopy(basename, "%d", "{index}");

    // Now replace all the instances of constant variables
    int mpi_rank = 0;
    #ifdef HAVE_MPI
        try {
            mpi_rank = MPI::COMM_WORLD.Get_rank();
        } catch(MPI::Exception e){
            std::ostringstream msg;
            msg << "Error getting MPI rank. " << std::endl
                << e.Get_error_code() << ": " << e.Get_error_string() << std::endl;
            LOG(L_ERROR, msg.str());
            throw;
        }
    #endif
    number_str.str(""); number_str << mpi_rank;
    replaceAll(newname, "{mpi_rank}", number_str.str());
    
    if(newname.find("{index}") == std::string::npos){
        // We cannot insert the file index anywhere, so just test if the file
        // does not exist
        f = fopen(newname.c_str(), "r");
        if(f){
            std::ostringstream msg;
            msg << "A report is trying to overwrite an existing file. "
                << std::endl << "'" << newname << "'" << std::endl;
            LOG(L_ERROR, msg.str());
            fclose(f);
            throw std::invalid_argument("Invalid file name pattern");
        }

        file(newname);
        return;
    }

    while(true){
        number_str.str(""); number_str << i;
        replaceAll(newname, "{index}", number_str.str());

        f = fopen(newname.c_str(), "r");
        if(!f)
            break;
        fclose(f);
        i++;
    }

    file(newname);
}

}}  // namespace
