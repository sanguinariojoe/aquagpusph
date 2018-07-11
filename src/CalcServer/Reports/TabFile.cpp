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
 * @brief Runtime output file.
 * (See Aqua::CalcServer::Reports::TabFile for details)
 */

#include <AuxiliarMethods.h>
#include <InputOutput/Logger.h>
#include <CalcServer/Reports/TabFile.h>

namespace Aqua{ namespace CalcServer{ namespace Reports{

TabFile::TabFile(const std::string tool_name,
                 const std::string fields,
                 const std::string output_file)
    : Report(tool_name, fields)
    , _output_file(output_file)
{
}

TabFile::~TabFile()
{
    if(_f.is_open()) _f.close();
}

void TabFile::setup()
{
    unsigned int i;

    std::ostringstream msg;
    msg << "Loading the report \"" << name() << "\"..." << std::endl;
    LOG(L_INFO, msg.str());

    _f.open(_output_file.c_str(), std::ios::out);

    Report::setup();

    // Write the header
    _f << "# ";
    std::vector<InputOutput::Variable*> vars = variables();
    for(i = 0; i < vars.size(); i++){
        _f << vars.at(i)->name() << " ";
    }
    _f << std::endl;
    _f.flush();
}

void TabFile::_execute()
{
    // Change break lines by spaces
    std::string out = replaceAllCopy(data(false, false), "\n", " ");
    _f << out << std::endl;
    _f.flush();
}

}}} // namespace
