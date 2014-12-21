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

#include <CalcServer/Reports/TabFile.h>
#include <ScreenManager.h>

namespace Aqua{ namespace CalcServer{ namespace Reports{

TabFile::TabFile(const char* tool_name,
                 const char* fields,
                 const char* output_file)
    : Report(tool_name, fields)
    , _output_file(NULL)
    , _f(NULL)
{
    _output_file = new char[strlen(output_file) + 1];
    strcpy(_output_file, output_file);
}

TabFile::~TabFile()
{
    if(_output_file) delete[] _output_file; _output_file=NULL;
    if(_f) fclose(_f); _f = NULL;
}

bool TabFile::setup()
{
    unsigned int i;
    char msg[1024];
    InputOutput::ScreenManager *S = InputOutput::ScreenManager::singleton();

    sprintf(msg,
            "Loading the report \"%s\"...\n",
            name());
    S->addMessageF(1, msg);

    // Open the output file
    _f = fopen(_output_file, "w");
    if(!_f){
        sprintf(msg,
                "The file \"%s\" cannot be written\n",
                _output_file);
        S->addMessageF(3, msg);
        return true;
    }

    if(Report::setup()){
        return true;
    }

    // Write the header
    fprintf(_f, "# ");
    std::deque<InputOutput::Variable*> vars = variables();
    for(i = 0; i < vars.size(); i++){
        fprintf(_f, "%s ", vars.at(i)->name());
    }
    fprintf(_f, "\n");

    return false;
}

bool TabFile::_execute()
{
    char* out = (char*)data(false, false);
    // Change break lines by spaces
    while(strchr(out, '\n')){
        strchr(out, '\n')[0] = ' ';
    }

    fprintf(_f, "%s\n", out);
    fflush(_f);
    return false;
}

}}} // namespace
