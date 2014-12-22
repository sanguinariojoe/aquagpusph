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
 * @brief Performance report.
 * (See Aqua::CalcServer::Reports::Performance for details)
 */

#include <CalcServer/Reports/Performance.h>
#include <CalcServer.h>
#include <ScreenManager.h>

namespace Aqua{ namespace CalcServer{ namespace Reports{

Performance::Performance(const char* tool_name,
                         const char* color,
                         bool bold)
    : Report(tool_name, "dummy_fields_string")
    , _color(NULL)
    , _bold(bold)
{
    _color = new char[strlen(color) + 1];
    strcpy(_color, color);
}

Performance::~Performance()
{
    if(_color) delete[] _color; _color=NULL;
}

bool Performance::setup()
{
    unsigned int i;
    char msg[1024];
    InputOutput::ScreenManager *S = InputOutput::ScreenManager::singleton();

    sprintf(msg,
            "Loading the report \"%s\"...\n",
            name());
    S->addMessageF(1, msg);

    // Set the color in lowercase
    for(i = 0; i < strlen(_color); i++){
        _color[i] = tolower(_color[i]);
    }

    return false;
}

size_t Performance::computeAllocatedMemory(){
    unsigned int i;
    size_t allocated_mem = 0;
    CalcServer *C = CalcServer::singleton();
    InputOutput::ScreenManager *S = InputOutput::ScreenManager::singleton();

    // Get the allocated memory in the variables
    InputOutput::Variables* vars = C->variables();
    allocated_mem += vars->allocatedMemory();

    // Gwet the additionally allocated memory in the tools
    std::deque<Tool*> tools = C->tools();
    for(i = 0; i < tools.size(); i++){
        allocated_mem += tools.at(i)->allocatedMemory();
    }

    return allocated_mem;
}

bool Performance::_execute()
{
    unsigned int i;
    CalcServer *C = CalcServer::singleton();
    InputOutput::ScreenManager *S = InputOutput::ScreenManager::singleton();
    char data[4096];

    size_t allocated_MB = computeAllocatedMemory() / (1024 * 1024);
    sprintf(data, "Performance:\nMemory=%16luMB\n", allocated_MB);

    // Add the tools time elapsed
    std::deque<Tool*> tools = C->tools();
    float elapsed = 0.f;
    float elapsed_var = 0.f;
    for(i = 0; i < tools.size(); i++){
        elapsed += tools.at(i)->elapsedTime();
        elapsed_var += tools.at(i)->elapsedTimeVariance();
    }
    sprintf(data, "%sElapsed=%16gs (+-%16gs)\n", data, elapsed, elapsed_var);


    // Replace the trailing space by a line break
    if(data[strlen(data) - 1] == ' ')
        data[strlen(data) - 1] = '\n';

    S->writeReport(data, _color, _bold);
    return false;
}

}}} // namespace
