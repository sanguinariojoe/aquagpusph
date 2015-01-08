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

#include <CalcServer/Tool.h>
#include <CalcServer.h>
#include <sys/time.h>

namespace Aqua{ namespace CalcServer{

Tool::Tool(const char* tool_name)
    : _name(NULL)
    , _allocated_memory(0)
    , _n_iters(0)
    , _elapsed_time(0.f)
    , _average_elapsed_time(0.f)
    , _squared_elapsed_time(0.f)
{
    name(tool_name);
}

Tool::~Tool()
{
    if(_name) delete[] _name; _name=NULL;
}

void Tool::name(const char* tool_name)
{
    if(_name) delete[] _name; _name=NULL;
    _name = new char[strlen(tool_name) + 1];
    strcpy(_name, tool_name);
}

bool Tool::execute()
{
    timeval tic, tac;
    gettimeofday(&tic, NULL);

    bool err = _execute();

    gettimeofday(&tac, NULL);

    float elapsed_seconds;
    elapsed_seconds = (float)(tac.tv_sec - tic.tv_sec);
    elapsed_seconds += (float)(tac.tv_usec - tic.tv_usec) * 1E-6f;

    addElapsedTime(elapsed_seconds);

    return err;
}

void Tool::addElapsedTime(float elapsed_time)
{
    _elapsed_time = elapsed_time;
    // Invert the average computation
    _average_elapsed_time *= _n_iters;
    _squared_elapsed_time *= _n_iters;
    // Add the new data
    _average_elapsed_time += elapsed_time;
    _squared_elapsed_time += elapsed_time * elapsed_time;
    // And average it again
    _n_iters++;
    _average_elapsed_time /= _n_iters;
    _squared_elapsed_time /= _n_iters;
}

}}  // namespace
