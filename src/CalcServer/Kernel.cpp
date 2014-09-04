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
 * @brief Aqua::CalcServer tools base class.
 * (See Aqua::CalcServer::Kernel for details)
 */

#include <CalcServer/Kernel.h>
#include <CalcServer.h>

namespace Aqua{ namespace CalcServer{

Kernel::Kernel(const char* kernel_name)
    : _name(0)
    #ifdef HAVE_GPUPROFILE
        , _time(0)
    #endif
{
    name(kernel_name);
}

Kernel::~Kernel()
{
    if(_name) delete[] _name; _name=0;
}

void Kernel::name(const char* kernel_name)
{
    if(_name) delete[] _name; _name=0;
    _name = new char[strlen(kernel_name)+1];
    strcpy(_name, kernel_name);
}

size_t Kernel::localWorkSize(unsigned int n, cl_command_queue queue)
{
    if(!n)
        n = CalcServer::CalcServer::singleton()->N;
    if(!queue)
        queue = CalcServer::CalcServer::singleton()->command_queue;
    return getLocalWorkSize(n,queue);
}

size_t Kernel::globalWorkSize(size_t size, unsigned int n)
{
    if(!n)
        n = CalcServer::CalcServer::singleton()->N;
    return getGlobalWorkSize(n, size);
}

}}  // namespace
