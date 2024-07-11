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
 * @brief Just a small demo/template to create installable tools
 * (see Aqua::CalcServer::InstallableDemo for details)
 */

#include "tool.hpp"
#include <aquagpusph/InputOutput/Logger.hpp>

extern "C" Aqua::CalcServer::InstallableDemo* create_object(
    const std::string name, bool once)
{
    return new Aqua::CalcServer::InstallableDemo(name, once);
}

namespace Aqua{ namespace CalcServer{

InstallableDemo::InstallableDemo(const std::string name, bool once)
    : Tool(name, once)
{
}

InstallableDemo::~InstallableDemo()
{
}

cl_event
InstallableDemo::_execute(const std::vector<cl_event> events)
{
    LOG(L_INFO, "Executing Installed demo tool!\n");
    return NULL;
}

}}  // namespaces
 
