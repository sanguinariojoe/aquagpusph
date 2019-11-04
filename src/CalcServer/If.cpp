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
 * @brief Check a condition to enable/disable all the tools in its scope
 * (See Aqua::CalcServer::If and Aqua::CalcServer::EndIf for details)
 */

#include <math.h>

#include <AuxiliarMethods.h>
#include <InputOutput/Logger.h>
#include <CalcServer.h>
#include <CalcServer/If.h>

namespace Aqua{ namespace CalcServer{

If::If(const std::string name, const std::string condition, bool once)
    : Tool(name, once)
    , _condition(condition)
    , _result(true)
{
}

If::~If()
{
}

void If::setup()
{
    std::ostringstream msg;
    msg << "Loading the tool \"" << name() << "\"..." << std::endl;
    LOG(L_INFO, msg.str());
}

cl_event If::_execute(const std::vector<cl_event> events)
{
    int result;
    InputOutput::Variables *vars = CalcServer::singleton()->variables();

    void *data = malloc(sizeof(int));
    if(!data){
        std::stringstream msg;
        msg << "Failure allocating memory for the integer result" << std::endl;
        LOG(L_ERROR, msg.str());
        throw std::bad_alloc();
    }

    vars->solve("int", _condition, data, "if_result");

    // Check the result
    memcpy(&result, data, sizeof(int));
    free(data);
    _result = result != 0;

    return NULL;
}

EndIf::EndIf(const std::string name, bool once)
    : Tool(name, once)
{
}

EndIf::~EndIf()
{
}

void EndIf::setup()
{
    std::ostringstream msg;
    msg << "Loading the tool \"" << name() << "\"..." << std::endl;
    LOG(L_INFO, msg.str());
}

}}  // namespaces
