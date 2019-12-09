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
 * @brief Check that a condition holds true, or trhow a fatal error otherwise.
 * (See Aqua::CalcServer::Assert for details)
 */

#include <math.h>

#include <AuxiliarMethods.h>
#include <InputOutput/Logger.h>
#include <CalcServer.h>
#include <CalcServer/Assert.h>

namespace Aqua{ namespace CalcServer{

Assert::Assert(const std::string name, const std::string condition, bool once)
    : Tool(name, once)
    , _condition(condition)
{
}

Assert::~Assert()
{
}

void Assert::setup()
{
    std::ostringstream msg;
    msg << "Loading the tool \"" << name() << "\"..." << std::endl;
    LOG(L_INFO, msg.str());
    Tool::setup();
}

cl_event Assert::_execute(const std::vector<cl_event> events)
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

    vars->solve("int", _condition, data, "assert_result");

    // Check the result
    memcpy(&result, data, sizeof(int));
    free(data);
    if(result == 0){
        std::stringstream msg;
        msg << "Assertion error. The expression \"" <<
               std::string(_condition) << "\" is false" << std::endl;
        LOG(L_ERROR, msg.str());
        throw std::runtime_error("Assertion error");
    }

    return NULL;
}

}}  // namespaces
