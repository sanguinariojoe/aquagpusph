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
 * @brief Set a scalar variable.
 * (See Aqua::CalcServer::SetScalar for details)
 */

#include <InputOutput/Logger.h>
#include <CalcServer/SetScalar.h>
#include <CalcServer.h>

namespace Aqua{ namespace CalcServer{

SetScalar::SetScalar(const std::string name,
                     const std::string var_name,
                     const std::string value,
                     bool once)
    : Tool(name, once)
    , _var_name(var_name)
    , _value(value)
    , _var(NULL)
{
}

SetScalar::~SetScalar()
{
}

void SetScalar::setup()
{
    std::ostringstream msg;
    msg << "Loading the tool \"" << name() << "\"..." << std::endl;
    LOG(L_INFO, msg.str());

    Tool::setup();
    variable();
}


cl_event SetScalar::_execute(const std::vector<cl_event> events)
{
    InputOutput::Variables *vars = CalcServer::singleton()->variables();

    void *data = malloc(_var->typesize());
    if(!data){
        std::stringstream msg;
        msg << "Failure allocating " << _var->typesize()
            << " bytes for the variable \""
            << _var->name() << "\"." << std::endl;
        LOG(L_ERROR, msg.str());
        throw std::bad_alloc();
    }

    try {
        vars->solve(_var->type(), _value, data, _var->name());
    } catch(...) {
        free(data);
        throw;
    }

    _var->set(data);
    free(data);
    // Ensure that the variable is populated
    vars->populate(_var);

    return NULL;
}

void SetScalar::variable()
{
    CalcServer *C = CalcServer::singleton();
    InputOutput::Variables *vars = C->variables();
    if(!vars->get(_var_name)){
        std::stringstream msg;
        msg << "The tool \"" << name()
            << "\" is asking the undeclared variable \""
            << _var_name << "\"." << std::endl;
        LOG(L_ERROR, msg.str());
        throw std::runtime_error("Invalid variable");
    }
    if(vars->get(_var_name)->type().find('*') != std::string::npos){
        std::stringstream msg;
        msg << "The tool \"" << name()
            << "\" is asking the variable \"" << _var_name
            << "\", which is an array." << std::endl;
        LOG(L_ERROR, msg.str());
        throw std::runtime_error("Invalid variable type");
    }
    _var = vars->get(_var_name);
}

}}  // namespaces
