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

#include <stdlib.h>
#include <math.h>
#include <vector>

#include <AuxiliarMethods.h>
#include <ProblemSetup.h>
#include <ScreenManager.h>
#include <CalcServer/SetScalar.h>
#include <CalcServer.h>

namespace Aqua{ namespace CalcServer{

SetScalar::SetScalar(const char *name, const char *var_name, const char *value)
    : Tool(name)
    , _var_name(NULL)
    , _value(NULL)
    , _var(NULL)
{
    _var_name = new char[strlen(var_name) + 1];
    strcpy(_var_name, var_name);
    _value = new char[strlen(value) + 1];
    strcpy(_value, value);
}

SetScalar::~SetScalar()
{
    if(_var_name) delete[] _var_name; _var_name=NULL;
    if(_value) delete[] _value; _value=NULL;
}

bool SetScalar::setup()
{
    char msg[1024];
    InputOutput::ScreenManager *S = InputOutput::ScreenManager::singleton();

    sprintf(msg,
            "Loading the tool \"%s\"...\n",
            name());
    S->addMessageF(L_INFO, msg);

    if(variable()){
        return true;
    }

    return false;
}


bool SetScalar::_execute()
{
    char msg[1024];
    InputOutput::ScreenManager *S = InputOutput::ScreenManager::singleton();
    CalcServer *C = CalcServer::singleton();
    InputOutput::Variables *vars = C->variables();

    void *data = malloc(_var->typesize());
    if(!data){
        sprintf(msg,
                "Failure allocating \"%lu\" bytes for variable \"%s\".\n",
                _var->typesize(),
                _var->name());
        S->addMessageF(L_ERROR, msg);
        free(data);
        return true;
    }

    try {
        vars->solve(_var->type(), _value, data, _var->name());
    } catch(...) {
        free(data);
        return true;
    }

    _var->set(data);
    free(data);
    // Ensure that the variable is populated
    try {
        vars->populate(_var);
    } catch (...) {
        return true;
    }

    return false;
}

bool SetScalar::variable()
{
    char msg[1024];
    InputOutput::ScreenManager *S = InputOutput::ScreenManager::singleton();
    CalcServer *C = CalcServer::singleton();
    InputOutput::Variables *vars = C->variables();
    if(!vars->get(_var_name)){
        sprintf(msg,
                "The tool \"%s\" is using the undeclared variable \"%s\".\n",
                name(),
                _var_name);
        S->addMessageF(L_ERROR, msg);
        return true;
    }
    if(vars->get(_var_name)->type().find('*') != std::string::npos){
        sprintf(msg,
                "The tool \"%s\" has received the array variable \"%s\".\n",
                name(),
                _var_name);
        S->addMessageF(L_ERROR, msg);
        return true;
    }
    _var = vars->get(_var_name);
    return false;
}

}}  // namespaces
