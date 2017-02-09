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
 * (See Aqua::CalcServer::Assert for details)
 */

#include <stdlib.h>
#include <math.h>
#include <vector>

#include <AuxiliarMethods.h>
#include <ProblemSetup.h>
#include <ScreenManager.h>
#include <CalcServer/Assert.h>
#include <CalcServer.h>

namespace Aqua{ namespace CalcServer{

Assert::Assert(const char *name, const char *condition)
    : Tool(name)
    , _condition(NULL)
{
    _condition = new char[strlen(condition) + 1];
    strcpy(_condition, condition);
}

Assert::~Assert()
{
    if(_condition) delete[] _condition; _condition=NULL;
}

bool Assert::setup()
{
    char msg[1024];
    InputOutput::ScreenManager *S = InputOutput::ScreenManager::singleton();

    sprintf(msg,
            "Loading the tool \"%s\"...\n",
            name());
    S->addMessageF(1, msg);

    return false;
}


bool Assert::_execute()
{
    std::stringstream msg;
    int result;
    InputOutput::ScreenManager *S = InputOutput::ScreenManager::singleton();
    CalcServer *C = CalcServer::singleton();
    InputOutput::Variables *vars = C->variables();

    void *data = malloc(sizeof(int));
    if(!data){
        msg << "Failure allocating memory for the integer result" << std::endl;
        S->addMessageF(3, msg.str().c_str());
        return true;
    }

    if(vars->solve("int", _condition, data, "assert_result")){
        return true;
    }

    // Check the result
    memcpy(&result, data, sizeof(int));
    free(data);
    if(result == 0){
        msg << "Assertion error. The expression \"" <<
               std::string(_condition) << "\" is false" << std::endl;
        S->addMessageF(3, msg.str().c_str());
        return true;
    }

    return false;
}

}}  // namespaces
