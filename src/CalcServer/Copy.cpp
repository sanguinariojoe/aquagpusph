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

#include <stdlib.h>
#include <math.h>
#include <vector>

#include <ProblemSetup.h>
#include <ScreenManager.h>
#include <CalcServer/Copy.h>
#include <CalcServer.h>

namespace Aqua{ namespace CalcServer{

Copy::Copy(const char *name, const char *input_name, const char *output_name)
    : Tool(name)
    , _input_name(NULL)
    , _output_name(NULL)
    , _input_var(NULL)
    , _output_var(NULL)
{
    _input_name = new char[strlen(input_name) + 1];
    strcpy(_input_name, input_name);
    _output_name = new char[strlen(output_name) + 1];
    strcpy(_output_name, output_name);
}

Copy::~Copy()
{
    if(_input_name) delete[] _input_name; _input_name=NULL;
    if(_output_name) delete[] _output_name; _output_name=NULL;
}

bool Copy::setup()
{
    char msg[1024];
    InputOutput::ScreenManager *S = InputOutput::ScreenManager::singleton();
    CalcServer *C = CalcServer::singleton();
    InputOutput::Variables *vars = C->variables();

    sprintf(msg,
            "Loading the tool \"%s\"...\n",
            name());
    S->addMessageF(1, msg);

    if(variables()){
        return true;
    }

    return false;
}


bool Copy::execute()
{
    unsigned int i;
    cl_int err_code;
    char msg[1024];
    InputOutput::ScreenManager *S = InputOutput::ScreenManager::singleton();
    CalcServer *C = CalcServer::singleton();

    err_code = clEnqueueCopyBuffer(C->command_queue(),
                                   *(cl_mem*)_input_var->get(),
                                   *(cl_mem*)_output_var->get(),
                                   0,
                                   0,
                                   _output_var->size(),
                                   0,
                                   NULL,
                                   NULL);
    if(err_code != CL_SUCCESS){
        sprintf(msg,
                "Failure during the execution of the tool \"%s\".\n",
                name());
        S->addMessageF(3, msg);
        S->printOpenCLError(err_code);
        return true;
    }

    return false;
}

bool Copy::variables()
{
    char msg[1024];
    InputOutput::ScreenManager *S = InputOutput::ScreenManager::singleton();
    CalcServer *C = CalcServer::singleton();
    InputOutput::Variables *vars = C->variables();
    if(!vars->get(_input_name)){
        sprintf(msg,
                "The tool \"%s\" has received undeclared variable \"%s\" as input.\n",
                name(),
                _input_name);
        S->addMessageF(3, msg);
        return true;
    }
    if(!strchr(vars->get(_input_name)->type(), '*')){
        sprintf(msg,
                "The tool \"%s\" has received the scalar variable \"%s\" as input.\n",
                name(),
                _input_name);
        S->addMessageF(3, msg);
        return true;
    }
    _input_var = (InputOutput::ArrayVariable *)vars->get(_input_name);
    size_t n_in = _input_var->size() / vars->typeToBytes(_input_var->type());
    if(!vars->get(_output_name)){
        sprintf(msg,
                "The tool \"%s\" has received undeclared variable \"%s\" as output.\n",
                name(),
                _output_name);
        S->addMessageF(3, msg);
        return true;
    }
    if(!strchr(vars->get(_output_name)->type(), '*')){
        sprintf(msg,
                "The tool \"%s\" has received the scalar variable \"%s\" as output.\n",
                name(),
                _output_name);
        S->addMessageF(3, msg);
        return true;
    }
    _output_var = (InputOutput::ArrayVariable *)vars->get(_output_name);
    size_t n_out = _input_var->size() / vars->typeToBytes(_input_var->type());
    if(!vars->isSameType(_input_var->type(), _output_var->type())){
        sprintf(msg,
                "The input and output types mismatch for the tool \"%s\".\n",
                name());
        S->addMessageF(3, msg);
        sprintf(msg,
                "\tInput variable \"%s\" is of type \"%s\".\n",
                _input_var->name(),
                _input_var->type());
        S->addMessageF(0, msg);
        sprintf(msg,
                "\tOutput variable \"%s\" is of type \"%s\".\n",
                _output_var->name(),
                _output_var->type());
        S->addMessageF(0, msg);
        return true;
    }
    if(n_in != n_out){
        sprintf(msg,
                "The input and output lengths mismatch for the tool \"%s\".\n",
                name());
        S->addMessageF(3, msg);
        sprintf(msg,
                "\tInput variable \"%s\" has a length n=%lu.\n",
                _input_var->name(),
                n_in);
        S->addMessageF(0, msg);
        sprintf(msg,
                "\tOutput variable \"%s\" has a length n=%lu.\n",
                _output_var->name(),
                n_out);
        S->addMessageF(0, msg);
        return true;
    }

    return false;
}

}}  // namespaces
