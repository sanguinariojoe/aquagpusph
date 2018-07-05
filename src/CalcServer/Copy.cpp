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
 * @brief Array copy tool.
 * (see Aqua::CalcServer::Copy for details)
 */

#include <stdlib.h>
#include <math.h>
#include <vector>

#include <ProblemSetup.h>
#include <ScreenManager.h>
#include <CalcServer/Copy.h>
#include <CalcServer.h>

namespace Aqua{ namespace CalcServer{

Copy::Copy(const std::string name,
           const std::string input_name,
           const std::string output_name)
    : Tool(name)
    , _input_name(input_name)
    , _output_name(output_name)
    , _input_var(NULL)
    , _output_var(NULL)
{
}

Copy::~Copy()
{
}

void Copy::setup()
{
    std::ostringstream msg;
    msg << "Loading the tool \"" << name() << "\"..." << std::endl;
    LOG(L_INFO, msg.str());

    variables();
}


void Copy::_execute()
{
    unsigned int i;
    cl_int err_code;
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
        std::stringstream msg;
        msg << "Failure executing the tool \"" <<
               name() << "\"." << std::endl;
        LOG(L_ERROR, msg.str());
        InputOutput::ScreenManager::singleton()->printOpenCLError(err_code);
        throw std::runtime_error("OpenCL execution error");
    }
}

void Copy::variables()
{
    CalcServer *C = CalcServer::singleton();
    InputOutput::Variables vars = C->variables();
    if(!vars.get(_input_name)){
        std::stringstream msg;
        msg << "The tool \"" << name()
            << "\" is asking the undeclared variable \""
            << _input_name << "\"." << std::endl;
        LOG(L_ERROR, msg.str());
        throw std::runtime_error("Invalid variable");
    }
    if(vars.get(_input_name)->type().find('*') == std::string::npos){
        std::stringstream msg;
        msg << "The tool \"" << name()
            << "\" may not use a scalar variable (\""
            << _input_name << "\")." << std::endl;
        LOG(L_ERROR, msg.str());
        throw std::runtime_error("Invalid variable type");
    }
    _input_var = (InputOutput::ArrayVariable *)vars.get(_input_name);
    size_t n_in = _input_var->size() / vars.typeToBytes(_input_var->type());
    if(!vars.get(_output_name)){
        std::stringstream msg;
        msg << "The tool \"" << name()
            << "\" is asking the undeclared variable \""
            << _output_name << "\"." << std::endl;
        LOG(L_ERROR, msg.str());
        throw std::runtime_error("Invalid variable");
    }
    if(vars.get(_output_name)->type().find('*') == std::string::npos){
        std::stringstream msg;
        msg << "The tool \"" << name()
            << "\" may not use a scalar variable (\""
            << _output_name << "\")." << std::endl;
        LOG(L_ERROR, msg.str());
        throw std::runtime_error("Invalid variable type");
    }
    _output_var = (InputOutput::ArrayVariable *)vars.get(_output_name);
    size_t n_out = _input_var->size() / vars.typeToBytes(_input_var->type());
    if(!vars.isSameType(_input_var->type(), _output_var->type())){
        std::stringstream msg;
        msg << "The input and output types mismatch for the tool \""
            << name() << "\"." << std::endl;
        LOG(L_ERROR, msg.str());
        msg.str("");
        msg << "\tInput variable \"" << _input_var->name()
            << "\" is of type \"" << _input_var->type() << "\"" << std::endl;
        LOG0(L_DEBUG, msg.str());
        msg << "\tOutput variable \"" << _output_var->name()
            << "\" is of type \"" << _output_var->type() << "\"" << std::endl;
        LOG0(L_DEBUG, msg.str());
        throw std::runtime_error("Incompatible types");
    }
    if(n_in != n_out){
        std::stringstream msg;
        msg << "Input and output lengths mismatch for the tool \""
            << name() << "\"." << std::endl;
        LOG(L_ERROR, msg.str());
        msg.str("");
        msg << "\tInput variable \"" << _input_var->name()
            << "\" has length " << n_in << std::endl;
        LOG0(L_DEBUG, msg.str());
        msg << "\tOutput variable \"" << _output_var->name()
            << "\" has length " << n_out << std::endl;
        LOG0(L_DEBUG, msg.str());
        throw std::runtime_error("Incompatible lenghts");
    }
}

}}  // namespaces
