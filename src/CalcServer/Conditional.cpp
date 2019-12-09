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
 * (See Aqua::CalcServer::Conditional and Aqua::CalcServer::End for details)
 */

#include <math.h>

#include <AuxiliarMethods.h>
#include <InputOutput/Logger.h>
#include <CalcServer.h>
#include <CalcServer/Conditional.h>

namespace Aqua{ namespace CalcServer{

Conditional::Conditional(const std::string name, const std::string condition, bool once)
    : Tool(name, once)
    , _condition(condition)
    , _ending_tool(NULL)
    , _result(true)
{
}

Conditional::~Conditional()
{
}

void Conditional::setup()
{
    std::vector<Tool*> tools = CalcServer::singleton()->tools();

    // Get the next tool in case the condition is not fulfilled
    int i = id_in_pipeline();
    if(i < 0){
        std::ostringstream msg;
        msg << "Tool \"" << name() << "\" cannot be found in the pipeline"
            << std::endl;
        LOG(L_ERROR, msg.str());
        throw std::runtime_error("Invalid tool");        
    }

    int scope = 1;
    while(i < tools.size() - 1) {
        i++;
        scope += tools.at(i)->scope_modifier();
        if(!scope)
            break;
    }
    if(scope > 0) {
        std::ostringstream msg;
        msg << "Unclosed scope for Tool \"" << name()
            << "\". Add and End tool" << std::endl;
        LOG(L_ERROR, msg.str());
        throw std::runtime_error("Unbalanced scope");
    }
    // We cannot use next_tool() attribute since that tool has not been setup
    // yet
    if(i == tools.size() - 1)
        _ending_tool = NULL;
    else
        _ending_tool = tools.at(i + 1);

    // Get the next tool in case the condition is fulfilled
    Tool::setup();
}

Tool* Conditional::next_tool()
{
    if(_result)
        return Tool::next_tool();
    return _ending_tool;
}

cl_event Conditional::_execute(const std::vector<cl_event> events)
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

While::While(const std::string name, const std::string condition, bool once)
    : Conditional(name, condition, once)
{
}

While::~While()
{
}

void While::setup()
{
    std::ostringstream msg;
    msg << "Loading the tool \"" << name() << "\"..." << std::endl;
    LOG(L_INFO, msg.str());
    Conditional::setup();
}

If::If(const std::string name, const std::string condition, bool once)
    : Conditional(name, condition, once)
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
    Conditional::setup();
}

Tool* If::next_tool()
{
    Tool* next_tool = Conditional::next_tool();
    _result = false;
    return next_tool;
}

End::End(const std::string name, bool once)
    : Tool(name, once)
{
}

End::~End()
{
}

void End::setup()
{
    std::ostringstream msg;
    msg << "Loading the tool \"" << name() << "\"..." << std::endl;
    LOG(L_INFO, msg.str());

    std::vector<Tool*> tools = CalcServer::singleton()->tools();

    // Locate the opening scope tool (which will be always the next tool on the
    // pipeline)
    int i = id_in_pipeline();
    if(i < 0){
        std::ostringstream msg;
        msg << "Tool \"" << name() << "\" cannot be found in the pipeline"
            << std::endl;
        LOG(L_ERROR, msg.str());
        throw std::runtime_error("Invalid tool");        
    }

    int scope = 1;
    while(i > 0) {
        i--;
        scope -= tools.at(i)->scope_modifier();
        if(!scope)
            break;
    }
    if(scope > 0) {
        std::ostringstream msg;
        msg << "Tool \"" << name()
            << "\" cannot be associated to any scope opening tool (If/While)"
            << std::endl;
        LOG(L_ERROR, msg.str());
        throw std::runtime_error("Unbalanced scope");
    }
    next_tool(tools.at(i));
}

}}  // namespaces
