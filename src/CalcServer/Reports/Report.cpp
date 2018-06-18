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
 * @brief Runtime output base class.
 * (See Aqua::CalcServer::Reports::Report for details)
 */

#include <CalcServer/Reports/Report.h>
#include <CalcServer.h>
#include <ScreenManager.h>

namespace Aqua{ namespace CalcServer{ namespace Reports{

Report::Report(const char* tool_name,
               const char* fields,
               unsigned int ipf,
               float fps)
    : Tool(tool_name)
    , _fields(NULL)
    , _ipf(ipf)
    , _fps(fps)
    , _iter(0)
    , _t(0.f)
    , _data(NULL)
{
    _fields = new char[strlen(fields) + 1];
    strcpy(_fields, fields);
}

Report::~Report()
{
    if(_fields) delete[] _fields; _fields=NULL;
    _vars_per_line.clear();
    _vars.clear();
    if(_data) delete[] _data; _data = NULL;
}

bool Report::setup()
{
    if(processFields(_fields)){
        return true;
    }

    return false;
}

const char* Report::data(bool with_title, bool with_names)
{
    unsigned int i, j, var_id=0;

    if(_data) delete[] _data; _data = NULL;
    _data = new char[dataLength(with_title)];
    strcpy(_data, "");

    // Create the title
    if(with_title){
        sprintf(_data, "%s:\n", name());
    }

    // Set the variable per lines
    for(i = 0; i < _vars_per_line.size(); i++){
        for(j = 0; j < _vars_per_line.at(i); j++){
            InputOutput::Variable* var = _vars.at(var_id);
            if(with_names){
                strcat(_data, var->name());
                strcat(_data, "=");
            }
            strcat(_data, var->asString());
            strcat(_data, " ");
            var_id++;
        }
        // Replace the trailing space by a line break
        _data[strlen(_data) - 1] = '\n';
    }

    return (const char*)_data;
}

bool Report::processFields(const char* input)
{
    unsigned int i;
    CalcServer *C = CalcServer::singleton();
    InputOutput::ScreenManager *S = InputOutput::ScreenManager::singleton();
    char fields[strlen(input) + 1];
    strcpy(fields, input);

    // Check if line breaks have been requested
    if(strchr(fields, ';')){
        // Now we need to store all the tokens before recursively calling this
        // function, otherwise the strtok(fields, " ,") will broke this one
        char *tok = strtok(fields, ";");
        std::deque<char *> tokens;
        while(tok){
            char *line = new char[strlen(tok) + 1];
            strcpy(line, tok);
            tokens.push_back(line);
            tok = strtok(NULL, ";");
        }
        // And now we can parse each line
        for(i = 0; i < tokens.size(); i++){
            tok = tokens.at(i);
            if(processFields(tok)){
                return true;
            }
            delete[] tok;
        }
        tokens.clear();
        return false;
    }

    // Now we now taht the fields should be process as a single line
    InputOutput::Variables* vars = C->variables();
    unsigned int vars_in_line = 0;
    char *field = strtok(fields, " ,");
    while(field){
        InputOutput::Variable *var = vars->get(field);
        if(!var){
            char msg[strlen(field) + 64];
            sprintf(msg,
                    "\"%s\" variable cannot be found.\n",
                    field);
            S->addMessageF(L_ERROR, msg);
            return true;
        }
        vars_in_line++;
        _vars.push_back(var);
        field = strtok(NULL, " ,");
    }
    _vars_per_line.push_back(vars_in_line);

    return false;
}

size_t Report::dataLength(bool with_title, bool with_names)
{
    unsigned int i;

    size_t len = 1;  // end of string
    // Title
    if(with_title){
        len += strlen(name()) + 2;
    }
    // Variables
    for(i = 0; i < _vars.size(); i++){
        InputOutput::Variable *var = _vars.at(i);
        if(with_names){
            len += strlen(var->name()) + 1;
        }
        len += strlen(var->asString()) + 1;
    }

    return len;
}

bool Report::mustUpdate()
{
    CalcServer *C = CalcServer::singleton();
    InputOutput::Variables* vars = C->variables();

    InputOutput::UIntVariable *iter_var =
        (InputOutput::UIntVariable*)vars->get("iter");
    unsigned int iter = *(unsigned int*)iter_var->get();
    InputOutput::FloatVariable *time_var =
        (InputOutput::FloatVariable*)vars->get("t");
    float t = *(float*)time_var->get();

    if(_ipf > 0){
        if(iter - _iter >= _ipf){
            _iter = iter;
            _t = t;
            return true;
        }
    }
    if(_fps > 0.f){
        if(t - _t >= 1.f / _fps){
            _iter = iter;
            _t = t;
            return true;
        }
    }
    return false;
}

}}} // namespace
