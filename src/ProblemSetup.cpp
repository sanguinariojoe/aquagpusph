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
 * @brief Simulation configuration data structures.
 * (See Aqua::InputOutput::ProblemSetup for details)
 */

#include <limits>

#include <ProblemSetup.h>
#include <ScreenManager.h>
#include <AuxiliarMethods.h>

namespace Aqua{ namespace InputOutput{

ProblemSetup::ProblemSetup()
{
    settings.init();

    time_opts.sim_end_mode = __NO_OUTPUT_MODE__;
    time_opts.sim_end_time = 0.f;
    time_opts.sim_end_step = 0;
    time_opts.sim_end_frame = 0;
    time_opts.log_mode = __NO_OUTPUT_MODE__;
    time_opts.log_fps = 0.f;
    time_opts.log_ipf = 0;
    time_opts.output_mode = __NO_OUTPUT_MODE__;
    time_opts.output_fps = 0.f;
    time_opts.output_ipf = 0;
}

ProblemSetup::~ProblemSetup()
{
    unsigned int i;
    settings.destroy();
    variables.destroy();
    definitions.destroy();
    for(i = 0; i < tools.size(); i++){
        delete tools.at(i);
    }
    tools.clear();
    for(i = 0; i < reports.size(); i++){
        delete reports.at(i);
    }
    reports.clear();
    for(i=0;i<sets.size();i++){
        delete sets.at(i);
    }
    sets.clear();
}

bool ProblemSetup::perform()
{
    unsigned int i;
    InputOutput::ScreenManager *S = InputOutput::ScreenManager::singleton();
    char msg[512];
    strcpy(msg, "");
    // Check for errors
    if(sets.size() == 0) {
        sprintf(msg, "No sets of particles have been added.\n");
        S->addMessageF(3, msg);
        return true;
    }
    return false;
}

void ProblemSetup::sphSettings::init()
{
    verbose_level = 1;
    platform_id = 0;
    device_id = 0;
    device_type = CL_DEVICE_TYPE_ALL;
}

void ProblemSetup::sphSettings::destroy()
{
}

void ProblemSetup::sphVariables::registerVariable(const char* name,
                                                  const char* type,
                                                  const char* length,
                                                  const char* value)
{
    size_t len;

    char *aux=NULL;
    len = strlen(name) + 1;
    aux = new char[len];
    strcpy(aux, name);
    names.push_back(aux);
    len = strlen(type) + 1;
    aux = new char[len];
    strcpy(aux, type);
    types.push_back(aux);
    len = strlen(length) + 1;
    aux = new char[len];
    strcpy(aux, length);
    lengths.push_back(aux);
    len = strlen(value) + 1;
    aux = new char[len];
    strcpy(aux, value);
    values.push_back(aux);
}

void ProblemSetup::sphVariables::destroy()
{
    unsigned int i;
    for(i=0;i<names.size();i++){
        delete[] names.at(i);
        delete[] types.at(i);
        delete[] lengths.at(i);
    }
    names.clear();
    types.clear();
    lengths.clear();
}

void ProblemSetup::sphDefinitions::registerDefinition(const char* name,
                                                      const char* value,
                                                      const bool evaluate)
{
    size_t len;

    char *aux=NULL;
    len = strlen(name) + 1;
    aux = new char[len];
    strcpy(aux, name);
    names.push_back(aux);
    len = strlen(value) + 1;
    aux = new char[len];
    strcpy(aux, value);
    values.push_back(aux);
    evaluations.push_back(evaluate);
}

bool ProblemSetup::sphDefinitions::isDefined(const char* name)
{
    unsigned int i;
    for(i = 0; i < names.size(); i++){
        if(!strcmp(name, names.at(i))){
            return true;
        }
    }

    return false;
}

void ProblemSetup::sphDefinitions::destroy()
{
    unsigned int i;
    for(i=0;i<names.size();i++){
        delete[] names.at(i);
        delete[] values.at(i);
    }
    names.clear();
    values.clear();
    evaluations.clear();
}

ProblemSetup::sphTool::sphTool()
{
}

ProblemSetup::sphTool::~sphTool()
{
    _data.clear();
}

void ProblemSetup::sphTool::set(const char* name,
                                const char* value)
{
    if(has(name)){
        _data[name] = value;
        return;
    }
    _data.insert(std::make_pair(name, value));
}

const char* ProblemSetup::sphTool::get(const char* name)
{
    if(!has(name)){
        return NULL;
    }
    return _data[name].c_str();
}

const char* ProblemSetup::sphTool::get(unsigned int index)
{
    unsigned int i = 0;
    for(std::map<std::string,std::string>::iterator it=_data.begin();
        it != _data.end();
        ++it){
        if(i == index){
            return it->second.c_str();
        }
        i++;
    }
    return NULL;
}

const char* ProblemSetup::sphTool::getName(unsigned int index)
{
    unsigned int i = 0;
    for(std::map<std::string,std::string>::iterator it=_data.begin();
        it != _data.end();
        ++it){
        if(i == index){
            return it->first.c_str();
        }
        i++;
    }
    return NULL;
}

bool ProblemSetup::sphTool::has(const char* name)
{
    std::map<std::string, std::string>::iterator var = _data.find(name);
    if(var != _data.end())
        return true;
    return false;
}

ProblemSetup::sphParticlesSet::sphParticlesSet()
    : _n(0)
    , _in_path(NULL)
    , _in_format(NULL)
    , _out_path(NULL)
    , _out_format(NULL)
{
}

ProblemSetup::sphParticlesSet::~sphParticlesSet()
{
    unsigned int i;
    for(i = 0; i < _snames.size(); i++){
        delete[] _snames.at(i); _snames.at(i) = NULL;
        delete[] _svalues.at(i); _svalues.at(i) = NULL;
    }
    _snames.clear();
    _svalues.clear();
    delete[] _in_path; _in_path=NULL;
    delete[] _in_format; _in_format=NULL;
    delete[] _out_path; _out_path=NULL;
    delete[] _out_format; _out_format=NULL;
    for(i = 0; i < _in_fields.size(); i++){
        delete[] _in_fields.at(i); _in_fields.at(i) = NULL;
    }
    _in_fields.clear();
    for(i = 0; i < _out_fields.size(); i++){
        delete[] _out_fields.at(i); _out_fields.at(i) = NULL;
    }
    _out_fields.clear();
}

void ProblemSetup::sphParticlesSet::addScalar(const char* name,
                                              const char* value)
{
    char* sname = new char[strlen(name) + 1];
    strcpy(sname, name);
    char* svalue = new char[strlen(value) + 1];
    strcpy(svalue, value);
    _snames.push_back(sname);
    _svalues.push_back(svalue);
}

void ProblemSetup::sphParticlesSet::input(const char* path,
                                          const char* format,
                                          const char* fields)
{
    unsigned int i;
    if(_in_path) delete[] _in_path;
    if(_in_format) delete[] _in_format;
    for(i = 0; i < _in_fields.size(); i++){
        delete[] _in_fields.at(i); _in_fields.at(i) = NULL;
    }
    _in_fields.clear();

    _in_path = new char[strlen(path) + 1];
    strcpy(_in_path, path);
    _in_format = new char[strlen(format) + 1];
    strcpy(_in_format, format);

    char auxfields[strlen(fields) + 1];
    strcpy(auxfields, fields);
    char *field = strtok(auxfields, " ,");
    while(field){
        char *aux = new char[strlen(field) + 1];
        strcpy(aux, field);
        _in_fields.push_back(aux);
        field = strtok (NULL, " ,");
    }
}

void ProblemSetup::sphParticlesSet::output(const char* path,
                                           const char* format,
                                           const char* fields)
{
    unsigned int i;
    if(_out_path) delete[] _out_path;
    if(_out_format) delete[] _out_format;
    for(i = 0; i < _out_fields.size(); i++){
        delete[] _out_fields.at(i); _out_fields.at(i) = NULL;
    }
    _out_fields.clear();

    _out_path = new char[strlen(path) + 1];
    strcpy(_out_path, path);
    _out_format = new char[strlen(format) + 1];
    strcpy(_out_format, format);

    char auxfields[strlen(fields) + 1];
    strcpy(auxfields, fields);
    char *field = strtok(auxfields, " ,");
    while(field){
        char *aux = new char[strlen(field) + 1];
        strcpy(aux, field);
        _out_fields.push_back(aux);
        field = strtok (NULL, " ,");
    }
}

}}  // namespace
