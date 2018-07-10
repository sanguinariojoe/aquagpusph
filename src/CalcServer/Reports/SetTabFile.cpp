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
 * @brief Runtime output file.
 * (See Aqua::CalcServer::Reports::SetTabFile for details)
 */

#include <CalcServer/Reports/SetTabFile.h>
#include <InputOutput/Logger.h>
#include <CalcServer.h>
#include <Variable.h>

namespace Aqua{ namespace CalcServer{ namespace Reports{

SetTabFile::SetTabFile(const std::string tool_name,
                       const std::string fields,
                       unsigned int first,
                       unsigned int n,
                       const std::string output_file,
                       unsigned int ipf,
                       float fps)
    : Report(tool_name, fields, ipf, fps)
    , _output_file(output_file)
{
    _bounds.x = first;
    _bounds.y = first + n;
}

SetTabFile::~SetTabFile()
{
    if(_f.is_open()) _f.close();
}

void SetTabFile::setup()
{
    unsigned int i, j;

    std::ostringstream msg;
    msg << "Loading the report \"" << name() << "\"..." << std::endl;
    LOG(L_INFO, msg.str());

    _f.open(_output_file.c_str(), std::ios::out);

    Report::setup();

    // Write the header
    _f << "# Time ";
    std::vector<InputOutput::Variable*> vars = variables();
    for(i = _bounds.x; i < _bounds.y; i++){
        for(j = 0; j < vars.size(); j++){
            _f << vars.at(j)->name() << "_" << i;
        }
    }
    _f << std::endl;
    _f.flush();
}

void SetTabFile::_execute()
{
    if(!mustUpdate()){
        return;
    }

    unsigned int i, j;
    cl_int err_code;
    CalcServer *C = CalcServer::singleton();

    // Print the time instant
    _f << C->variables()->get("t")->asString() << " ";

    // Get the data to be printed
    std::vector<InputOutput::Variable*> vars = variables();
    for(i = 0; i < vars.size(); i++){
        if(vars.at(i)->type().find('*') == std::string::npos){
            std::stringstream msg;
            msg << "The report \"" << name()
                << "\" may not save scalar variables (\""
                << vars.at(i)->name() << "\")." << std::endl;
            LOG(L_ERROR, msg.str());
            throw std::runtime_error("Invalid variable type");
        }
        InputOutput::ArrayVariable *var = (InputOutput::ArrayVariable*)vars.at(i);
        size_t typesize = C->variables()->typeToBytes(var->type());
        size_t len = var->size() / typesize;
        if(len < bounds().y){
            std::stringstream msg;
            msg << "The report \"" << name()
                << "\" may not save field \""
                << vars.at(i)->name() << "\" because is not long enough."
                << std::endl;
            LOG(L_ERROR, msg.str());
            throw std::runtime_error("Invalid variable type");
        }
    }
    std::vector<void*> data = download(vars);

    // Print the data
    for(i = 0; i < bounds().y - bounds().x; i++){
        for(j = 0; j < vars.size(); j++){
            InputOutput::ArrayVariable *var = (InputOutput::ArrayVariable*)vars.at(j);
            const std::string type_name = var->type();
            if(!type_name.compare("int*")){
                int* v = (int*)data.at(j);
                _f << v[i] << " ";
            }
            else if(!type_name.compare("unsigned int*")){
                unsigned int* v = (unsigned int*)data.at(j);
                _f << v[i] << " ";
            }
            else if(!type_name.compare("float*")){
                float* v = (float*)data.at(j);
                _f << v[i] << " ";
            }
            else if(!type_name.compare("ivec*")){
                #ifdef HAVE_3D
                    ivec* v = (ivec*)data.at(j);
                    _f << "(" << v[i].x << ","
                              << v[i].y << ","
                              << v[i].z << ","
                              << v[i].w << ") ";
                #else
                    ivec* v = (ivec*)data.at(j);
                    _f << "(" << v[i].x << ","
                              << v[i].y << ") ";
                #endif // HAVE_3D
            }
            else if(!type_name.compare("ivec2*")){
                ivec2* v = (ivec2*)data.at(j);
                _f << "(" << v[i].x << ","
                          << v[i].y << ") ";
            }
            else if(!type_name.compare("ivec3*")){
                ivec3* v = (ivec3*)data.at(j);
                _f << "(" << v[i].x << ","
                          << v[i].y << ","
                          << v[i].z << ") ";
            }
            else if(!type_name.compare("ivec4*")){
                ivec4* v = (ivec4*)data.at(j);
                _f << "(" << v[i].x << ","
                          << v[i].y << ","
                          << v[i].z << ","
                          << v[i].w << ") ";
            }
            else if(!type_name.compare("uivec*")){
                #ifdef HAVE_3D
                    uivec* v = (uivec*)data.at(j);
                    _f << "(" << v[i].x << ","
                              << v[i].y << ","
                              << v[i].z << ","
                              << v[i].w << ") ";
                #else
                    uivec* v = (uivec*)data.at(j);
                    _f << "(" << v[i].x << ","
                              << v[i].y << ") ";
                #endif // HAVE_3D
            }
            else if(!type_name.compare("uivec2*")){
                uivec2* v = (uivec2*)data.at(j);
                _f << "(" << v[i].x << ","
                          << v[i].y << ") ";
            }
            else if(!type_name.compare("uivec3*")){
                uivec3* v = (uivec3*)data.at(j);
                _f << "(" << v[i].x << ","
                          << v[i].y << ","
                          << v[i].z << ") ";
            }
            else if(!type_name.compare("uivec4*")){
                uivec4* v = (uivec4*)data.at(j);
                _f << "(" << v[i].x << ","
                          << v[i].y << ","
                          << v[i].z << ","
                          << v[i].w << ") ";
            }
            else if(!type_name.compare("vec*")){
                #ifdef HAVE_3D
                    vec* v = (vec*)data.at(j);
                    _f << "(" << v[i].x << ","
                              << v[i].y << ","
                              << v[i].z << ","
                              << v[i].w << ") ";
                #else
                    vec* v = (vec*)data.at(j);
                    _f << "(" << v[i].x << ","
                              << v[i].y << ") ";
                #endif // HAVE_3D
            }
            else if(!type_name.compare("vec2*")){
                vec2* v = (vec2*)data.at(j);
                _f << "(" << v[i].x << ","
                          << v[i].y << ") ";
            }
            else if(!type_name.compare("vec3*")){
                vec3* v = (vec3*)data.at(j);
                _f << "(" << v[i].x << ","
                          << v[i].y << ","
                          << v[i].z << ") ";
            }
            else if(!type_name.compare("vec4*")){
                vec4* v = (vec4*)data.at(j);
                    _f << "(" << v[i].x << ","
                              << v[i].y << ","
                              << v[i].z << ","
                              << v[i].w << ") ";
            }
        }
    }
    _f << std::endl;
    _f.flush();

    for(i = 0; i < vars.size(); i++){
        free(data.at(i)); data.at(i) = NULL;
    }
    data.clear();
}

std::vector<void*> SetTabFile::download(std::vector<InputOutput::Variable*> vars)
{
    std::vector<void*> data;
    std::vector<cl_event> events;  // vector storage is continuous memory
    size_t typesize, len;
    unsigned int i, j;
    cl_int err_code;
    CalcServer *C = CalcServer::singleton();

    for(i = 0; i < vars.size(); i++){
        InputOutput::ArrayVariable *var = (InputOutput::ArrayVariable*)vars.at(i);
        typesize = C->variables()->typeToBytes(var->type());
        len = var->size() / typesize;
        if(len < bounds().y){
            std::stringstream msg;
            msg << "The report \"" << name()
                << "\" may not save field \""
                << vars.at(i)->name() << "\" because is not long enough."
                << std::endl;
            LOG(L_ERROR, msg.str());
            clearList(&data);
            throw std::runtime_error("Invalid variable type");
        }
        void *store = malloc(typesize * (bounds().y - bounds().x));
        if(!store){
            std::stringstream msg;
            msg << "Failure allocating " << typesize * (bounds().y - bounds().x)
                << " bytes for the field \"" << vars.at(i)->name()
                << "\"." << std::endl;
            clearList(&data);
            throw std::bad_alloc();
        }
        data.push_back(store);

        cl_event event;
        try {
            event = C->getUnsortedMem(var->name().c_str(),
                                      typesize * bounds().x,
                                      typesize * (bounds().y - bounds().x),
                                      store);
        } catch(...) {
            clearList(&data);
            throw;
        }

        events.push_back(event);
    }

    // Wait until all the data has been downloaded
    err_code = clWaitForEvents(events.size(),
                               events.data());
    if(err_code != CL_SUCCESS){
        LOG(L_ERROR, "Failure waiting for the variables download.\n");
        InputOutput::Logger::singleton()->printOpenCLError(err_code);
        clearList(&data);
        throw std::runtime_error("OpenCL error");
    }

    // Destroy the events
    for(i = 0; i < events.size(); i++){
        err_code = clReleaseEvent(events.at(i));
        if(err_code != CL_SUCCESS){
            LOG(L_ERROR, "Failure releasing the events.\n");
            InputOutput::Logger::singleton()->printOpenCLError(err_code);
            clearList(&data);
            throw std::runtime_error("OpenCL error");
        }
    }

    return data;
}

void SetTabFile::clearList(std::vector<void*> *data)
{
    unsigned int i;
    for(i = 0; i < data->size(); i++){
        if(data->at(i)) free(data->at(i)); data->at(i)=NULL;
    }
    data->clear();
}

}}} // namespace
