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
#include <ScreenManager.h>
#include <CalcServer.h>
#include <Variable.h>

namespace Aqua{ namespace CalcServer{ namespace Reports{

SetTabFile::SetTabFile(const char* tool_name,
                       const char* fields,
                       unsigned int first,
                       unsigned int n,
                       const char* output_file)
    : Report(tool_name, fields)
    , _output_file(NULL)
    , _f(NULL)
{
    _bounds.x = first;
    _bounds.y = first + n;

    _output_file = new char[strlen(output_file) + 1];
    strcpy(_output_file, output_file);
}

SetTabFile::~SetTabFile()
{
    if(_output_file) delete[] _output_file; _output_file=NULL;
    if(_f) fclose(_f); _f = NULL;
}

bool SetTabFile::setup()
{
    unsigned int i, j;
    char msg[1024];
    InputOutput::ScreenManager *S = InputOutput::ScreenManager::singleton();

    sprintf(msg,
            "Loading the report \"%s\"...\n",
            name());
    S->addMessageF(1, msg);

    // Open the output file
    _f = fopen(_output_file, "w");
    if(!_f){
        sprintf(msg,
                "The file \"%s\" cannot be written\n",
                _output_file);
        S->addMessageF(3, msg);
        return true;
    }

    if(Report::setup()){
        return true;
    }

    // Write the header
    fprintf(_f, "# Time ");
    std::deque<InputOutput::Variable*> vars = variables();
    for(i = _bounds.x; i < _bounds.y; i++){
        for(j = 0; j < vars.size(); j++){
            fprintf(_f, "%s_%u ", vars.at(j)->name(), i);
        }
    }
    fprintf(_f, "\n");
    fflush(_f);

    return false;
}

bool SetTabFile::execute()
{
    unsigned int i, j;
    cl_int err_code;
    char msg[256];
    InputOutput::ScreenManager *S = InputOutput::ScreenManager::singleton();
	CalcServer *C = CalcServer::singleton();

    // Print the time instant
    fprintf(_f, "%s ", C->variables()->get("t")->asString());

    // Get the data to be printed
    std::deque<InputOutput::Variable*> vars = variables();
    for(i = 0; i < vars.size(); i++){
        if(!strchr(vars.at(i)->type(), '*')){
            sprintf(msg,
                    "\"%s\" field has been set to be saved, but it was declared as an scalar.\n",
                    vars.at(i)->name());
            S->addMessageF(3, msg);
            return true;
        }
        InputOutput::ArrayVariable *var = (InputOutput::ArrayVariable*)vars.at(i);
        size_t typesize = C->variables()->typeToBytes(var->type());
        size_t len = var->size() / typesize;
        if(len < bounds().y){
            sprintf(msg,
                    "Failure saving \"%s\" field, which has not length enough.\n",
                    vars.at(i)->name());
            S->addMessageF(3, msg);
            return true;
        }
    }
    std::deque<void*> data = download(vars);
    if(!data.size()){
        return true;
    }

    // Print the data
    for(i = 0; i < bounds().y - bounds().x; i++){
        for(j = 0; j < vars.size(); j++){
            InputOutput::ArrayVariable *var = (InputOutput::ArrayVariable*)vars.at(j);
            const char* type_name = var->type();
            if(!strcmp(type_name, "int*")){
                int* v = (int*)data.at(j);
                fprintf(_f, "%16d ", v[i]);
            }
            else if(!strcmp(type_name, "unsigned int*")){
                unsigned int* v = (unsigned int*)data.at(j);
                fprintf(_f, "%16u ", v[i]);
            }
            else if(!strcmp(type_name, "float*")){
                float* v = (float*)data.at(j);
                fprintf(_f, "%16g ", v[i]);
            }
            else if(!strcmp(type_name, "ivec*")){
                #ifdef HAVE_3D
                    ivec* v = (ivec*)data.at(j);
                    fprintf(_f, "(%16d,%16d,%16d,%16d) ",
                            v[i].x, v[i].y, v[i].z, v[i].w);
                #else
                    ivec* v = (ivec*)data.at(j);
                    fprintf(_f, "(%16d,%16d) ", v[i].x, v[i].y);
                #endif // HAVE_3D
            }
            else if(!strcmp(type_name, "ivec2*")){
                ivec2* v = (ivec2*)data.at(j);
                fprintf(_f, "(%16d,%16d) ", v[i].x, v[i].y);
            }
            else if(!strcmp(type_name, "ivec3*")){
                ivec3* v = (ivec3*)data.at(j);
                fprintf(_f, "(%16d,%16d,%16d) ", v[i].x, v[i].y, v[i].z);
            }
            else if(!strcmp(type_name, "ivec4*")){
                ivec4* v = (ivec4*)data.at(j);
                fprintf(_f, "(%16d,%16d,%16d,%16d) ",
                        v[i].x, v[i].y, v[i].z, v[i].w);
            }
            else if(!strcmp(type_name, "uivec*")){
                #ifdef HAVE_3D
                    uivec* v = (uivec*)data.at(j);
                    fprintf(_f, "(%16u,%16u,%16u,%16u) ",
                            v[i].x, v[i].y, v[i].z, v[i].w);
                #else
                    uivec* v = (uivec*)data.at(j);
                    fprintf(_f, "(%16u,%16u) ", v[i].x, v[i].y);
                #endif // HAVE_3D
            }
            else if(!strcmp(type_name, "uivec2*")){
                uivec2* v = (uivec2*)data.at(j);
                fprintf(_f, "(%16u,%16u) ", v[i].x, v[i].y);
            }
            else if(!strcmp(type_name, "uivec3*")){
                uivec3* v = (uivec3*)data.at(j);
                fprintf(_f, "(%16u,%16u,%16u) ", v[i].x, v[i].y, v[i].z);
            }
            else if(!strcmp(type_name, "uivec4*")){
                uivec4* v = (uivec4*)data.at(j);
                fprintf(_f, "(%16u,%16u,%16u,%16u) ",
                        v[i].x, v[i].y, v[i].z, v[i].w);
            }
            else if(!strcmp(type_name, "vec*")){
                #ifdef HAVE_3D
                    vec* v = (vec*)data.at(j);
                    fprintf(_f, "(%16g,%16g,%16g,%16g) ",
                            v[i].x, v[i].y, v[i].z, v[i].w);
                #else
                    vec* v = (vec*)data.at(j);
                    fprintf(_f, "(%16g,%16g) ", v[i].x, v[i].y);
                #endif // HAVE_3D
            }
            else if(!strcmp(type_name, "vec2*")){
                vec2* v = (vec2*)data.at(j);
                fprintf(_f, "(%16g,%16g) ", v[i].x, v[i].y);
            }
            else if(!strcmp(type_name, "vec3*")){
                vec3* v = (vec3*)data.at(j);
                fprintf(_f, "(%16u,%16u,%16u) ", v[i].x, v[i].y, v[i].z);
            }
            else if(!strcmp(type_name, "vec4*")){
                vec4* v = (vec4*)data.at(j);
                fprintf(_f, "(%16u,%16u,%16u,%16u) ",
                        v[i].x, v[i].y, v[i].z, v[i].w);
            }
        }
        fprintf(_f, "\n");
        fflush(_f);
    }

    for(i = 0; i < vars.size(); i++){
        free(data.at(i)); data.at(i) = NULL;
    }
    data.clear();

    fflush(_f);
    return false;
}

std::deque<void*> SetTabFile::download(std::deque<InputOutput::Variable*> vars)
{
    std::deque<void*> data;
    std::vector<cl_event> events;  // vector storage is continuous memory
    size_t typesize, len;
    unsigned int i, j;
    cl_int err_code;
    char msg[256];
	InputOutput::ScreenManager *S = InputOutput::ScreenManager::singleton();
	CalcServer *C = CalcServer::singleton();

    for(i = 0; i < vars.size(); i++){
        InputOutput::ArrayVariable *var = (InputOutput::ArrayVariable*)vars.at(i);
        typesize = C->variables()->typeToBytes(var->type());
        len = var->size() / typesize;
        if(len < bounds().y){
            sprintf(msg,
                    "Failure saving \"%s\" field, which has not length enough.\n",
                    vars.at(i)->name());
            S->addMessageF(3, msg);
            sprintf(msg,
                    "length = %u was required, but just %lu was found.\n",
                    bounds().y,
                    len);
            S->addMessage(0, msg);
            clearList(&data);
            return data;
        }
        void *store = malloc(typesize * (bounds().y - bounds().x));
        if(!store){
            sprintf(msg,
                    "Failure allocating memory for \"%s\" field.\n",
                    vars.at(i)->name());
            S->addMessageF(3, msg);
            clearList(&data);
            return data;
        }
        data.push_back(store);

        cl_event event = C->getUnsortedMem(var->name(),
                                           typesize * bounds().x,
                                           typesize * (bounds().y - bounds().x),
                                           store);
        if(!event){
            clearList(&data);
            return data;
        }

        events.push_back(event);
    }

    // Wait until all the data has been downloaded
    err_code = clWaitForEvents(events.size(),
                               events.data());
    if(err_code != CL_SUCCESS){
        S->addMessageF(3, "Failure waiting for the variables download.\n");
        S->printOpenCLError(err_code);
        clearList(&data);
        return data;
    }

    // Destroy the events
    for(i = 0; i < events.size(); i++){
        err_code = clReleaseEvent(events.at(i));
        if(err_code != CL_SUCCESS){
            S->addMessageF(3, "Failure releasing the events.\n");
            S->printOpenCLError(err_code);
            clearList(&data);
            return data;
        }
    }

    return data;
}

void SetTabFile::clearList(std::deque<void*> *data)
{
    unsigned int i;
    for(i = 0; i < data->size(); i++){
        if(data->at(i)) free(data->at(i)); data->at(i)=NULL;
    }
    data->clear();
}

}}} // namespace
