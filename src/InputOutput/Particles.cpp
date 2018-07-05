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
 * @brief Particles files manager.
 * (See Aqua::InputOutput::Particles for details)
 */

#include <stdlib.h>
#include <string.h>
#include <vector>

#include <InputOutput/Particles.h>
#include <ScreenManager.h>
#include <CalcServer.h>
#include <AuxiliarMethods.h>

namespace Aqua{ namespace InputOutput{

Particles::Particles(ProblemSetup sim_data,
                     unsigned int first,
                     unsigned int n,
                     unsigned int iset)
    : _sim_data(sim_data)
    , _iset(iset)
    , _output_file(NULL)
{
    _bounds.x = first;
    _bounds.y = first + n;
}

Particles::~Particles()
{
    if(_output_file)
        delete[] _output_file;
    _output_file = NULL;
}

bool Particles::loadDefault()
{
    unsigned int i;
    cl_int err_code;
    ArrayVariable * var;
    cl_mem mem;
    ScreenManager *S = ScreenManager::singleton();
    CalcServer::CalcServer *C = CalcServer::CalcServer::singleton();

    unsigned int n = bounds().y - bounds().x;
    unsigned int *iset = new unsigned int[n];
    unsigned int *id = new unsigned int[n];

    if(!iset || !id){
        S->addMessageF(L_ERROR, "Failure allocating memory.\n");
    }

    for(i = 0; i < n; i++){
        iset[i] = setId();
        id[i] = bounds().x + i;
    }

    Variables vars = C->variables();
    var = (ArrayVariable*)vars.get("iset");
    mem = *(cl_mem*)var->get();
    err_code = clEnqueueWriteBuffer(C->command_queue(),
                                    mem,
                                    CL_TRUE,
                                    sizeof(unsigned int) * bounds().x,
                                    sizeof(unsigned int) * n,
                                    iset,
                                    0,
                                    NULL,
                                    NULL);
    if(err_code != CL_SUCCESS){
        S->addMessageF(L_ERROR, "Failure sending variable \"iset\" to the server.\n");
        S->printOpenCLError(err_code);
    }
    var = (ArrayVariable*)vars.get("id");
    mem = *(cl_mem*)var->get();
    err_code = clEnqueueWriteBuffer(C->command_queue(),
                                    mem,
                                    CL_TRUE,
                                    sizeof(unsigned int) * bounds().x,
                                    sizeof(unsigned int) * n,
                                    id,
                                    0,
                                    NULL,
                                    NULL);
    if(err_code != CL_SUCCESS){
        S->addMessageF(L_ERROR, "Failure sending variable \"id\" to the server.\n");
        S->printOpenCLError(err_code);
    }
    var = (ArrayVariable*)vars.get("id_sorted");
    mem = *(cl_mem*)var->get();
    err_code = clEnqueueWriteBuffer(C->command_queue(),
                                    mem,
                                    CL_TRUE,
                                    sizeof(unsigned int) * bounds().x,
                                    sizeof(unsigned int) * n,
                                    id,
                                    0,
                                    NULL,
                                    NULL);
    if(err_code != CL_SUCCESS){
        S->addMessageF(L_ERROR, "Failure sending variable \"id_sorted\" to the server.\n");
        S->printOpenCLError(err_code);
    }
    var = (ArrayVariable*)vars.get("id_unsorted");
    mem = *(cl_mem*)var->get();
    err_code = clEnqueueWriteBuffer(C->command_queue(),
                                    mem,
                                    CL_TRUE,
                                    sizeof(unsigned int) * bounds().x,
                                    sizeof(unsigned int) * n,
                                    id,
                                    0,
                                    NULL,
                                    NULL);
    if(err_code != CL_SUCCESS){
        S->addMessageF(L_ERROR, "Failure sending variable \"id_unsorted\" to the server.\n");
        S->printOpenCLError(err_code);
    }

    delete[] iset; iset = NULL;
    delete[] id; id = NULL;
}

void Particles::file(const char* filename)
{
    size_t len;

    if(_output_file)
        delete[] _output_file;
    _output_file = NULL;

    if(!filename)
        return;

    len = strlen(filename) + 1;
    _output_file = new char[len];
    strcpy(_output_file, filename);
}

unsigned int Particles::file(const char* basename,
                             unsigned int startindex,
                             unsigned int digits)
{
    FILE *f;
    char *newname = NULL, *orig_pos, *dest_pos;
    size_t len;
    unsigned int i = startindex, j;

    if(!basename)
        return 0;

    if(!strstr(basename, "%d")){
        // We cannot replace nothing in the basename, just test if the file
        // does not exist
        f = fopen(basename, "r");
        if(f){
            // The fail already exist, so we cannot operate
            fclose(f);
            return 0;
        }

        file(basename);
        return 1;
    }

    while(true){
        if(newname)
            delete[] newname;
        newname = NULL;

        len = strlen(basename) - 1 + max(numberOfDigits(i), digits);
        newname = new char[len];

        // Copy all the string
        strcpy(newname, basename);
        // Replace the number
        dest_pos = strstr(newname, "%d");
        for(j = 0; j < digits - numberOfDigits(i); j++){
            strcpy(dest_pos, "0");
            dest_pos += 1;
        }
        sprintf(dest_pos, "%u", i);
        // Copy the rest of the original string after the inserted number
        dest_pos += numberOfDigits(i);
        orig_pos = (char*)strstr(basename, "%d") + 2;
        strcpy(dest_pos, orig_pos);

        f = fopen(newname, "r");
        if(!f){
            // We found an available slot
            break;
        }
        fclose(f);

        i++;
    }

    file(newname);
    delete[] newname;
    return i + 1;
}

std::deque<void*> Particles::download(std::vector<std::string> fields)
{
    std::deque<void*> data;
    std::vector<cl_event> events;  // vector storage is continuous memory
    size_t typesize, len;
    unsigned int i, j;
    cl_int err_code;
    char msg[256];
	ScreenManager *S = ScreenManager::singleton();
	CalcServer::CalcServer *C = CalcServer::CalcServer::singleton();
    Variables vars = C->variables();

    for(i = 0; i < fields.size(); i++){
        if(!vars.get(fields.at(i).c_str())){
            sprintf(msg,
                    "Undeclared \"%s\" field cannot be downloaded.\n",
                    fields.at(i).c_str());
            S->addMessageF(L_ERROR, msg);
            clearList(&data);
            return data;
        }
        if(vars.get(fields.at(i).c_str())->type().find('*') == std::string::npos){
            sprintf(msg,
                    "\"%s\" field is a scalar.\n",
                    fields.at(i).c_str());
            S->addMessageF(L_ERROR, msg);
            clearList(&data);
            return data;
        }
        ArrayVariable *var = (ArrayVariable*)vars.get(fields.at(i).c_str());
        typesize = vars.typeToBytes(var->type());
        len = var->size() / typesize;
        if(len < bounds().y){
            sprintf(msg,
                    "Failure saving \"%s\" field, which has not length enough.\n",
                    fields.at(i).c_str());
            S->addMessageF(L_ERROR, msg);
            sprintf(msg,
                    "length = %u was required, but just %lu was found.\n",
                    bounds().y,
                    len);
            S->addMessage(L_DEBUG, msg);
            clearList(&data);
            return data;
        }
        void *store = malloc(typesize * (bounds().y - bounds().x));
        if(!store){
            sprintf(msg,
                    "Failure allocating memory for \"%s\" field.\n",
                    fields.at(i).c_str());
            S->addMessageF(L_ERROR, msg);
            clearList(&data);
            return data;
        }
        data.push_back(store);

        cl_event event = C->getUnsortedMem(var->name().c_str(),
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
        S->addMessageF(L_ERROR, "Failure waiting for the variables download.\n");
        S->printOpenCLError(err_code);
        clearList(&data);
        return data;
    }

    // Destroy the events
    for(i = 0; i < events.size(); i++){
        err_code = clReleaseEvent(events.at(i));
        if(err_code != CL_SUCCESS){
            S->addMessageF(L_ERROR, "Failure releasing the events.\n");
            S->printOpenCLError(err_code);
            clearList(&data);
            return data;
        }
    }

    return data;
}

void Particles::clearList(std::deque<void*> *data)
{
    unsigned int i;
    for(i = 0; i < data->size(); i++){
        if(data->at(i)) free(data->at(i)); data->at(i)=NULL;
    }
    data->clear();
}

}}  // namespace
