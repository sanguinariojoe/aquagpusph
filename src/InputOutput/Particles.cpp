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

#include <string>

#include <InputOutput/Particles.h>
#include <InputOutput/Logger.h>
#include <CalcServer.h>
#include <AuxiliarMethods.h>

namespace Aqua{ namespace InputOutput{

Particles::Particles(ProblemSetup& sim_data,
                     unsigned int first,
                     unsigned int n,
                     unsigned int iset)
    : _sim_data(sim_data)
    , _iset(iset)
{
    _bounds.x = first;
    _bounds.y = first + n;
}

Particles::~Particles()
{
}

void Particles::loadDefault()
{
    unsigned int i;
    cl_int err_code;
    ArrayVariable *var;
    cl_mem mem;
    CalcServer::CalcServer *C = CalcServer::CalcServer::singleton();

    unsigned int n = bounds().y - bounds().x;
    unsigned int *iset = new unsigned int[n];
    unsigned int *id = new unsigned int[n];

    if(!iset || !id){
        LOG(L_ERROR, "Failure allocating memory.\n");
    }

    for(i = 0; i < n; i++){
        iset[i] = setId();
        id[i] = bounds().x + i;
    }

    Variables *vars = C->variables();
    var = (ArrayVariable*)vars->get("iset");
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
        LOG(L_ERROR, "Failure sending variable \"iset\" to the server.\n");
        Logger::singleton()->printOpenCLError(err_code);
        throw std::runtime_error("OpenCL error");
    }
    var = (ArrayVariable*)vars->get("id");
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
        LOG(L_ERROR, "Failure sending variable \"id\" to the server.\n");
        Logger::singleton()->printOpenCLError(err_code);
        throw std::runtime_error("OpenCL error");
    }
    var = (ArrayVariable*)vars->get("id_sorted");
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
        LOG(L_ERROR, "Failure sending variable \"id_sorted\" to the server.\n");
        Logger::singleton()->printOpenCLError(err_code);
        throw std::runtime_error("OpenCL error");
    }
    var = (ArrayVariable*)vars->get("id_unsorted");
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
        LOG(L_ERROR, "Failure sending variable \"id_unsorted\" to the server.\n");
        Logger::singleton()->printOpenCLError(err_code);
        throw std::runtime_error("OpenCL error");
    }

    delete[] iset; iset = NULL;
    delete[] id; id = NULL;
}

unsigned int Particles::file(const std::string basename,
                             unsigned int startindex,
                             unsigned int digits)
{
    unsigned int i = startindex;
    try {
        file(newFilePath(basename, i, digits));
    } catch(std::invalid_argument e) {
        std::ostringstream msg;
        msg << "It is forbidden to overwrite particles output file '"
            << setStrConstantsCopy(basename) << "'" << std::endl;
        LOG(L_ERROR, msg.str());
        throw;
    }
    return i;
}

std::vector<void*> Particles::download(std::vector<std::string> fields)
{
    std::vector<void*> data;
    std::vector<cl_event> events;
    size_t typesize, len;
    cl_int err_code;
    CalcServer::CalcServer *C = CalcServer::CalcServer::singleton();
    Variables *vars = C->variables();

    for(auto field : fields){
        if(!vars->get(field)){
            std::ostringstream msg;
            msg << "Can't download undeclared variable \"" << field
                << "\"." << std::endl;
            LOG(L_ERROR, msg.str());
            clearList(&data);
            throw std::runtime_error("Invalid variable");
        }
        if(vars->get(field)->type().find('*') == std::string::npos){
            std::ostringstream msg;
            msg << "Variable \"" << field << "\" is a scalar." << std::endl;
            LOG(L_ERROR, msg.str());
            clearList(&data);
            throw std::runtime_error("Invalid variable type");
        }
        ArrayVariable *var = (ArrayVariable*)vars->get(field);
        typesize = vars->typeToBytes(var->type());
        len = var->size() / typesize;
        if(len < bounds().y){
            std::ostringstream msg;
            msg << "Variable \"" << field << "\" is not long enough." << std::endl;
            LOG(L_ERROR, msg.str());
            msg.str("");
            msg << "length = " << bounds().y << "is required, but just "
                << len << " components are available." << std::endl;
            LOG0(L_DEBUG, msg.str());
            clearList(&data);
            throw std::runtime_error("Invalid variable length");
        }
        void *store = malloc(typesize * (bounds().y - bounds().x));
        if(!store){
            std::ostringstream msg;
            msg << "Failure allocating " << typesize * (bounds().y - bounds().x)
                << "bytes for variable \"" << field << "\"." << std::endl;
            LOG(L_ERROR, msg.str());
            clearList(&data);
            throw std::bad_alloc();
        }
        data.push_back(store);

        cl_event event;
        try {
            event = C->getUnsortedMem(
                var->name().c_str(),
                typesize * bounds().x,
                typesize * (bounds().y - bounds().x),
                store);
        } catch (...) {
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
        Logger::singleton()->printOpenCLError(err_code);
        clearList(&data);
        throw std::runtime_error("OpenCL error");
    }

    // Destroy the events
    for(auto event : events){
        err_code = clReleaseEvent(event);
        if(err_code != CL_SUCCESS){
            LOG(L_ERROR, "Failure releasing the events.\n");
            Logger::singleton()->printOpenCLError(err_code);
            clearList(&data);
            throw std::runtime_error("OpenCL error");
        }
    }

    return data;
}

void Particles::clearList(std::vector<void*> *data)
{
    for(auto d : *data){
        if(d)
            free(d);
    }
    data->clear();
}

}}  // namespace
