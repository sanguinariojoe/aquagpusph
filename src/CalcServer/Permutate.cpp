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
 * @brief Particles sorting/unsorting tool.
 * (See Aqua::CalcServer::Permutate for details)
 */

#include <ProblemSetup.h>
#include <ScreenManager.h>
#include <CalcServer/Permutate.h>
#include <CalcServer.h>

namespace Aqua{ namespace CalcServer{

Permutate::Permutate()
    : Kernel("Permutate")
    , _path(0)
    , _program(0)
    , _kernel(0)
    , _global_work_size(0)
    , _local_work_size(0)
{
    InputOutput::ScreenManager *S = InputOutput::ScreenManager::singleton();
    InputOutput::ProblemSetup *P = InputOutput::ProblemSetup::singleton();
    unsigned int i, str_len = strlen(P->OpenCL_kernels.permutate);
    if(str_len <= 0) {
        S->addMessageF(3, "The path of the kernel is empty.\n");
        exit(EXIT_FAILURE);
    }
    _path = new char[str_len+4];
    if(!_path) {
        S->addMessageF(3, "Memory cannot be allocated for the path.\n");
        exit(EXIT_FAILURE);
    }
    strcpy(_path, P->OpenCL_kernels.permutate);
    strcat(_path, ".cl");

    _local_work_size  = localWorkSize();
    if(!_local_work_size){
        S->addMessageF(3, "I cannot get a valid local work size for the required computation tool.\n");
        exit(EXIT_FAILURE);
    }
    _global_work_size = globalWorkSize(_local_work_size);
    if(setupOpenCL()) {
        exit(EXIT_FAILURE);
    }
    S->addMessageF(1, "Permutate ready to work!\n");
}

Permutate::~Permutate()
{
    unsigned int i;
    if(_kernel)clReleaseKernel(_kernel); _kernel=0;
    if(_program)clReleaseProgram(_program); _program=0;
    if(_path) delete[] _path; _path=0;
    for(i = 0; i < _backup.size(); i++){
        if(_backup.at(i))
            clReleaseMemObject(_backup.at(i));
        _backup.at(i)=NULL;
    }
}

bool Permutate::sort()
{
    CalcServer *C = CalcServer::singleton();
    return execute(C->permutation_inverse);
}

bool Permutate::unsort()
{
    CalcServer *C = CalcServer::singleton();
    return execute(C->permutation);
}

bool Permutate::execute(cl_mem permutations)
{
    unsigned int i;
    InputOutput::ScreenManager *S = InputOutput::ScreenManager::singleton();
    CalcServer *C = CalcServer::singleton();
    cl_int err_code=0;

    // Duplicate the data arrays
    cl_event events[_mems.size()];
    for(i = 0; i < _mems.size(); i++){
        events[i] = NULL;
        err_code = clEnqueueCopyBuffer(C->command_queue,
                                       _mems.at(i),
                                       _backup.at(i),
                                       0,
                                       0,
                                       _mems_size.at(i),
                                       0,
                                       NULL,
                                       &(events[i]));
        if(err_code != CL_SUCCESS){
            S->addMessageF(3, "Failure duplicating memory buffer.\n");
            return true;
        }
    }

    // Set the input and output arguments
    for(i = 0; i < _mems.size(); i++){
        err_code |= sendArgument(_kernel,
                                 2 * i,
                                 sizeof(cl_mem),
                                 (void*)&(_mems.at(i)));
        err_code |= sendArgument(_kernel,
                                 2 * i + 1,
                                 sizeof(cl_mem),
                                 (void*)&(_backup.at(i)));
    }
    // Set the used permutations
    err_code |= sendArgument(_kernel,
                             2 * _mems.size(),
                             sizeof(cl_mem),
                             (void*)&permutations);
    err_code |= sendArgument(_kernel,
                             2 * _mems.size() + 1,
                             sizeof(cl_uint),
                             (void*)&(C->N));
    if(err_code != CL_SUCCESS){
        S->addMessageF(3, "Failure setting kernel arguments.\n");
        return true;
    }

    // Call the algorithm
    #ifdef HAVE_GPUPROFILE
        cl_event event;
        err_code = clEnqueueNDRangeKernel(C->command_queue,
                                          _kernel,
                                          1,
                                          NULL,
                                          &_global_work_size,
                                          &_local_work_size,
                                          _mems.size(),
                                          events,
                                          &event);
    #else
        err_code = clEnqueueNDRangeKernel(C->command_queue,
                                          _kernel,
                                          1,
                                          NULL,
                                          &_global_work_size,
                                          &_local_work_size,
                                          _mems.size(),
                                          events,
                                          NULL);
    #endif
    if(err_code != CL_SUCCESS) {
        S->addMessageF(3, "I cannot execute the kernel.\n");
        S->printOpenCLError(err_code);
        return true;
    }
    for(i = 0; i < _mems.size(); i++){
        err_code = clReleaseEvent(events[i]);
        if(err_code != CL_SUCCESS) {
            S->addMessageF(3, "Failure releasing an event.\n");
            S->printOpenCLError(err_code);
            return true;
        }
    }
    #ifdef HAVE_GPUPROFILE
        err_code = clWaitForEvents(1, &event);
        if(err_code != CL_SUCCESS) {
            S->addMessage(3, "Impossible to wait for the kernels end.\n");
            S->printOpenCLError(err_code);
            return true;
        }
        err_code |= clGetEventProfilingInfo(event,
                                            CL_PROFILING_COMMAND_END,
                                            sizeof(cl_ulong),
                                            &end,
                                            0);
        if(err_code != CL_SUCCESS) {
            S->addMessage(3, "I cannot profile the kernel execution.\n");
            S->printOpenCLError(err_code);
            return true;
        }
        err_code |= clGetEventProfilingInfo(event,
                                            CL_PROFILING_COMMAND_START,
                                            sizeof(cl_ulong),
                                            &start,
                                            0);
        if(err_code != CL_SUCCESS) {
            S->addMessage(3, "I cannot profile the kernel execution.\n");
            S->printOpenCLError(err_code);
            return true;
        }
        profileTime(profileTime() + (end - start)/1000.f);  // 10^-3 ms
    #endif

    return false;
}

bool Permutate::setupMems()
{
    unsigned int i = 0;
    char msg[256];
    InputOutput::ScreenManager *S = InputOutput::ScreenManager::singleton();
    CalcServer *C = CalcServer::singleton();
    _mems.push_back(C->ifluid);
    _mems.push_back(C->imove);
    for(i = i; i < _mems.size(); i++){
        _mems_size.push_back(C->N * sizeof(cl_int));
        _backup.push_back(C->allocMemory(C->N * sizeof(cl_int)));
        if(!_backup.at(i)) {
            sprintf(msg,
                    "Fail allocating memory (%lu bytes).\n",
                    C->N * sizeof(cl_int));
            S->addMessageF(3, msg);
            return true;
        }
    }
    // 3 vectorials
    _mems.push_back(C->pos);
    _mems.push_back(C->normal);
    _mems.push_back(C->v);
    _mems.push_back(C->f);
    for(i = i; i < _mems.size(); i++){
        _mems_size.push_back(C->N * sizeof(vec));
        _backup.push_back(C->allocMemory(C->N * sizeof(vec)));
        if(!_backup.at(i)) {
            sprintf(msg,
                    "Fail allocating memory (%lu bytes).\n",
                    C->N * sizeof(vec));
            S->addMessageF(3, msg);
            return true;
        }
    }
    // 6 floats
    _mems.push_back(C->dens);
    _mems.push_back(C->press);
    _mems.push_back(C->mass);
    _mems.push_back(C->drdt);
    _mems.push_back(C->drdt_F);
    _mems.push_back(C->shepard);
    for(i = i; i < _mems.size(); i++){
        _mems_size.push_back(C->N * sizeof(cl_float));
        _backup.push_back(C->allocMemory(C->N * sizeof(cl_float)));
        if(!_backup.at(i)) {
            sprintf(msg,
                    "Fail allocating memory (%lu bytes).\n",
                    C->N * sizeof(cl_float));
            S->addMessageF(3, msg);
            return true;
        }
    }
    sprintf(msg, "\tAllocated memory = %lu bytes\n",
            C->allocated_mem);
    S->addMessage(1, msg);
    return false;
}

bool Permutate::setupOpenCL()
{
    InputOutput::ScreenManager *S = InputOutput::ScreenManager::singleton();
    CalcServer *C = CalcServer::singleton();
    char msg[1024];
    cl_int err_code;
    cl_device_id device;
    cl_ulong local_mem, required_local_mem;

    if(setupMems())
        return true;

    err_code = clGetCommandQueueInfo(C->command_queue,
                                     CL_QUEUE_DEVICE,
                                     sizeof(cl_device_id),
                                     &device,
                                     NULL);
    if(err_code != CL_SUCCESS) {
        S->addMessageF(3, "I cannot query the device from the command queue.\n");
        S->printOpenCLError(err_code);
        return true;
    }
    err_code = clGetDeviceInfo(device,
                               CL_DEVICE_LOCAL_MEM_SIZE,
                               sizeof(local_mem),
                               &local_mem,
                               NULL);
    if(err_code != CL_SUCCESS) {
        S->addMessageF(3, "Failure getting the local memory available on the device.\n");
        S->printOpenCLError(err_code);
        return true;
    }
    if(!loadKernelFromFile(&_kernel,
                           &_program,
                           C->context,
                           C->device,
                           _path,
                           "Permutate",
                           ""))
        return true;
    if(_program) clReleaseProgram(_program); _program=0;
    if(!loadKernelFromFile(&_kernel,
                           &_program,
                           C->context,
                           C->device,
                           _path,
                           "Permutate",
                           ""))
        return true;
    if(_program)clReleaseProgram(_program); _program=0;
    // Test if there are enough local memory
    err_code = clGetKernelWorkGroupInfo(_kernel,
                                        device,
                                        CL_KERNEL_LOCAL_MEM_SIZE,
                                        sizeof(cl_ulong),
                                        &required_local_mem,
                                        NULL);
    if(err_code != CL_SUCCESS) {
        S->addMessageF(3, "Failure getting kernel memory usage.\n");
        S->printOpenCLError(err_code);
        return true;
    }
    if(local_mem < required_local_mem){
        S->addMessageF(3, "Not enough local memory for the execution.\n");
        sprintf(msg, "\tNeeds %lu bytes, but only %lu bytes are available.\n",
               required_local_mem, local_mem);
        S->addMessage(0, msg);
        return true;
    }

    // Test if the local work group size must be modified
    size_t local_work_size=0;
    err_code = clGetKernelWorkGroupInfo(_kernel,
                                        device,
                                        CL_KERNEL_WORK_GROUP_SIZE,
                                        sizeof(size_t),
                                        &local_work_size,
                                        NULL);
    if(err_code != CL_SUCCESS) {
        S->addMessageF(3, "A valid maximum local work group size cannot be found.\n");
        S->printOpenCLError(err_code);
        return true;
    }
    if(local_work_size < _local_work_size)
        _local_work_size  = local_work_size;
    // Look for a better local work group size
    err_code = clGetKernelWorkGroupInfo(_kernel,
                                        device,
                                        CL_KERNEL_PREFERRED_WORK_GROUP_SIZE_MULTIPLE,
                                        sizeof(size_t),
                                        &local_work_size,
                                        NULL);
    if(err_code != CL_SUCCESS) {
        S->addMessageF(3, "Preferred local work group size cannot be queried.\n");
        S->printOpenCLError(err_code);
        return true;
    }
    _local_work_size  = (_local_work_size/local_work_size) * local_work_size;
    _global_work_size = globalWorkSize(_local_work_size);
    return false;
}

}}  // namespace
