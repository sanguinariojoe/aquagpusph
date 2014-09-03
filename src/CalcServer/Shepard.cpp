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

#include <ProblemSetup.h>
#include <ScreenManager.h>
#include <CalcServer/Shepard.h>
#include <CalcServer.h>

namespace Aqua{ namespace CalcServer{

Shepard::Shepard()
    : Kernel("Shepard")
    , _path(0)
    , _program(0)
    , _kernel(0)
{
    InputOutput::ScreenManager *S = InputOutput::ScreenManager::singleton();
    InputOutput::ProblemSetup *P = InputOutput::ProblemSetup::singleton();
    if(!P->SPH_opts.has_shepard)
        return;

    int str_len = strlen(P->OpenCL_kernels.shepard);
    if(str_len <= 0) {
        S->addMessageF(3, "The path of the kernel is empty.\n");
        exit(EXIT_FAILURE);
    }
    _path = new char[str_len+4];
    if(!_path) {
        S->addMessageF(3, "Memory cannot be allocated for the path.\n");
        exit(EXIT_FAILURE);
    }
    strcpy(_path, P->OpenCL_kernels.shepard);
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
    S->addMessageF(1, "Shepard ready to work!\n");
}

Shepard::~Shepard()
{
    if(_kernel)clReleaseKernel(_kernel); _kernel=0;
    if(_program)clReleaseProgram(_program); _program=0;
    if(_path) delete[] _path; _path=0;
}

bool Shepard::execute()
{
    InputOutput::ProblemSetup *P = InputOutput::ProblemSetup::singleton();
    if(!P->SPH_opts.has_shepard)
        return false;
    InputOutput::ScreenManager *S = InputOutput::ScreenManager::singleton();
    CalcServer *C = CalcServer::singleton();
    cl_int err_code=0;

    err_code |= sendArgument(_kernel,
                             0,
                             sizeof(cl_mem),
                             (void*)&(C->imove));
    err_code |= sendArgument(_kernel,
                             1,
                             sizeof(cl_mem),
                             (void*)&(C->f));
    err_code |= sendArgument(_kernel,
                             2,
                             sizeof(cl_mem),
                             (void*)&(C->drdt));
    err_code |= sendArgument(_kernel,
                             3,
                             sizeof(cl_mem),
                             (void*)&(C->drdt_F));
    err_code |= sendArgument(_kernel,
                             4,
                             sizeof(cl_mem),
                             (void*)&(C->shepard));
    err_code |= sendArgument(_kernel,
                             5,
                             sizeof(cl_uint),
                             (void*)&(C->n));
    if(err_code != CL_SUCCESS) {
        S->addMessageF(3, "Failure sending variables to Shepard computation kernel.\n");
        return true;
    }
    //! Execute the kernel
    size_t globalWorkSize = getGlobalWorkSize(C->n, _local_work_size);
    #ifdef HAVE_GPUPROFILE
        cl_event event;
        cl_ulong end, start;
        err_code = clEnqueueNDRangeKernel(C->command_queue,
                                          _kernel,
                                          1,
                                          NULL,
                                          &globalWorkSize,
                                          NULL,
                                          0,
                                          NULL,
                                          &event);
    #else
        err_code = clEnqueueNDRangeKernel(C->command_queue,
                                          _kernel,
                                          1,
                                          NULL,
                                          &globalWorkSize,
                                          NULL,
                                          0,
                                          NULL,
                                          NULL);
    #endif
    if(err_code != CL_SUCCESS) {
        S->addMessageF(3, "I cannot execute the kernel.\n");
        S->printOpenCLError(err_code);
        return true;
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

bool Shepard::setupOpenCL()
{
    InputOutput::ScreenManager *S = InputOutput::ScreenManager::singleton();
    InputOutput::ProblemSetup *P  = InputOutput::ProblemSetup::singleton();
    CalcServer *C = CalcServer::singleton();
    cl_int err_code;

    char args[256];
    strcpy(args, "");
    if(P->SPH_opts.has_shepard & 1){
        strcat(args, "-D__FORCE_CORRECTION__ ");
    }
    if(P->SPH_opts.has_shepard & 2){
        strcat(args, "-D__DENS_CORRECTION__ ");
    }
    if(!loadKernelFromFile(&_kernel,
                           &_program,
                           C->context,
                           C->device,
                           _path,
                           "Shepard",
                           args))
        return true;
    if(_program)clReleaseProgram(_program); _program=0;

    cl_device_id device;
    size_t local_work_size=0;
    err_code = clGetCommandQueueInfo(C->command_queue,
                                     CL_QUEUE_DEVICE,
                                     sizeof(cl_device_id),
                                     &device,
                                     NULL);
    if(err_code != CL_SUCCESS) {
        S->addMessageF(3, "I Cannot get the device from the command queue.\n");
        S->printOpenCLError(err_code);
        return true;
    }
    err_code = clGetKernelWorkGroupInfo(_kernel,
                                        device,
                                        CL_KERNEL_WORK_GROUP_SIZE,
                                        sizeof(size_t),
                                        &local_work_size,
                                        NULL);
    if(err_code != CL_SUCCESS) {
        S->addMessageF(3, "Failure retrieving the maximum local work size.\n");
        S->printOpenCLError(err_code);
        return true;
    }
    if(local_work_size < _local_work_size)
        _local_work_size  = local_work_size;
    _global_work_size = globalWorkSize(_local_work_size);
    return false;
}

}}  // namespaces
