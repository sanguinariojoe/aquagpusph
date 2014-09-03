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

#include <math.h>

#include <ProblemSetup.h>
#include <ScreenManager.h>
#include <CalcServer/Boundary/ElasticBounce.h>
#include <CalcServer.h>

namespace Aqua{ namespace CalcServer{ namespace Boundary{

ElasticBounce::ElasticBounce()
    : Kernel("ElasticBounce")
    , _path(NULL)
    , _program(NULL)
    , _kernel(NULL)
{
    InputOutput::ScreenManager *S = InputOutput::ScreenManager::singleton();
    InputOutput::ProblemSetup *P = InputOutput::ProblemSetup::singleton();
    if(P->SPH_opts.elastic_dist <= 0.f){
        // Elastic bounce is disabled
        return;
    }
    if(   P->SPH_opts.boundary_type!=0  // ElasticBounce boundary
       && P->SPH_opts.boundary_type!=1  // Fixed particles boundary
       && P->SPH_opts.boundary_type!=2  // DeLeffe boundary
       )
        return;

    int str_len = strlen(P->OpenCL_kernels.elastic_bounce);
    if(str_len <= 0) {
        S->addMessageF(3, "The path of the kernel is empty.\n");
        exit(EXIT_FAILURE);
    }
    _path = new char[str_len+4];
    if(!_path) {
        S->addMessageF(3, "Memory cannot be allocated for the path.\n");
        exit(EXIT_FAILURE);
    }
    strcpy(_path, P->OpenCL_kernels.elastic_bounce);
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

    vec deltar = P->SPH_opts.deltar;
    #ifdef HAVE_3D
        _r = sqrt(2.f) * (deltar.x + deltar.y + deltar.z) / 3.f;
    #else
        _r = sqrt(2.f) * (deltar.x + deltar.y) / 2.f;
    #endif
    S->addMessageF(1, "ElasticBounce boundary condition ready to work!\n");
}

ElasticBounce::~ElasticBounce()
{
    if(_kernel)clReleaseKernel(_kernel); _kernel=0;
    if(_program)clReleaseProgram(_program); _program=0;
    if(_path) delete[] _path; _path=0;
}

bool ElasticBounce::execute()
{
    InputOutput::ProblemSetup *P = InputOutput::ProblemSetup::singleton();
    if(P->SPH_opts.elastic_dist <= 0.f){
        // Elastic bounce is disabled
        return false;
    }
    if(   P->SPH_opts.boundary_type!=0  // ElasticBounce boundary
       && P->SPH_opts.boundary_type!=1  // Fixed particles boundary
       && P->SPH_opts.boundary_type!=2) // DeLeffe boundary
        return false;
    InputOutput::ScreenManager *S = InputOutput::ScreenManager::singleton();
    CalcServer *C = CalcServer::singleton();
    cl_int err_code=0;

    err_code  = sendArgument(_kernel, 0, sizeof(cl_mem),
                             (void*)&(C->imove));
    err_code |= sendArgument(_kernel, 1, sizeof(cl_mem),
                             (void*)&(C->pos));
    err_code |= sendArgument(_kernel, 2, sizeof(cl_mem),
                             (void*)&(C->normal));
    err_code |= sendArgument(_kernel, 3, sizeof(cl_mem),
                             (void*)&(C->v));
    err_code |= sendArgument(_kernel, 4, sizeof(cl_mem),
                             (void*)&(C->f));
    err_code |= sendArgument(_kernel, 5, sizeof(cl_mem),
                             (void*)&(C->icell));
    err_code |= sendArgument(_kernel, 6, sizeof(cl_mem),
                             (void*)&(C->ihoc));
    err_code |= sendArgument(_kernel, 7, sizeof(cl_uint),
                             (void*)&(C->N));
    err_code |= sendArgument(_kernel, 8, sizeof(cl_float),
                             (void*)&(C->dt));
    err_code |= sendArgument(_kernel, 9, sizeof(uivec),
                             (void*)&(C->num_cells_vec));
    err_code |= sendArgument(_kernel, 10, sizeof(vec),
                             (void*)&(C->g));
    err_code |= sendArgument(_kernel, 11, sizeof(cl_float),
                             (void*)&_r);
    if(err_code != CL_SUCCESS) {
        S->addMessageF(3, "Failure sending variables to boundary computation kernel.\n");
        return true;
    }

    #ifdef HAVE_GPUPROFILE
        cl_event event;
        cl_ulong end, start;
        profileTime(0.f);
        err_code = clEnqueueNDRangeKernel(C->command_queue,
                                          _kernel,
                                          1,
                                          NULL,
                                          &_global_work_size,
                                          NULL,
                                          0,
                                          NULL,
                                          &event);
    #else
        err_code = clEnqueueNDRangeKernel(C->command_queue,
                                          _kernel,
                                          1,
                                          NULL,
                                          &_global_work_size,
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

bool ElasticBounce::setupOpenCL()
{
    InputOutput::ScreenManager *S = InputOutput::ScreenManager::singleton();
    InputOutput::ProblemSetup *P = InputOutput::ProblemSetup::singleton();
    CalcServer *C = CalcServer::singleton();
    char msg[1024];
    cl_int err_code;
    cl_device_id device;
    cl_ulong local_mem, required_local_mem;
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
    char flags[512];
    sprintf(flags, "-D__ELASTIC_FACTOR__=%ff", P->SPH_opts.elastic_factor);
    sprintf(flags, "%s -D__MIN_BOUND_DIST__=%ff", flags, fabs(P->SPH_opts.elastic_dist));
    if(P->SPH_opts.elastic_dist < 0.f){
        sprintf(flags, "%s -D__FORCE_MIN_BOUND_DIST__", flags);
    }
    if(!loadKernelFromFile(&_kernel,
                           &_program,
                           C->context,
                           C->device,
                           _path,
                           "Boundary",
                           flags))
        return true;
    if(_program)clReleaseProgram(_program); _program=0;

    err_code = clGetKernelWorkGroupInfo(_kernel,
                                        device,
                                        CL_KERNEL_LOCAL_MEM_SIZE,
                                        sizeof(cl_ulong),
                                        &required_local_mem,
                                        NULL);
    if(err_code != CL_SUCCESS) {
        S->addMessageF(3, "Error retrieving the used local memory.\n");
        S->printOpenCLError(err_code);
        return true;
    }
    if(local_mem < required_local_mem){
        S->addMessageF(3, "There are not enough local memory in the device.\n");
        sprintf(msg, "\tNeeds %lu bytes, but only %lu bytes are available.\n",
               required_local_mem, local_mem);
        S->addMessage(0, msg);
        return true;
    }

    size_t local_work_size=0;
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

    err_code = clGetKernelWorkGroupInfo(_kernel,
                                        device,
                                        CL_KERNEL_PREFERRED_WORK_GROUP_SIZE_MULTIPLE,
                                        sizeof(size_t),
                                        &local_work_size,
                                        NULL);
    if(err_code != CL_SUCCESS) {
        S->addMessageF(3, "I cannot query the preferred local work size");
        S->printOpenCLError(err_code);
        return true;
    }
    _local_work_size = (_local_work_size/local_work_size) * local_work_size;
    _global_work_size = globalWorkSize(_local_work_size);
    return false;
}

}}}  // namespaces
