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
#include <TimeManager.h>
#include <ScreenManager.h>
#include <CalcServer/Predictor.h>
#include <CalcServer.h>

namespace Aqua{ namespace CalcServer{

Predictor::Predictor()
	: Kernel("Predictor")
	, _path(0)
	, _program(0)
	, _kernel(0)
{
	InputOutput::ProblemSetup *P  = InputOutput::ProblemSetup::singleton();
	InputOutput::ScreenManager *S = InputOutput::ScreenManager::singleton();
	unsigned int str_len = strlen(P->OpenCL_kernels.predictor);
	if(str_len <= 0) {
	    S->addMessageF(3, "The path of the predictor kernel is empty.\n");
	    exit(EXIT_FAILURE);
	}
	_path = new char[str_len+4];
	if(!_path) {
	    S->addMessageF(3, "Memory cannot be allocated for the path.\n");
	    exit(EXIT_FAILURE);
	}
	strcpy(_path, P->OpenCL_kernels.predictor);
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
	S->addMessageF(1, "Predictor ready to work!\n");
}

Predictor::~Predictor()
{
	if(_kernel)clReleaseKernel(_kernel); _kernel=0;
	if(_program)clReleaseProgram(_program); _program=0;
	if(_path)delete[] _path; _path=0;
}

bool Predictor::execute()
{
	InputOutput::ProblemSetup *P = InputOutput::ProblemSetup::singleton();
	InputOutput::TimeManager *T = InputOutput::TimeManager::singleton();
	InputOutput::ScreenManager *S = InputOutput::ScreenManager::singleton();
	CalcServer *C = CalcServer::singleton();
	cl_int err_code=0;

	float t = T->time();
	err_code |= sendArgument(_kernel, 17, sizeof(cl_uint),
                             (void*)&(C->N));
	err_code |= sendArgument(_kernel, 18, sizeof(cl_float),
                             (void*)&t);
	err_code |= sendArgument(_kernel, 19, sizeof(cl_float),
                             (void*)&(C->dt));
	err_code |= sendArgument(_kernel, 20, sizeof(cl_float),
                             (void*)&(C->cs));
	err_code |= sendArgument(_kernel, 21, sizeof(vec),
                             (void*)&(C->g));
	err_code |= sendArgument(_kernel, 22, sizeof(cl_float),
                             (void*)&(P->SPH_opts.rho_min));
	err_code |= sendArgument(_kernel, 23, sizeof(cl_float),
                             (void*)&(P->SPH_opts.rho_max));
	if(err_code != CL_SUCCESS) {
		S->addMessageF(3, "I cannot send a variable to the kernel.\n");
	    return true;
	}
	//! 2nd.- Execute the kernel
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

bool Predictor::setupOpenCL()
{
	InputOutput::ScreenManager *S = InputOutput::ScreenManager::singleton();
	CalcServer *C = CalcServer::singleton();
	cl_int err_code=0;
	if(!loadKernelFromFile(&_kernel, &_program, C->context, C->device, _path,
                           "Predictor", ""))
	    return true;
	err_code |= sendArgument(_kernel, 0, sizeof(cl_mem),
                             (void*)&(C->imove));
	err_code |= sendArgument(_kernel, 1, sizeof(cl_mem),
                             (void*)&(C->ifluid));
	err_code |= sendArgument(_kernel, 2, sizeof(cl_mem),
                             (void*)&(C->pos));
	err_code |= sendArgument(_kernel, 3, sizeof(cl_mem),
                             (void*)&(C->v));
	err_code |= sendArgument(_kernel, 4, sizeof(cl_mem),
                             (void*)&(C->f));
	err_code |= sendArgument(_kernel, 5, sizeof(cl_mem),
                             (void*)&(C->dens));
	err_code |= sendArgument(_kernel, 6, sizeof(cl_mem),
                             (void*)&(C->mass));
	err_code |= sendArgument(_kernel, 7, sizeof(cl_mem),
                             (void*)&(C->drdt));
	err_code |= sendArgument(_kernel, 8, sizeof(cl_mem),
                             (void*)&(C->posin));
	err_code |= sendArgument(_kernel, 9, sizeof(cl_mem),
                             (void*)&(C->vin));
	err_code |= sendArgument(_kernel, 10, sizeof(cl_mem),
                             (void*)&(C->fin));
	err_code |= sendArgument(_kernel, 11, sizeof(cl_mem),
                             (void*)&(C->densin));
	err_code |= sendArgument(_kernel, 12, sizeof(cl_mem),
                             (void*)&(C->massin));
	err_code |= sendArgument(_kernel, 13, sizeof(cl_mem),
                             (void*)&(C->drdtin));
	err_code |= sendArgument(_kernel, 14, sizeof(cl_mem),
                             (void*)&(C->press));
	err_code |= sendArgument(_kernel, 15, sizeof(cl_mem),
                             (void*)&(C->refd));
	err_code |= sendArgument(_kernel, 16, sizeof(cl_mem),
                             (void*)&(C->gamma));
	if(err_code)
	    return true;

	cl_device_id device;
	size_t local_work_size=0;
	err_code = clGetCommandQueueInfo(C->command_queue,
                                     CL_QUEUE_DEVICE,
	                                 sizeof(cl_device_id),
                                     &device,
                                     NULL);
	if(err_code != CL_SUCCESS) {
		S->addMessageF(3, "I Cannot get the device from the command queue.\n");
	    return true;
	}
	err_code = clGetKernelWorkGroupInfo(_kernel,device,
                                        CL_KERNEL_WORK_GROUP_SIZE,
	                                    sizeof(size_t),
                                        &local_work_size,
                                        NULL);
	if(err_code != CL_SUCCESS) {
		S->addMessageF(3, "Failure retrieving the maximum local work size.\n");
	    return true;
	}
	if(local_work_size < _local_work_size)
	    _local_work_size  = local_work_size;
	_global_work_size = globalWorkSize(_local_work_size);
	return false;
}

}}
