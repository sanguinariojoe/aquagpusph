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

// ----------------------------------------------------------------------------
// Include the Problem setup manager header
// ----------------------------------------------------------------------------
#include <ProblemSetup.h>

// ----------------------------------------------------------------------------
// Include the Time manager header
// ----------------------------------------------------------------------------
#include <TimeManager.h>

// ----------------------------------------------------------------------------
// Include the Problem setup manager header
// ----------------------------------------------------------------------------
#include <ScreenManager.h>

// ----------------------------------------------------------------------------
// Include the main header
// ----------------------------------------------------------------------------
#include <CalcServer/Corrector.h>

// ----------------------------------------------------------------------------
// Include the calculation server
// ----------------------------------------------------------------------------
#include <CalcServer.h>

namespace Aqua{ namespace CalcServer{

Corrector::Corrector()
	: Kernel("Corrector")
	, mPath(0)
	, program(0)
	, kernel(0)
	, clClampVKernel(0)
{
	//! 1st.- Get data
	InputOutput::ScreenManager *S = InputOutput::ScreenManager::singleton();
	InputOutput::ProblemSetup *P  = InputOutput::ProblemSetup::singleton();
	unsigned int nChar = strlen(P->OpenCL_kernels.corrector);
	if(nChar <= 0) {
	    S->addMessage(3, "(Corrector::Corrector): Path of corrector kernel is empty.\n");
	    exit(EXIT_FAILURE);
	}
	mPath = new char[nChar+4];
	if(!mPath) {
	    S->addMessage(3, "(Corrector::Corrector): Can't allocate memory for path.\n");
	    exit(EXIT_FAILURE);
	}
	strcpy(mPath, P->OpenCL_kernels.corrector);
	strcat(mPath, ".cl");
	//! 2nd.- Setup the kernel
	local_work_size  = localWorkSize();
	if(!local_work_size){
	    S->addMessage(3, "(Corrector::Corrector): No valid local work size for required computation.\n");
	    exit(EXIT_FAILURE);
	}
	global_work_size = globalWorkSize(local_work_size);
	if(setupOpenCL()) {
	    exit(EXIT_FAILURE);
	}
	S->addMessage(1, "(Corrector::Corrector): Corrector ready to work!\n");
}

Corrector::~Corrector()
{
	if(kernel)clReleaseKernel(kernel); kernel=0;
	if(clClampVKernel)clReleaseKernel(clClampVKernel); clClampVKernel=0;
	if(program)clReleaseProgram(program); program=0;
	if(mPath)delete[] mPath; mPath=0;
}

bool Corrector::execute()
{
	InputOutput::ProblemSetup *P  = InputOutput::ProblemSetup::singleton();
	InputOutput::TimeManager *T   = InputOutput::TimeManager::singleton();
	InputOutput::ScreenManager *S = InputOutput::ScreenManager::singleton();
	CalcServer *C = CalcServer::singleton();
	cl_int err_code;
	float t = T->time();
	/** Execute the velocity clamping if has been requested. If minimum time step is
	 * not able to do it (\f$ dt_{min} \le 0 \mbox{s} \f$) we skip this step
	 */
	if( (P->time_opts.velocity_clamp) && (P->time_opts.dt_min > 0.f) && (t >= 0.f)){
	    err_code  = sendArgument(clClampVKernel, 5, sizeof(cl_float), (void*)&(P->time_opts.dt_min));
	    err_code |= sendArgument(clClampVKernel, 6, sizeof(vec),      (void*)&(C->g));
	    err_code |= sendArgument(clClampVKernel, 7, sizeof(cl_uint),  (void*)&(C->N));
	    if(err_code != CL_SUCCESS) {
	        S->addMessage(3, "(Corrector::execute): Can't send variable to velocity clamping kernel.\n");
	        return true;
	    }
	    #ifdef HAVE_GPUPROFILE
	        cl_event event;
	        cl_ulong end, start;
	        profileTime(0.f);
	        err_code = clEnqueueNDRangeKernel(C->command_queue, clClampVKernel, 1, NULL, &global_work_size, NULL, 0, NULL, &event);
	    #else
	        err_code = clEnqueueNDRangeKernel(C->command_queue, clClampVKernel, 1, NULL, &global_work_size, NULL, 0, NULL, NULL);
	    #endif
	    if(err_code != CL_SUCCESS) {
	        S->addMessage(3, "(Corrector::execute): Can't execute the velocity clmaping kernel.\n");
	        if(err_code == CL_INVALID_WORK_GROUP_SIZE) {
	            S->addMessage(0, "\tInvalid local work group size.\n");
	        }
	        else if(err_code == CL_OUT_OF_RESOURCES) {
	            S->addMessage(0, "\tDevice out of resources.\n");
	        }
	        else if(err_code == CL_MEM_OBJECT_ALLOCATION_FAILURE) {
	            S->addMessage(0, "\tAllocation error at device.\n");
	        }
	        else if(err_code == CL_OUT_OF_HOST_MEMORY) {
	            S->addMessage(0, "\tfailure to allocate resources required by the OpenCL implementation on the host.\n");
	        }
	        return true;
	    }
	    #ifdef HAVE_GPUPROFILE
	        err_code = clWaitForEvents(1, &event);
	        if(err_code != CL_SUCCESS) {
	            S->addMessage(3, "(Corrector::execute): Can't wait to velocity clamping kernel end.\n");
	            return true;
	        }
	        err_code |= clGetEventProfilingInfo(event, CL_PROFILING_COMMAND_END, sizeof(cl_ulong), &end, 0);
	        err_code |= clGetEventProfilingInfo(event, CL_PROFILING_COMMAND_START, sizeof(cl_ulong), &start, 0);
	        if(err_code != CL_SUCCESS) {
	            S->addMessage(3, "(Corrector::execute): Can't profile velocity clamping kernel execution.\n");
	            return true;
	        }
	        profileTime(profileTime() + (end - start)/1000.f);  // 10^-3 ms
	    #endif
	}
	//! Execute the corrector
	err_code  = sendArgument(kernel, 16, sizeof(cl_uint), (void*)&(C->N));
	err_code |= sendArgument(kernel, 17, sizeof(cl_float), (void*)&t);
	err_code |= sendArgument(kernel, 18, sizeof(cl_float), (void*)&(C->dt));
	if(err_code != CL_SUCCESS) {
		S->addMessage(3, "(Corrector::execute): Can't send variable to kernel.\n");
	    return true;
	}
	#ifdef HAVE_GPUPROFILE
	    cl_event event;
	    cl_ulong end, start;
	    profileTime(0.f);
	    err_code = clEnqueueNDRangeKernel(C->command_queue, kernel, 1, NULL, &global_work_size, NULL, 0, NULL, &event);
	#else
	    err_code = clEnqueueNDRangeKernel(C->command_queue, kernel, 1, NULL, &global_work_size, NULL, 0, NULL, NULL);
	#endif
	if(err_code != CL_SUCCESS) {
		S->addMessage(3, "(Corrector::execute): Can't execute the kernel.\n");
	    if(err_code == CL_INVALID_WORK_GROUP_SIZE) {
	        S->addMessage(0, "\tInvalid local work group size.\n");
	    }
	    else if(err_code == CL_OUT_OF_RESOURCES) {
	        S->addMessage(0, "\tDevice out of resources.\n");
	    }
	    else if(err_code == CL_MEM_OBJECT_ALLOCATION_FAILURE) {
	        S->addMessage(0, "\tAllocation error at device.\n");
	    }
	    else if(err_code == CL_OUT_OF_HOST_MEMORY) {
	        S->addMessage(0, "\tfailure to allocate resources required by the OpenCL implementation on the host.\n");
	    }
	    return true;
	}
	#ifdef HAVE_GPUPROFILE
	    err_code = clWaitForEvents(1, &event);
	    if(err_code != CL_SUCCESS) {
	        S->addMessage(3, "(Corrector::execute): Can't wait to kernels end.\n");
	        return true;
	    }
	    err_code |= clGetEventProfilingInfo(event, CL_PROFILING_COMMAND_END, sizeof(cl_ulong), &end, 0);
	    err_code |= clGetEventProfilingInfo(event, CL_PROFILING_COMMAND_START, sizeof(cl_ulong), &start, 0);
	    if(err_code != CL_SUCCESS) {
	        S->addMessage(3, "(Corrector::execute): Can't profile kernel execution.\n");
	        return true;
	    }
	    profileTime(profileTime() + (end - start)/1000.f);  // 10^-3 ms
	#endif
	return false;
}

bool Corrector::setupOpenCL()
{
	InputOutput::ProblemSetup *P = InputOutput::ProblemSetup::singleton();
	InputOutput::ScreenManager *S = InputOutput::ScreenManager::singleton();
	CalcServer *C = CalcServer::singleton();
	cl_int err_code;
	if(!loadKernelFromFile(&kernel, &program, C->context, C->device, mPath, "Corrector", ""))
	    return true;
	err_code  = sendArgument(kernel,  0, sizeof(cl_mem), (void*)&(C->imove));
	err_code |= sendArgument(kernel,  1, sizeof(cl_mem), (void*)&(C->pos));
	err_code |= sendArgument(kernel,  2, sizeof(cl_mem), (void*)&(C->v));
	err_code |= sendArgument(kernel,  3, sizeof(cl_mem), (void*)&(C->f));
	err_code |= sendArgument(kernel,  4, sizeof(cl_mem), (void*)&(C->dens));
	err_code |= sendArgument(kernel,  5, sizeof(cl_mem), (void*)&(C->mass));
	err_code |= sendArgument(kernel,  6, sizeof(cl_mem), (void*)&(C->drdt));
	err_code |= sendArgument(kernel,  7, sizeof(cl_mem), (void*)&(C->hp));
	err_code |= sendArgument(kernel,  8, sizeof(cl_mem), (void*)&(C->posin));
	err_code |= sendArgument(kernel,  9, sizeof(cl_mem), (void*)&(C->vin));
	err_code |= sendArgument(kernel, 10, sizeof(cl_mem), (void*)&(C->fin));
	err_code |= sendArgument(kernel, 11, sizeof(cl_mem), (void*)&(C->densin));
	err_code |= sendArgument(kernel, 12, sizeof(cl_mem), (void*)&(C->massin));
	err_code |= sendArgument(kernel, 13, sizeof(cl_mem), (void*)&(C->drdtin));
	err_code |= sendArgument(kernel, 14, sizeof(cl_mem), (void*)&(C->hpin));
	err_code |= sendArgument(kernel, 15, sizeof(cl_mem), (void*)&(C->sigma));
	if(err_code)
	    return true;
	if( (P->time_opts.velocity_clamp) && (P->time_opts.dt_min > 0) ){
	    if(program)clReleaseProgram(program); program=0;
	    if(!loadKernelFromFile(&clClampVKernel, &program, C->context, C->device, mPath, "ClampVel",""))
	        return true;
	    err_code  = sendArgument(clClampVKernel, 0, sizeof(cl_mem), (void*)&(C->imove));
	    err_code |= sendArgument(clClampVKernel, 1, sizeof(cl_mem), (void*)&(C->v));
	    err_code |= sendArgument(clClampVKernel, 2, sizeof(cl_mem), (void*)&(C->f));
	    err_code |= sendArgument(clClampVKernel, 3, sizeof(cl_mem), (void*)&(C->hp));
	    err_code |= sendArgument(clClampVKernel, 4, sizeof(cl_mem), (void*)&(C->fin));
	    if(err_code)
	        return true;
	}
	//! Test for right work group size
	cl_device_id device;
	size_t localWorkGroupSize=0;
	err_code |= clGetCommandQueueInfo(C->command_queue,CL_QUEUE_DEVICE,
	                                sizeof(cl_device_id),&device, NULL);
	if(err_code != CL_SUCCESS) {
		S->addMessage(3, "(Corrector::setupOpenCL): Can't get device from command queue.\n");
	    return true;
	}
	err_code |= clGetKernelWorkGroupInfo(kernel,device,CL_KERNEL_WORK_GROUP_SIZE,
	                                   sizeof(size_t), &localWorkGroupSize, NULL);
	if(err_code != CL_SUCCESS) {
		S->addMessage(3, "(Corrector::setupOpenCL): Can't get maximum local work group size.\n");
	    return true;
	}
	if(localWorkGroupSize < local_work_size)
	    local_work_size  = localWorkGroupSize;
	global_work_size = globalWorkSize(local_work_size);
	return false;
}

}}  // namespace
