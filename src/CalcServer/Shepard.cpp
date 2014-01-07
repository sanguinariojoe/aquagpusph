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
// Include the Problem setup manager header
// ----------------------------------------------------------------------------
#include <ScreenManager.h>

// ----------------------------------------------------------------------------
// Include the main header
// ----------------------------------------------------------------------------
#include <CalcServer/Shepard.h>

// ----------------------------------------------------------------------------
// Include the calculation server
// ----------------------------------------------------------------------------
#include <CalcServer.h>

namespace Aqua{ namespace CalcServer{

Shepard::Shepard()
	: Kernel("Shepard")
	, _path(0)
	, program(0)
	, kernel(0)
{
	InputOutput::ScreenManager *S = InputOutput::ScreenManager::singleton();
	InputOutput::ProblemSetup *P = InputOutput::ProblemSetup::singleton();
	if(!P->SPH_opts.has_shepard)  // Shepard correction has not been selected
	    return;
	//! 1st.- Get data
	int nChar = strlen(P->OpenCL_kernels.shepard);
	if(nChar <= 0) {
	    S->addMessage(3, "(Shepard::Shepard): _path of Shepard kernel is empty.\n");
	    exit(EXIT_FAILURE);
	}
	_path = new char[nChar+4];
	if(!_path) {
	    S->addMessage(3, "(Shepard::Shepard): Can't allocate memory for path.\n");
	    exit(EXIT_FAILURE);
	}
	strcpy(_path, P->OpenCL_kernels.shepard);
	strcat(_path, ".cl");
	//! 2nd.- Setup the kernel
	_local_work_size  = localWorkSize();
	if(!_local_work_size){
	    S->addMessage(3, "(Shepard::Shepard): I cannot get a valid local work size for the required computation tool.\n");
	    exit(EXIT_FAILURE);
	}
	_global_work_size = globalWorkSize(_local_work_size);
	if(setupOpenCL()) {
	    exit(EXIT_FAILURE);
	}
	S->addMessage(1, "(Shepard::Shepard): Shepard ready to work!\n");
}

Shepard::~Shepard()
{
	if(kernel)clReleaseKernel(kernel); kernel=0;
	if(program)clReleaseProgram(program); program=0;
	if(_path) delete[] _path; _path=0;
}

bool Shepard::execute()
{
	InputOutput::ProblemSetup *P = InputOutput::ProblemSetup::singleton();
	if(!P->SPH_opts.has_shepard)  // Shepard condition has not been selected
	    return false;
	InputOutput::ScreenManager *S = InputOutput::ScreenManager::singleton();
	CalcServer *C = CalcServer::singleton();
	cl_int err_code=0;
	//! Send all variables to the server
	err_code |= sendArgument(kernel,  0, sizeof(cl_mem  ), (void*)&(C->imove));
	err_code |= sendArgument(kernel,  1, sizeof(cl_mem  ), (void*)&(C->f));
	err_code |= sendArgument(kernel,  2, sizeof(cl_mem  ), (void*)&(C->drdt));
	err_code |= sendArgument(kernel,  3, sizeof(cl_mem  ), (void*)&(C->shepard));
	err_code |= sendArgument(kernel,  4, sizeof(cl_uint ), (void*)&(C->n));
	if(err_code != CL_SUCCESS) {
		S->addMessage(3, "(Shepard::Shepard): Can't send arguments to Shepard computation kernel.\n");
	    return true;
	}
	//! Execute the kernel
	size_t globalWorkSize = getGlobalWorkSize(C->n, _local_work_size);
	#ifdef HAVE_GPUPROFILE
	    cl_event event;
	    cl_ulong end, start;
	    err_code = clEnqueueNDRangeKernel(C->command_queue, kernel, 1, NULL, &globalWorkSize, NULL, 0, NULL, &event);
	#else
	    err_code = clEnqueueNDRangeKernel(C->command_queue, kernel, 1, NULL, &globalWorkSize, NULL, 0, NULL, NULL);
	#endif
	if(err_code != CL_SUCCESS) {
		S->addMessage(3, "(Shepard::Shepard): Can't execute the kernel.\n");
	    if(err_code == CL_INVALID_WORK_GROUP_SIZE)
	        S->addMessage(0, "\tInvalid local work group size.\n");
	    else if(err_code == CL_OUT_OF_RESOURCES)
	        S->addMessage(0, "\tDevice out of resources.\n");
	    else if(err_code == CL_MEM_OBJECT_ALLOCATION_FAILURE)
	        S->addMessage(0, "\tAllocation error at device.\n");
	    else if(err_code == CL_OUT_OF_HOST_MEMORY)
	        S->addMessage(0, "\tfailure to allocate resources required by the OpenCL implementation on the host.\n");
	    return true;
	}
	//! Profile the kernel execution
	#ifdef HAVE_GPUPROFILE
	    err_code = clWaitForEvents(1, &event);
	    if(err_code != CL_SUCCESS) {
	        S->addMessage(3, "(Shepard::Shepard): Can't wait to kernels end.\n");
	        return true;
	    }
	    err_code |= clGetEventProfilingInfo(event, CL_PROFILING_COMMAND_END, sizeof(cl_ulong), &end, 0);
	    err_code |= clGetEventProfilingInfo(event, CL_PROFILING_COMMAND_START, sizeof(cl_ulong), &start, 0);
	    if(err_code != CL_SUCCESS) {
	        S->addMessage(3, "(Shepard::Shepard): Can't profile kernel execution.\n");
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
	// Create arguments
	char args[256];
	strcpy(args, "");
	if(P->SPH_opts.has_shepard & 1){
	    strcat(args, "-D__FORCE_CORRECTION__ ");
	}
	if(P->SPH_opts.has_shepard & 2){
	    strcat(args, "-D__DENS_CORRECTION__ ");
	}
	if(!loadKernelFromFile(&kernel, &program, C->context, C->device, _path, "Shepard", args))
	    return true;
	if(program)clReleaseProgram(program); program=0;
	//! Test for right work group size
	cl_device_id device;
	size_t localWorkGroupSize=0;
	err_code |= clGetCommandQueueInfo(C->command_queue,CL_QUEUE_DEVICE,
	                                sizeof(cl_device_id),&device, NULL);
	if(err_code != CL_SUCCESS) {
		S->addMessage(3, "(Predictor::setupOpenCL): I Cannot get the device from the command queue.\n");
	    return true;
	}
	err_code |= clGetKernelWorkGroupInfo(kernel,device,CL_KERNEL_WORK_GROUP_SIZE,
	                                   sizeof(size_t), &localWorkGroupSize, NULL);
	if(err_code != CL_SUCCESS) {
		S->addMessage(3, "(Predictor::setupOpenCL): Can't get maximum local work group size.\n");
	    return true;
	}
	if(localWorkGroupSize < _local_work_size)
	    _local_work_size  = localWorkGroupSize;
	_global_work_size = globalWorkSize(_local_work_size);
	return false;
}

}}  // namespaces
