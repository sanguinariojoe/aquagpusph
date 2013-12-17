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
#include <CalcServer/Domain.h>

// ----------------------------------------------------------------------------
// Include the calculation server
// ----------------------------------------------------------------------------
#include <CalcServer.h>

namespace Aqua{ namespace CalcServer{

Domain::Domain()
	: Kernel("Domain")
	, mPath(0)
	, clProgram(0)
	, clKernel(0)
{
	//! 1st.- Get data
	InputOutput::ScreenManager *S = InputOutput::ScreenManager::singleton();
	InputOutput::ProblemSetup *P  = InputOutput::ProblemSetup::singleton();
	if(!P->SPHParameters.hasDomain)
	    return;
	unsigned int nChar = strlen(P->OpenCLKernels.Domain);
	if(nChar <= 0) {
	    S->addMessage(3, "(Domain::Domain): Path of Domain kernel is empty.\n");
	    exit(EXIT_FAILURE);
	}
	mPath = new char[nChar+4];
	if(!mPath) {
	    S->addMessage(3, "(Domain::Domain): Can't allocate memory for path.\n");
	    exit(EXIT_FAILURE);
	}
	strcpy(mPath, P->OpenCLKernels.Domain);
	strcat(mPath, ".cl");
	//! 2nd.- Setup the kernel
	clLocalWorkSize  = localWorkSize();
	if(!clLocalWorkSize){
	    S->addMessage(3, "(Domain::Domain): No valid local work size for required computation.\n");
	    exit(EXIT_FAILURE);
	}
	clGlobalWorkSize = globalWorkSize(clLocalWorkSize);
	if(setupOpenCL()) {
	    exit(EXIT_FAILURE);
	}
	S->addMessage(1, "(Domain::Domain): Domain ready to work!\n");
}

Domain::~Domain()
{
	InputOutput::ProblemSetup *P = InputOutput::ProblemSetup::singleton();
	if(!P->SPHParameters.hasDomain)
	    return;
	if(clKernel)clReleaseKernel(clKernel); clKernel=0;
	if(clProgram)clReleaseProgram(clProgram); clProgram=0;
	if(mPath)delete[] mPath; mPath=0;
}

bool Domain::execute()
{
	InputOutput::ProblemSetup *P = InputOutput::ProblemSetup::singleton();
	if(!P->SPHParameters.hasDomain)
	    return false;
	InputOutput::ScreenManager *S = InputOutput::ScreenManager::singleton();
	CalcServer *C = CalcServer::singleton();
	cl_int clFlag=0;
	//! 1st.- Send arguments
	clFlag |= sendArgument(clKernel, 5, sizeof(cl_uint), (void*)&(C->n));
	clFlag |= sendArgument(clKernel, 6, sizeof(vec    ), (void*)&(P->SPHParameters.minDomain));
	clFlag |= sendArgument(clKernel, 7, sizeof(vec    ), (void*)&(P->SPHParameters.maxDomain));
	if(clFlag != CL_SUCCESS) {
		S->addMessage(3, "(Domain::execute): Can't send variable to kernel.\n");
	    return true;
	}
	//! 2nd.- Execute the kernel
	#ifdef HAVE_GPUPROFILE
	    cl_event event;
	    cl_ulong end, start;
	    profileTime(0.f);
	    clFlag = clEnqueueNDRangeKernel(C->clComQueue, clKernel, 1, NULL, &clGlobalWorkSize, NULL, 0, NULL, &event);
	#else
	    clFlag = clEnqueueNDRangeKernel(C->clComQueue, clKernel, 1, NULL, &clGlobalWorkSize, NULL, 0, NULL, NULL);
	#endif
	if(clFlag != CL_SUCCESS) {
		S->addMessage(3, "(Domain::execute): Can't execute the kernel.\n");
	    if(clFlag == CL_INVALID_WORK_GROUP_SIZE)
	        S->addMessage(0, "\tInvalid local work group size.\n");
	    else if(clFlag == CL_OUT_OF_RESOURCES)
	        S->addMessage(0, "\tDevice out of resources.\n");
	    else if(clFlag == CL_MEM_OBJECT_ALLOCATION_FAILURE)
	        S->addMessage(0, "\tAllocation error at device.\n");
	    else if(clFlag == CL_OUT_OF_HOST_MEMORY)
	        S->addMessage(0, "\tfailure to allocate resources required by the OpenCL implementation on the host.\n");
	    return true;
	}
	#ifdef HAVE_GPUPROFILE
	    clFlag = clWaitForEvents(1, &event);
	    if(clFlag != CL_SUCCESS) {
	        S->addMessage(3, "(Domain::execute): Can't wait to kernels end.\n");
	        return true;
	    }
	    clFlag |= clGetEventProfilingInfo(event, CL_PROFILING_COMMAND_END, sizeof(cl_ulong), &end, 0);
	    clFlag |= clGetEventProfilingInfo(event, CL_PROFILING_COMMAND_START, sizeof(cl_ulong), &start, 0);
	    if(clFlag != CL_SUCCESS) {
	        S->addMessage(3, "(Domain::execute): Can't profile kernel execution.\n");
	        return true;
	    }
	    profileTime(profileTime() + (end - start)/1000.f);  // 10^-3 ms
	#endif
	return false;
}

bool Domain::setupOpenCL()
{
	InputOutput::ScreenManager *S = InputOutput::ScreenManager::singleton();
	CalcServer *C = CalcServer::singleton();
	cl_int clFlag=0;
	if(!loadKernelFromFile(&clKernel, &clProgram, C->clContext, C->clDevice, mPath, "Domain", ""))
	    return true;
	clFlag |= sendArgument(clKernel, 0, sizeof(cl_mem ), (void*)&(C->imove));
	clFlag |= sendArgument(clKernel, 1, sizeof(cl_mem ), (void*)&(C->posin));
	clFlag |= sendArgument(clKernel, 2, sizeof(cl_mem ), (void*)&(C->vin));
	clFlag |= sendArgument(clKernel, 3, sizeof(cl_mem ), (void*)&(C->fin));
	clFlag |= sendArgument(clKernel, 4, sizeof(cl_mem ), (void*)&(C->mass));
	clFlag |= sendArgument(clKernel, 5, sizeof(cl_uint), (void*)&(C->n));
	if(clFlag)
	    return true;
	//! Test for right work group size
	cl_device_id device;
	size_t localWorkGroupSize=0;
	clFlag |= clGetCommandQueueInfo(C->clComQueue,CL_QUEUE_DEVICE,
	                                sizeof(cl_device_id),&device, NULL);
	if(clFlag != CL_SUCCESS) {
		S->addMessage(3, "(Domain::setupOpenCL): Can't get device from command queue.\n");
	    return true;
	}
	clFlag |= clGetKernelWorkGroupInfo(clKernel,device,CL_KERNEL_WORK_GROUP_SIZE,
	                                   sizeof(size_t), &localWorkGroupSize, NULL);
	if(clFlag != CL_SUCCESS) {
		S->addMessage(3, "(Domain::setupOpenCL): Can't get maximum local work group size.\n");
	    return true;
	}
	if(localWorkGroupSize < clLocalWorkSize)
	    clLocalWorkSize  = localWorkGroupSize;
	clGlobalWorkSize = globalWorkSize(clLocalWorkSize);
	return false;
}

}}
