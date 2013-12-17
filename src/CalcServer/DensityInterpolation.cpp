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
#include <CalcServer/DensityInterpolation.h>

// ----------------------------------------------------------------------------
// Include the calculation server
// ----------------------------------------------------------------------------
#include <CalcServer.h>

namespace Aqua{ namespace CalcServer{

DensityInterpolation::DensityInterpolation()
	: Kernel("DensityInterpolation")
	, mPath(0)
	, clProgram(0)
	, clKernel(0)
	, isLocalMemory(true)
{
	//! 1st.- Get data
	InputOutput::ScreenManager *S = InputOutput::ScreenManager::singleton();
	InputOutput::ProblemSetup *P  = InputOutput::ProblemSetup::singleton();
	if(!P->SPHParameters.DensSteps)  // Density interpolation disabled
	    return;
	unsigned int nChar = strlen(P->OpenCL_kernels.dens_int);
	if(nChar <= 0) {
	    S->addMessage(3, "(DensityInterpolation::DensityInterpolation): Path of kernel is empty.\n");
	    exit(EXIT_FAILURE);
	}
	mPath = new char[nChar+4];
	if(!mPath) {
	    S->addMessage(3, "(DensityInterpolation::DensityInterpolation): Can't allocate memory for path.\n");
	    exit(EXIT_FAILURE);
	}
	strcpy(mPath, P->OpenCL_kernels.dens_int);
	strcat(mPath, ".cl");
	//! 2nd.- Setup the kernel
	clLocalWorkSize  = localWorkSize();
	if(!clLocalWorkSize){
	    S->addMessage(3, "(DensityInterpolation::DensityInterpolation): No valid local work size for required computation.\n");
	    exit(EXIT_FAILURE);
	}
	clGlobalWorkSize = globalWorkSize(clLocalWorkSize);
	if(setupOpenCL()) {
	    exit(EXIT_FAILURE);
	}
	S->addMessage(1, "(DensityInterpolation::DensityInterpolation): DensityInterpolation ready to work!\n");
}

DensityInterpolation::~DensityInterpolation()
{
	if(clKernel)clReleaseKernel(clKernel); clKernel=0;
	if(clProgram)clReleaseProgram(clProgram); clProgram=0;
	if(mPath)delete[] mPath; mPath=0;
}

bool DensityInterpolation::execute()
{
	InputOutput::ProblemSetup *P = InputOutput::ProblemSetup::singleton();
	if(!P->SPHParameters.DensSteps)  // Density interpolation disabled
	    return false;
	InputOutput::ScreenManager *S = InputOutput::ScreenManager::singleton();
	CalcServer *C = CalcServer::singleton();
	cl_int clFlag=0;
	//! Send data to device
	clFlag |= sendArgument(clKernel,  0, sizeof(cl_mem  ), (void*)&(C->dens));
	clFlag |= sendArgument(clKernel,  1, sizeof(cl_mem  ), (void*)&(C->imovein));
	clFlag |= sendArgument(clKernel,  2, sizeof(cl_mem  ), (void*)&(C->posin));
	clFlag |= sendArgument(clKernel,  3, sizeof(cl_mem  ), (void*)&(C->hpin));
	clFlag |= sendArgument(clKernel,  4, sizeof(cl_mem  ), (void*)&(C->massin));
	clFlag |= sendArgument(clKernel,  5, sizeof(cl_mem  ), (void*)&(C->shepard));
	clFlag |= sendArgument(clKernel,  6, sizeof(cl_mem  ), (void*)&(C->lcell));
	clFlag |= sendArgument(clKernel,  7, sizeof(cl_mem  ), (void*)&(C->ihoc));
	clFlag |= sendArgument(clKernel,  8, sizeof(cl_mem  ), (void*)&(C->permutation));
	clFlag |= sendArgument(clKernel,  9, sizeof(cl_mem  ), (void*)&(C->reversePermutation));
	clFlag |= sendArgument(clKernel, 10, sizeof(cl_uint ), (void*)&(C->N));
	clFlag |= sendArgument(clKernel, 11, sizeof(cl_float), (void*)&(C->hfac));
	clFlag |= sendArgument(clKernel, 12, sizeof(uivec   ), (void*)&(C->lvec));
	if(isLocalMemory){
	    clFlag |= sendArgument(clKernel, 13, clLocalWorkSize*sizeof(cl_float), NULL);
	}
	if(clFlag != CL_SUCCESS) {
		S->addMessage(3, "(DensityInterpolation::execute): Can't send arguments to kernel.\n");
	    return true;
	}
	//! Call to execute the kernel
	#ifdef HAVE_GPUPROFILE
	    cl_event event;
	    cl_ulong end, start;
	    profileTime(0.f);
	    clFlag = clEnqueueNDRangeKernel(C->clComQueue, clKernel, 1, NULL, &clGlobalWorkSize, NULL, 0, NULL, &event);
	#else
	    clFlag = clEnqueueNDRangeKernel(C->clComQueue, clKernel, 1, NULL, &clGlobalWorkSize, NULL, 0, NULL, NULL);
	#endif
	if(clFlag != CL_SUCCESS) {
		S->addMessage(3, "(DensityInterpolation::execute): Can't execute the kernel.\n");
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
	//! Profile the kernel execution
	#ifdef HAVE_GPUPROFILE
	    clFlag = clWaitForEvents(1, &event);
	    if(clFlag != CL_SUCCESS) {
	        S->addMessage(3, "(DensityInterpolation::execute): Can't wait to kernels end.\n");
	        return true;
	    }
	    clFlag |= clGetEventProfilingInfo(event, CL_PROFILING_COMMAND_END, sizeof(cl_ulong), &end, 0);
	    clFlag |= clGetEventProfilingInfo(event, CL_PROFILING_COMMAND_START, sizeof(cl_ulong), &start, 0);
	    if(clFlag != CL_SUCCESS) {
	        S->addMessage(3, "(DensityInterpolation::execute): Can't profile kernel execution.\n");
	        return true;
	    }
	    profileTime(profileTime() + (end - start)/1000.f);  // 10^-3 ms
	#endif

	return false;
}

bool DensityInterpolation::setupOpenCL()
{
	InputOutput::ScreenManager *S = InputOutput::ScreenManager::singleton();
	InputOutput::ProblemSetup *P  = InputOutput::ProblemSetup::singleton();
	CalcServer *C = CalcServer::singleton();
	char msg[1024];
	cl_int clFlag;
	cl_device_id device;
	cl_ulong localMem, reqLocalMem;
	clFlag |= clGetCommandQueueInfo(C->clComQueue,CL_QUEUE_DEVICE,
	                                sizeof(cl_device_id),&device, NULL);
	if(clFlag != CL_SUCCESS) {
		S->addMessage(3, "(DensityInterpolation::setupOpenCL): Can't get device from command queue.\n");
	    return true;
	}
	clFlag |= clGetDeviceInfo(device, CL_DEVICE_LOCAL_MEM_SIZE, sizeof(localMem), &localMem, NULL);
	if(clFlag != CL_SUCCESS) {
		S->addMessage(3, "(DensityInterpolation::setupOpenCL): Can't get local memory available on device.\n");
	    return true;
	}
	if(!loadKernelFromFile(&clKernel, &clProgram, C->clContext, C->clDevice, mPath, "DensityInterpolation", ""))
	    return true;
	if(clProgram)clReleaseProgram(clProgram); clProgram=0;
	//! Test if there are enough local memory
	clFlag |= clGetKernelWorkGroupInfo(clKernel,device,CL_KERNEL_LOCAL_MEM_SIZE,
	                                   sizeof(cl_ulong), &reqLocalMem, NULL);
	if(clFlag != CL_SUCCESS) {
		S->addMessage(3, "(DensityInterpolation::setupOpenCL): Can't get kernel memory usage.\n");
	    return true;
	}
	if(localMem < reqLocalMem){
		S->addMessage(3, "(DensityInterpolation::setupOpenCL): Not enough local memory for execution.\n");
	    sprintf(msg, "\tNeeds %lu bytes, but only %lu bytes are available.\n",
	           reqLocalMem, localMem);
	    S->addMessage(0, msg);
	    return true;
	}
	//! Test if local work gorup size must be modified
	size_t localWorkGroupSize=0;
	clFlag |= clGetKernelWorkGroupInfo(clKernel,device,CL_KERNEL_WORK_GROUP_SIZE,
	                                   sizeof(size_t), &localWorkGroupSize, NULL);
	if(clFlag != CL_SUCCESS) {
		S->addMessage(3, "(DensityInterpolation::setupOpenCL): Can't get maximum local work group size.\n");
	    return true;
	}
	if(localWorkGroupSize < clLocalWorkSize)
	    clLocalWorkSize  = localWorkGroupSize;
	//! Look for better local work group size
	clFlag |= clGetKernelWorkGroupInfo(clKernel,device,CL_KERNEL_PREFERRED_WORK_GROUP_SIZE_MULTIPLE,
	                                   sizeof(size_t), &localWorkGroupSize, NULL);
	if(clFlag != CL_SUCCESS) {
		S->addMessage(3, "(DensityInterpolation::setupOpenCL): Can't get preferred local work group size.\n");
	    return true;
	}
	clLocalWorkSize  = (clLocalWorkSize/localWorkGroupSize) * localWorkGroupSize;
	clGlobalWorkSize = globalWorkSize(clLocalWorkSize);
	//! Test if computation can be accelerated with local memory
	reqLocalMem += clLocalWorkSize*(sizeof(cl_float));
	if(localMem < reqLocalMem){
		S->addMessage(2, "(DensityInterpolation::setupOpenCL): Not enough local memory.\n");
	    sprintf(msg, "\tNeeds %lu bytes, but only %lu bytes are available.\n",
	           reqLocalMem, localMem);
	    S->addMessage(0, msg);
	    S->addMessage(0, "\tLocal memory usage will be avoided therefore.\n");
	    isLocalMemory = false;
	    char options[19]; strcpy(options,"-D__NO_LOCAL_MEM__");
	    if(clKernel)clReleaseKernel(clKernel); clKernel=0;
	    if(!loadKernelFromFile(&clKernel, &clProgram, C->clContext, C->clDevice, mPath, "DensityInterpolation", options))
	        return true;
	    if(clProgram)clReleaseProgram(clProgram); clProgram=0;
	}
	return false;
}

}}  // namespace
