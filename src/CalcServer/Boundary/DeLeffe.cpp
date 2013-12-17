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
#include <CalcServer/Boundary/DeLeffe.h>

// ----------------------------------------------------------------------------
// Include the calculation server
// ----------------------------------------------------------------------------
#include <CalcServer.h>

namespace Aqua{ namespace CalcServer{ namespace Boundary{

DeLeffe::DeLeffe()
	: Kernel("DeLeffe")
	, mPath(0)
	, clProgram(0)
	, clVerticesKernel(0)
	, clBoundaryKernel(0)
	, isLocalMemory(true)
{
	InputOutput::ScreenManager *S = InputOutput::ScreenManager::singleton();
	InputOutput::ProblemSetup *P = InputOutput::ProblemSetup::singleton();
	if(P->SPHParameters.Boundary!=2)  // DeLeffe condition has not been selected
	    return;
	//! 1st.- Get data
	int nChar = strlen(P->OpenCL_kernels.de_Leffe);
	if(nChar <= 0) {
	    S->addMessage(3, "(DeLeffe::DeLeffe): mPath of DeLeffe kernel is empty.\n");
	    exit(EXIT_FAILURE);
	}
	mPath = new char[nChar+4];
	if(!mPath) {
	    S->addMessage(3, "(DeLeffe::DeLeffe): Can't allocate memory for path.\n");
	    exit(EXIT_FAILURE);
	}
	strcpy(mPath, P->OpenCL_kernels.de_Leffe);
	strcat(mPath, ".cl");
	//! 2nd.- Setup the kernel
	clLocalWorkSize  = localWorkSize();
	if(!clLocalWorkSize){
	    S->addMessage(3, "(DeLeffe::DeLeffe): No valid local work size for required computation.\n");
	    exit(EXIT_FAILURE);
	}
	clGlobalWorkSize = globalWorkSize(clLocalWorkSize);
	if(setupOpenCL()) {
	    exit(EXIT_FAILURE);
	}
	S->addMessage(1, "(DeLeffe::DeLeffe): DeLeffe ready to work!\n");
}

DeLeffe::~DeLeffe()
{
	if(clBoundaryKernel)clReleaseKernel(clBoundaryKernel); clBoundaryKernel=0;
	if(clVerticesKernel)clReleaseKernel(clVerticesKernel); clVerticesKernel=0;
	if(clProgram)clReleaseProgram(clProgram); clProgram=0;
	if(mPath) delete[] mPath; mPath=0;
}

bool DeLeffe::execute()
{
	InputOutput::ProblemSetup *P = InputOutput::ProblemSetup::singleton();
	if(P->SPHParameters.Boundary!=2)  // DeLeffe condition has not been selected
	    return false;
	if(vertices())
	    return true;
	if(boundary())
	    return true;
	return false;
}

bool DeLeffe::vertices()
{
	InputOutput::ScreenManager *S = InputOutput::ScreenManager::singleton();
	CalcServer *C = CalcServer::singleton();
	cl_int clFlag=0;
	//! Send all variables to the server
	clFlag |= sendArgument(clVerticesKernel,  0, sizeof(cl_mem  ), (void*)&(C->imove));
	clFlag |= sendArgument(clVerticesKernel,  1, sizeof(cl_mem  ), (void*)&(C->f));
	clFlag |= sendArgument(clVerticesKernel,  2, sizeof(cl_mem  ), (void*)&(C->drdt));
	clFlag |= sendArgument(clVerticesKernel,  3, sizeof(cl_mem  ), (void*)&(C->press));
	clFlag |= sendArgument(clVerticesKernel,  4, sizeof(cl_mem  ), (void*)&(C->pressin));
	clFlag |= sendArgument(clVerticesKernel,  5, sizeof(cl_mem  ), (void*)&(C->dens));
	clFlag |= sendArgument(clVerticesKernel,  6, sizeof(cl_mem  ), (void*)&(C->densin));
	clFlag |= sendArgument(clVerticesKernel,  7, sizeof(cl_mem  ), (void*)&(C->refd));
	clFlag |= sendArgument(clVerticesKernel,  8, sizeof(cl_mem  ), (void*)&(C->ifluid));
	clFlag |= sendArgument(clVerticesKernel,  9, sizeof(cl_mem  ), (void*)&(C->gamma));
	clFlag |= sendArgument(clVerticesKernel, 10, sizeof(cl_mem  ), (void*)&(C->normal));
	clFlag |= sendArgument(clVerticesKernel, 11, sizeof(cl_mem  ), (void*)&(C->normalin));
	clFlag |= sendArgument(clVerticesKernel, 12, sizeof(cl_mem  ), (void*)&(C->shepard));
	clFlag |= sendArgument(clVerticesKernel, 13, sizeof(cl_mem  ), (void*)&(C->reversePermutation));
	clFlag |= sendArgument(clVerticesKernel, 14, sizeof(cl_uint ), (void*)&(C->n));
	clFlag |= sendArgument(clVerticesKernel, 15, sizeof(cl_float), (void*)&(C->cs));
	if(clFlag != CL_SUCCESS) {
		S->addMessage(3, "(DeLeffe::vertices): Can't send arguments to Vertex set kernel.\n");
	    return true;
	}
	//! Execute the kernel
	#ifdef HAVE_GPUPROFILE
	    cl_event event;
	    cl_ulong end, start;
	    clFlag = clEnqueueNDRangeKernel(C->clComQueue, clVerticesKernel, 1, NULL, &clGlobalWorkSize, NULL, 0, NULL, &event);
	#else
	    clFlag = clEnqueueNDRangeKernel(C->clComQueue, clVerticesKernel, 1, NULL, &clGlobalWorkSize, NULL, 0, NULL, NULL);
	#endif
	if(clFlag != CL_SUCCESS) {
		S->addMessage(3, "(DeLeffe::vertices): Can't execute the kernel.\n");
	    if(clFlag == CL_INVALID_KERNEL_ARGS)
	        S->addMessage(0, "\tInvalid kernel arguments.\n");
	    else if(clFlag == CL_INVALID_WORK_GROUP_SIZE)
	        S->addMessage(0, "\tInvalid local work group size.\n");
	    else if(clFlag == CL_INVALID_WORK_ITEM_SIZE)
	        S->addMessage(0, "\tInvalid local work group size (greather than maximum allowed value).\n");
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
	        S->addMessage(3, "(DeLeffe::vertices): Can't wait to kernels end.\n");
	        return true;
	    }
	    clFlag |= clGetEventProfilingInfo(event, CL_PROFILING_COMMAND_END, sizeof(cl_ulong), &end, 0);
	    clFlag |= clGetEventProfilingInfo(event, CL_PROFILING_COMMAND_START, sizeof(cl_ulong), &start, 0);
	    if(clFlag != CL_SUCCESS) {
	        S->addMessage(3, "(DeLeffe::vertices): Can't profile kernel execution.\n");
	        return true;
	    }
	    profileTime(profileTime() + (end - start)/1000.f);  // 10^-3 ms
	#endif

	return false;
}

bool DeLeffe::boundary()
{
	InputOutput::ScreenManager *S = InputOutput::ScreenManager::singleton();
	CalcServer *C = CalcServer::singleton();
	cl_int clFlag=0;
	//! Send all variables to the server
	clFlag |= sendArgument(clBoundaryKernel,  0, sizeof(cl_mem  ), (void*)&(C->ifluidin));
	clFlag |= sendArgument(clBoundaryKernel,  1, sizeof(cl_mem  ), (void*)&(C->imovein));
	clFlag |= sendArgument(clBoundaryKernel,  2, sizeof(cl_mem  ), (void*)&(C->posin));
	clFlag |= sendArgument(clBoundaryKernel,  3, sizeof(cl_mem  ), (void*)&(C->normalin));
	clFlag |= sendArgument(clBoundaryKernel,  4, sizeof(cl_mem  ), (void*)&(C->vin));
	clFlag |= sendArgument(clBoundaryKernel,  5, sizeof(cl_mem  ), (void*)&(C->densin));
	clFlag |= sendArgument(clBoundaryKernel,  6, sizeof(cl_mem  ), (void*)&(C->hpin));
	clFlag |= sendArgument(clBoundaryKernel,  7, sizeof(cl_mem  ), (void*)&(C->pressin));
	clFlag |= sendArgument(clBoundaryKernel,  8, sizeof(cl_mem  ), (void*)&(C->massin));
	clFlag |= sendArgument(clBoundaryKernel,  9, sizeof(cl_mem  ), (void*)&(C->Viscdyn));
	clFlag |= sendArgument(clBoundaryKernel, 10, sizeof(cl_mem  ), (void*)&(C->f));
	clFlag |= sendArgument(clBoundaryKernel, 11, sizeof(cl_mem  ), (void*)&(C->drdt));
	clFlag |= sendArgument(clBoundaryKernel, 12, sizeof(cl_mem  ), (void*)&(C->gradShepard));
	clFlag |= sendArgument(clBoundaryKernel, 13, sizeof(cl_mem  ), (void*)&(C->lcell));
	clFlag |= sendArgument(clBoundaryKernel, 14, sizeof(cl_mem  ), (void*)&(C->ihoc));
	clFlag |= sendArgument(clBoundaryKernel, 15, sizeof(cl_mem  ), (void*)&(C->permutation));
	clFlag |= sendArgument(clBoundaryKernel, 16, sizeof(cl_mem  ), (void*)&(C->reversePermutation));
	clFlag |= sendArgument(clBoundaryKernel, 17, sizeof(cl_uint ), (void*)&(C->N));
	clFlag |= sendArgument(clBoundaryKernel, 18, sizeof(cl_float), (void*)&(C->hfac));
	clFlag |= sendArgument(clBoundaryKernel, 19, sizeof(uivec   ), (void*)&(C->lvec));
	if(isLocalMemory){
	    clFlag |= sendArgument(clBoundaryKernel, 20, clLocalWorkSize*sizeof(vec     ), NULL);
	    clFlag |= sendArgument(clBoundaryKernel, 21, clLocalWorkSize*sizeof(cl_float), NULL);
	    clFlag |= sendArgument(clBoundaryKernel, 22, clLocalWorkSize*sizeof(vec     ), NULL);
	}
	if(clFlag != CL_SUCCESS) {
		S->addMessage(3, "(DeLeffe::boundary): Can't send arguments to boundary computation kernel.\n");
	    return true;
	}
	//! Execute the kernel
	#ifdef HAVE_GPUPROFILE
	    cl_event event;
	    cl_ulong end, start;
	    profileTime(0.f);
	    clFlag = clEnqueueNDRangeKernel(C->clComQueue, clBoundaryKernel, 1, NULL, &clGlobalWorkSize, &clLocalWorkSize, 0, NULL, &event);
	#else
	    clFlag = clEnqueueNDRangeKernel(C->clComQueue, clBoundaryKernel, 1, NULL, &clGlobalWorkSize, &clLocalWorkSize, 0, NULL, NULL);
	#endif
	if(clFlag != CL_SUCCESS) {
		S->addMessage(3, "(DeLeffe::boundary): Can't execute the kernel.\n");
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
	        S->addMessage(3, "(DeLeffe::boundary): Can't wait to kernels end.\n");
	        return true;
	    }
	    clFlag |= clGetEventProfilingInfo(event, CL_PROFILING_COMMAND_END, sizeof(cl_ulong), &end, 0);
	    clFlag |= clGetEventProfilingInfo(event, CL_PROFILING_COMMAND_START, sizeof(cl_ulong), &start, 0);
	    if(clFlag != CL_SUCCESS) {
	        S->addMessage(3, "(DeLeffe::boundary): Can't profile kernel execution.\n");
	        return true;
	    }
	    profileTime(profileTime() + (end - start)/1000.f);  // 10^-3 ms
	#endif
	return false;
}

bool DeLeffe::setupOpenCL()
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
		S->addMessage(3, "(DeLeffe::setupOpenCL): Can't get device from command queue.\n");
	    return true;
	}
	clFlag |= clGetDeviceInfo(device, CL_DEVICE_LOCAL_MEM_SIZE, sizeof(localMem), &localMem, NULL);
	if(clFlag != CL_SUCCESS) {
		S->addMessage(3, "(DeLeffe::setupOpenCL): Can't get local memory available on device.\n");
	    return true;
	}
	if(!loadKernelFromFile(&clVerticesKernel, &clProgram, C->clContext, C->clDevice, mPath, "Vertices", ""))
	    return true;
	if(clProgram)clReleaseProgram(clProgram); clProgram=0;
	if(!loadKernelFromFile(&clBoundaryKernel, &clProgram, C->clContext, C->clDevice, mPath, "Boundary", ""))
	    return true;
	if(clProgram)clReleaseProgram(clProgram); clProgram=0;
	//! Test if there are enough local memory
	clFlag |= clGetKernelWorkGroupInfo(clVerticesKernel,device,CL_KERNEL_LOCAL_MEM_SIZE,
	                                   sizeof(cl_ulong), &reqLocalMem, NULL);
	if(clFlag != CL_SUCCESS) {
		S->addMessage(3, "(DeLeffe::setupOpenCL): Can't get vertices kernel memory usage.\n");
	    return true;
	}
	if(localMem < reqLocalMem){
		S->addMessage(3, "(DeLeffe::setupOpenCL): Not enough local memory for vertices execution.\n");
	    sprintf(msg, "\tNeeds %lu bytes, but only %lu bytes are available.\n",
	           reqLocalMem, localMem);
	    S->addMessage(0, msg);
	    return true;
	}
	clFlag |= clGetKernelWorkGroupInfo(clBoundaryKernel,device,CL_KERNEL_LOCAL_MEM_SIZE,
	                                   sizeof(cl_ulong), &reqLocalMem, NULL);
	if(clFlag != CL_SUCCESS) {
		S->addMessage(3, "(DeLeffe::setupOpenCL): Can't get boundary kernel memory usage.\n");
	    return true;
	}
	if(localMem < reqLocalMem){
		S->addMessage(3, "(DeLeffe::setupOpenCL): Not enough local memory for boundary execution.\n");
	    sprintf(msg, "\tNeeds %lu bytes, but only %lu bytes are available.\n",
	           reqLocalMem, localMem);
	    S->addMessage(0, msg);
	    return true;
	}
	//! Test if local work gorup size must be modified
	size_t localWorkGroupSize=0;
	clFlag |= clGetKernelWorkGroupInfo(clVerticesKernel,device,CL_KERNEL_WORK_GROUP_SIZE,
	                                   sizeof(size_t), &localWorkGroupSize, NULL);
	if(clFlag != CL_SUCCESS) {
		S->addMessage(3, "(DeLeffe::setupOpenCL): Can't get vertices maximum local work group size.\n");
	    return true;
	}
	if(localWorkGroupSize < clLocalWorkSize)
	    clLocalWorkSize  = localWorkGroupSize;
	clFlag |= clGetKernelWorkGroupInfo(clBoundaryKernel,device,CL_KERNEL_WORK_GROUP_SIZE,
	                                   sizeof(size_t), &localWorkGroupSize, NULL);
	if(clFlag != CL_SUCCESS) {
		S->addMessage(3, "(DeLeffe::setupOpenCL): Can't get boundary maximum local work group size.\n");
	    return true;
	}
	if(localWorkGroupSize < clLocalWorkSize)
	    clLocalWorkSize  = localWorkGroupSize;
	//! Look for better local work group size
	clFlag |= clGetKernelWorkGroupInfo(clBoundaryKernel,device,CL_KERNEL_PREFERRED_WORK_GROUP_SIZE_MULTIPLE,
	                                   sizeof(size_t), &localWorkGroupSize, NULL);
	if(clFlag != CL_SUCCESS) {
		S->addMessage(3, "(DeLeffe::setupOpenCL): Can't get boundary preferred local work group size.\n");
	    return true;
	}
	clLocalWorkSize  = (clLocalWorkSize/localWorkGroupSize) * localWorkGroupSize;
	clGlobalWorkSize = globalWorkSize(clLocalWorkSize);
	//! Test if computation can be accelerated with local memory
	reqLocalMem += clLocalWorkSize*(  sizeof(vec     )
	                                + sizeof(cl_float)
	                                + sizeof(vec     ));
	if(localMem < reqLocalMem){
		S->addMessage(2, "(DeLeffe::setupOpenCL): Not enough local memory for boundary.\n");
	    sprintf(msg, "\tNeeds %lu bytes, but only %lu bytes are available.\n",
	           reqLocalMem, localMem);
	    S->addMessage(0, msg);
	    S->addMessage(0, "\tLocal memory usage will be avoided therefore.\n");
	    isLocalMemory = false;
	    char options[19]; strcpy(options,"-D__NO_LOCAL_MEM__");
	    if(clBoundaryKernel)clReleaseKernel(clBoundaryKernel); clBoundaryKernel=0;
	    if(!loadKernelFromFile(&clBoundaryKernel, &clProgram, C->clContext, C->clDevice, mPath, "Boundary", options))
	        return true;
	    if(clProgram)clReleaseProgram(clProgram); clProgram=0;
	}
	return false;
}

}}}  // namespaces
