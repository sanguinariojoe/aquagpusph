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
#include <CalcServer/LinkList.h>

// ----------------------------------------------------------------------------
// Include the calculation server
// ----------------------------------------------------------------------------
#include <CalcServer.h>

namespace Aqua{ namespace CalcServer{

LinkList::LinkList()
	: Kernel("LinkList")
	, mPath(0)
	, clProgram(0)
	, clLcellKernel(0)
	, clIhocKernel(0)
	, clLLKernel(0)
	, mRadixSort(0)

{
	InputOutput::ProblemSetup *P  = InputOutput::ProblemSetup::singleton();
	InputOutput::ScreenManager *S = InputOutput::ScreenManager::singleton();
	CalcServer *C = CalcServer::singleton();
	int nChar;
	//! 1st.- Get data
	nChar = strlen(P->OpenCL_kernels.link_list);
	if(nChar <= 0) {
	    printf("ERROR (LinkList::Init): Path of LinkList kernel is empty.\n");
	    exit(EXIT_FAILURE);
	}
	mPath = new char[nChar+4];
	if(!mPath) {
	    printf("ERROR (LinkList::Init): Can't allocate memory for path.\n");
	    exit(EXIT_FAILURE);
	}
	strcpy(mPath, P->OpenCL_kernels.link_list);
	strcat(mPath, ".cl");
	//! 2nd.- Setup the kernels
	clLocalWorkSize  = localWorkSize();
	if(!clLocalWorkSize){
	    printf("ERROR (LinkList::Init): No valid local work size for required computation.\n");
	    exit(EXIT_FAILURE);
	}
	clGlobalWorkSize = globalWorkSize(clLocalWorkSize);
	if(setupOpenCL()) {
	    exit(EXIT_FAILURE);
	}
	//! 3rd.- Built radix sort
	mRadixSort = new RadixSort();
	printf("\tINFO (LinkList::Init): LinkList ready to work!\n");
}

LinkList::~LinkList()
{
	InputOutput::ScreenManager *S = InputOutput::ScreenManager::singleton();
	if(clLcellKernel)clReleaseKernel(clLcellKernel); clLcellKernel=0;
	if(clIhocKernel)clReleaseKernel(clIhocKernel); clIhocKernel=0;
	if(clLLKernel)clReleaseKernel(clLLKernel); clLLKernel=0;
	if(clProgram)clReleaseProgram(clProgram); clProgram=0;
	if(mPath)delete[] mPath; mPath=0;
	S->addMessage(1, "(LinkList::~LinkList): Destroying radix sort processor...\n");
	if(mRadixSort)delete mRadixSort; mRadixSort=0;
}

bool LinkList::execute()
{
	InputOutput::ScreenManager *S = InputOutput::ScreenManager::singleton();
	CalcServer *C = CalcServer::singleton();
	cl_int clFlag=0;
	//! Allocate memory if needed
	if(allocLinkList())
	     return true;
	//! Set initial values for cells
	clFlag |= sendArgument(clIhocKernel, 0, sizeof(cl_mem), (void*)&(C->ihoc));
	clFlag |= sendArgument(clIhocKernel, 1, sizeof(cl_mem), (void*)&(C->isValidCell));
	clFlag |= sendArgument(clIhocKernel, 2, sizeof(cl_uint), (void*)&(C->lxydim));
	clFlag |= sendArgument(clIhocKernel, 3, sizeof(cl_uint), (void*)&(C->N));
	if(clFlag != CL_SUCCESS) {
		S->addMessage(3, (char *)"(LinkList::Execute): Can't send arguments to to ihoc clearer kernel.\n");
	    return true;
	}
	clLocalWorkSize  = localWorkSize(C->lxy);
	if(!clLocalWorkSize){
	    S->addMessage(3, (char*)"(LinkList::Execute): No valid local work size for ihoc clearer.\n");
	    return true;
	}
	clGlobalWorkSize = globalWorkSize(clLocalWorkSize, C->lxy);
	#ifdef HAVE_GPUPROFILE
	    cl_event event;
	    cl_ulong end, start;
	    profileTime(0.f);
	    clFlag = clEnqueueNDRangeKernel(C->clComQueue, clIhocKernel, 1, NULL, &clGlobalWorkSize, NULL, 0, NULL, &event);
	#else
	    clFlag = clEnqueueNDRangeKernel(C->clComQueue, clIhocKernel, 1, NULL, &clGlobalWorkSize, NULL, 0, NULL, NULL);
	#endif
	if(clFlag != CL_SUCCESS) {
		S->addMessage(3, "(LinkList::Execute): Can't execute the kernel \"InitIhoc\".\n");
	    if(clFlag == CL_INVALID_PROGRAM_EXECUTABLE)
	        S->addMessage(0, "\tInvalid program (Compile errors maybe?).\n");
	    else if(clFlag == CL_INVALID_COMMAND_QUEUE)
	        S->addMessage(0, "\tInvalid command queue.\n");
	    else if(clFlag == CL_INVALID_KERNEL)
	        S->addMessage(0, "\tKernel is not a valid object.\n");
	    else if(clFlag == CL_INVALID_CONTEXT)
	        S->addMessage(0, "\tContext associated to command queue don't match qith kernel context.\n");
	    else if(clFlag == CL_INVALID_KERNEL_ARGS)
	        S->addMessage(0, "\tOne or more arguments are invalid (maybe don't specified).\n");
	    else if(clFlag == CL_INVALID_WORK_DIMENSION)
	        S->addMessage(0, "\tDimension must be a value between 1 and 3.\n");
	    else if(clFlag == CL_INVALID_WORK_GROUP_SIZE)
	        S->addMessage(0, "\tInvalid local work group size.\n");
	    else if(clFlag == CL_INVALID_WORK_ITEM_SIZE)
	        S->addMessage(0, "\tLocal work group size is out of bounds.\n");
	    else if(clFlag == CL_INVALID_GLOBAL_OFFSET)
	        S->addMessage(0, "\tGlobal offset must be NULL.\n");
	    else if(clFlag == CL_OUT_OF_RESOURCES)
	        S->addMessage(0, "\tDevice out of resources.\n");
	    else if(clFlag == CL_MEM_OBJECT_ALLOCATION_FAILURE)
	        S->addMessage(0, "\tAllocation error at device.\n");
	    else if(clFlag == CL_INVALID_EVENT_WAIT_LIST)
	        S->addMessage(0, "\tInvalid event wait instruction.\n");
	    else if(clFlag == CL_OUT_OF_HOST_MEMORY)
	        S->addMessage(0, "\tfailure to allocate resources required by the OpenCL implementation on the host.\n");
	    return true;
	}
	#ifdef HAVE_GPUPROFILE
	    clFlag = clWaitForEvents(1, &event);
	    if(clFlag != CL_SUCCESS) {
	        S->addMessage(3, "(LinkList::Execute): Can't wait to kernels end.\n");
	        return true;
	    }
	    clFlag |= clGetEventProfilingInfo(event, CL_PROFILING_COMMAND_END, sizeof(cl_ulong), &end, 0);
	    clFlag |= clGetEventProfilingInfo(event, CL_PROFILING_COMMAND_START, sizeof(cl_ulong), &start, 0);
	    if(clFlag != CL_SUCCESS) {
	        S->addMessage(3, "(LinkList::Execute): Can't profile kernel execution.\n");
	        return true;
	    }
	    profileTime(profileTime() + (end - start)/1000.f);  // 10^-3 ms
	#endif
	//! Find each particle cell
	clFlag |= sendArgument(clLcellKernel, 0, sizeof(cl_mem  ), (void*)&C->lcell);
	clFlag |= sendArgument(clLcellKernel, 1, sizeof(cl_mem  ), (void*)&C->isValidCell);
	clFlag |= sendArgument(clLcellKernel, 2, sizeof(cl_mem  ), (void*)&C->imove);
	clFlag |= sendArgument(clLcellKernel, 3, sizeof(cl_mem  ), (void*)&C->pos);
	clFlag |= sendArgument(clLcellKernel, 4, sizeof(cl_uint ), (void*)&C->N);
	clFlag |= sendArgument(clLcellKernel, 5, sizeof(cl_uint ), (void*)&C->nLcell);
	clFlag |= sendArgument(clLcellKernel, 6, sizeof(vec     ), (void*)&C->posmin);
	clFlag |= sendArgument(clLcellKernel, 7, sizeof(cl_float), (void*)&C->rdist);
	clFlag |= sendArgument(clLcellKernel, 8, sizeof(cl_uint ), (void*)&C->lxy);
	clFlag |= sendArgument(clLcellKernel, 9, sizeof(uivec   ), (void*)&C->lvec);
	if(clFlag != CL_SUCCESS) {
		S->addMessage(3, "(LinkList::Execute): Can't send arguments to cell allocation kernel.\n");
	    return true;
	}
	clLocalWorkSize  = localWorkSize(C->nLcell);
	if(!clLocalWorkSize){
	    S->addMessage(3, "(LinkList::Execute): No valid local work size for cell allocation.\n");
	    return true;
	}
	clGlobalWorkSize = globalWorkSize(clLocalWorkSize, C->nLcell);
	#ifdef HAVE_GPUPROFILE
	    clFlag = clEnqueueNDRangeKernel(C->clComQueue, clLcellKernel, 1, NULL, &clGlobalWorkSize, &clLocalWorkSize, 0, NULL, &event);
	#else
	    clFlag = clEnqueueNDRangeKernel(C->clComQueue, clLcellKernel, 1, NULL, &clGlobalWorkSize, &clLocalWorkSize, 0, NULL, NULL);
	#endif
	if(clFlag != CL_SUCCESS) {
		S->addMessage(3, "(LinkList::Execute): Can't execute the kernel \"LCell\".\n");
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
	        S->addMessage(3, "(LinkList::Execute): Can't wait to kernels end.\n");
	        return true;
	    }
	    clFlag |= clGetEventProfilingInfo(event, CL_PROFILING_COMMAND_END, sizeof(cl_ulong), &end, 0);
	    clFlag |= clGetEventProfilingInfo(event, CL_PROFILING_COMMAND_START, sizeof(cl_ulong), &start, 0);
	    if(clFlag != CL_SUCCESS) {
	        S->addMessage(3, "(LinkList::Execute): Can't profile kernel execution.\n");
	        return true;
	    }
	    profileTime(profileTime() + (end - start)/1000.f);  // 10^-3 ms
	#endif
	//! Sort particle
	if(mRadixSort->sort())
	    return true;
	#ifdef HAVE_GPUPROFILE
	    profileTime(profileTime() + mRadixSort->profileTime());  // 10^-3 ms
	#endif
	//! Set heads of chain
	unsigned int nMinOne = C->N-1;
	clFlag |= sendArgument(clLLKernel, 0, sizeof(cl_mem), (void*)&(C->lcell));
	clFlag |= sendArgument(clLLKernel, 1, sizeof(cl_mem), (void*)&(C->ihoc));
	clFlag |= sendArgument(clLLKernel, 2, sizeof(cl_uint), (void*)&nMinOne);
	if(clFlag != CL_SUCCESS) {
		S->addMessage(3, "(LinkList::Execute): Can't send arguments to ihoc perform kernel.\n");
	    return true;
	}
	clLocalWorkSize  = localWorkSize(nMinOne);
	if(!clLocalWorkSize){
	    S->addMessage(3, "(LinkList::Execute): No valid local work size for ihoc computation.\n");
	    return true;
	}
	clGlobalWorkSize = globalWorkSize(clLocalWorkSize, nMinOne);
	#ifdef HAVE_GPUPROFILE
	    clFlag = clEnqueueNDRangeKernel(C->clComQueue, clLLKernel, 1, NULL, &clGlobalWorkSize, NULL, 0, NULL, &event);
	#else
	    clFlag = clEnqueueNDRangeKernel(C->clComQueue, clLLKernel, 1, NULL, &clGlobalWorkSize, NULL, 0, NULL, NULL);
	#endif
	if(clFlag != CL_SUCCESS) {
		S->addMessage(3, "(LinkList::Execute): Can't execute the kernel \"LinkList\".\n");
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
	        S->addMessage(3, "(LinkList::Execute): Can't wait to kernels end.\n");
	        return true;
	    }
	    clFlag |= clGetEventProfilingInfo(event, CL_PROFILING_COMMAND_END, sizeof(cl_ulong), &end, 0);
	    clFlag |= clGetEventProfilingInfo(event, CL_PROFILING_COMMAND_START, sizeof(cl_ulong), &start, 0);
	    if(clFlag != CL_SUCCESS) {
	        S->addMessage(3, "(LinkList::Execute): Can't profile kernel execution.\n");
	        return true;
	    }
	    profileTime(profileTime() + (end - start)/1000.f);  // 10^-3 ms
	#endif
	return false;
}

bool LinkList::setupOpenCL()
{
	InputOutput::ScreenManager *S = InputOutput::ScreenManager::singleton();
	CalcServer *C = CalcServer::singleton();
	int clFlag=0;
	if(!loadKernelFromFile(&clLcellKernel, &clProgram, C->clContext, C->clDevice, mPath, "LCell", ""))
	    return true;
	if(clProgram)clReleaseProgram(clProgram); clProgram=0;
	if(!loadKernelFromFile(&clIhocKernel, &clProgram, C->clContext, C->clDevice, mPath, "InitIhoc", ""))
	    return true;
	if(clProgram)clReleaseProgram(clProgram); clProgram=0;
	if(!loadKernelFromFile(&clLLKernel, &clProgram, C->clContext, C->clDevice, mPath, "LinkList", ""))
	    return true;
	if(clProgram)clReleaseProgram(clProgram); clProgram=0;
	clFlag |= sendArgument(clLcellKernel,  0, sizeof(cl_mem), (void*)&(C->pos));
	if(clFlag)
	    return true;
	//! Test for right work group size
	cl_device_id device;
	size_t localWorkGroupSize=0;
	clFlag |= clGetCommandQueueInfo(C->clComQueue,CL_QUEUE_DEVICE,
	                                sizeof(cl_device_id),&device, NULL);
	if(clFlag != CL_SUCCESS) {
		S->addMessage(3, "(LinkList::setupOpenCL): Can't get device from command queue.\n");
	    return true;
	}
	clFlag |= clGetKernelWorkGroupInfo(clLcellKernel,device,CL_KERNEL_WORK_GROUP_SIZE,
	                                   sizeof(size_t), &localWorkGroupSize, NULL);
	if(clFlag != CL_SUCCESS) {
		S->addMessage(3, "(LinkList::setupOpenCL): Can't get Lcell maximum local work group size.\n");
	    return true;
	}
	if(localWorkGroupSize < clLocalWorkSize)
	    clLocalWorkSize  = localWorkGroupSize;
	clFlag |= clGetKernelWorkGroupInfo(clIhocKernel,device,CL_KERNEL_WORK_GROUP_SIZE,
	                                   sizeof(size_t), &localWorkGroupSize, NULL);
	if(clFlag != CL_SUCCESS) {
		S->addMessage(3, "(LinkList::setupOpenCL): Can't get ihoc maximum local work group size.\n");
	    return true;
	}
	if(localWorkGroupSize < clLocalWorkSize)
	    clLocalWorkSize  = localWorkGroupSize;
	clFlag |= clGetKernelWorkGroupInfo(clLLKernel,device,CL_KERNEL_WORK_GROUP_SIZE,
	                                   sizeof(size_t), &localWorkGroupSize, NULL);
	if(clFlag != CL_SUCCESS) {
		S->addMessage(3, "(LinkList::setupOpenCL): Can't get ll maximum local work group size.\n");
	    return true;
	}
	if(localWorkGroupSize < clLocalWorkSize)
	    clLocalWorkSize  = localWorkGroupSize;
	clGlobalWorkSize = globalWorkSize(clLocalWorkSize);
	return false;
}

bool LinkList::allocLinkList()
{
	InputOutput::ScreenManager *S = InputOutput::ScreenManager::singleton();
	CalcServer *C = CalcServer::singleton();
	char Log[256];
	if(C->lxy > C->lxydim)
	{
	    sprintf(Log, "(LinkList::sphAllocLinkList): Number of cells increased [%d -> %d]\n",C->lxydim,C->lxy);
		S->addMessage(2, Log);
	    if(C->lxydim > 0) {
	        if(C->ihoc)clReleaseMemObject(C->ihoc); C->ihoc=0;
	        if(C->isValidCell)clReleaseMemObject(C->isValidCell); C->isValidCell=0;
			C->AllocatedMem -= C->lxydim * sizeof( cl_uint );
			C->AllocatedMem -= C->lxydim * sizeof( cl_short );
	    }
		if(C->allocMemory(&(C->ihoc), C->lxy * sizeof( cl_uint ))) {
	        sprintf(Log, "(LinkList::sphAllocLinkList): Fail allocating memory for ihoc (%u bytes).\n", (unsigned int)(C->lxy * sizeof( cl_uint )) );
		    S->addMessage(3, Log);
		    return true;
		}
		if(C->allocMemory(&(C->isValidCell), C->lxy * sizeof( cl_short ))) {
	        sprintf(Log, "(LinkList::sphAllocLinkList): Fail allocating memory for isValidCell (%u bytes).\n", (unsigned int)(C->lxy * sizeof( cl_short )) );
		    S->addMessage(3, Log);
		    return true;
		}
		sprintf(Log, "\tAllocated memory = %u bytes\n", (unsigned int)C->AllocatedMem);
		S->addMessage(1, Log);
	    C->lxydim = C->lxy;
	    clLocalWorkSize = 256;
	    clGlobalWorkSize = globalWorkSize(clLocalWorkSize);
	}
	return false;
}

}}  // namespace
