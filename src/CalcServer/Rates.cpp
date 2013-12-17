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
#include <CalcServer/Rates.h>

// ----------------------------------------------------------------------------
// Include the calculation server
// ----------------------------------------------------------------------------
#include <CalcServer.h>

namespace Aqua{ namespace CalcServer{

Rates::Rates()
	: Kernel("Rates")
	, mPath(0)
	, clProgram(0)
	, clKernel(0)
	, clSortKernel(0)
	, clGlobalWorkSize(0)
	, clLocalWorkSize(0)
	, isDelta(false)
	, isLocalMemory(true)
{
	InputOutput::ScreenManager *S = InputOutput::ScreenManager::singleton();
	InputOutput::ProblemSetup *P  = InputOutput::ProblemSetup::singleton();
	unsigned int i, nChar = strlen(P->OpenCL_kernels.rates);
	if(nChar <= 0) {
	    S->addMessage(3, "(Rates::Rates): Path of rates kernel is empty.\n");
	    exit(EXIT_FAILURE);
	}
	mPath = new char[nChar+4];
	if(!mPath) {
	    S->addMessage(3, "(Rates::Rates): Can't allocate memory for path.\n");
	    exit(EXIT_FAILURE);
	}
	strcpy(mPath, P->OpenCL_kernels.rates);
	strcat(mPath, ".cl");
    for(i=0;i<P->n_fluids;i++){
        if(P->fluids[i].delta > 0.f){
            isDelta = true;
            break;
        }
    }

	clLocalWorkSize  = localWorkSize();
	if(!clLocalWorkSize){
	    S->addMessage(3, "(Rates::Rates): No valid local work size for required computation.\n");
	    exit(EXIT_FAILURE);
	}
	clGlobalWorkSize = globalWorkSize(clLocalWorkSize);
	if(setupOpenCL()) {
	    exit(EXIT_FAILURE);
	}
	S->addMessage(1, "(Rates::Rates): Rates ready to work!\n");
}

Rates::~Rates()
{
	if(clKernel)clReleaseKernel(clKernel); clKernel=0;
	if(clSortKernel)clReleaseKernel(clSortKernel); clSortKernel=0;
	if(clProgram)clReleaseProgram(clProgram); clProgram=0;
	if(mPath)delete[] mPath; mPath=0;
}

bool Rates::execute()
{
	InputOutput::ScreenManager *S = InputOutput::ScreenManager::singleton();
	CalcServer *C = CalcServer::singleton();
	cl_int clFlag=0;
	// Sort the rates kernel input data
	clFlag |= sendArgument(clSortKernel,  0, sizeof(cl_mem), (void*)&(C->ifluid));
	clFlag |= sendArgument(clSortKernel,  1, sizeof(cl_mem), (void*)&(C->ifluidin));
	clFlag |= sendArgument(clSortKernel,  2, sizeof(cl_mem), (void*)&(C->imove));
	clFlag |= sendArgument(clSortKernel,  3, sizeof(cl_mem), (void*)&(C->imovein));
	clFlag |= sendArgument(clSortKernel,  4, sizeof(cl_mem), (void*)&(C->pos));
	clFlag |= sendArgument(clSortKernel,  5, sizeof(cl_mem), (void*)&(C->v));
	clFlag |= sendArgument(clSortKernel,  6, sizeof(cl_mem), (void*)&(C->hp));
	clFlag |= sendArgument(clSortKernel,  7, sizeof(cl_mem), (void*)&(C->dens));
	clFlag |= sendArgument(clSortKernel,  8, sizeof(cl_mem), (void*)&(C->press));
	clFlag |= sendArgument(clSortKernel,  9, sizeof(cl_mem), (void*)&(C->mass));
	clFlag |= sendArgument(clSortKernel, 10, sizeof(cl_mem), (void*)&(C->posin));
	clFlag |= sendArgument(clSortKernel, 11, sizeof(cl_mem), (void*)&(C->vin));
	clFlag |= sendArgument(clSortKernel, 12, sizeof(cl_mem), (void*)&(C->hpin));
	clFlag |= sendArgument(clSortKernel, 13, sizeof(cl_mem), (void*)&(C->densin));
	clFlag |= sendArgument(clSortKernel, 14, sizeof(cl_mem), (void*)&(C->pressin));
	clFlag |= sendArgument(clSortKernel, 15, sizeof(cl_mem), (void*)&(C->massin));
	clFlag |= sendArgument(clSortKernel, 16, sizeof(cl_mem), (void*)&(C->permutation));
	clFlag |= sendArgument(clSortKernel, 17, sizeof(cl_mem), (void*)&(C->reversePermutation));
	clFlag |= sendArgument(clSortKernel, 18, sizeof(cl_uint),  (void*)&(C->N));
	if(clFlag != CL_SUCCESS) {
		S->addMessage(3, "(Rates::execute): Can't send arguments to sort kernel.\n");
	    return true;
	}
	#ifdef HAVE_GPUPROFILE
	    cl_event event;
	    cl_ulong end, start;
	    profileTime(0.f);
	    clFlag = clEnqueueNDRangeKernel(C->clComQueue, clSortKernel, 1, NULL, &clGlobalWorkSize, NULL, 0, NULL, &event);
	#else
	    clFlag = clEnqueueNDRangeKernel(C->clComQueue, clSortKernel, 1, NULL, &clGlobalWorkSize, NULL, 0, NULL, NULL);
	#endif
	if(clFlag != CL_SUCCESS) {
		S->addMessage(3, "(Rates::execute): Can't execute the sorting kernel.\n");
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
	        S->addMessage(3, "(Rates::execute): Can't wait to sorting kernel ends.\n");
	        return true;
	    }
	    clFlag |= clGetEventProfilingInfo(event, CL_PROFILING_COMMAND_END, sizeof(cl_ulong), &end, 0);
	    clFlag |= clGetEventProfilingInfo(event, CL_PROFILING_COMMAND_START, sizeof(cl_ulong), &start, 0);
	    if(clFlag != CL_SUCCESS) {
	        S->addMessage(3, "(Rates::execute): Can't profile sorting kernel execution.\n");
	        return true;
	    }
	    profileTime(profileTime() + (end - start)/1000.f);  // 10^-3 ms
	#endif
	// Compute the variation rates
	clFlag |= sendArgument(clKernel,  0, sizeof(cl_mem  ), (void*)&(C->ifluidin));
	clFlag |= sendArgument(clKernel,  1, sizeof(cl_mem  ), (void*)&(C->imovein));
	clFlag |= sendArgument(clKernel,  2, sizeof(cl_mem  ), (void*)&(C->posin));
	clFlag |= sendArgument(clKernel,  3, sizeof(cl_mem  ), (void*)&(C->vin));
	clFlag |= sendArgument(clKernel,  4, sizeof(cl_mem  ), (void*)&(C->densin));
	clFlag |= sendArgument(clKernel,  5, sizeof(cl_mem  ), (void*)&(C->hpin));
	clFlag |= sendArgument(clKernel,  6, sizeof(cl_mem  ), (void*)&(C->massin));
	clFlag |= sendArgument(clKernel,  7, sizeof(cl_mem  ), (void*)&(C->pressin));
	clFlag |= sendArgument(clKernel,  8, sizeof(cl_mem  ), (void*)&(C->visc_kin));
	clFlag |= sendArgument(clKernel,  9, sizeof(cl_mem  ), (void*)&(C->visc_dyn_corrected));
	clFlag |= sendArgument(clKernel, 10, sizeof(cl_mem  ), (void*)&(C->f));
	clFlag |= sendArgument(clKernel, 11, sizeof(cl_mem  ), (void*)&(C->drdt));
	clFlag |= sendArgument(clKernel, 12, sizeof(cl_mem  ), (void*)&(C->drdt_F));
	clFlag |= sendArgument(clKernel, 13, sizeof(cl_mem  ), (void*)&(C->sigma));
	clFlag |= sendArgument(clKernel, 14, sizeof(cl_mem  ), (void*)&(C->shepard));
	clFlag |= sendArgument(clKernel, 15, sizeof(cl_mem  ), (void*)&(C->gradShepard));
	clFlag |= sendArgument(clKernel, 16, sizeof(cl_mem  ), (void*)&(C->lcell));
	clFlag |= sendArgument(clKernel, 17, sizeof(cl_mem  ), (void*)&(C->ihoc));
	clFlag |= sendArgument(clKernel, 18, sizeof(cl_mem  ), (void*)&(C->isValidCell));
	clFlag |= sendArgument(clKernel, 19, sizeof(cl_mem  ), (void*)&(C->permutation));
	clFlag |= sendArgument(clKernel, 20, sizeof(cl_mem  ), (void*)&(C->reversePermutation));
	clFlag |= sendArgument(clKernel, 21, sizeof(cl_mem  ), (void*)&(C->sensorMode));
	clFlag |= sendArgument(clKernel, 22, sizeof(cl_uint ), (void*)&(C->n));
	clFlag |= sendArgument(clKernel, 23, sizeof(cl_uint ), (void*)&(C->N));
	clFlag |= sendArgument(clKernel, 24, sizeof(uivec   ), (void*)&(C->lvec));
	clFlag |= sendArgument(clKernel, 25, sizeof(vec     ), (void*)&(C->g));
	unsigned int nAddedArgs = 0;
	if(isDelta) {
	    clFlag |= sendArgument(clKernel, 26, sizeof(cl_mem), (void*)&(C->refd));
	    clFlag |= sendArgument(clKernel, 27, sizeof(cl_mem), (void*)&(C->delta));
	    clFlag |= sendArgument(clKernel, 28, sizeof(cl_float), (void*)&(C->dt));
	    clFlag |= sendArgument(clKernel, 29, sizeof(cl_float), (void*)&(C->cs));
        nAddedArgs = 4;
	}
	if(isLocalMemory) {
	    clFlag |= sendArgument(clKernel, 26+nAddedArgs, clLocalWorkSize*sizeof(cl_float), NULL);
	    clFlag |= sendArgument(clKernel, 27+nAddedArgs, clLocalWorkSize*sizeof(vec     ), NULL);
	    clFlag |= sendArgument(clKernel, 28+nAddedArgs, clLocalWorkSize*sizeof(cl_float), NULL);
	    clFlag |= sendArgument(clKernel, 29+nAddedArgs, clLocalWorkSize*sizeof(cl_float), NULL);
	    clFlag |= sendArgument(clKernel, 30+nAddedArgs, clLocalWorkSize*sizeof(cl_float), NULL);
	    clFlag |= sendArgument(clKernel, 31+nAddedArgs, clLocalWorkSize*sizeof(vec     ), NULL);
	}
	if(clFlag != CL_SUCCESS) {
		S->addMessage(3, "(Rates::execute): Can't send arguments to kernel.\n");
	    return true;
	}
	#ifdef HAVE_GPUPROFILE
	    clFlag = clEnqueueNDRangeKernel(C->clComQueue, clKernel, 1, NULL, &clGlobalWorkSize, &clLocalWorkSize, 0, NULL, &event);
	#else
	    clFlag = clEnqueueNDRangeKernel(C->clComQueue, clKernel, 1, NULL, &clGlobalWorkSize, &clLocalWorkSize, 0, NULL, NULL);
	#endif
	if(clFlag != CL_SUCCESS) {
		S->addMessage(3, "(Rates::execute): Can't execute the kernel.\n");
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
	        S->addMessage(3, "(Rates::execute): Can't wait to kernels end.\n");
	        return true;
	    }
	    clFlag |= clGetEventProfilingInfo(event, CL_PROFILING_COMMAND_END, sizeof(cl_ulong), &end, 0);
	    clFlag |= clGetEventProfilingInfo(event, CL_PROFILING_COMMAND_START, sizeof(cl_ulong), &start, 0);
	    if(clFlag != CL_SUCCESS) {
	        S->addMessage(3, "(Rates::execute): Can't profile kernel execution.\n");
	        return true;
	    }
	    profileTime(profileTime() + (end - start)/1000.f);  // 10^-3 ms
	#endif
	return false;
}

bool Rates::setupOpenCL()
{
	InputOutput::ScreenManager *S = InputOutput::ScreenManager::singleton();
	CalcServer *C = CalcServer::singleton();
	char msg[1024];
	cl_int clFlag;
	cl_device_id device;
	cl_ulong localMem, reqLocalMem;
	clFlag |= clGetCommandQueueInfo(C->clComQueue,CL_QUEUE_DEVICE,
	                                sizeof(cl_device_id),&device, NULL);
	if(clFlag != CL_SUCCESS) {
		S->addMessage(3, "(Rates::setupOpenCL): Can't get device from command queue.\n");
	    return true;
	}
	clFlag |= clGetDeviceInfo(device, CL_DEVICE_LOCAL_MEM_SIZE, sizeof(localMem), &localMem, NULL);
	if(clFlag != CL_SUCCESS) {
		S->addMessage(3, "(Rates::setupOpenCL): Can't get local memory available on device.\n");
	    return true;
	}
	if(!loadKernelFromFile(&clSortKernel, &clProgram, C->clContext, C->clDevice, mPath, "SortData", ""))
	    return true;
	if(clProgram)clReleaseProgram(clProgram); clProgram=0;
	char args[32]; strcpy(args, "");
	if(isDelta)
        strcat(args, "-D__DELTA_SPH__");
	if(!loadKernelFromFile(&clKernel, &clProgram, C->clContext, C->clDevice, mPath, "Rates", args))
	    return true;
	if(clProgram)clReleaseProgram(clProgram); clProgram=0;
	// Test if there are enough local memory
	clFlag |= clGetKernelWorkGroupInfo(clSortKernel,device,CL_KERNEL_LOCAL_MEM_SIZE,
	                                   sizeof(cl_ulong), &reqLocalMem, NULL);
	if(clFlag != CL_SUCCESS) {
		S->addMessage(3, "(Rates::setupOpenCL): Can't get sort kernel memory usage.\n");
	    return true;
	}
	if(localMem < reqLocalMem){
		S->addMessage(3, "(Rates::setupOpenCL): Not enough local memory for sort execution.\n");
	    sprintf(msg, "\tNeeds %lu bytes, but only %lu bytes are available.\n",
	           reqLocalMem, localMem);
	    S->addMessage(0, msg);
	    return true;
	}
	clFlag |= clGetKernelWorkGroupInfo(clKernel,device,CL_KERNEL_LOCAL_MEM_SIZE,
	                                   sizeof(cl_ulong), &reqLocalMem, NULL);
	if(clFlag != CL_SUCCESS) {
		S->addMessage(3, "(Rates::setupOpenCL): Can't get rates kernel memory usage.\n");
	    return true;
	}
	if(localMem < reqLocalMem){
		S->addMessage(3, "(Rates::setupOpenCL): Not enough local memory for rates execution.\n");
	    sprintf(msg, "\tNeeds %lu bytes, but only %lu bytes are available.\n",
	           reqLocalMem, localMem);
	    S->addMessage(0, msg);
	    return true;
	}
	// Test if local work group size must be modified
	size_t localWorkGroupSize=0;
	clFlag |= clGetKernelWorkGroupInfo(clSortKernel,device,CL_KERNEL_WORK_GROUP_SIZE,
	                                   sizeof(size_t), &localWorkGroupSize, NULL);
	if(clFlag != CL_SUCCESS) {
		S->addMessage(3, "(Rates::setupOpenCL): Can't get sort maximum local work group size.\n");
	    return true;
	}
	if(localWorkGroupSize < clLocalWorkSize)
	    clLocalWorkSize  = localWorkGroupSize;
	clFlag |= clGetKernelWorkGroupInfo(clKernel,device,CL_KERNEL_WORK_GROUP_SIZE,
	                                   sizeof(size_t), &localWorkGroupSize, NULL);
	if(clFlag != CL_SUCCESS) {
		S->addMessage(3, "(Rates::setupOpenCL): Can't get rates maximum local work group size.\n");
	    return true;
	}
	if(localWorkGroupSize < clLocalWorkSize)
	    clLocalWorkSize  = localWorkGroupSize;
	// Look for a better local work group size
	clFlag |= clGetKernelWorkGroupInfo(clKernel,device,CL_KERNEL_PREFERRED_WORK_GROUP_SIZE_MULTIPLE,
	                                   sizeof(size_t), &localWorkGroupSize, NULL);
	if(clFlag != CL_SUCCESS) {
		S->addMessage(3, "(Rates::setupOpenCL): Can't get rates preferred local work group size.\n");
	    return true;
	}
	clLocalWorkSize  = (clLocalWorkSize/localWorkGroupSize) * localWorkGroupSize;
	clGlobalWorkSize = globalWorkSize(clLocalWorkSize);
	// Test if the computation can be accelerated with local memory
	reqLocalMem += clLocalWorkSize*(  sizeof(cl_float)
	                                + sizeof(vec     )
	                                + sizeof(cl_float)
	                                + sizeof(cl_float)
	                                + sizeof(cl_float)
	                                + sizeof(vec     ));
	if(localMem < reqLocalMem){
		S->addMessage(2, "(Rates::setupOpenCL): Not enough local memory for rates.\n");
	    sprintf(msg, "\tNeeds %lu bytes, but only %lu bytes are available.\n",
	           reqLocalMem, localMem);
	    S->addMessage(0, msg);
	    S->addMessage(0, "\tLocal memory usage will be avoided therefore.\n");
	    isLocalMemory = false;
	    char options[19]; strcpy(options,"-D__NO_LOCAL_MEM__");
	    if(clKernel)clReleaseKernel(clKernel); clKernel=0;
	    if(!loadKernelFromFile(&clKernel, &clProgram, C->clContext, C->clDevice, mPath, "Rates", options))
	        return true;
	    if(clProgram)clReleaseProgram(clProgram); clProgram=0;
	}
	return false;
}

}}  // namespace
