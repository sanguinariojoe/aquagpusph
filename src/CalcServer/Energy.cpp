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

#include <TimeManager.h>

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
#include <CalcServer/Energy.h>

// ----------------------------------------------------------------------------
// Include the calculation server
// ----------------------------------------------------------------------------
#include <CalcServer.h>

namespace Aqua{ namespace CalcServer{

Energy::Energy()
	: Kernel("Energy")
	, mDevEnergy(NULL)
	, mTime(0.f)
	, mPath(NULL)
	, clProgram(NULL)
	, clKernel(NULL)
	, clGlobalWorkSize(0)
	, clLocalWorkSize(0)
	, mReduction(NULL)
{
	InputOutput::ScreenManager *S = InputOutput::ScreenManager::singleton();
	InputOutput::ProblemSetup *P  = InputOutput::ProblemSetup::singleton();
	unsigned int nChar = strlen(P->OpenCLKernels.Energy);
	if(nChar <= 0) {
		S->addMessage(3, "(Energy::Energy): Path of the kernels (2D) is empty.\n");
		exit(EXIT_FAILURE);
	}
	mPath = new char[nChar+4];
	if(!mPath) {
		S->addMessage(3, "(Energy::Energy): Can't allocate memory for path.\n");
		exit(EXIT_FAILURE);
	}
	strcpy(mPath, P->OpenCLKernels.Energy);
	strcat(mPath, ".cl");

	clLocalWorkSize  = localWorkSize();
	if(!clLocalWorkSize){
	    S->addMessage(3, "(Energy::Energy): No valid local work size for required computation.\n");
	    exit(EXIT_FAILURE);
	}
	clGlobalWorkSize = globalWorkSize(clLocalWorkSize);
	if(setupEnergy()) {
		exit(EXIT_FAILURE);
	}
	if(setupReduction()) {
		exit(EXIT_FAILURE);
	}
	mEnergy.x = 0.f;
	mEnergy.y = 0.f;
	mEnergy.z = 0.f;
	mEnergy.w = 0.f;
	S->addMessage(1, "(Energy::Energy): Energy ready to work!\n");
}

Energy::~Energy()
{
	InputOutput::ScreenManager *S = InputOutput::ScreenManager::singleton();
	S->addMessage(1, "(Energy::~Energy): Destroying reduction processor...\n");
	if(mReduction) delete mReduction;
	if(mDevEnergy)clReleaseMemObject(mDevEnergy); mDevEnergy=0;
	if(clKernel)clReleaseKernel(clKernel); clKernel=0;
	if(clProgram)clReleaseProgram(clProgram); clProgram=0;
	if(mPath)delete[] mPath; mPath=0;
}

bool Energy::execute()
{
	InputOutput::ScreenManager *S = InputOutput::ScreenManager::singleton();
	InputOutput::TimeManager *T   = InputOutput::TimeManager::singleton();
	CalcServer *C = CalcServer::singleton();
	int clFlag=0;
	clFlag |= sendArgument(clKernel, 13, sizeof(cl_float), (void*)&(C->cs));
	clFlag |= sendArgument(clKernel, 14, sizeof(vec     ), (void*)&(C->g));
	if(clFlag != CL_SUCCESS) {
		S->addMessage(3, "(Energy::Execute): Fail senting the kernel parameters.\n");
		S->printOpenCLError(clFlag);
	}
	#ifdef HAVE_GPUPROFILE
		cl_event event;
		cl_ulong end, start;
		profileTime(0.f);
		clFlag = clEnqueueNDRangeKernel(C->clComQueue, clKernel, 1, NULL, &clGlobalWorkSize, NULL, 0, NULL, &event);
	#else
		clFlag = clEnqueueNDRangeKernel(C->clComQueue, clKernel, 1, NULL, &clGlobalWorkSize, NULL, 0, NULL, NULL);
	#endif
	if(clFlag != CL_SUCCESS) {
		S->addMessage(3, "(Energy::Execute): Can't execute the energy calculation kernel.\n");
		if(clFlag == CL_INVALID_WORK_GROUP_SIZE) {
			S->addMessage(0, "\tInvalid local work group size.\n");
		}
		else if(clFlag == CL_OUT_OF_RESOURCES) {
			S->addMessage(0, "\tDevice out of resources.\n");
		}
		else if(clFlag == CL_MEM_OBJECT_ALLOCATION_FAILURE) {
			S->addMessage(0, "\tAllocation error at device.\n");
		}
		else if(clFlag == CL_OUT_OF_HOST_MEMORY) {
			S->addMessage(0, "\tFailure allocating resources required by the OpenCL implementation on the host.\n");
		}
		return true;
	}
	#ifdef HAVE_GPUPROFILE
		clFlag = clWaitForEvents(1, &event);
		if(clFlag != CL_SUCCESS) {
			S->addMessage(3, "(Energy::Execute): Can't wait to kernels end.\n");
			return true;
		}
		clFlag |= clGetEventProfilingInfo(event, CL_PROFILING_COMMAND_END, sizeof(cl_ulong), &end, 0);
		clFlag |= clGetEventProfilingInfo(event, CL_PROFILING_COMMAND_START, sizeof(cl_ulong), &start, 0);
		if(clFlag != CL_SUCCESS) {
			S->addMessage(3, "(Energy::Execute): Can't profile kernel execution.\n");
			return true;
		}
		profileTime(profileTime() + (end - start)/1000.f);  // 10^-3 ms
	#endif
    cl_mem output;
    output = mReduction->execute();
    if(!output)
        return true;
    // Get the computed compnents
	if(C->getData((void *)&mDEnergyDT, output, sizeof(vec4)))
		return true;
    // Integrate the elastic energy and add it to the total one
    float dt = T->time() - mTime;
    mTime = T->time();
    mEnergy.x += mDEnergyDT.x * dt;
    mEnergy.y += mDEnergyDT.y * dt;
    mEnergy.z += mDEnergyDT.z;
    mEnergy.w += mDEnergyDT.w;
	return false;
}

bool Energy::setupEnergy()
{
	InputOutput::ScreenManager *S = InputOutput::ScreenManager::singleton();
	CalcServer *C = CalcServer::singleton();
	char msg[1024];
	if(!loadKernelFromFile(&clKernel, &clProgram, C->clContext, C->clDevice, mPath, "Energy", "-Dvec4=float4"))
		return true;
	cl_int clFlag = C->allocMemory(&mDevEnergy, C->n * sizeof( vec4 ));
	if(clFlag)
		return true;
	sprintf(msg, "\tAllocated memory = %u bytes\n", (unsigned int)C->AllocatedMem);
	S->addMessage(0, msg);
	clFlag  = sendArgument(clKernel,  0, sizeof(cl_mem  ), (void*)&mDevEnergy);
	clFlag |= sendArgument(clKernel,  1, sizeof(cl_mem  ), (void*)&(C->imove));
	clFlag |= sendArgument(clKernel,  2, sizeof(cl_mem  ), (void*)&(C->ifluid));
	clFlag |= sendArgument(clKernel,  3, sizeof(cl_mem  ), (void*)&(C->pos));
	clFlag |= sendArgument(clKernel,  4, sizeof(cl_mem  ), (void*)&(C->v));
	clFlag |= sendArgument(clKernel,  5, sizeof(cl_mem  ), (void*)&(C->mass));
	clFlag |= sendArgument(clKernel,  6, sizeof(cl_mem  ), (void*)&(C->dens));
	clFlag |= sendArgument(clKernel,  7, sizeof(cl_mem  ), (void*)&(C->press));
	clFlag |= sendArgument(clKernel,  8, sizeof(cl_mem  ), (void*)&(C->drdtin));
	clFlag |= sendArgument(clKernel,  9, sizeof(cl_mem  ), (void*)&(C->drdt_F));
	clFlag |= sendArgument(clKernel, 10, sizeof(cl_mem  ), (void*)&(C->fin));
	clFlag |= sendArgument(clKernel, 11, sizeof(cl_mem  ), (void*)&(C->refd));
	clFlag |= sendArgument(clKernel, 12, sizeof(cl_mem  ), (void*)&(C->gamma));
	clFlag |= sendArgument(clKernel, 13, sizeof(cl_float), (void*)&(C->cs));
	clFlag |= sendArgument(clKernel, 14, sizeof(vec     ), (void*)&(C->g));
	clFlag |= sendArgument(clKernel, 15, sizeof(cl_uint ), (void*)&(C->n));
	if(clFlag)
		return true;
	//! Test for right work group size
	cl_device_id device;
	size_t localWorkGroupSize=0;
	clFlag |= clGetCommandQueueInfo(C->clComQueue,CL_QUEUE_DEVICE,
	                                sizeof(cl_device_id),&device, NULL);
	if(clFlag != CL_SUCCESS) {
		S->addMessage(3, "(Energy::setupEnergy): Can't get device from command queue.\n");
	    return true;
	}
	clFlag |= clGetKernelWorkGroupInfo(clKernel,device,CL_KERNEL_WORK_GROUP_SIZE,
	                                   sizeof(size_t), &localWorkGroupSize, NULL);
	if(clFlag != CL_SUCCESS) {
		S->addMessage(3, "(Energy::setupEnergy): Can't get maximum local work group size.\n");
	    return true;
	}
	if(localWorkGroupSize < clLocalWorkSize)
	    clLocalWorkSize  = localWorkGroupSize;
	clGlobalWorkSize = globalWorkSize(clLocalWorkSize);
	return false;
}

bool Energy::setupReduction()
{
	CalcServer *C = CalcServer::singleton();
    mReduction = new Reduction(mDevEnergy, C->n, "float4", "(float4)(0.f,0.f,0.f,0.f)", "c = a + b;");
	return false;
}

}}  // namespace
