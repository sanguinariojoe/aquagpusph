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
#include <CalcServer/Torque.h>

// ----------------------------------------------------------------------------
// Include the calculation server
// ----------------------------------------------------------------------------
#include <CalcServer.h>

namespace Aqua{ namespace CalcServer{

Torque::Torque()
	: Kernel("Torque")
	, mDevTorque(NULL)
	, mDevForce(NULL)
	, clProgram(NULL)
	, mPath(NULL)
	, clKernel(NULL)
	, clGlobalWorkSize(0)
	, clLocalWorkSize(0)
	, torqueReduction(NULL)
	, forceReduction(NULL)
{
	InputOutput::ScreenManager *S = InputOutput::ScreenManager::singleton();
	InputOutput::ProblemSetup *P  = InputOutput::ProblemSetup::singleton();
	unsigned int nChar = strlen(P->OpenCL_kernels.torque);
	if(nChar <= 0) {
		S->addMessage(3, "(Torque::Torque): Path of reduction kernels (2D) is empty.\n");
		exit(EXIT_FAILURE);
	}
	mPath = new char[nChar+4];
	if(!mPath) {
		S->addMessage(3, "(Torque::Torque): Can't allocate memory for path.\n");
		exit(EXIT_FAILURE);
	}
	strcpy(mPath, P->OpenCL_kernels.torque);
	strcat(mPath, ".cl");

	clLocalWorkSize  = localWorkSize();
	if(!clLocalWorkSize){
	    S->addMessage(3, "(Torque::Torque): No valid local work size for required computation.\n");
	    exit(EXIT_FAILURE);
	}
	clGlobalWorkSize = globalWorkSize(clLocalWorkSize);
	if(setupTorque()) {
		exit(EXIT_FAILURE);
	}
	if(setupReduction()) {
		exit(EXIT_FAILURE);
	}
	mTorque.x = 0.f;
	mTorque.y = 0.f;
	mForce.x  = 0.f;
	mForce.y  = 0.f;
	#ifdef HAVE_3D
		mTorque.z = 0.f;
		mForce.z  = 0.f;
	#endif
	S->addMessage(1, "(Torque::Torque): Torque ready to work!\n");
}

Torque::~Torque()
{
	InputOutput::ScreenManager *S = InputOutput::ScreenManager::singleton();
	S->addMessage(1, "(Torque::~Torque): Destroying torque reduction processor...\n");
	if(torqueReduction) delete torqueReduction;
	S->addMessage(1, "(Torque::~Torque): Destroying force reduction processor...\n");
	if(forceReduction)  delete forceReduction;
	if(mDevTorque)clReleaseMemObject(mDevTorque); mDevTorque=0;
	if(mDevForce)clReleaseMemObject(mDevForce); mDevForce=0;
	if(clKernel)clReleaseKernel(clKernel); clKernel=0;
	if(clProgram)clReleaseProgram(clProgram); clProgram=0;
	if(mPath)delete[] mPath; mPath=0;
}

bool Torque::execute()
{
	InputOutput::ScreenManager *S = InputOutput::ScreenManager::singleton();
	CalcServer *C = CalcServer::singleton();
	unsigned int i;
	int clFlag=0;

	clFlag |= sendArgument(clKernel, 10, sizeof(vec), (void*)&(C->g));
	clFlag |= sendArgument(clKernel, 11, sizeof(vec), (void*)&mCOR);
	if(clFlag != CL_SUCCESS) {
		S->addMessage(3, "(Torque::Execute): Can't send variable to kernel.\n");
		return true;
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
		S->addMessage(3, "(Torque::Execute): Can't execute torque calculation kernel.\n");
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
			S->addMessage(0, "\tfailure to allocate resources required by the OpenCL implementation on the host.\n");
		}
		return true;
	}
	#ifdef HAVE_GPUPROFILE
		clFlag = clWaitForEvents(1, &event);
		if(clFlag != CL_SUCCESS) {
			S->addMessage(3, "(Torque::Execute): Can't wait to kernels end.\n");
			return true;
		}
		clFlag |= clGetEventProfilingInfo(event, CL_PROFILING_COMMAND_END, sizeof(cl_ulong), &end, 0);
		clFlag |= clGetEventProfilingInfo(event, CL_PROFILING_COMMAND_START, sizeof(cl_ulong), &start, 0);
		if(clFlag != CL_SUCCESS) {
			S->addMessage(3, "(Torque::Execute): Can't profile kernel execution.\n");
			return true;
		}
		profileTime(profileTime() + (end - start)/1000.f);  // 10^-3 ms
	#endif
    cl_mem output;
    output = torqueReduction->execute();
    if(!output)
        return true;
	if(C->getData((void *)&mTorque, output, sizeof(vec)))
		return true;
    output = forceReduction->execute();
    if(!output)
        return true;
	if(C->getData((void *)&mForce, output, sizeof(vec)))
		return true;
	return false;
}

bool Torque::setupTorque()
{
	InputOutput::ScreenManager *S = InputOutput::ScreenManager::singleton();
	CalcServer *C = CalcServer::singleton();
	char msg[1024];
	cl_int clFlag=0;
	if(!loadKernelFromFile(&clKernel, &clProgram, C->clContext, C->clDevice, mPath, "Torque", ""))
		return true;
	mDevTorque = C->allocMemory(C->n * sizeof( vec ));
	mDevForce = C->allocMemory(C->n * sizeof( vec ));
	if(!mDevTorque || !mDevForce)
		return true;
	sprintf(msg, "\tAllocated memory = %u bytes\n", (unsigned int)C->AllocatedMem);
	S->addMessage(0, msg);
	clFlag  = sendArgument(clKernel,  0, sizeof(cl_mem ), (void*)&mDevTorque);
	clFlag |= sendArgument(clKernel,  1, sizeof(cl_mem ), (void*)&mDevForce);
	clFlag |= sendArgument(clKernel,  2, sizeof(cl_mem ), (void*)&(C->imove));
	clFlag |= sendArgument(clKernel,  3, sizeof(cl_mem ), (void*)&(C->ifluid));
	clFlag |= sendArgument(clKernel,  4, sizeof(cl_mem ), (void*)&(C->pos));
	clFlag |= sendArgument(clKernel,  5, sizeof(cl_mem ), (void*)&(C->f));
	clFlag |= sendArgument(clKernel,  6, sizeof(cl_mem ), (void*)&(C->mass));
	clFlag |= sendArgument(clKernel,  7, sizeof(cl_mem ), (void*)&(C->dens));
	clFlag |= sendArgument(clKernel,  8, sizeof(cl_mem ), (void*)&(C->refd));
	clFlag |= sendArgument(clKernel,  9, sizeof(cl_uint), (void*)&(C->n));
	clFlag |= sendArgument(clKernel, 10, sizeof(vec	   ), (void*)&(C->g));
	if(clFlag)
		return true;
	//! Test for right work group size
	cl_device_id device;
	size_t localWorkGroupSize=0;
	clFlag |= clGetCommandQueueInfo(C->clComQueue,CL_QUEUE_DEVICE,
	                                sizeof(cl_device_id),&device, NULL);
	if(clFlag != CL_SUCCESS) {
		S->addMessage(3, "(Torque::setupTorque): Can't get device from command queue.\n");
	    return true;
	}
	clFlag |= clGetKernelWorkGroupInfo(clKernel,device,CL_KERNEL_WORK_GROUP_SIZE,
	                                   sizeof(size_t), &localWorkGroupSize, NULL);
	if(clFlag != CL_SUCCESS) {
		S->addMessage(3, "(Torque::setupTorque): Can't get maximum local work group size.\n");
	    return true;
	}
	if(localWorkGroupSize < clLocalWorkSize)
	    clLocalWorkSize  = localWorkGroupSize;
	clGlobalWorkSize = globalWorkSize(clLocalWorkSize);
	return false;
}

bool Torque::setupReduction()
{
	CalcServer *C = CalcServer::singleton();
	#ifdef HAVE_3D
        torqueReduction = new Reduction(mDevTorque, C->n, "vec", "(vec)(0.f,0.f,0.f,0.f)", "c = a + b;");
        forceReduction = new Reduction(mDevForce, C->n, "vec", "(vec)(0.f,0.f,0.f,0.f)", "c = a + b;");
    #else
        torqueReduction = new Reduction(mDevTorque, C->n, "vec", "(vec)(0.f,0.f)", "c = a + b;");
        forceReduction = new Reduction(mDevForce, C->n, "vec", "(vec)(0.f,0.f)", "c = a + b;");
    #endif
	return false;
}

}}  // namespace
