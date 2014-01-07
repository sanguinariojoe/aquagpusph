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
	, program(NULL)
	, _path(NULL)
	, kernel(NULL)
	, _global_work_size(0)
	, _local_work_size(0)
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
	_path = new char[nChar+4];
	if(!_path) {
		S->addMessage(3, "(Torque::Torque): Can't allocate memory for path.\n");
		exit(EXIT_FAILURE);
	}
	strcpy(_path, P->OpenCL_kernels.torque);
	strcat(_path, ".cl");

	_local_work_size  = localWorkSize();
	if(!_local_work_size){
	    S->addMessage(3, "(Torque::Torque): I cannot get a valid local work size for the required computation tool.\n");
	    exit(EXIT_FAILURE);
	}
	_global_work_size = globalWorkSize(_local_work_size);
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
	if(kernel)clReleaseKernel(kernel); kernel=0;
	if(program)clReleaseProgram(program); program=0;
	if(_path)delete[] _path; _path=0;
}

bool Torque::execute()
{
	InputOutput::ScreenManager *S = InputOutput::ScreenManager::singleton();
	CalcServer *C = CalcServer::singleton();
	unsigned int i;
	int err_code=0;

	err_code |= sendArgument(kernel, 10, sizeof(vec), (void*)&(C->g));
	err_code |= sendArgument(kernel, 11, sizeof(vec), (void*)&mCOR);
	if(err_code != CL_SUCCESS) {
		S->addMessage(3, "(Torque::Execute): Can't send variable to kernel.\n");
		return true;
	}
	#ifdef HAVE_GPUPROFILE
		cl_event event;
		cl_ulong end, start;
		profileTime(0.f);
		err_code = clEnqueueNDRangeKernel(C->command_queue, kernel, 1, NULL, &_global_work_size, NULL, 0, NULL, &event);
	#else
		err_code = clEnqueueNDRangeKernel(C->command_queue, kernel, 1, NULL, &_global_work_size, NULL, 0, NULL, NULL);
	#endif
	if(err_code != CL_SUCCESS) {
		S->addMessage(3, "(Torque::Execute): Can't execute torque calculation kernel.\n");
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
			S->addMessage(3, "(Torque::Execute): Can't wait to kernels end.\n");
			return true;
		}
		err_code |= clGetEventProfilingInfo(event, CL_PROFILING_COMMAND_END, sizeof(cl_ulong), &end, 0);
		err_code |= clGetEventProfilingInfo(event, CL_PROFILING_COMMAND_START, sizeof(cl_ulong), &start, 0);
		if(err_code != CL_SUCCESS) {
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
	cl_int err_code=0;
	if(!loadKernelFromFile(&kernel, &program, C->context, C->device, _path, "Torque", ""))
		return true;
	mDevTorque = C->allocMemory(C->n * sizeof( vec ));
	mDevForce = C->allocMemory(C->n * sizeof( vec ));
	if(!mDevTorque || !mDevForce)
		return true;
	sprintf(msg, "\tAllocated memory = %u bytes\n", (unsigned int)C->allocated_mem);
	S->addMessage(0, msg);
	err_code  = sendArgument(kernel,  0, sizeof(cl_mem ), (void*)&mDevTorque);
	err_code |= sendArgument(kernel,  1, sizeof(cl_mem ), (void*)&mDevForce);
	err_code |= sendArgument(kernel,  2, sizeof(cl_mem ), (void*)&(C->imove));
	err_code |= sendArgument(kernel,  3, sizeof(cl_mem ), (void*)&(C->ifluid));
	err_code |= sendArgument(kernel,  4, sizeof(cl_mem ), (void*)&(C->pos));
	err_code |= sendArgument(kernel,  5, sizeof(cl_mem ), (void*)&(C->f));
	err_code |= sendArgument(kernel,  6, sizeof(cl_mem ), (void*)&(C->mass));
	err_code |= sendArgument(kernel,  7, sizeof(cl_mem ), (void*)&(C->dens));
	err_code |= sendArgument(kernel,  8, sizeof(cl_mem ), (void*)&(C->refd));
	err_code |= sendArgument(kernel,  9, sizeof(cl_uint), (void*)&(C->n));
	err_code |= sendArgument(kernel, 10, sizeof(vec	   ), (void*)&(C->g));
	if(err_code)
		return true;
	//! Test for right work group size
	cl_device_id device;
	size_t localWorkGroupSize=0;
	err_code |= clGetCommandQueueInfo(C->command_queue,CL_QUEUE_DEVICE,
	                                sizeof(cl_device_id),&device, NULL);
	if(err_code != CL_SUCCESS) {
		S->addMessage(3, "(Torque::setupTorque): I Cannot get the device from the command queue.\n");
	    return true;
	}
	err_code |= clGetKernelWorkGroupInfo(kernel,device,CL_KERNEL_WORK_GROUP_SIZE,
	                                   sizeof(size_t), &localWorkGroupSize, NULL);
	if(err_code != CL_SUCCESS) {
		S->addMessage(3, "(Torque::setupTorque): Can't get maximum local work group size.\n");
	    return true;
	}
	if(localWorkGroupSize < _local_work_size)
	    _local_work_size  = localWorkGroupSize;
	_global_work_size = globalWorkSize(_local_work_size);
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
