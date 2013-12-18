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
#include <CalcServer/Bounds.h>

// ----------------------------------------------------------------------------
// Include the calculation server
// ----------------------------------------------------------------------------
#include <CalcServer.h>

namespace Aqua{ namespace CalcServer{

Bounds::Bounds()
	: Kernel("Bounds")
	, mDevMem(NULL)
	, program(NULL)
	, mPath(NULL)
	, clMaxCoordsKernel(NULL)
	, clMinCoordsKernel(NULL)
	, clMaxVelKernel(NULL)
	, clMinVelKernel(NULL)
	, global_work_size(0)
	, local_work_size(0)
	, maxCoordsReduction(NULL)
	, minCoordsReduction(NULL)
	, maxVelReduction(NULL)
	, minVelReduction(NULL)
{
	InputOutput::ScreenManager *S = InputOutput::ScreenManager::singleton();
	InputOutput::ProblemSetup *P  = InputOutput::ProblemSetup::singleton();
	unsigned int nChar = strlen(P->OpenCL_kernels.bounds);
	if(nChar <= 0) {
		S->addMessage(3, "(Bounds::Bounds): Path of bounds kernels is empty.\n");
		exit(EXIT_FAILURE);
	}
	mPath = new char[nChar+4];
	if(!mPath) {
		S->addMessage(3, "(Bounds::Bounds): Can't allocate memory for path.\n");
		exit(EXIT_FAILURE);
	}
	strcpy(mPath, P->OpenCL_kernels.bounds);
	strcat(mPath, ".cl");

	local_work_size  = localWorkSize();
	if(!local_work_size){
	    S->addMessage(3, "(Bounds::Bounds): No valid local work size for required computation.\n");
	    exit(EXIT_FAILURE);
	}
	global_work_size = globalWorkSize(local_work_size);
	if(setupBounds()) {
		exit(EXIT_FAILURE);
	}
	if(setupReduction()) {
		exit(EXIT_FAILURE);
	}
	mMaxCoords.x = 0.f;
	mMaxCoords.y = 0.f;
	mMinCoords.x  = 0.f;
	mMinCoords.y  = 0.f;
	mMaxVel.x = 0.f;
	mMaxVel.y = 0.f;
	mMinVel.x  = 0.f;
	mMinVel.y  = 0.f;
	#ifdef HAVE_3D
		mMaxCoords.z = 0.f;
		mMinCoords.z  = 0.f;
		mMaxVel.z = 0.f;
		mMinVel.z  = 0.f;
	#endif
	S->addMessage(1, "(Bounds::Bounds): Bounds ready to work!\n");
}

Bounds::~Bounds()
{
	InputOutput::ScreenManager *S = InputOutput::ScreenManager::singleton();
	S->addMessage(1, "(Bounds::~Bounds): Destroying maximum coordinates reduction processor...\n");
	if(maxCoordsReduction) delete maxCoordsReduction;
	S->addMessage(1, "(Bounds::~Bounds): Destroying minimum coordinates reduction processor...\n");
	if(minCoordsReduction)  delete minCoordsReduction;
	S->addMessage(1, "(Bounds::~Bounds): Destroying maximum velocity reduction processor...\n");
	if(maxVelReduction) delete maxVelReduction;
	S->addMessage(1, "(Bounds::~Bounds): Destroying minimum velocity reduction processor...\n");
	if(minVelReduction)  delete minVelReduction;
	if(mDevMem)clReleaseMemObject(mDevMem); mDevMem=0;
	if(clMaxCoordsKernel)clReleaseKernel(clMaxCoordsKernel); clMaxCoordsKernel=0;
	if(clMinCoordsKernel)clReleaseKernel(clMinCoordsKernel); clMinCoordsKernel=0;
	if(clMaxVelKernel)clReleaseKernel(clMaxVelKernel); clMaxVelKernel=0;
	if(clMinVelKernel)clReleaseKernel(clMinVelKernel); clMinVelKernel=0;
	if(program)clReleaseProgram(program); program=0;
	if(mPath)delete[] mPath; mPath=0;
}

bool Bounds::execute()
{
    InputOutput::ScreenManager *S = InputOutput::ScreenManager::singleton();
	CalcServer *C = CalcServer::singleton();

    if(execute(&mMaxCoords, __BOUNDS_COORDS_MAX_OP__)){
        S->addMessage(0, "\tDuring the maximum coordinates computation.\n");
        return true;
    }
    if(execute(&mMinCoords, __BOUNDS_COORDS_MIN_OP__)){
        S->addMessage(0, "\tDuring the minimum coordinates computation.\n");
        return true;
    }
    if(execute(&mMaxVel, __BOUNDS_VEL_MAX_OP__)){
        S->addMessage(0, "\tDuring the maximum velocity computation.\n");
        return true;
    }
    if(execute(&mMinVel, __BOUNDS_VEL_MIN_OP__)){
        S->addMessage(0, "\tDuring the minimum velocity computation.\n");
        return true;
    }
    return false;
}

bool Bounds::execute(vec *output, int op)
{
	InputOutput::ScreenManager *S = InputOutput::ScreenManager::singleton();
	CalcServer *C = CalcServer::singleton();
	int err_code=0;
    // Get the desired tools
	cl_kernel kernel = NULL;
	Reduction *reduction = NULL;
	if(op == __BOUNDS_COORDS_MAX_OP__){
        kernel = clMaxCoordsKernel;
        reduction = maxCoordsReduction;
	}
	if(op == __BOUNDS_COORDS_MIN_OP__){
        kernel = clMinCoordsKernel;
        reduction = minCoordsReduction;
	}
	if(op == __BOUNDS_VEL_MAX_OP__){
        kernel = clMaxVelKernel;
        reduction = maxVelReduction;
	}
	if(op == __BOUNDS_VEL_MIN_OP__){
        kernel = clMinVelKernel;
        reduction = minVelReduction;
	}
    // Prepare the data array, discarding fixed particles
	#ifdef HAVE_GPUPROFILE
		cl_event event;
		cl_ulong end, start;
		profileTime(0.f);
		err_code = clEnqueueNDRangeKernel(C->command_queue, kernel, 1, NULL, &global_work_size, NULL, 0, NULL, &event);
	#else
		err_code = clEnqueueNDRangeKernel(C->command_queue, kernel, 1, NULL, &global_work_size, NULL, 0, NULL, NULL);
	#endif
	if(err_code != CL_SUCCESS) {
		S->addMessage(3, "(Bounds::Execute): Can't execute bounds data setup kernel.\n");
        S->printOpenCLError(err_code);
		return true;
	}
	#ifdef HAVE_GPUPROFILE
		err_code = clWaitForEvents(1, &event);
		if(err_code != CL_SUCCESS) {
			S->addMessage(3, "(Bounds::Execute): Can't wait to the kernel end.\n");
            S->printOpenCLError(err_code);
			return true;
		}
		err_code |= clGetEventProfilingInfo(event, CL_PROFILING_COMMAND_END, sizeof(cl_ulong), &end, 0);
		err_code |= clGetEventProfilingInfo(event, CL_PROFILING_COMMAND_START, sizeof(cl_ulong), &start, 0);
		if(err_code != CL_SUCCESS) {
			S->addMessage(3, "(Bounds::Execute): Can't profile the kernel execution.\n");
            S->printOpenCLError(err_code);
			return true;
		}
		profileTime(profileTime() + (end - start)/1000.f);  // 10^-3 ms
	#endif
	// Compute the maximum/minimum
    cl_mem devOutput = reduction->execute();
    if(!devOutput)
        return true;
	if(C->getData((void *)output, devOutput, sizeof(vec)))
		return true;
	return false;
}

bool Bounds::setupBounds()
{
	InputOutput::ScreenManager *S = InputOutput::ScreenManager::singleton();
	CalcServer *C = CalcServer::singleton();
	cl_int err_code = 0;
	char msg[1024];
	mDevMem = C->allocMemory(C->n * sizeof( vec ));
	if(!mDevMem)
		return true;
	sprintf(msg, "\tAllocated memory = %u bytes\n", (unsigned int)C->allocated_mem);
	S->addMessage(0, msg);
	if(!loadKernelFromFile(&clMaxCoordsKernel, &program, C->context, C->device, mPath, "MaximumCoordsFilter", ""))
		return true;
	err_code  = sendArgument(clMaxCoordsKernel,  0, sizeof(cl_mem ), (void*)&mDevMem);
	err_code |= sendArgument(clMaxCoordsKernel,  1, sizeof(cl_mem ), (void*)&(C->imove));
	err_code |= sendArgument(clMaxCoordsKernel,  2, sizeof(cl_mem ), (void*)&(C->pos));
	err_code |= sendArgument(clMaxCoordsKernel,  3, sizeof(cl_uint), (void*)&(C->n));
	if(err_code)
		return true;
	if(!loadKernelFromFile(&clMinCoordsKernel, &program, C->context, C->device, mPath, "MinimumCoordsFilter", ""))
		return true;
	err_code  = sendArgument(clMinCoordsKernel,  0, sizeof(cl_mem ), (void*)&mDevMem);
	err_code |= sendArgument(clMinCoordsKernel,  1, sizeof(cl_mem ), (void*)&(C->imove));
	err_code |= sendArgument(clMinCoordsKernel,  2, sizeof(cl_mem ), (void*)&(C->pos));
	err_code |= sendArgument(clMinCoordsKernel,  3, sizeof(cl_uint), (void*)&(C->n));
	if(err_code)
		return true;
	if(!loadKernelFromFile(&clMaxVelKernel, &program, C->context, C->device, mPath, "MaximumVelFilter", ""))
		return true;
	err_code  = sendArgument(clMaxVelKernel,  0, sizeof(cl_mem ), (void*)&mDevMem);
	err_code |= sendArgument(clMaxVelKernel,  1, sizeof(cl_mem ), (void*)&(C->imove));
	err_code |= sendArgument(clMaxVelKernel,  2, sizeof(cl_mem ), (void*)&(C->v));
	err_code |= sendArgument(clMaxVelKernel,  3, sizeof(cl_uint), (void*)&(C->n));
	if(err_code)
		return true;
	if(!loadKernelFromFile(&clMinVelKernel, &program, C->context, C->device, mPath, "MinimumVelFilter", ""))
		return true;
	err_code  = sendArgument(clMinVelKernel,  0, sizeof(cl_mem ), (void*)&mDevMem);
	err_code |= sendArgument(clMinVelKernel,  1, sizeof(cl_mem ), (void*)&(C->imove));
	err_code |= sendArgument(clMinVelKernel,  2, sizeof(cl_mem ), (void*)&(C->v));
	err_code |= sendArgument(clMinVelKernel,  3, sizeof(cl_uint), (void*)&(C->n));
	if(err_code)
		return true;
	//! Test for right work group size
	cl_device_id device;
	size_t localWorkGroupSize=0;
	err_code |= clGetCommandQueueInfo(C->command_queue,CL_QUEUE_DEVICE,
	                                sizeof(cl_device_id),&device, NULL);
	if(err_code != CL_SUCCESS) {
		S->addMessage(3, "(Bounds::setupBounds): Can't get device from command queue.\n");
	    S->printOpenCLError(err_code);
	    return true;
	}
	err_code |= clGetKernelWorkGroupInfo(clMaxCoordsKernel,device,CL_KERNEL_WORK_GROUP_SIZE,
	                                   sizeof(size_t), &localWorkGroupSize, NULL);
	if(err_code != CL_SUCCESS) {
		S->addMessage(3, "(Bounds::setupBounds): Can't get maximum local work group size.\n");
	    S->printOpenCLError(err_code);
	    return true;
	}
	if(localWorkGroupSize < local_work_size)
	    local_work_size  = localWorkGroupSize;
	global_work_size = globalWorkSize(local_work_size);
	return false;
}

bool Bounds::setupReduction()
{
	CalcServer *C = CalcServer::singleton();
    char operation[512];
    strcpy(operation, "");
    strcat(operation, "c.x = (a.x > b.x) ? a.x : b.x;\n");
    strcat(operation, "\tc.y = (a.y > b.y) ? a.y : b.y;\n");
    #ifdef HAVE_3D
        strcat(operation, "\tc.z = (a.z > b.z) ? a.z : b.z;\n");
        strcat(operation, "\tc.w = 0.f;\n");
        maxCoordsReduction = new Reduction(mDevMem, C->n, "vec", "(vec)(-INFINITY,-INFINITY,-INFINITY,0.f)", operation);
    #else
        maxCoordsReduction = new Reduction(mDevMem, C->n, "vec", "(vec)(-INFINITY,-INFINITY)", operation);
    #endif
    strcpy(operation, "");
    strcat(operation, "c.x = (a.x < b.x) ? a.x : b.x;\n");
    strcat(operation, "\tc.y = (a.y < b.y) ? a.y : b.y;\n");
    #ifdef HAVE_3D
        strcat(operation, "\tc.z = (a.z < b.z) ? a.z : b.z;\n");
        strcat(operation, "\tc.w = 0.f;\n");
        minCoordsReduction = new Reduction(mDevMem, C->n, "vec", "(vec)(INFINITY,INFINITY,INFINITY,0.f)", operation);
    #else
        minCoordsReduction = new Reduction(mDevMem, C->n, "vec", "(vec)(INFINITY,INFINITY)", operation);
    #endif
    maxVelReduction = new Reduction(mDevMem, C->n, "vec", "VEC_ZERO", "c = (dot(a,a) > dot(b,b)) ? a : b;\n");
    #ifdef HAVE_3D
        minVelReduction = new Reduction(mDevMem, C->n, "vec", "(vec)(INFINITY,INFINITY,INFINITY,0.f)", "c = (dot(a,a) < dot(b,b)) ? a : b;\n");
    #else
        minVelReduction = new Reduction(mDevMem, C->n, "vec", "(vec)(INFINITY,INFINITY)", "c = (dot(a,a) < dot(b,b)) ? a : b;\n");
    #endif

	return false;
}

}}  // namespace
