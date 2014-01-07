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
#include <CalcServer/Reduction.h>

// ----------------------------------------------------------------------------
// Include the calculation server
// ----------------------------------------------------------------------------
#include <CalcServer.h>

namespace Aqua{ namespace CalcServer{

Reduction::Reduction(cl_mem input, unsigned int N, const char* type, const char* identity, const char* operation)
	: Kernel("Reduction")
	, _path(0)
	, _program(0)
	, kernels(0)
{
	InputOutput::ScreenManager *S = InputOutput::ScreenManager::singleton();
	InputOutput::ProblemSetup *P = InputOutput::ProblemSetup::singleton();

	int nChar = strlen(P->OpenCL_kernels.reduction);
	if(nChar <= 0) {
	    S->addMessage(3, "(Reduction::Reduction): _path of Reduction kernel is empty.\n");
	    exit(EXIT_FAILURE);
	}
	_path = new char[nChar+4];
	if(!_path) {
	    S->addMessage(3, "(Reduction::Reduction): Can't allocate memory for path.\n");
	    exit(EXIT_FAILURE);
	}
	strcpy(_path, P->OpenCL_kernels.reduction);
	strcat(_path, ".cl");

    mInput  = input;
    mMems.push_back(input);
    mN.push_back(N);
    if(setupOpenCL(type, identity, operation))
        exit(EXIT_FAILURE);
	S->addMessage(1, "(Reduction::Reduction): Reduction ready to work!\n");
}

Reduction::~Reduction()
{
    unsigned int i;
    for(i=1;i<mMems.size();i++){ // Don't try to remove input memory object, strat with i=1
        if(mMems.at(i))clReleaseMemObject(mMems.at(i)); mMems.at(i)=NULL;
    }
    for(i=0;i<kernels.size();i++){
        if(kernels.at(i))clReleaseKernel(kernels.at(i)); kernels.at(i)=NULL;
    }
    kernels.clear();
	if(_program)clReleaseProgram(_program); _program=0;
	if(_path) delete[] _path; _path=0;
	mGSize.clear();
	mLSize.clear();
}

cl_mem Reduction::execute()
{
	InputOutput::ScreenManager *S = InputOutput::ScreenManager::singleton();
	InputOutput::ProblemSetup *P = InputOutput::ProblemSetup::singleton();
	CalcServer *C = CalcServer::singleton();
    cl_int flag;
    unsigned int i;
	for(i=0;i<kernels.size();i++){
        size_t _global_work_size = mGSize.at(i);
        size_t _local_work_size  = mLSize.at(i);
		#ifdef HAVE_GPUPROFILE
			cl_event event;
			cl_ulong end, start;
			flag = clEnqueueNDRangeKernel(C->command_queue, kernels.at(i), 1, NULL, &_global_work_size, &_local_work_size, 0, NULL, &event);
		#else
			flag = clEnqueueNDRangeKernel(C->command_queue, kernels.at(i), 1, NULL, &_global_work_size, &_local_work_size, 0, NULL, NULL);
		#endif
		if(flag != CL_SUCCESS) {
			S->addMessage(3, "(Reduction::execute): I cannot execute the kernel.\n");
			if(flag == CL_INVALID_WORK_GROUP_SIZE)
				S->addMessage(0, "\tInvalid local work group size.\n");
			else if(flag == CL_OUT_OF_RESOURCES)
				S->addMessage(0, "\tDevice out of resources.\n");
			else if(flag == CL_MEM_OBJECT_ALLOCATION_FAILURE)
				S->addMessage(0, "\tAllocation error at device.\n");
			else if(flag == CL_OUT_OF_HOST_MEMORY)
				S->addMessage(0, "\tfailure to allocate resources required by the OpenCL implementation on the host.\n");
			return NULL;
		}
		#ifdef HAVE_GPUPROFILE
			flag = clWaitForEvents(1, &event);
			if(flag != CL_SUCCESS) {
				S->addMessage(3, "(Reduction::Execute): Impossible to wait for the kernels end.\n");
				return NULL;
			}
			flag |= clGetEventProfilingInfo(event, CL_PROFILING_COMMAND_END, sizeof(cl_ulong), &end, 0);
			flag |= clGetEventProfilingInfo(event, CL_PROFILING_COMMAND_START, sizeof(cl_ulong), &start, 0);
			if(flag != CL_SUCCESS) {
				S->addMessage(3, "(Reduction::Execute): I cannot profile the kernel execution.\n");
				return NULL;
			}
			profileTime(profileTime() + (end - start)/1000.f);  // 10^-3 ms
		#endif
	}
	return mMems.at(mMems.size()-1);
}

bool Reduction::setInput(cl_mem input)
{
    mMems.at(0) = input;
    cl_int flag = sendArgument(kernels.at(0), 0, sizeof(cl_mem ), (void*)&(input));
    if(flag != CL_SUCCESS)
        return true;
    return false;
}

bool Reduction::setupOpenCL(const char* type, const char* identity, const char* operation)
{
	InputOutput::ScreenManager *S = InputOutput::ScreenManager::singleton();
	CalcServer *C = CalcServer::singleton();
    size_t inputSize, dataSize;
    cl_int flag;
    char msg[1024]; strcpy(msg, "");
    // Get elements data size to can allocate local memory later
    flag = clGetMemObjectInfo(mInput, CL_MEM_SIZE, sizeof(size_t), &inputSize, NULL);
    if(flag != CL_SUCCESS)
        return true;
    dataSize = inputSize / mN.at(0);
    // Create a header for the source code where the operation will be placed
    char header[512];
    strcpy(header, "T reduce(T a, T b) \n");
    strcat(header, "{ \n");
    strcat(header, "\tT c; \n");
    strcat(header, "\t");
    strcat(header, operation);
    strcat(header, "\n");
    strcat(header, "\treturn c;\n");
    strcat(header, "} \n");
    // Starts a dummy kernel in order to study the local size that can be used
    size_t lsize = __CL_MAX_LOCALSIZE__;
    char args[512];
    sprintf(args, "-DT=%s -DIDENTITY=%s -DLOCAL_WORK_SIZE=%luu", type, identity, lsize);
    cl_kernel kernel;
    size_t maxlsize = loadKernelFromFile(&kernel, &_program, C->context, C->device, _path, "Reduction", args, header);
    if(maxlsize < __CL_MIN_LOCALSIZE__){
        S->addMessage(3, "(Reduction::Reduction): Reduction can't be performed due to insufficient local memory\n");
        sprintf(msg, "\t%lu elements can be executed, but __CL_MIN_LOCALSIZE__=%lu\n", maxlsize, __CL_MIN_LOCALSIZE__);
        S->addMessage(0, msg);
        return true;
    }
    lsize = maxlsize;
    if(!isPowerOf2(lsize)){
        lsize = nextPowerOf2(lsize) / 2;
    }
    if(lsize < __CL_MIN_LOCALSIZE__){
        S->addMessage(3, "(Reduction::Reduction): Reduction can't be performed due to insufficient local memory\n");
        sprintf(msg, "\t%lu elements may be executed, but __CL_MIN_LOCALSIZE__=%lu\n", lsize, __CL_MIN_LOCALSIZE__);
        S->addMessage(0, msg);
        return true;
    }
    sprintf(args, "-DT=%s -DIDENTITY=%s -DLOCAL_WORK_SIZE=%luu", type, identity, lsize);
	if(kernel)clReleaseKernel(kernel); kernel=NULL;
	if(_program)clReleaseProgram(_program); _program=NULL;
    // Now we can start a loop while the amount of resulting data was upper than one
    unsigned int N = mN.at(0);
    mN.clear();
    unsigned int i=0;
    while(N > 1){
        // Get work sizes
        mN.push_back(N);                                 // Input data size
        mLSize.push_back(lsize);                         // Feasible local size
        mGSize.push_back(roundUp(N, lsize));             // Global size
        mNGroups.push_back(mGSize.at(i) / mLSize.at(i)); // Number of work groups (and amount of output data)
        // Build the output memory object
        cl_mem output = NULL;
        output = C->allocMemory(mNGroups.at(i) * dataSize);
        if(!output){
            S->addMessage(3, "(Reduction::Reduction): Can't create output array.\n");
            return true;
        }
        mMems.push_back(output);
        // Build the kernel
        if(!loadKernelFromFile(&kernel, &_program, C->context, C->device, _path, "Reduction", args, header))
            return true;
        kernels.push_back(kernel);
        if(_program)clReleaseProgram(_program); _program=NULL;

        flag  = CL_SUCCESS;
        flag |= sendArgument(kernel, 0, sizeof(cl_mem ), (void*)&(mMems.at(i)));
        flag |= sendArgument(kernel, 1, sizeof(cl_mem ), (void*)&(mMems.at(i+1)));
        flag |= sendArgument(kernel, 2, sizeof(cl_uint), (void*)&(N));
        flag |= sendArgument(kernel, 3, lsize*dataSize , NULL);
        if(flag != CL_SUCCESS){
            S->addMessage(3, "(Reduction::Reduction): Arguments set failed\n");
            return true;
        }
        // Setup next step
        sprintf(msg, "\tStep %u, %u elements reduced to %u\n", i, N, mNGroups.at(i));
        S->addMessage(0, msg);
        N = mNGroups.at(i);
        i++;
    }
	return false;
}

}}  // namespaces
