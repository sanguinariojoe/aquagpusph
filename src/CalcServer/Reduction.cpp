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

#include <stdlib.h>
#include <math.h>
#include <vector>

#include <ProblemSetup.h>
#include <ScreenManager.h>
#include <CalcServer/Reduction.h>
#include <CalcServer.h>

namespace Aqua{ namespace CalcServer{

Reduction::Reduction(cl_mem input,
                     unsigned int n,
                     const char* type,
                     const char* null_val,
                     const char* operation)
	: Kernel("Reduction")
	, _path(0)
	, _program(0)
	, _kernels(0)
{
	InputOutput::ScreenManager *S = InputOutput::ScreenManager::singleton();
	InputOutput::ProblemSetup *P = InputOutput::ProblemSetup::singleton();

	int str_len = strlen(P->OpenCL_kernels.reduction);
	if(str_len <= 0) {
	    S->addMessageF(3, "The path of Reduction kernel is empty.\n");
	    exit(EXIT_FAILURE);
	}
	_path = new char[str_len+4];
	if(!_path) {
	    S->addMessageF(3, "Memory cannot be allocated for the path.\n");
	    exit(EXIT_FAILURE);
	}
	strcpy(_path, P->OpenCL_kernels.reduction);
	strcat(_path, ".cl");

    _input  = input;
    _mems.push_back(input);
    _n.push_back(n);
    if(setupOpenCL(type, null_val, operation))
        exit(EXIT_FAILURE);
	S->addMessageF(1, "Reduction ready to work!\n");
}

Reduction::~Reduction()
{
    unsigned int i;
    for(i=1;i<_mems.size();i++){
        if(_mems.at(i))
            clReleaseMemObject(_mems.at(i));
            _mems.at(i)=NULL;
    }
    for(i=0;i<_kernels.size();i++){
        if(_kernels.at(i))
            clReleaseKernel(_kernels.at(i));
        _kernels.at(i)=NULL;
    }
    _kernels.clear();
	if(_program) clReleaseProgram(_program); _program=0;
	if(_path) delete[] _path; _path=0;
	_global_work_sizes.clear();
	_local_work_sizes.clear();
}

cl_mem Reduction::execute()
{
	InputOutput::ScreenManager *S = InputOutput::ScreenManager::singleton();
	InputOutput::ProblemSetup *P = InputOutput::ProblemSetup::singleton();
	CalcServer *C = CalcServer::singleton();
    cl_int err_code;
    unsigned int i;
	for(i=0;i<_kernels.size();i++){
        size_t _global_work_size = _global_work_sizes.at(i);
        size_t _local_work_size  = _local_work_sizes.at(i);
		#ifdef HAVE_GPUPROFILE
			cl_event event;
			cl_ulong end, start;
			err_code = clEnqueueNDRangeKernel(C->command_queue,
                                              _kernels.at(i),
                                              1,
                                              NULL,
                                              &_global_work_size,
                                              &_local_work_size,
                                              0,
                                              NULL,
                                              &event);
		#else
			err_code = clEnqueueNDRangeKernel(C->command_queue,
                                              _kernels.at(i),
                                              1,
                                              NULL,
                                              &_global_work_size,
                                              &_local_work_size,
                                              0,
                                              NULL,
                                              NULL);
		#endif
		if(err_code != CL_SUCCESS) {
			S->addMessageF(3, "I cannot execute the kernel.\n");
            S->printOpenCLError(err_code);
            return NULL;
        }
        #ifdef HAVE_GPUPROFILE
            err_code = clWaitForEvents(1, &event);
            if(err_code != CL_SUCCESS) {
                S->addMessage(3, "Impossible to wait for the kernels end.\n");
                S->printOpenCLError(err_code);
                return NULL;
            }
            err_code |= clGetEventProfilingInfo(event,
                                                CL_PROFILING_COMMAND_END,
                                                sizeof(cl_ulong),
                                                &end,
                                                0);
            if(err_code != CL_SUCCESS) {
                S->addMessage(3, "I cannot profile the kernel execution.\n");
                S->printOpenCLError(err_code);
                return NULL;
            }
            err_code |= clGetEventProfilingInfo(event,
                                                CL_PROFILING_COMMAND_START,
                                                sizeof(cl_ulong),
                                                &start,
                                                0);
            if(err_code != CL_SUCCESS) {
                S->addMessage(3, "I cannot profile the kernel execution.\n");
                S->printOpenCLError(err_code);
                return NULL;
            }
			profileTime(profileTime() + (end - start)/1000.f);  // 10^-3 ms
		#endif
	}
	return _mems.at(_mems.size()-1);
}

bool Reduction::setInput(cl_mem input)
{
    _mems.at(0) = input;
    cl_int err_code = sendArgument(_kernels.at(0),
                               0,
                               sizeof(cl_mem),
                               (void*)&(input));
    if(err_code != CL_SUCCESS)
        return true;
    return false;
}

bool Reduction::setupOpenCL(const char* type,
                            const char* null_val,
                            const char* operation)
{
	InputOutput::ScreenManager *S = InputOutput::ScreenManager::singleton();
	CalcServer *C = CalcServer::singleton();
    size_t input_size, data_size;
    cl_int err_code;
    char msg[1024]; strcpy(msg, "");
    // Get the elements data size to can allocate local memory later
    err_code = clGetMemObjectInfo(_input,
                                  CL_MEM_SIZE,
                                  sizeof(size_t),
                                  &input_size,
                                  NULL);
    if(err_code != CL_SUCCESS)
        return true;
    data_size = input_size / _n.at(0);
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
    sprintf(args,
            "-DT=%s -DIDENTITY=%s -DLOCAL_WORK_SIZE=%luu",
            type,
            null_val,
            lsize);
    cl_kernel kernel;
    size_t maxlsize = loadKernelFromFile(&kernel,
                                         &_program,
                                         C->context,
                                         C->device,
                                         _path,
                                         "Reduction",
                                         args,
                                         header);
    if(maxlsize < __CL_MIN_LOCALSIZE__){
        S->addMessageF(3, "Reduction can't be performed due to insufficient local memory\n");
        sprintf(msg,
                "\t%lu elements can be executed, but __CL_MIN_LOCALSIZE__=%lu\n",
                maxlsize,
                __CL_MIN_LOCALSIZE__);
        S->addMessage(0, msg);
        return true;
    }
    lsize = maxlsize;
    if(!isPowerOf2(lsize)){
        lsize = nextPowerOf2(lsize) / 2;
    }
    if(lsize < __CL_MIN_LOCALSIZE__){
        S->addMessageF(3, "Reduction can't be performed due to insufficient local memory\n");
        sprintf(msg,
                "\t%lu elements may be executed, but __CL_MIN_LOCALSIZE__=%lu\n",
                lsize,
                __CL_MIN_LOCALSIZE__);
        S->addMessage(0, msg);
        return true;
    }
    sprintf(args,
            "-DT=%s -DIDENTITY=%s -DLOCAL_WORK_SIZE=%luu",
            type,
            null_val,
            lsize);
	if(kernel)clReleaseKernel(kernel); kernel=NULL;
	if(_program)clReleaseProgram(_program); _program=NULL;
    // Now we can start a loop while the amount of output data is greater than
    // one
    unsigned int n = _n.at(0);
    _n.clear();
    unsigned int i=0;
    while(n > 1){
        // Get work sizes
        _n.push_back(n);
        _local_work_sizes.push_back(lsize);
        _global_work_sizes.push_back(roundUp(n, lsize));
        _number_groups.push_back(
            _global_work_sizes.at(i) / _local_work_sizes.at(i)
        );
        // Build the output memory object
        cl_mem output = NULL;
        output = C->allocMemory(_number_groups.at(i) * data_size);
        if(!output){
            S->addMessageF(3, "Can't create an output array.\n");
            return true;
        }
        _mems.push_back(output);
        // Build the kernel
        if(!loadKernelFromFile(&kernel,
                               &_program,
                               C->context,
                               C->device,
                               _path,
                               "Reduction",
                               args,
                               header))
            return true;
        _kernels.push_back(kernel);
        if(_program)clReleaseProgram(_program); _program=NULL;

        err_code = CL_SUCCESS;
        err_code |= sendArgument(kernel,
                                 0,
                                 sizeof(cl_mem),
                                 (void*)&(_mems.at(i)));
        err_code |= sendArgument(kernel,
                                 1,
                                 sizeof(cl_mem),
                                 (void*)&(_mems.at(i+1)));
        err_code |= sendArgument(kernel,
                                 2,
                                 sizeof(cl_uint),
                                 (void*)&(n));
        err_code |= sendArgument(kernel,
                                 3,
                                 lsize*data_size ,
                                 NULL);
        if(err_code != CL_SUCCESS){
            S->addMessageF(3, "Arguments setting failed\n");
            return true;
        }
        // Setup next step
        sprintf(msg,
                "\tStep %u, %u elements reduced to %u\n",
                i,
                n,
                _number_groups.at(i));
        S->addMessage(0, msg);
        n = _number_groups.at(i);
        i++;
    }
	return false;
}

}}  // namespaces
