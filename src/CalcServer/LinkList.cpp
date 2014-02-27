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

#include <ProblemSetup.h>
#include <ScreenManager.h>
#include <CalcServer/LinkList.h>
#include <CalcServer.h>

namespace Aqua{ namespace CalcServer{

LinkList::LinkList()
	: Kernel("LinkList")
	, _path(0)
	, _program(0)
	, _icell_kernel(0)
	, _ihoc_kernel(0)
	, _ll_kernel(0)
	, _radix_sort(0)

{
	InputOutput::ProblemSetup *P  = InputOutput::ProblemSetup::singleton();
	InputOutput::ScreenManager *S = InputOutput::ScreenManager::singleton();
	CalcServer *C = CalcServer::singleton();
	int str_len;

	str_len = strlen(P->OpenCL_kernels.link_list);
	if(str_len <= 0) {
		S->addMessageF(3, "Path of Link-List kernel is empty.\n");
		exit(EXIT_FAILURE);
	}
	_path = new char[str_len+4];
	if(!_path) {
		S->addMessageF(3, "Memory cannot be allocated for the path.\n");
		exit(EXIT_FAILURE);
	}
	strcpy(_path, P->OpenCL_kernels.link_list);
	strcat(_path, ".cl");
	//! 2nd.- Setup the kernels
	_local_work_size  = localWorkSize();
	if(!_local_work_size){
	    S->addMessageF(3, "I cannot get a valid local work size for the required computation tool.\n");
	    exit(EXIT_FAILURE);
	}
	_global_work_size = globalWorkSize(_local_work_size);
	if(setupOpenCL()) {
	    exit(EXIT_FAILURE);
	}
	//! 3rd.- Built radix sort
	_radix_sort = new RadixSort();
	S->addMessageF(1, "LinkList ready to work!\n");
}

LinkList::~LinkList()
{
	InputOutput::ScreenManager *S = InputOutput::ScreenManager::singleton();
	if(_icell_kernel)clReleaseKernel(_icell_kernel); _icell_kernel=0;
	if(_ihoc_kernel)clReleaseKernel(_ihoc_kernel); _ihoc_kernel=0;
	if(_ll_kernel)clReleaseKernel(_ll_kernel); _ll_kernel=0;
	if(_program)clReleaseProgram(_program); _program=0;
	if(_path)delete[] _path; _path=0;
	S->addMessageF(1, "Destroying radix sort processor...\n");
	if(_radix_sort)delete _radix_sort; _radix_sort=0;
}

bool LinkList::execute()
{
	InputOutput::ScreenManager *S = InputOutput::ScreenManager::singleton();
	CalcServer *C = CalcServer::singleton();
	cl_int err_code=0;

	if(allocLinkList())
	     return true;

    // Clear the previous ihoc
    // =======================
	err_code |= sendArgument(_ihoc_kernel,
                             0,
                             sizeof(cl_mem),
                             (void*)&(C->ihoc));
	err_code |= sendArgument(_ihoc_kernel,
                             1,
                             sizeof(cl_uint),
                             (void*)&(C->num_cells_allocated));
	err_code |= sendArgument(_ihoc_kernel,
                             2,
                             sizeof(cl_uint),
                             (void*)&(C->N));
	if(err_code != CL_SUCCESS) {
		S->addMessageF(3, "Failure sending variables to the ihoc clearer kernel.\n");
	    return true;
	}
	_local_work_size  = localWorkSize(C->num_cells);
	if(!_local_work_size){
	    S->addMessageF(3, "No valid local work size has been found for the ihoc clearer.\n");
	    return true;
	}
	_global_work_size = globalWorkSize(_local_work_size, C->num_cells);
	#ifdef HAVE_GPUPROFILE
	    cl_event event;
	    cl_ulong end, start;
	    profileTime(0.f);
	    err_code = clEnqueueNDRangeKernel(C->command_queue,
                                          _ihoc_kernel,
                                          1,
                                          NULL,
                                          &_global_work_size,
                                          NULL,
                                          0,
                                          NULL,
                                          &event);
	#else
	    err_code = clEnqueueNDRangeKernel(C->command_queue,
                                          _ihoc_kernel,
                                          1,
                                          NULL,
                                          &_global_work_size,
                                          NULL,
                                          0,
                                          NULL,
                                          NULL);
	#endif
	if(err_code != CL_SUCCESS) {
		S->addMessageF(3, "I cannot execute the kernel \"InitIhoc\".\n");
        S->printOpenCLError(err_code);
	    return true;
	}
	#ifdef HAVE_GPUPROFILE
	    err_code = clWaitForEvents(1, &event);
		if(err_code != CL_SUCCESS) {
			S->addMessageF(3, "Impossible to wait for the ihoc clearer kernel end.\n");
            S->printOpenCLError(err_code);
			return true;
		}
	    err_code = clGetEventProfilingInfo(event,
                                           CL_PROFILING_COMMAND_END,
                                           sizeof(cl_ulong),
                                           &end,
                                           0);
		if(err_code != CL_SUCCESS) {
			S->addMessageF(3, "I cannot profile the ihoc clearer kernel execution.\n");
            S->printOpenCLError(err_code);
			return true;
		}
	    err_code = clGetEventProfilingInfo(event,
                                           CL_PROFILING_COMMAND_START,
                                           sizeof(cl_ulong),
                                           &start, 0);
		if(err_code != CL_SUCCESS) {
			S->addMessageF(3, "I cannot profile the ihoc clearer kernel execution.\n");
            S->printOpenCLError(err_code);
			return true;
		}
	    profileTime(profileTime() + (end - start)/1000.f);  // 10^-3 ms
	#endif

    // Compute the particles cells
    // ===========================
	err_code |= sendArgument(_icell_kernel,
                             0,
                             sizeof(cl_mem),
                             (void*)&C->icell);
	err_code |= sendArgument(_icell_kernel,
                             1,
                             sizeof(cl_mem),
                             (void*)&C->pos);
	err_code |= sendArgument(_icell_kernel,
                             2,
                             sizeof(cl_uint),
                             (void*)&C->N);
	err_code |= sendArgument(_icell_kernel,
                             3,
                             sizeof(cl_uint),
                             (void*)&C->num_icell);
	err_code |= sendArgument(_icell_kernel,
                             4,
                             sizeof(vec),
                             (void*)&C->pos_min);
	err_code |= sendArgument(_icell_kernel,
                             5,
                             sizeof(cl_float),
                             (void*)&C->cell_length);
	err_code |= sendArgument(_icell_kernel,
                             6,
                             sizeof(cl_uint),
                             (void*)&C->num_cells);
	err_code |= sendArgument(_icell_kernel,
                             7,
                             sizeof(uivec),
                             (void*)&C->num_cells_vec);
	if(err_code != CL_SUCCESS) {
		S->addMessageF(3, "Failure sending variables to the cells localization kernel.\n");
	    return true;
	}
	_local_work_size  = localWorkSize(C->num_icell);
	if(!_local_work_size){
	    S->addMessageF(3, "No valid local work size has been found for the cells localization.\n");
	    return true;
	}
	_global_work_size = globalWorkSize(_local_work_size, C->num_icell);
	#ifdef HAVE_GPUPROFILE
	    err_code = clEnqueueNDRangeKernel(C->command_queue,
                                          _icell_kernel,
                                          1,
                                          NULL,
                                          &_global_work_size,
                                          &_local_work_size,
                                          0,
                                          NULL,
                                          &event);
	#else
	    err_code = clEnqueueNDRangeKernel(C->command_queue,
                                          _icell_kernel,
                                          1,
                                          NULL,
                                          &_global_work_size,
                                          &_local_work_size,
                                          0,
                                          NULL,
                                          NULL);
	#endif
	if(err_code != CL_SUCCESS) {
		S->addMessageF(3, "I cannot execute the kernel \"LCell\".\n");
        S->printOpenCLError(err_code);
	    return true;
	}
	#ifdef HAVE_GPUPROFILE
	    err_code = clWaitForEvents(1, &event);
		if(err_code != CL_SUCCESS) {
			S->addMessageF(3, "Impossible to wait for the cells localization kernel end.\n");
            S->printOpenCLError(err_code);
			return true;
		}
	    err_code = clGetEventProfilingInfo(event,
                                           CL_PROFILING_COMMAND_END,
                                           sizeof(cl_ulong),
                                           &end,
                                           0);
		if(err_code != CL_SUCCESS) {
			S->addMessageF(3, "I cannot profile the cells localization kernel execution.\n");
            S->printOpenCLError(err_code);
			return true;
		}
	    err_code = clGetEventProfilingInfo(event,
                                           CL_PROFILING_COMMAND_START,
                                           sizeof(cl_ulong),
                                           &start, 0);
		if(err_code != CL_SUCCESS) {
			S->addMessageF(3, "I cannot profile the cells localization kernel execution.\n");
            S->printOpenCLError(err_code);
			return true;
		}
	    profileTime(profileTime() + (end - start)/1000.f);  // 10^-3 ms
	#endif

    // Sort the particles by their cells
    // =================================
	if(_radix_sort->sort())
	    return true;
	#ifdef HAVE_GPUPROFILE
	    profileTime(profileTime() + _radix_sort->profileTime());  // 10^-3 ms
	#endif

    // Compute the heads of chain for the sorted particles
    // ===================================================
	unsigned int nMinOne = C->N-1;
	err_code |= sendArgument(_ll_kernel,
                             0,
                             sizeof(cl_mem),
                             (void*)&(C->icell));
	err_code |= sendArgument(_ll_kernel,
                             1,
                             sizeof(cl_mem),
                             (void*)&(C->ihoc));
	err_code |= sendArgument(_ll_kernel,
                             2,
                             sizeof(cl_uint),
                             (void*)&nMinOne);
	if(err_code != CL_SUCCESS) {
		S->addMessageF(3, "Failure sending variables to the Link-List kernel.\n");
	    return true;
	}
	_local_work_size  = localWorkSize(nMinOne);
	if(!_local_work_size){
	    S->addMessageF(3, "No valid local work size has been found for the Link-List.\n");
	    return true;
	}
	_global_work_size = globalWorkSize(_local_work_size, nMinOne);
	#ifdef HAVE_GPUPROFILE
	    err_code = clEnqueueNDRangeKernel(C->command_queue,
                                          _ll_kernel,
                                          1,
                                          NULL,
                                          &_global_work_size,
                                          NULL,
                                          0,
                                          NULL,
                                          &event);
	#else
	    err_code = clEnqueueNDRangeKernel(C->command_queue,
                                          _ll_kernel,
                                          1,
                                          NULL,
                                          &_global_work_size,
                                          NULL,
                                          0,
                                          NULL,
                                          NULL);
	#endif
	if(err_code != CL_SUCCESS) {
		S->addMessageF(3, "I cannot execute the kernel \"LinkList\".\n");
        S->printOpenCLError(err_code);
	    return true;
	}
	#ifdef HAVE_GPUPROFILE
	    err_code = clWaitForEvents(1, &event);
		if(err_code != CL_SUCCESS) {
			S->addMessageF(3, "Impossible to wait for the Link-List kernel end.\n");
            S->printOpenCLError(err_code);
			return true;
		}
	    err_code = clGetEventProfilingInfo(event,
                                           CL_PROFILING_COMMAND_END,
                                           sizeof(cl_ulong),
                                           &end,
                                           0);
		if(err_code != CL_SUCCESS) {
			S->addMessageF(3, "I cannot profile the Link-List kernel execution.\n");
            S->printOpenCLError(err_code);
			return true;
		}
	    err_code = clGetEventProfilingInfo(event,
                                           CL_PROFILING_COMMAND_START,
                                           sizeof(cl_ulong),
                                           &start, 0);
		if(err_code != CL_SUCCESS) {
			S->addMessageF(3, "I cannot profile the Link-List kernel execution.\n");
            S->printOpenCLError(err_code);
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
	int err_code=0;
	if(!loadKernelFromFile(&_icell_kernel,
                           &_program,
                           C->context,
                           C->device,
                           _path,
                           "LCell",
                           ""))
	    return true;
	if(_program)clReleaseProgram(_program); _program=0;
	if(!loadKernelFromFile(&_ihoc_kernel,
                           &_program,
                           C->context,
                           C->device,
                           _path,
                           "InitIhoc",
                           ""))
	    return true;
	if(_program)clReleaseProgram(_program); _program=0;
	if(!loadKernelFromFile(&_ll_kernel,
                           &_program,
                           C->context,
                           C->device,
                           _path,
                           "LinkList",
                           ""))
	    return true;
	if(_program)clReleaseProgram(_program); _program=0;
	err_code = sendArgument(_icell_kernel,
                            0,
                            sizeof(cl_mem),
                            (void*)&(C->pos));
	if(err_code)
	    return true;

	cl_device_id device;
	size_t local_work_size=0;
	err_code |= clGetCommandQueueInfo(C->command_queue,
                                      CL_QUEUE_DEVICE,
	                                  sizeof(cl_device_id),
                                      &device,
                                      NULL);
	if(err_code != CL_SUCCESS) {
		S->addMessageF(3, "I Cannot get the device from the command queue.\n");
        S->printOpenCLError(err_code);
	    return true;
	}
	err_code |= clGetKernelWorkGroupInfo(_icell_kernel,
                                         device,
                                         CL_KERNEL_WORK_GROUP_SIZE,
	                                     sizeof(size_t),
                                         &local_work_size,
                                         NULL);
	if(err_code != CL_SUCCESS) {
		S->addMessageF(3, "Failure computing the maximum local work size for Lcell.\n");
        S->printOpenCLError(err_code);
	    return true;
	}
	if(local_work_size < _local_work_size)
	    _local_work_size  = local_work_size;
	err_code |= clGetKernelWorkGroupInfo(_ihoc_kernel,
                                         device,
                                         CL_KERNEL_WORK_GROUP_SIZE,
	                                     sizeof(size_t),
                                         &local_work_size,
                                         NULL);
	if(err_code != CL_SUCCESS) {
		S->addMessageF(3, "Failure computing the maximum local work size for ihoc clearer.\n");
        S->printOpenCLError(err_code);
	    return true;
	}
	if(local_work_size < _local_work_size)
	    _local_work_size  = local_work_size;
	err_code |= clGetKernelWorkGroupInfo(_ll_kernel,
                                         device,
                                         CL_KERNEL_WORK_GROUP_SIZE,
	                                     sizeof(size_t),
                                         &local_work_size, NULL);
	if(err_code != CL_SUCCESS) {
		S->addMessageF(3, "Failure computing the maximum local work size for Link-List.\n");
        S->printOpenCLError(err_code);
	    return true;
	}
	if(local_work_size < _local_work_size)
	    _local_work_size  = local_work_size;
	_global_work_size = globalWorkSize(_local_work_size);
	return false;
}

bool LinkList::allocLinkList()
{
	InputOutput::ScreenManager *S = InputOutput::ScreenManager::singleton();
	CalcServer *C = CalcServer::singleton();
	char msg[256];
	if(C->num_cells > C->num_cells_allocated)
	{
	    sprintf(msg,
                "Number of cells increased [%d -> %d]\n",
                C->num_cells_allocated,
                C->num_cells);
		S->addMessageF(2, msg);
	    if(C->num_cells_allocated > 0) {
	        if(C->ihoc)
                clReleaseMemObject(C->ihoc);
            C->ihoc=NULL;
			C->allocated_mem -= C->num_cells_allocated * sizeof( cl_uint );
	    }
	    C->ihoc = C->allocMemory(C->num_cells * sizeof( cl_uint ));
		if(!C->ihoc) {
	        sprintf(msg,
                    "Fail allocating memory for ihoc (%u bytes).\n",
                    (unsigned int)(C->num_cells * sizeof( cl_uint )) );
		    S->addMessageF(3, msg);
		    return true;
		}
		sprintf(msg, "\tAllocated memory = %u bytes\n",
                (unsigned int)C->allocated_mem);
		S->addMessage(1, msg);
	    C->num_cells_allocated = C->num_cells;
	    _local_work_size = 256;
	    _global_work_size = globalWorkSize(_local_work_size);
	}
	return false;
}

}}  // namespace
