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
	, _path(0)
	, program(0)
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
	_path = new char[nChar+4];
	if(!_path) {
	    printf("ERROR (LinkList::Init): Can't allocate memory for path.\n");
	    exit(EXIT_FAILURE);
	}
	strcpy(_path, P->OpenCL_kernels.link_list);
	strcat(_path, ".cl");
	//! 2nd.- Setup the kernels
	_local_work_size  = localWorkSize();
	if(!_local_work_size){
	    printf("ERROR (LinkList::Init): I cannot get a valid local work size for the required computation tool.\n");
	    exit(EXIT_FAILURE);
	}
	_global_work_size = globalWorkSize(_local_work_size);
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
	if(program)clReleaseProgram(program); program=0;
	if(_path)delete[] _path; _path=0;
	S->addMessage(1, "(LinkList::~LinkList): Destroying radix sort processor...\n");
	if(mRadixSort)delete mRadixSort; mRadixSort=0;
}

bool LinkList::execute()
{
	InputOutput::ScreenManager *S = InputOutput::ScreenManager::singleton();
	CalcServer *C = CalcServer::singleton();
	cl_int err_code=0;
	//! Allocate memory if needed
	if(allocLinkList())
	     return true;
	//! Set initial values for cells
	err_code |= sendArgument(clIhocKernel, 0, sizeof(cl_mem), (void*)&(C->ihoc));
	err_code |= sendArgument(clIhocKernel, 1, sizeof(cl_mem), (void*)&(C->cell_has_particles));
	err_code |= sendArgument(clIhocKernel, 2, sizeof(cl_uint), (void*)&(C->num_cells_allocated));
	err_code |= sendArgument(clIhocKernel, 3, sizeof(cl_uint), (void*)&(C->N));
	if(err_code != CL_SUCCESS) {
		S->addMessage(3, (char *)"(LinkList::Execute): Can't send arguments to to ihoc clearer kernel.\n");
	    return true;
	}
	_local_work_size  = localWorkSize(C->num_cells);
	if(!_local_work_size){
	    S->addMessage(3, (char*)"(LinkList::Execute): No valid local work size for ihoc clearer.\n");
	    return true;
	}
	_global_work_size = globalWorkSize(_local_work_size, C->num_cells);
	#ifdef HAVE_GPUPROFILE
	    cl_event event;
	    cl_ulong end, start;
	    profileTime(0.f);
	    err_code = clEnqueueNDRangeKernel(C->command_queue, clIhocKernel, 1, NULL, &_global_work_size, NULL, 0, NULL, &event);
	#else
	    err_code = clEnqueueNDRangeKernel(C->command_queue, clIhocKernel, 1, NULL, &_global_work_size, NULL, 0, NULL, NULL);
	#endif
	if(err_code != CL_SUCCESS) {
		S->addMessage(3, "(LinkList::Execute): Can't execute the kernel \"InitIhoc\".\n");
	    if(err_code == CL_INVALID_PROGRAM_EXECUTABLE)
	        S->addMessage(0, "\tInvalid program (Compile errors maybe?).\n");
	    else if(err_code == CL_INVALID_COMMAND_QUEUE)
	        S->addMessage(0, "\tInvalid command queue.\n");
	    else if(err_code == CL_INVALID_KERNEL)
	        S->addMessage(0, "\tKernel is not a valid object.\n");
	    else if(err_code == CL_INVALID_CONTEXT)
	        S->addMessage(0, "\tContext associated to command queue don't match qith kernel context.\n");
	    else if(err_code == CL_INVALID_KERNEL_ARGS)
	        S->addMessage(0, "\tOne or more arguments are invalid (maybe don't specified).\n");
	    else if(err_code == CL_INVALID_WORK_DIMENSION)
	        S->addMessage(0, "\tDimension must be a value between 1 and 3.\n");
	    else if(err_code == CL_INVALID_WORK_GROUP_SIZE)
	        S->addMessage(0, "\tInvalid local work group size.\n");
	    else if(err_code == CL_INVALID_WORK_ITEM_SIZE)
	        S->addMessage(0, "\tLocal work group size is out of bounds.\n");
	    else if(err_code == CL_INVALID_GLOBAL_OFFSET)
	        S->addMessage(0, "\tGlobal offset must be NULL.\n");
	    else if(err_code == CL_OUT_OF_RESOURCES)
	        S->addMessage(0, "\tDevice out of resources.\n");
	    else if(err_code == CL_MEM_OBJECT_ALLOCATION_FAILURE)
	        S->addMessage(0, "\tAllocation error at device.\n");
	    else if(err_code == CL_INVALID_EVENT_WAIT_LIST)
	        S->addMessage(0, "\tInvalid event wait instruction.\n");
	    else if(err_code == CL_OUT_OF_HOST_MEMORY)
	        S->addMessage(0, "\tfailure to allocate resources required by the OpenCL implementation on the host.\n");
	    return true;
	}
	#ifdef HAVE_GPUPROFILE
	    err_code = clWaitForEvents(1, &event);
	    if(err_code != CL_SUCCESS) {
	        S->addMessage(3, "(LinkList::Execute): Can't wait to kernels end.\n");
	        return true;
	    }
	    err_code |= clGetEventProfilingInfo(event, CL_PROFILING_COMMAND_END, sizeof(cl_ulong), &end, 0);
	    err_code |= clGetEventProfilingInfo(event, CL_PROFILING_COMMAND_START, sizeof(cl_ulong), &start, 0);
	    if(err_code != CL_SUCCESS) {
	        S->addMessage(3, "(LinkList::Execute): Can't profile kernel execution.\n");
	        return true;
	    }
	    profileTime(profileTime() + (end - start)/1000.f);  // 10^-3 ms
	#endif
	//! Find each particle cell
	err_code |= sendArgument(clLcellKernel, 0, sizeof(cl_mem  ), (void*)&C->icell);
	err_code |= sendArgument(clLcellKernel, 1, sizeof(cl_mem  ), (void*)&C->cell_has_particles);
	err_code |= sendArgument(clLcellKernel, 2, sizeof(cl_mem  ), (void*)&C->imove);
	err_code |= sendArgument(clLcellKernel, 3, sizeof(cl_mem  ), (void*)&C->pos);
	err_code |= sendArgument(clLcellKernel, 4, sizeof(cl_uint ), (void*)&C->N);
	err_code |= sendArgument(clLcellKernel, 5, sizeof(cl_uint ), (void*)&C->num_icell);
	err_code |= sendArgument(clLcellKernel, 6, sizeof(vec     ), (void*)&C->pos_min);
	err_code |= sendArgument(clLcellKernel, 7, sizeof(cl_float), (void*)&C->cell_length);
	err_code |= sendArgument(clLcellKernel, 8, sizeof(cl_uint ), (void*)&C->num_cells);
	err_code |= sendArgument(clLcellKernel, 9, sizeof(uivec   ), (void*)&C->num_cells_vec);
	if(err_code != CL_SUCCESS) {
		S->addMessage(3, "(LinkList::Execute): Can't send arguments to cell allocation kernel.\n");
	    return true;
	}
	_local_work_size  = localWorkSize(C->num_icell);
	if(!_local_work_size){
	    S->addMessage(3, "(LinkList::Execute): No valid local work size for cell allocation.\n");
	    return true;
	}
	_global_work_size = globalWorkSize(_local_work_size, C->num_icell);
	#ifdef HAVE_GPUPROFILE
	    err_code = clEnqueueNDRangeKernel(C->command_queue, clLcellKernel, 1, NULL, &_global_work_size, &_local_work_size, 0, NULL, &event);
	#else
	    err_code = clEnqueueNDRangeKernel(C->command_queue, clLcellKernel, 1, NULL, &_global_work_size, &_local_work_size, 0, NULL, NULL);
	#endif
	if(err_code != CL_SUCCESS) {
		S->addMessage(3, "(LinkList::Execute): Can't execute the kernel \"LCell\".\n");
	    if(err_code == CL_INVALID_WORK_GROUP_SIZE)
	        S->addMessage(0, "\tInvalid local work group size.\n");
	    else if(err_code == CL_OUT_OF_RESOURCES)
	        S->addMessage(0, "\tDevice out of resources.\n");
	    else if(err_code == CL_MEM_OBJECT_ALLOCATION_FAILURE)
	        S->addMessage(0, "\tAllocation error at device.\n");
	    else if(err_code == CL_OUT_OF_HOST_MEMORY)
	        S->addMessage(0, "\tfailure to allocate resources required by the OpenCL implementation on the host.\n");
	    return true;
	}
	#ifdef HAVE_GPUPROFILE
	    err_code = clWaitForEvents(1, &event);
	    if(err_code != CL_SUCCESS) {
	        S->addMessage(3, "(LinkList::Execute): Can't wait to kernels end.\n");
	        return true;
	    }
	    err_code |= clGetEventProfilingInfo(event, CL_PROFILING_COMMAND_END, sizeof(cl_ulong), &end, 0);
	    err_code |= clGetEventProfilingInfo(event, CL_PROFILING_COMMAND_START, sizeof(cl_ulong), &start, 0);
	    if(err_code != CL_SUCCESS) {
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
	err_code |= sendArgument(clLLKernel, 0, sizeof(cl_mem), (void*)&(C->icell));
	err_code |= sendArgument(clLLKernel, 1, sizeof(cl_mem), (void*)&(C->ihoc));
	err_code |= sendArgument(clLLKernel, 2, sizeof(cl_uint), (void*)&nMinOne);
	if(err_code != CL_SUCCESS) {
		S->addMessage(3, "(LinkList::Execute): Can't send arguments to ihoc perform kernel.\n");
	    return true;
	}
	_local_work_size  = localWorkSize(nMinOne);
	if(!_local_work_size){
	    S->addMessage(3, "(LinkList::Execute): No valid local work size for ihoc computation.\n");
	    return true;
	}
	_global_work_size = globalWorkSize(_local_work_size, nMinOne);
	#ifdef HAVE_GPUPROFILE
	    err_code = clEnqueueNDRangeKernel(C->command_queue, clLLKernel, 1, NULL, &_global_work_size, NULL, 0, NULL, &event);
	#else
	    err_code = clEnqueueNDRangeKernel(C->command_queue, clLLKernel, 1, NULL, &_global_work_size, NULL, 0, NULL, NULL);
	#endif
	if(err_code != CL_SUCCESS) {
		S->addMessage(3, "(LinkList::Execute): Can't execute the kernel \"LinkList\".\n");
	    if(err_code == CL_INVALID_WORK_GROUP_SIZE)
	        S->addMessage(0, "\tInvalid local work group size.\n");
	    else if(err_code == CL_OUT_OF_RESOURCES)
	        S->addMessage(0, "\tDevice out of resources.\n");
	    else if(err_code == CL_MEM_OBJECT_ALLOCATION_FAILURE)
	        S->addMessage(0, "\tAllocation error at device.\n");
	    else if(err_code == CL_OUT_OF_HOST_MEMORY)
	        S->addMessage(0, "\tfailure to allocate resources required by the OpenCL implementation on the host.\n");
	    return true;
	}
	#ifdef HAVE_GPUPROFILE
	    err_code = clWaitForEvents(1, &event);
	    if(err_code != CL_SUCCESS) {
	        S->addMessage(3, "(LinkList::Execute): Can't wait to kernels end.\n");
	        return true;
	    }
	    err_code |= clGetEventProfilingInfo(event, CL_PROFILING_COMMAND_END, sizeof(cl_ulong), &end, 0);
	    err_code |= clGetEventProfilingInfo(event, CL_PROFILING_COMMAND_START, sizeof(cl_ulong), &start, 0);
	    if(err_code != CL_SUCCESS) {
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
	int err_code=0;
	if(!loadKernelFromFile(&clLcellKernel, &program, C->context, C->device, _path, "LCell", ""))
	    return true;
	if(program)clReleaseProgram(program); program=0;
	if(!loadKernelFromFile(&clIhocKernel, &program, C->context, C->device, _path, "InitIhoc", ""))
	    return true;
	if(program)clReleaseProgram(program); program=0;
	if(!loadKernelFromFile(&clLLKernel, &program, C->context, C->device, _path, "LinkList", ""))
	    return true;
	if(program)clReleaseProgram(program); program=0;
	err_code |= sendArgument(clLcellKernel,  0, sizeof(cl_mem), (void*)&(C->pos));
	if(err_code)
	    return true;
	//! Test for right work group size
	cl_device_id device;
	size_t localWorkGroupSize=0;
	err_code |= clGetCommandQueueInfo(C->command_queue,CL_QUEUE_DEVICE,
	                                sizeof(cl_device_id),&device, NULL);
	if(err_code != CL_SUCCESS) {
		S->addMessage(3, "(LinkList::setupOpenCL): I Cannot get the device from the command queue.\n");
	    return true;
	}
	err_code |= clGetKernelWorkGroupInfo(clLcellKernel,device,CL_KERNEL_WORK_GROUP_SIZE,
	                                   sizeof(size_t), &localWorkGroupSize, NULL);
	if(err_code != CL_SUCCESS) {
		S->addMessage(3, "(LinkList::setupOpenCL): Can't get Lcell maximum local work group size.\n");
	    return true;
	}
	if(localWorkGroupSize < _local_work_size)
	    _local_work_size  = localWorkGroupSize;
	err_code |= clGetKernelWorkGroupInfo(clIhocKernel,device,CL_KERNEL_WORK_GROUP_SIZE,
	                                   sizeof(size_t), &localWorkGroupSize, NULL);
	if(err_code != CL_SUCCESS) {
		S->addMessage(3, "(LinkList::setupOpenCL): Can't get ihoc maximum local work group size.\n");
	    return true;
	}
	if(localWorkGroupSize < _local_work_size)
	    _local_work_size  = localWorkGroupSize;
	err_code |= clGetKernelWorkGroupInfo(clLLKernel,device,CL_KERNEL_WORK_GROUP_SIZE,
	                                   sizeof(size_t), &localWorkGroupSize, NULL);
	if(err_code != CL_SUCCESS) {
		S->addMessage(3, "(LinkList::setupOpenCL): Can't get ll maximum local work group size.\n");
	    return true;
	}
	if(localWorkGroupSize < _local_work_size)
	    _local_work_size  = localWorkGroupSize;
	_global_work_size = globalWorkSize(_local_work_size);
	return false;
}

bool LinkList::allocLinkList()
{
	InputOutput::ScreenManager *S = InputOutput::ScreenManager::singleton();
	CalcServer *C = CalcServer::singleton();
	char Log[256];
	if(C->num_cells > C->num_cells_allocated)
	{
	    sprintf(Log, "Number of cells increased [%d -> %d]\n",
                C->num_cells_allocated,C->num_cells);
		S->addMessageF(2, Log);
	    if(C->num_cells_allocated > 0) {
	        if(C->ihoc)clReleaseMemObject(C->ihoc); C->ihoc=NULL;
	        if(C->cell_has_particles)clReleaseMemObject(C->cell_has_particles); C->cell_has_particles=NULL;
			C->allocated_mem -= C->num_cells_allocated * sizeof( cl_uint );
			C->allocated_mem -= C->num_cells_allocated * sizeof( cl_short );
	    }
	    C->ihoc = C->allocMemory(C->num_cells * sizeof( cl_uint ));
		if(!C->ihoc) {
	        sprintf(Log, "Fail allocating memory for ihoc (%u bytes).\n",
                    (unsigned int)(C->num_cells * sizeof( cl_uint )) );
		    S->addMessageF(3, Log);
		    return true;
		}
		C->cell_has_particles = C->allocMemory(C->num_cells * sizeof( cl_short ));
		if(!C->cell_has_particles) {
	        sprintf(Log, "Fail allocating memory for cell_has_particles (%u bytes).\n",
                    (unsigned int)(C->num_cells * sizeof( cl_short )) );
		    S->addMessageF(3, Log);
		    return true;
		}
		sprintf(Log, "\tAllocated memory = %u bytes\n",
                (unsigned int)C->allocated_mem);
		S->addMessage(1, Log);
	    C->num_cells_allocated = C->num_cells;
	    _local_work_size = 256;
	    _global_work_size = globalWorkSize(_local_work_size);
	}
	return false;
}

}}  // namespace
