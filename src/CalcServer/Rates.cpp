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
	, _path(0)
	, program(0)
	, kernel(0)
	, clSortKernel(0)
	, _global_work_size(0)
	, _local_work_size(0)
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
	_path = new char[nChar+4];
	if(!_path) {
	    S->addMessage(3, "(Rates::Rates): Can't allocate memory for path.\n");
	    exit(EXIT_FAILURE);
	}
	strcpy(_path, P->OpenCL_kernels.rates);
	strcat(_path, ".cl");
    for(i=0;i<P->n_fluids;i++){
        if(P->fluids[i].delta > 0.f){
            isDelta = true;
            break;
        }
    }

	_local_work_size  = localWorkSize();
	if(!_local_work_size){
	    S->addMessage(3, "(Rates::Rates): I cannot get a valid local work size for the required computation tool.\n");
	    exit(EXIT_FAILURE);
	}
	_global_work_size = globalWorkSize(_local_work_size);
	if(setupOpenCL()) {
	    exit(EXIT_FAILURE);
	}
	S->addMessage(1, "(Rates::Rates): Rates ready to work!\n");
}

Rates::~Rates()
{
	if(kernel)clReleaseKernel(kernel); kernel=0;
	if(clSortKernel)clReleaseKernel(clSortKernel); clSortKernel=0;
	if(program)clReleaseProgram(program); program=0;
	if(_path)delete[] _path; _path=0;
}

bool Rates::execute()
{
	InputOutput::ScreenManager *S = InputOutput::ScreenManager::singleton();
	CalcServer *C = CalcServer::singleton();
	cl_int err_code=0;
	// Sort the rates kernel input data
	err_code |= sendArgument(clSortKernel,  0, sizeof(cl_mem), (void*)&(C->ifluid));
	err_code |= sendArgument(clSortKernel,  1, sizeof(cl_mem), (void*)&(C->ifluidin));
	err_code |= sendArgument(clSortKernel,  2, sizeof(cl_mem), (void*)&(C->imove));
	err_code |= sendArgument(clSortKernel,  3, sizeof(cl_mem), (void*)&(C->imovein));
	err_code |= sendArgument(clSortKernel,  4, sizeof(cl_mem), (void*)&(C->pos));
	err_code |= sendArgument(clSortKernel,  5, sizeof(cl_mem), (void*)&(C->v));
	err_code |= sendArgument(clSortKernel,  6, sizeof(cl_mem), (void*)&(C->hp));
	err_code |= sendArgument(clSortKernel,  7, sizeof(cl_mem), (void*)&(C->dens));
	err_code |= sendArgument(clSortKernel,  8, sizeof(cl_mem), (void*)&(C->press));
	err_code |= sendArgument(clSortKernel,  9, sizeof(cl_mem), (void*)&(C->mass));
	err_code |= sendArgument(clSortKernel, 10, sizeof(cl_mem), (void*)&(C->posin));
	err_code |= sendArgument(clSortKernel, 11, sizeof(cl_mem), (void*)&(C->vin));
	err_code |= sendArgument(clSortKernel, 12, sizeof(cl_mem), (void*)&(C->hpin));
	err_code |= sendArgument(clSortKernel, 13, sizeof(cl_mem), (void*)&(C->densin));
	err_code |= sendArgument(clSortKernel, 14, sizeof(cl_mem), (void*)&(C->pressin));
	err_code |= sendArgument(clSortKernel, 15, sizeof(cl_mem), (void*)&(C->massin));
	err_code |= sendArgument(clSortKernel, 16, sizeof(cl_mem), (void*)&(C->permutation));
	err_code |= sendArgument(clSortKernel, 17, sizeof(cl_mem), (void*)&(C->permutation_inverse));
	err_code |= sendArgument(clSortKernel, 18, sizeof(cl_uint),  (void*)&(C->N));
	if(err_code != CL_SUCCESS) {
		S->addMessage(3, "(Rates::execute): Can't send arguments to sort kernel.\n");
	    return true;
	}
	#ifdef HAVE_GPUPROFILE
	    cl_event event;
	    cl_ulong end, start;
	    profileTime(0.f);
	    err_code = clEnqueueNDRangeKernel(C->command_queue, clSortKernel, 1, NULL, &_global_work_size, NULL, 0, NULL, &event);
	#else
	    err_code = clEnqueueNDRangeKernel(C->command_queue, clSortKernel, 1, NULL, &_global_work_size, NULL, 0, NULL, NULL);
	#endif
	if(err_code != CL_SUCCESS) {
		S->addMessage(3, "(Rates::execute): Can't execute the sorting kernel.\n");
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
	        S->addMessage(3, "(Rates::execute): Can't wait to sorting kernel ends.\n");
	        return true;
	    }
	    err_code |= clGetEventProfilingInfo(event, CL_PROFILING_COMMAND_END, sizeof(cl_ulong), &end, 0);
	    err_code |= clGetEventProfilingInfo(event, CL_PROFILING_COMMAND_START, sizeof(cl_ulong), &start, 0);
	    if(err_code != CL_SUCCESS) {
	        S->addMessage(3, "(Rates::execute): Can't profile sorting kernel execution.\n");
	        return true;
	    }
	    profileTime(profileTime() + (end - start)/1000.f);  // 10^-3 ms
	#endif
	// Compute the variation rates
	err_code |= sendArgument(kernel,  0, sizeof(cl_mem  ), (void*)&(C->ifluidin));
	err_code |= sendArgument(kernel,  1, sizeof(cl_mem  ), (void*)&(C->imovein));
	err_code |= sendArgument(kernel,  2, sizeof(cl_mem  ), (void*)&(C->posin));
	err_code |= sendArgument(kernel,  3, sizeof(cl_mem  ), (void*)&(C->vin));
	err_code |= sendArgument(kernel,  4, sizeof(cl_mem  ), (void*)&(C->densin));
	err_code |= sendArgument(kernel,  5, sizeof(cl_mem  ), (void*)&(C->hpin));
	err_code |= sendArgument(kernel,  6, sizeof(cl_mem  ), (void*)&(C->massin));
	err_code |= sendArgument(kernel,  7, sizeof(cl_mem  ), (void*)&(C->pressin));
	err_code |= sendArgument(kernel,  8, sizeof(cl_mem  ), (void*)&(C->visc_kin));
	err_code |= sendArgument(kernel,  9, sizeof(cl_mem  ), (void*)&(C->visc_dyn_corrected));
	err_code |= sendArgument(kernel, 10, sizeof(cl_mem  ), (void*)&(C->f));
	err_code |= sendArgument(kernel, 11, sizeof(cl_mem  ), (void*)&(C->drdt));
	err_code |= sendArgument(kernel, 12, sizeof(cl_mem  ), (void*)&(C->drdt_F));
	err_code |= sendArgument(kernel, 13, sizeof(cl_mem  ), (void*)&(C->sigma));
	err_code |= sendArgument(kernel, 14, sizeof(cl_mem  ), (void*)&(C->shepard));
	err_code |= sendArgument(kernel, 15, sizeof(cl_mem  ), (void*)&(C->shepard_gradient));
	err_code |= sendArgument(kernel, 16, sizeof(cl_mem  ), (void*)&(C->icell));
	err_code |= sendArgument(kernel, 17, sizeof(cl_mem  ), (void*)&(C->ihoc));
	err_code |= sendArgument(kernel, 18, sizeof(cl_mem  ), (void*)&(C->cell_has_particles));
	err_code |= sendArgument(kernel, 19, sizeof(cl_mem  ), (void*)&(C->permutation));
	err_code |= sendArgument(kernel, 20, sizeof(cl_mem  ), (void*)&(C->permutation_inverse));
	err_code |= sendArgument(kernel, 21, sizeof(cl_mem  ), (void*)&(C->sensor_mode));
	err_code |= sendArgument(kernel, 22, sizeof(cl_uint ), (void*)&(C->n));
	err_code |= sendArgument(kernel, 23, sizeof(cl_uint ), (void*)&(C->N));
	err_code |= sendArgument(kernel, 24, sizeof(uivec   ), (void*)&(C->num_cells_vec));
	err_code |= sendArgument(kernel, 25, sizeof(vec     ), (void*)&(C->g));
	unsigned int nAddedArgs = 0;
	if(isDelta) {
	    err_code |= sendArgument(kernel, 26, sizeof(cl_mem), (void*)&(C->refd));
	    err_code |= sendArgument(kernel, 27, sizeof(cl_mem), (void*)&(C->delta));
	    err_code |= sendArgument(kernel, 28, sizeof(cl_float), (void*)&(C->dt));
	    err_code |= sendArgument(kernel, 29, sizeof(cl_float), (void*)&(C->cs));
        nAddedArgs = 4;
	}
	if(isLocalMemory) {
	    err_code |= sendArgument(kernel, 26+nAddedArgs, _local_work_size*sizeof(cl_float), NULL);
	    err_code |= sendArgument(kernel, 27+nAddedArgs, _local_work_size*sizeof(vec     ), NULL);
	    err_code |= sendArgument(kernel, 28+nAddedArgs, _local_work_size*sizeof(cl_float), NULL);
	    err_code |= sendArgument(kernel, 29+nAddedArgs, _local_work_size*sizeof(cl_float), NULL);
	    err_code |= sendArgument(kernel, 30+nAddedArgs, _local_work_size*sizeof(cl_float), NULL);
	    err_code |= sendArgument(kernel, 31+nAddedArgs, _local_work_size*sizeof(vec     ), NULL);
	}
	if(err_code != CL_SUCCESS) {
		S->addMessage(3, "(Rates::execute): Can't send arguments to kernel.\n");
	    return true;
	}
	#ifdef HAVE_GPUPROFILE
	    err_code = clEnqueueNDRangeKernel(C->command_queue, kernel, 1, NULL, &_global_work_size, &_local_work_size, 0, NULL, &event);
	#else
	    err_code = clEnqueueNDRangeKernel(C->command_queue, kernel, 1, NULL, &_global_work_size, &_local_work_size, 0, NULL, NULL);
	#endif
	if(err_code != CL_SUCCESS) {
		S->addMessage(3, "(Rates::execute): Can't execute the kernel.\n");
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
	        S->addMessage(3, "(Rates::execute): Can't wait to kernels end.\n");
	        return true;
	    }
	    err_code |= clGetEventProfilingInfo(event, CL_PROFILING_COMMAND_END, sizeof(cl_ulong), &end, 0);
	    err_code |= clGetEventProfilingInfo(event, CL_PROFILING_COMMAND_START, sizeof(cl_ulong), &start, 0);
	    if(err_code != CL_SUCCESS) {
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
	cl_int err_code;
	cl_device_id device;
	cl_ulong localMem, reqLocalMem;
	err_code |= clGetCommandQueueInfo(C->command_queue,CL_QUEUE_DEVICE,
	                                sizeof(cl_device_id),&device, NULL);
	if(err_code != CL_SUCCESS) {
		S->addMessage(3, "(Rates::setupOpenCL): I Cannot get the device from the command queue.\n");
	    return true;
	}
	err_code |= clGetDeviceInfo(device, CL_DEVICE_LOCAL_MEM_SIZE, sizeof(localMem), &localMem, NULL);
	if(err_code != CL_SUCCESS) {
		S->addMessage(3, "(Rates::setupOpenCL): Can't get local memory available on device.\n");
	    return true;
	}
	if(!loadKernelFromFile(&clSortKernel, &program, C->context, C->device, _path, "SortData", ""))
	    return true;
	if(program)clReleaseProgram(program); program=0;
	char args[32]; strcpy(args, "");
	if(isDelta)
        strcat(args, "-D__DELTA_SPH__");
	if(!loadKernelFromFile(&kernel, &program, C->context, C->device, _path, "Rates", args))
	    return true;
	if(program)clReleaseProgram(program); program=0;
	// Test if there are enough local memory
	err_code |= clGetKernelWorkGroupInfo(clSortKernel,device,CL_KERNEL_LOCAL_MEM_SIZE,
	                                   sizeof(cl_ulong), &reqLocalMem, NULL);
	if(err_code != CL_SUCCESS) {
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
	err_code |= clGetKernelWorkGroupInfo(kernel,device,CL_KERNEL_LOCAL_MEM_SIZE,
	                                   sizeof(cl_ulong), &reqLocalMem, NULL);
	if(err_code != CL_SUCCESS) {
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
	err_code |= clGetKernelWorkGroupInfo(clSortKernel,device,CL_KERNEL_WORK_GROUP_SIZE,
	                                   sizeof(size_t), &localWorkGroupSize, NULL);
	if(err_code != CL_SUCCESS) {
		S->addMessage(3, "(Rates::setupOpenCL): Can't get sort maximum local work group size.\n");
	    return true;
	}
	if(localWorkGroupSize < _local_work_size)
	    _local_work_size  = localWorkGroupSize;
	err_code |= clGetKernelWorkGroupInfo(kernel,device,CL_KERNEL_WORK_GROUP_SIZE,
	                                   sizeof(size_t), &localWorkGroupSize, NULL);
	if(err_code != CL_SUCCESS) {
		S->addMessage(3, "(Rates::setupOpenCL): Can't get rates maximum local work group size.\n");
	    return true;
	}
	if(localWorkGroupSize < _local_work_size)
	    _local_work_size  = localWorkGroupSize;
	// Look for a better local work group size
	err_code |= clGetKernelWorkGroupInfo(kernel,device,CL_KERNEL_PREFERRED_WORK_GROUP_SIZE_MULTIPLE,
	                                   sizeof(size_t), &localWorkGroupSize, NULL);
	if(err_code != CL_SUCCESS) {
		S->addMessage(3, "(Rates::setupOpenCL): Can't get rates preferred local work group size.\n");
	    return true;
	}
	_local_work_size  = (_local_work_size/localWorkGroupSize) * localWorkGroupSize;
	_global_work_size = globalWorkSize(_local_work_size);
	// Test if the computation can be accelerated with local memory
	reqLocalMem += _local_work_size*(  sizeof(cl_float)
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
	    if(kernel)clReleaseKernel(kernel); kernel=0;
	    if(!loadKernelFromFile(&kernel, &program, C->context, C->device, _path, "Rates", options))
	        return true;
	    if(program)clReleaseProgram(program); program=0;
	}
	return false;
}

}}  // namespace
