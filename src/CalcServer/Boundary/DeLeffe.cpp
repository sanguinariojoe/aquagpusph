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
#include <CalcServer/Boundary/DeLeffe.h>

// ----------------------------------------------------------------------------
// Include the calculation server
// ----------------------------------------------------------------------------
#include <CalcServer.h>

namespace Aqua{ namespace CalcServer{ namespace Boundary{

DeLeffe::DeLeffe()
	: Kernel("DeLeffe")
	, _path(0)
	, _program(0)
	, clVerticesKernel(0)
	, clBoundaryKernel(0)
	, _use_local_mem(true)
{
	InputOutput::ScreenManager *S = InputOutput::ScreenManager::singleton();
	InputOutput::ProblemSetup *P = InputOutput::ProblemSetup::singleton();
	if(P->SPH_opts.boundary_type!=2)  // DeLeffe condition has not been selected
	    return;
	//! 1st.- Get data
	int str_len = strlen(P->OpenCL_kernels.de_Leffe);
	if(str_len <= 0) {
	    S->addMessage(3, "(DeLeffe::DeLeffe): _path of DeLeffe kernel is empty.\n");
	    exit(EXIT_FAILURE);
	}
	_path = new char[str_len+4];
	if(!_path) {
	    S->addMessage(3, "(DeLeffe::DeLeffe): Memory cannot be allocated for the path.\n");
	    exit(EXIT_FAILURE);
	}
	strcpy(_path, P->OpenCL_kernels.de_Leffe);
	strcat(_path, ".cl");
	//! 2nd.- Setup the kernel
	_local_work_size  = localWorkSize();
	if(!_local_work_size){
	    S->addMessage(3, "(DeLeffe::DeLeffe): I cannot get a valid local work size for the required computation tool.\n");
	    exit(EXIT_FAILURE);
	}
	_global_work_size = globalWorkSize(_local_work_size);
	if(setupOpenCL()) {
	    exit(EXIT_FAILURE);
	}
	S->addMessage(1, "(DeLeffe::DeLeffe): DeLeffe ready to work!\n");
}

DeLeffe::~DeLeffe()
{
	if(clBoundaryKernel)clReleaseKernel(clBoundaryKernel); clBoundaryKernel=0;
	if(clVerticesKernel)clReleaseKernel(clVerticesKernel); clVerticesKernel=0;
	if(_program)clReleaseProgram(_program); _program=0;
	if(_path) delete[] _path; _path=0;
}

bool DeLeffe::execute()
{
	InputOutput::ProblemSetup *P = InputOutput::ProblemSetup::singleton();
	if(P->SPH_opts.boundary_type!=2)  // DeLeffe condition has not been selected
	    return false;
	if(vertices())
	    return true;
	if(boundary())
	    return true;
	return false;
}

bool DeLeffe::vertices()
{
	InputOutput::ScreenManager *S = InputOutput::ScreenManager::singleton();
	CalcServer *C = CalcServer::singleton();
	cl_int err_code=0;
	//! Send all variables to the server
	err_code |= sendArgument(clVerticesKernel,  0, sizeof(cl_mem  ), (void*)&(C->imove));
	err_code |= sendArgument(clVerticesKernel,  1, sizeof(cl_mem  ), (void*)&(C->f));
	err_code |= sendArgument(clVerticesKernel,  2, sizeof(cl_mem  ), (void*)&(C->drdt));
	err_code |= sendArgument(clVerticesKernel,  3, sizeof(cl_mem  ), (void*)&(C->press));
	err_code |= sendArgument(clVerticesKernel,  4, sizeof(cl_mem  ), (void*)&(C->pressin));
	err_code |= sendArgument(clVerticesKernel,  5, sizeof(cl_mem  ), (void*)&(C->dens));
	err_code |= sendArgument(clVerticesKernel,  6, sizeof(cl_mem  ), (void*)&(C->densin));
	err_code |= sendArgument(clVerticesKernel,  7, sizeof(cl_mem  ), (void*)&(C->refd));
	err_code |= sendArgument(clVerticesKernel,  8, sizeof(cl_mem  ), (void*)&(C->ifluid));
	err_code |= sendArgument(clVerticesKernel,  9, sizeof(cl_mem  ), (void*)&(C->gamma));
	err_code |= sendArgument(clVerticesKernel, 10, sizeof(cl_mem  ), (void*)&(C->normal));
	err_code |= sendArgument(clVerticesKernel, 11, sizeof(cl_mem  ), (void*)&(C->normalin));
	err_code |= sendArgument(clVerticesKernel, 12, sizeof(cl_mem  ), (void*)&(C->shepard));
	err_code |= sendArgument(clVerticesKernel, 13, sizeof(cl_mem  ), (void*)&(C->permutation_inverse));
	err_code |= sendArgument(clVerticesKernel, 14, sizeof(cl_uint ), (void*)&(C->n));
	err_code |= sendArgument(clVerticesKernel, 15, sizeof(cl_float), (void*)&(C->cs));
	if(err_code != CL_SUCCESS) {
		S->addMessage(3, "(DeLeffe::vertices): Can't send arguments to Vertex set kernel.\n");
	    return true;
	}
	//! Execute the kernel
	#ifdef HAVE_GPUPROFILE
	    cl_event event;
	    cl_ulong end, start;
	    err_code = clEnqueueNDRangeKernel(C->command_queue, clVerticesKernel, 1, NULL, &_global_work_size, NULL, 0, NULL, &event);
	#else
	    err_code = clEnqueueNDRangeKernel(C->command_queue, clVerticesKernel, 1, NULL, &_global_work_size, NULL, 0, NULL, NULL);
	#endif
	if(err_code != CL_SUCCESS) {
		S->addMessage(3, "(DeLeffe::vertices): I cannot execute the kernel.\n");
	    if(err_code == CL_INVALID_KERNEL_ARGS)
	        S->addMessage(0, "\tInvalid kernel arguments.\n");
	    else if(err_code == CL_INVALID_WORK_GROUP_SIZE)
	        S->addMessage(0, "\tInvalid local work group size.\n");
	    else if(err_code == CL_INVALID_WORK_ITEM_SIZE)
	        S->addMessage(0, "\tInvalid local work group size (greather than maximum allowed value).\n");
	    else if(err_code == CL_OUT_OF_RESOURCES)
	        S->addMessage(0, "\tDevice out of resources.\n");
	    else if(err_code == CL_MEM_OBJECT_ALLOCATION_FAILURE)
	        S->addMessage(0, "\tAllocation error at device.\n");
	    else if(err_code == CL_OUT_OF_HOST_MEMORY)
	        S->addMessage(0, "\tfailure to allocate resources required by the OpenCL implementation on the host.\n");
	    return true;
	}
	//! Profile the kernel execution
	#ifdef HAVE_GPUPROFILE
	    err_code = clWaitForEvents(1, &event);
	    if(err_code != CL_SUCCESS) {
	        S->addMessage(3, "(DeLeffe::vertices): Impossible to wait for the kernels end.\n");
	        return true;
	    }
	    err_code |= clGetEventProfilingInfo(event, CL_PROFILING_COMMAND_END, sizeof(cl_ulong), &end, 0);
	    err_code |= clGetEventProfilingInfo(event, CL_PROFILING_COMMAND_START, sizeof(cl_ulong), &start, 0);
	    if(err_code != CL_SUCCESS) {
	        S->addMessage(3, "(DeLeffe::vertices): I cannot profile the kernel execution.\n");
	        return true;
	    }
	    profileTime(profileTime() + (end - start)/1000.f);  // 10^-3 ms
	#endif

	return false;
}

bool DeLeffe::boundary()
{
	InputOutput::ScreenManager *S = InputOutput::ScreenManager::singleton();
	CalcServer *C = CalcServer::singleton();
	cl_int err_code=0;
	//! Send all variables to the server
	err_code |= sendArgument(clBoundaryKernel,  0, sizeof(cl_mem  ), (void*)&(C->ifluidin));
	err_code |= sendArgument(clBoundaryKernel,  1, sizeof(cl_mem  ), (void*)&(C->imovein));
	err_code |= sendArgument(clBoundaryKernel,  2, sizeof(cl_mem  ), (void*)&(C->posin));
	err_code |= sendArgument(clBoundaryKernel,  3, sizeof(cl_mem  ), (void*)&(C->normalin));
	err_code |= sendArgument(clBoundaryKernel,  4, sizeof(cl_mem  ), (void*)&(C->vin));
	err_code |= sendArgument(clBoundaryKernel,  5, sizeof(cl_mem  ), (void*)&(C->densin));
	err_code |= sendArgument(clBoundaryKernel,  6, sizeof(cl_mem  ), (void*)&(C->hpin));
	err_code |= sendArgument(clBoundaryKernel,  7, sizeof(cl_mem  ), (void*)&(C->pressin));
	err_code |= sendArgument(clBoundaryKernel,  8, sizeof(cl_mem  ), (void*)&(C->massin));
	err_code |= sendArgument(clBoundaryKernel,  9, sizeof(cl_mem  ), (void*)&(C->visc_dyn));
	err_code |= sendArgument(clBoundaryKernel, 10, sizeof(cl_mem  ), (void*)&(C->f));
	err_code |= sendArgument(clBoundaryKernel, 11, sizeof(cl_mem  ), (void*)&(C->drdt));
	err_code |= sendArgument(clBoundaryKernel, 12, sizeof(cl_mem  ), (void*)&(C->shepard_gradient));
	err_code |= sendArgument(clBoundaryKernel, 13, sizeof(cl_mem  ), (void*)&(C->icell));
	err_code |= sendArgument(clBoundaryKernel, 14, sizeof(cl_mem  ), (void*)&(C->ihoc));
	err_code |= sendArgument(clBoundaryKernel, 15, sizeof(cl_mem  ), (void*)&(C->permutation));
	err_code |= sendArgument(clBoundaryKernel, 16, sizeof(cl_mem  ), (void*)&(C->permutation_inverse));
	err_code |= sendArgument(clBoundaryKernel, 17, sizeof(cl_uint ), (void*)&(C->N));
	err_code |= sendArgument(clBoundaryKernel, 18, sizeof(cl_float), (void*)&(C->hfac));
	err_code |= sendArgument(clBoundaryKernel, 19, sizeof(uivec   ), (void*)&(C->num_cells_vec));
	if(_use_local_mem){
	    err_code |= sendArgument(clBoundaryKernel, 20, _local_work_size*sizeof(vec     ), NULL);
	    err_code |= sendArgument(clBoundaryKernel, 21, _local_work_size*sizeof(cl_float), NULL);
	    err_code |= sendArgument(clBoundaryKernel, 22, _local_work_size*sizeof(vec     ), NULL);
	}
	if(err_code != CL_SUCCESS) {
		S->addMessage(3, "(DeLeffe::boundary): Can't send arguments to boundary computation kernel.\n");
	    return true;
	}
	//! Execute the kernel
	#ifdef HAVE_GPUPROFILE
	    cl_event event;
	    cl_ulong end, start;
	    profileTime(0.f);
	    err_code = clEnqueueNDRangeKernel(C->command_queue, clBoundaryKernel, 1, NULL, &_global_work_size, &_local_work_size, 0, NULL, &event);
	#else
	    err_code = clEnqueueNDRangeKernel(C->command_queue, clBoundaryKernel, 1, NULL, &_global_work_size, &_local_work_size, 0, NULL, NULL);
	#endif
	if(err_code != CL_SUCCESS) {
		S->addMessage(3, "(DeLeffe::boundary): I cannot execute the kernel.\n");
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
	//! Profile the kernel execution
	#ifdef HAVE_GPUPROFILE
	    err_code = clWaitForEvents(1, &event);
	    if(err_code != CL_SUCCESS) {
	        S->addMessage(3, "(DeLeffe::boundary): Impossible to wait for the kernels end.\n");
	        return true;
	    }
	    err_code |= clGetEventProfilingInfo(event, CL_PROFILING_COMMAND_END, sizeof(cl_ulong), &end, 0);
	    err_code |= clGetEventProfilingInfo(event, CL_PROFILING_COMMAND_START, sizeof(cl_ulong), &start, 0);
	    if(err_code != CL_SUCCESS) {
	        S->addMessage(3, "(DeLeffe::boundary): I cannot profile the kernel execution.\n");
	        return true;
	    }
	    profileTime(profileTime() + (end - start)/1000.f);  // 10^-3 ms
	#endif
	return false;
}

bool DeLeffe::setupOpenCL()
{
	InputOutput::ScreenManager *S = InputOutput::ScreenManager::singleton();
	InputOutput::ProblemSetup *P  = InputOutput::ProblemSetup::singleton();
	CalcServer *C = CalcServer::singleton();
	char msg[1024];
	cl_int err_code;
	cl_device_id device;
	cl_ulong local_mem, required_local_mem;
	err_code |= clGetCommandQueueInfo(C->command_queue,CL_QUEUE_DEVICE,
	                                sizeof(cl_device_id),&device, NULL);
	if(err_code != CL_SUCCESS) {
		S->addMessage(3, "(DeLeffe::setupOpenCL): I Cannot get the device from the command queue.\n");
	    return true;
	}
	err_code |= clGetDeviceInfo(device, CL_DEVICE_LOCAL_MEM_SIZE, sizeof(local_mem), &local_mem, NULL);
	if(err_code != CL_SUCCESS) {
		S->addMessage(3, "(DeLeffe::setupOpenCL): Can't get local memory available on device.\n");
	    return true;
	}
	if(!loadKernelFromFile(&clVerticesKernel, &_program, C->context, C->device, _path, "Vertices", ""))
	    return true;
	if(_program)clReleaseProgram(_program); _program=0;
	if(!loadKernelFromFile(&clBoundaryKernel, &_program, C->context, C->device, _path, "Boundary", ""))
	    return true;
	if(_program)clReleaseProgram(_program); _program=0;
	//! Test if there are enough local memory
	err_code |= clGetKernelWorkGroupInfo(clVerticesKernel,device,CL_KERNEL_LOCAL_MEM_SIZE,
	                                   sizeof(cl_ulong), &required_local_mem, NULL);
	if(err_code != CL_SUCCESS) {
		S->addMessage(3, "(DeLeffe::setupOpenCL): Can't get vertices kernel memory usage.\n");
	    return true;
	}
	if(local_mem < required_local_mem){
		S->addMessage(3, "(DeLeffe::setupOpenCL): Not enough local memory for vertices execution.\n");
	    sprintf(msg, "\tNeeds %lu bytes, but only %lu bytes are available.\n",
	           required_local_mem, local_mem);
	    S->addMessage(0, msg);
	    return true;
	}
	err_code |= clGetKernelWorkGroupInfo(clBoundaryKernel,device,CL_KERNEL_LOCAL_MEM_SIZE,
	                                   sizeof(cl_ulong), &required_local_mem, NULL);
	if(err_code != CL_SUCCESS) {
		S->addMessage(3, "(DeLeffe::setupOpenCL): Can't get boundary kernel memory usage.\n");
	    return true;
	}
	if(local_mem < required_local_mem){
		S->addMessage(3, "(DeLeffe::setupOpenCL): Not enough local memory for boundary execution.\n");
	    sprintf(msg, "\tNeeds %lu bytes, but only %lu bytes are available.\n",
	           required_local_mem, local_mem);
	    S->addMessage(0, msg);
	    return true;
	}
	//! Test if local work gorup size must be modified
	size_t local_work_size=0;
	err_code |= clGetKernelWorkGroupInfo(clVerticesKernel,device,CL_KERNEL_WORK_GROUP_SIZE,
	                                   sizeof(size_t), &local_work_size, NULL);
	if(err_code != CL_SUCCESS) {
		S->addMessage(3, "(DeLeffe::setupOpenCL): Can't get vertices maximum local work group size.\n");
	    return true;
	}
	if(local_work_size < _local_work_size)
	    _local_work_size  = local_work_size;
	err_code |= clGetKernelWorkGroupInfo(clBoundaryKernel,device,CL_KERNEL_WORK_GROUP_SIZE,
	                                   sizeof(size_t), &local_work_size, NULL);
	if(err_code != CL_SUCCESS) {
		S->addMessage(3, "(DeLeffe::setupOpenCL): Can't get boundary maximum local work group size.\n");
	    return true;
	}
	if(local_work_size < _local_work_size)
	    _local_work_size  = local_work_size;
	//! Look for better local work group size
	err_code |= clGetKernelWorkGroupInfo(clBoundaryKernel,device,CL_KERNEL_PREFERRED_WORK_GROUP_SIZE_MULTIPLE,
	                                   sizeof(size_t), &local_work_size, NULL);
	if(err_code != CL_SUCCESS) {
		S->addMessage(3, "(DeLeffe::setupOpenCL): Can't get boundary preferred local work group size.\n");
	    return true;
	}
	_local_work_size  = (_local_work_size/local_work_size) * local_work_size;
	_global_work_size = globalWorkSize(_local_work_size);
	//! Test if computation can be accelerated with local memory
	required_local_mem += _local_work_size*(  sizeof(vec     )
	                                + sizeof(cl_float)
	                                + sizeof(vec     ));
	if(local_mem < required_local_mem){
		S->addMessage(2, "(DeLeffe::setupOpenCL): Not enough local memory for boundary.\n");
	    sprintf(msg, "\tNeeds %lu bytes, but only %lu bytes are available.\n",
	           required_local_mem, local_mem);
	    S->addMessage(0, msg);
	    S->addMessage(0, "\tLocal memory usage will be avoided therefore.\n");
	    _use_local_mem = false;
	    char options[19]; strcpy(options,"-D__NO_LOCAL_MEM__");
	    if(clBoundaryKernel)clReleaseKernel(clBoundaryKernel); clBoundaryKernel=0;
	    if(!loadKernelFromFile(&clBoundaryKernel, &_program, C->context, C->device, _path, "Boundary", options))
	        return true;
	    if(_program)clReleaseProgram(_program); _program=0;
	}
	return false;
}

}}}  // namespaces
