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
#include <CalcServer/Boundary/ElasticBounce.h>

// ----------------------------------------------------------------------------
// Include the calculation server
// ----------------------------------------------------------------------------
#include <CalcServer.h>

namespace Aqua{ namespace CalcServer{ namespace Boundary{

ElasticBounce::ElasticBounce()
	: Kernel("ElasticBounce")
	, _path(0)
	, _program(0)
	, _kernel(0)
{
	InputOutput::ScreenManager *S = InputOutput::ScreenManager::singleton();
	InputOutput::ProblemSetup *P  = InputOutput::ProblemSetup::singleton();
	if(   P->SPH_opts.boundary_type!=0  // ElasticBounce boundary
	   && P->SPH_opts.boundary_type!=1  // Fixed particles boundary
	   && P->SPH_opts.boundary_type!=2  // DeLeffe boundary
	   )
	    return;
	//! 1st.- Get data
	int str_len = strlen(P->OpenCL_kernels.elastic_bounce);
	if(str_len <= 0) {
	    S->addMessage(3, "(ElasticBounce::ElasticBounce): _path of ElasticBounce kernel is empty.\n");
	    exit(EXIT_FAILURE);
	}
	_path = new char[str_len+4];
	if(!_path) {
	    S->addMessage(3, "(ElasticBounce::ElasticBounce): Memory cannot be allocated for the path.\n");
	    exit(EXIT_FAILURE);
	}
	strcpy(_path, P->OpenCL_kernels.elastic_bounce);
	strcat(_path, ".cl");
	//! 2nd.- Setup the kernel
	_local_work_size  = localWorkSize();
	if(!_local_work_size){
	    S->addMessage(3, "(ElasticBounce::ElasticBounce): I cannot get a valid local work size for the required computation tool.\n");
	    exit(EXIT_FAILURE);
	}
	_global_work_size = globalWorkSize(_local_work_size);
	if(setupOpenCL()) {
	    exit(EXIT_FAILURE);
	}
	S->addMessage(1, "(ElasticBounce::ElasticBounce): ElasticBounce boundary condition ready to work!\n");
}

ElasticBounce::~ElasticBounce()
{
	if(_kernel)clReleaseKernel(_kernel); _kernel=0;
	if(_program)clReleaseProgram(_program); _program=0;
	if(_path) delete[] _path; _path=0;
}

bool ElasticBounce::execute()
{
	//! Test if boundary condition must be executed
	InputOutput::ProblemSetup *P = InputOutput::ProblemSetup::singleton();
	if(   P->SPH_opts.boundary_type!=0  // ElasticBounce boundary
	   && P->SPH_opts.boundary_type!=1  // Fixed particles boundary
	   && P->SPH_opts.boundary_type!=2  // DeLeffe boundary
	   )
	    return false;
	InputOutput::ScreenManager *S = InputOutput::ScreenManager::singleton();
	CalcServer *C = CalcServer::singleton();
	cl_int err_code=0;
	//! Send variables to server
	err_code  = sendArgument(_kernel,  0, sizeof(cl_mem  ), (void*)&(C->imovein));
	err_code |= sendArgument(_kernel,  1, sizeof(cl_mem  ), (void*)&(C->posin));
	err_code |= sendArgument(_kernel,  2, sizeof(cl_mem  ), (void*)&(C->v));
	err_code |= sendArgument(_kernel,  3, sizeof(cl_mem  ), (void*)&(C->f));
	err_code |= sendArgument(_kernel,  4, sizeof(cl_mem  ), (void*)&(C->fin));
	err_code |= sendArgument(_kernel,  5, sizeof(cl_mem  ), (void*)&(C->normal));
	err_code |= sendArgument(_kernel,  6, sizeof(cl_mem  ), (void*)&(C->hpin));
	err_code |= sendArgument(_kernel,  7, sizeof(cl_mem  ), (void*)&(C->pos));
	err_code |= sendArgument(_kernel,  8, sizeof(cl_mem  ), (void*)&(C->icell));
	err_code |= sendArgument(_kernel,  9, sizeof(cl_mem  ), (void*)&(C->ihoc));
	err_code |= sendArgument(_kernel, 10, sizeof(cl_mem  ), (void*)&(C->permutation));
	err_code |= sendArgument(_kernel, 11, sizeof(cl_mem  ), (void*)&(C->permutation_inverse));
	err_code |= sendArgument(_kernel, 12, sizeof(cl_uint ), (void*)&(C->N));
	err_code |= sendArgument(_kernel, 13, sizeof(cl_float), (void*)&(C->dt));
	err_code |= sendArgument(_kernel, 14, sizeof(uivec   ), (void*)&(C->num_cells_vec));
	err_code |= sendArgument(_kernel, 15, sizeof(vec     ), (void*)&(C->g));
	if(err_code != CL_SUCCESS) {
		S->addMessage(3, "(ElasticBounce::Boundary): Can't send arguments to boundary computation kernel.\n");
	    return true;
	}
	//! Execute the kernel
	#ifdef HAVE_GPUPROFILE
	    cl_event event;
	    cl_ulong end, start;
	    profileTime(0.f);
	    err_code = clEnqueueNDRangeKernel(C->command_queue, _kernel, 1, NULL, &_global_work_size, NULL, 0, NULL, &event);
	#else
	    err_code = clEnqueueNDRangeKernel(C->command_queue, _kernel, 1, NULL, &_global_work_size, NULL, 0, NULL, NULL);
	#endif
	if(err_code != CL_SUCCESS) {
		S->addMessage(3, "(ElasticBounce::Boundary): I cannot execute the kernel.\n");
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
	        S->addMessage(3, "(ElasticBounce::Boundary): Impossible to wait for the kernels end.\n");
	        return true;
	    }
	    err_code |= clGetEventProfilingInfo(event, CL_PROFILING_COMMAND_END, sizeof(cl_ulong), &end, 0);
	    err_code |= clGetEventProfilingInfo(event, CL_PROFILING_COMMAND_START, sizeof(cl_ulong), &start, 0);
	    if(err_code != CL_SUCCESS) {
	        S->addMessage(3, "(ElasticBounce::Boundary): I cannot profile the kernel execution.\n");
	        return true;
	    }
	    profileTime(profileTime() + (end - start)/1000.f);  // 10^-3 ms
	#endif
	return false;
}

bool ElasticBounce::setupOpenCL()
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
		S->addMessage(3, "(ElasticBounce::setupOpenCL): I Cannot get the device from the command queue.\n");
	    return true;
	}
	err_code |= clGetDeviceInfo(device, CL_DEVICE_LOCAL_MEM_SIZE, sizeof(local_mem), &local_mem, NULL);
	if(err_code != CL_SUCCESS) {
		S->addMessage(3, "(ElasticBounce::setupOpenCL): Can't get local memory available on device.\n");
	    return true;
	}
	char flags[512];
	sprintf(flags, "-D__ELASTIC_FACTOR__=%ff", P->SPH_opts.elastic_factor);
	sprintf(flags, "%s -D__MIN_BOUND_DIST__=%ff", flags, fabs(P->SPH_opts.elastic_dist));
    if(P->SPH_opts.elastic_dist < 0.f){
        sprintf(flags, "%s -D__FORCE_MIN_BOUND_DIST__", flags);
    }
	if(!loadKernelFromFile(&_kernel, &_program, C->context, C->device, _path, "Boundary", flags))
	    return true;
	if(_program)clReleaseProgram(_program); _program=0;
	//! Test if there are enough local memory
	err_code |= clGetKernelWorkGroupInfo(_kernel,device,CL_KERNEL_LOCAL_MEM_SIZE,
	                                   sizeof(cl_ulong), &required_local_mem, NULL);
	if(err_code != CL_SUCCESS) {
		S->addMessage(3, "(ElasticBounce::setupOpenCL): Error retrieving the used local memory.\n");
	    return true;
	}
	if(local_mem < required_local_mem){
		S->addMessage(3, "(ElasticBounce::setupOpenCL): There are not enough local memory in the device.\n");
	    sprintf(msg, "\tNeeds %lu bytes, but only %lu bytes are available.\n",
	           required_local_mem, local_mem);
	    S->addMessage(0, msg);
	    return true;
	}
	//! Test if local work gorup size must be modified
	size_t local_work_size=0;
	err_code |= clGetKernelWorkGroupInfo(_kernel,device,CL_KERNEL_WORK_GROUP_SIZE,
	                                   sizeof(size_t), &local_work_size, NULL);
	if(err_code != CL_SUCCESS) {
		S->addMessage(3, "(ElasticBounce::setupOpenCL): Failure retrieving the maximum local work size.\n");
	    return true;
	}
	if(local_work_size < _local_work_size)
	    _local_work_size  = local_work_size;
	//! Look for better local work group size
	err_code |= clGetKernelWorkGroupInfo(_kernel,device,CL_KERNEL_PREFERRED_WORK_GROUP_SIZE_MULTIPLE,
	                                   sizeof(size_t), &local_work_size, NULL);
	if(err_code != CL_SUCCESS) {
		S->addMessage(3, "(ElasticBounce::setupOpenCL): I cannot query the preferred local work size");
	    return true;
	}
	_local_work_size  = (_local_work_size/local_work_size) * local_work_size;
	_global_work_size = globalWorkSize(_local_work_size);
	return false;
}

}}}  // namespaces
