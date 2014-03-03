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
#include <CalcServer/Boundary/DeLeffe.h>
#include <CalcServer.h>

namespace Aqua{ namespace CalcServer{ namespace Boundary{

DeLeffe::DeLeffe()
	: Kernel("DeLeffe")
	, _path(0)
	, _program(0)
	, _setup_kernel(0)
	, _boundary_kernel(0)
{
	InputOutput::ScreenManager *S = InputOutput::ScreenManager::singleton();
	InputOutput::ProblemSetup *P = InputOutput::ProblemSetup::singleton();
	if(P->SPH_opts.boundary_type!=2)  // DeLeffe condition has not been selected
	    return;

	int str_len = strlen(P->OpenCL_kernels.de_Leffe);
	if(str_len <= 0) {
	    S->addMessageF(3, "The path of the kernel is empty.\n");
	    exit(EXIT_FAILURE);
	}
	_path = new char[str_len+4];
	if(!_path) {
	    S->addMessageF(3, "Memory cannot be allocated for the path.\n");
	    exit(EXIT_FAILURE);
	}
	strcpy(_path, P->OpenCL_kernels.de_Leffe);
	strcat(_path, ".cl");

	_local_work_size  = localWorkSize();
	if(!_local_work_size){
	    S->addMessageF(3, "I cannot get a valid local work size for the required computation tool.\n");
	    exit(EXIT_FAILURE);
	}
	_global_work_size = globalWorkSize(_local_work_size);
	if(setupOpenCL()) {
	    exit(EXIT_FAILURE);
	}
	S->addMessageF(1, "DeLeffe ready to work!\n");
}

DeLeffe::~DeLeffe()
{
	if(_boundary_kernel)clReleaseKernel(_boundary_kernel); _boundary_kernel=0;
	if(_setup_kernel)clReleaseKernel(_setup_kernel); _setup_kernel=0;
	if(_program)clReleaseProgram(_program); _program=0;
	if(_path) delete[] _path; _path=0;
}

bool DeLeffe::execute()
{
	InputOutput::ProblemSetup *P = InputOutput::ProblemSetup::singleton();
	if(P->SPH_opts.boundary_type!=2)
	    return false;
	if(elements())
	    return true;
	if(boundary())
	    return true;
	return false;
}

bool DeLeffe::elements()
{
	InputOutput::ScreenManager *S = InputOutput::ScreenManager::singleton();
	CalcServer *C = CalcServer::singleton();
	cl_int err_code=0;

	err_code |= sendArgument(_setup_kernel,
                             0, sizeof(cl_mem),
                             (void*)&(C->imove));
	err_code |= sendArgument(_setup_kernel,
                             1,
                             sizeof(cl_mem),
                             (void*)&(C->f));
	err_code |= sendArgument(_setup_kernel,
                             2,
                             sizeof(cl_mem),
                             (void*)&(C->drdt));
	err_code |= sendArgument(_setup_kernel,
                             3,
                             sizeof(cl_mem),
                             (void*)&(C->press));
	err_code |= sendArgument(_setup_kernel,
                             4,
                             sizeof(cl_mem),
                             (void*)&(C->dens));
	err_code |= sendArgument(_setup_kernel,
                             5,
                             sizeof(cl_mem),
                             (void*)&(C->refd));
	err_code |= sendArgument(_setup_kernel,
                             6,
                             sizeof(cl_mem),
                             (void*)&(C->ifluid));
	err_code |= sendArgument(_setup_kernel,
                             7,
                             sizeof(cl_mem),
                             (void*)&(C->gamma));
	err_code |= sendArgument(_setup_kernel,
                             8,
                             sizeof(cl_mem),
                             (void*)&(C->shepard));
	err_code |= sendArgument(_setup_kernel,
                             9,
                             sizeof(cl_uint),
                             (void*)&(C->n));
	err_code |= sendArgument(_setup_kernel,
                             10,
                             sizeof(cl_float),
                             (void*)&(C->cs));
	if(err_code != CL_SUCCESS) {
		S->addMessageF(3, "Failure sending variables to Vertex set kernel.\n");
	    return true;
	}

	#ifdef HAVE_GPUPROFILE
	    cl_event event;
	    cl_ulong end, start;
	    err_code = clEnqueueNDRangeKernel(C->command_queue,
                                          _setup_kernel,
                                          1,
                                          NULL,
                                          &_global_work_size,
                                          NULL,
                                          0,
                                          NULL,
                                          &event);
	#else
	    err_code = clEnqueueNDRangeKernel(C->command_queue,
                                          _setup_kernel,
                                          1,
                                          NULL,
                                          &_global_work_size,
                                          NULL,
                                          0,
                                          NULL,
                                          NULL);
	#endif
	if(err_code != CL_SUCCESS) {
		S->addMessageF(3, "I cannot execute the kernel.\n");
        S->printOpenCLError(err_code);
	    return true;
	}
	#ifdef HAVE_GPUPROFILE
	    err_code = clWaitForEvents(1, &event);
	    if(err_code != CL_SUCCESS) {
	        S->addMessage(3, "Impossible to wait for the kernels end.\n");
            S->printOpenCLError(err_code);
	        return true;
	    }
	    err_code |= clGetEventProfilingInfo(event,
                                            CL_PROFILING_COMMAND_END,
                                            sizeof(cl_ulong),
                                            &end,
                                            0);
	    if(err_code != CL_SUCCESS) {
	        S->addMessage(3, "I cannot profile the kernel execution.\n");
            S->printOpenCLError(err_code);
	        return true;
	    }
	    err_code |= clGetEventProfilingInfo(event,
                                            CL_PROFILING_COMMAND_START,
                                            sizeof(cl_ulong),
                                            &start,
                                            0);
	    if(err_code != CL_SUCCESS) {
	        S->addMessage(3, "I cannot profile the kernel execution.\n");
            S->printOpenCLError(err_code);
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

	err_code |= sendArgument(_boundary_kernel,
                             0,
                             sizeof(cl_mem),
                             (void*)&(C->ifluid));
	err_code |= sendArgument(_boundary_kernel,
                             1,
                             sizeof(cl_mem),
                             (void*)&(C->imove));
	err_code |= sendArgument(_boundary_kernel,
                             2,
                             sizeof(cl_mem),
                             (void*)&(C->pos));
	err_code |= sendArgument(_boundary_kernel,
                             3,
                             sizeof(cl_mem),
                             (void*)&(C->normal));
	err_code |= sendArgument(_boundary_kernel,
                             4,
                             sizeof(cl_mem),
                             (void*)&(C->v));
	err_code |= sendArgument(_boundary_kernel,
                             5,
                             sizeof(cl_mem),
                             (void*)&(C->dens));
	err_code |= sendArgument(_boundary_kernel,
                             6,
                             sizeof(cl_mem),
                             (void*)&(C->press));
	err_code |= sendArgument(_boundary_kernel,
                             7,
                             sizeof(cl_mem),
                             (void*)&(C->mass));
	err_code |= sendArgument(_boundary_kernel,
                             8,
                             sizeof(cl_mem),
                             (void*)&(C->visc_dyn));
	err_code |= sendArgument(_boundary_kernel,
                             9,
                             sizeof(cl_mem),
                             (void*)&(C->f));
	err_code |= sendArgument(_boundary_kernel,
                             10,
                             sizeof(cl_mem),
                             (void*)&(C->drdt));
	err_code |= sendArgument(_boundary_kernel,
                             11,
                             sizeof(cl_mem),
                             (void*)&(C->icell));
	err_code |= sendArgument(_boundary_kernel,
                             12,
                             sizeof(cl_mem),
                             (void*)&(C->ihoc));
	err_code |= sendArgument(_boundary_kernel,
                             13,
                             sizeof(cl_uint),
                             (void*)&(C->N));
	err_code |= sendArgument(_boundary_kernel,
                             14,
                             sizeof(cl_float),
                             (void*)&(C->hfac));
	err_code |= sendArgument(_boundary_kernel,
                             15, sizeof(uivec),
                             (void*)&(C->num_cells_vec));
	if(err_code != CL_SUCCESS) {
		S->addMessageF(3, "Failure sending variables to boundary computation kernel.\n");
	    return true;
	}

	#ifdef HAVE_GPUPROFILE
	    cl_event event;
	    cl_ulong end, start;
	    profileTime(0.f);
	    err_code = clEnqueueNDRangeKernel(C->command_queue,
                                          _boundary_kernel,
                                          1,
                                          NULL,
                                          &_global_work_size,
                                          &_local_work_size,
                                          0,
                                          NULL,
                                          &event);
	#else
	    err_code = clEnqueueNDRangeKernel(C->command_queue,
                                          _boundary_kernel,
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
	    return true;
	}
	#ifdef HAVE_GPUPROFILE
	    err_code = clWaitForEvents(1, &event);
	    if(err_code != CL_SUCCESS) {
	        S->addMessage(3, "Impossible to wait for the kernels end.\n");
            S->printOpenCLError(err_code);
	        return true;
	    }
	    err_code |= clGetEventProfilingInfo(event,
                                            CL_PROFILING_COMMAND_END,
                                            sizeof(cl_ulong),
                                            &end,
                                            0);
	    if(err_code != CL_SUCCESS) {
	        S->addMessage(3, "I cannot profile the kernel execution.\n");
            S->printOpenCLError(err_code);
	        return true;
	    }
	    err_code |= clGetEventProfilingInfo(event,
                                            CL_PROFILING_COMMAND_START,
                                            sizeof(cl_ulong),
                                            &start,
                                            0);
	    if(err_code != CL_SUCCESS) {
	        S->addMessage(3, "I cannot profile the kernel execution.\n");
            S->printOpenCLError(err_code);
	        return true;
	    }
	    profileTime(profileTime() + (end - start)/1000.f);  // 10^-3 ms
	#endif
	return false;
}

bool DeLeffe::setupOpenCL()
{
	InputOutput::ScreenManager *S = InputOutput::ScreenManager::singleton();
	InputOutput::ProblemSetup *P = InputOutput::ProblemSetup::singleton();
	CalcServer *C = CalcServer::singleton();
	char msg[1024];
	cl_int err_code;
	cl_device_id device;
	cl_ulong local_mem, required_local_mem;
	err_code |= clGetCommandQueueInfo(C->command_queue,CL_QUEUE_DEVICE,
	                                sizeof(cl_device_id),&device, NULL);
	if(err_code != CL_SUCCESS) {
		S->addMessageF(3, "I Cannot get the device from the command queue.\n");
        S->printOpenCLError(err_code);
	    return true;
	}
	err_code |= clGetDeviceInfo(device,
                                CL_DEVICE_LOCAL_MEM_SIZE,
                                sizeof(local_mem),
                                &local_mem,
                                NULL);
	if(err_code != CL_SUCCESS) {
		S->addMessageF(3, "Failure getting the local memory available on the device.\n");
        S->printOpenCLError(err_code);
	    return true;
	}
	if(!loadKernelFromFile(&_setup_kernel,
                           &_program,
                           C->context,
                           C->device,
                           _path,
                           "Vertices",
                           ""))
	    return true;
	if(_program)clReleaseProgram(_program); _program=0;
	if(!loadKernelFromFile(&_boundary_kernel,
                           &_program,
                           C->context,
                           C->device,
                           _path,
                           "Boundary",
                           ""))
	    return true;
	if(_program)clReleaseProgram(_program); _program=0;

	err_code |= clGetKernelWorkGroupInfo(_setup_kernel,
                                         device,
                                         CL_KERNEL_LOCAL_MEM_SIZE,
	                                     sizeof(cl_ulong),
                                         &required_local_mem,
                                         NULL);
	if(err_code != CL_SUCCESS) {
		S->addMessageF(3, "Failure getting the kernel memory usage (boundary elements setup).\n");
        S->printOpenCLError(err_code);
	    return true;
	}
	if(local_mem < required_local_mem){
		S->addMessageF(3, "Not enough local memory kernel (boundary elements setup).\n");
	    sprintf(msg, "\tNeeds %lu bytes, but only %lu bytes are available.\n",
	           required_local_mem, local_mem);
	    S->addMessage(0, msg);
	    return true;
	}
	err_code |= clGetKernelWorkGroupInfo(_boundary_kernel,device,CL_KERNEL_LOCAL_MEM_SIZE,
	                                   sizeof(cl_ulong), &required_local_mem, NULL);
	if(err_code != CL_SUCCESS) {
		S->addMessageF(3, "Failure getting the kernel memory usage (boundary condition).\n");
        S->printOpenCLError(err_code);
	    return true;
	}
	if(local_mem < required_local_mem){
		S->addMessageF(3, "Not enough local memory kernel (boundary condition).\n");
	    sprintf(msg, "\tNeeds %lu bytes, but only %lu bytes are available.\n",
	           required_local_mem, local_mem);
	    S->addMessage(0, msg);
	    return true;
	}

	size_t local_work_size=0;
	err_code |= clGetKernelWorkGroupInfo(_setup_kernel,
                                         device,
                                         CL_KERNEL_WORK_GROUP_SIZE,
	                                     sizeof(size_t),
                                         &local_work_size,
                                         NULL);
	if(err_code != CL_SUCCESS) {
		S->addMessageF(3, "I cannont get the maximum local work size (boundary elements setup).\n");
        S->printOpenCLError(err_code);
	    return true;
	}
	if(local_work_size < _local_work_size)
	    _local_work_size  = local_work_size;
	err_code |= clGetKernelWorkGroupInfo(_boundary_kernel,
                                         device,
                                         CL_KERNEL_WORK_GROUP_SIZE,
	                                     sizeof(size_t),
                                         &local_work_size, NULL);
	if(err_code != CL_SUCCESS) {
		S->addMessageF(3, "I cannont get the maximum local work size (boundary condition).\n");
        S->printOpenCLError(err_code);
	    return true;
	}
	if(local_work_size < _local_work_size)
	    _local_work_size  = local_work_size;

	err_code |= clGetKernelWorkGroupInfo(_boundary_kernel,
                                         device,
                                         CL_KERNEL_PREFERRED_WORK_GROUP_SIZE_MULTIPLE,
	                                     sizeof(size_t),
                                         &local_work_size,
                                         NULL);
	if(err_code != CL_SUCCESS) {
		S->addMessageF(3, "I cannot query the preferred local work size (boundary condition).\n");
        S->printOpenCLError(err_code);
	    return true;
	}
	_local_work_size  = (_local_work_size/local_work_size) * local_work_size;
	_global_work_size = globalWorkSize(_local_work_size);

	required_local_mem += _local_work_size*(
          sizeof(vec)
        + sizeof(cl_float));
    char options[256]; strcpy(options, "");
	if(local_mem < required_local_mem){
		S->addMessageF(2, "Not enough local memory for boundary.\n");
	    sprintf(msg, "\tNeeds %lu bytes, but only %lu bytes are available.\n",
	           required_local_mem, local_mem);
	    S->addMessage(0, msg);
	    S->addMessage(0, "\tLocal memory usage will be avoided therefore.\n");
	}
	else{
        sprintf(options, "-DLOCAL_MEM_SIZE=%lu", _local_work_size);
	}
    if(_boundary_kernel)clReleaseKernel(_boundary_kernel); _boundary_kernel=0;
    if(!loadKernelFromFile(&_boundary_kernel,
                           &_program,
                           C->context,
                           C->device,
                           _path,
                           "Boundary",
                           options))
        return true;
    if(_program)clReleaseProgram(_program); _program=0;
	return false;
}

}}}  // namespaces
