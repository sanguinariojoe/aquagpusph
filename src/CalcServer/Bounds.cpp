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
#include <CalcServer/Bounds.h>
#include <CalcServer.h>

namespace Aqua{ namespace CalcServer{

Bounds::Bounds()
	: Kernel("Bounds")
	, _device_mem(NULL)
	, _program(NULL)
	, _path(NULL)
	, _pos_max_kernel(NULL)
	, _pos_min_kernel(NULL)
	, _vel_max_kernel(NULL)
	, _vel_min_kernel(NULL)
	, _global_work_size(0)
	, _local_work_size(0)
	, _pos_max_reduction(NULL)
	, _pos_min_reduction(NULL)
	, _vel_max_reduction(NULL)
	, _vel_min_reduction(NULL)
{
	InputOutput::ScreenManager *S = InputOutput::ScreenManager::singleton();
	InputOutput::ProblemSetup *P  = InputOutput::ProblemSetup::singleton();
	unsigned int str_len = strlen(P->OpenCL_kernels.bounds);
	if(str_len <= 0) {
		S->addMessageF(3, "Path of bounds kernels is empty.\n");
		exit(EXIT_FAILURE);
	}
	_path = new char[str_len+4];
	if(!_path) {
		S->addMessageF(3, "Memory cannot be allocated for the path.\n");
		exit(EXIT_FAILURE);
	}
	strcpy(_path, P->OpenCL_kernels.bounds);
	strcat(_path, ".cl");

	_local_work_size  = localWorkSize();
	if(!_local_work_size){
	    S->addMessageF(3, "I cannot get a valid local work size for the required computation tool.\n");
	    exit(EXIT_FAILURE);
	}
	_global_work_size = globalWorkSize(_local_work_size);
	if(setupBounds()){
		exit(EXIT_FAILURE);
	}
	if(setupReduction()){
		exit(EXIT_FAILURE);
	}
	_pos_max.x = 0.f;
	_pos_max.y = 0.f;
	_pos_min.x = 0.f;
	_pos_min.y = 0.f;
	_vel_max.x = 0.f;
	_vel_max.y = 0.f;
	_vel_min.x = 0.f;
	_vel_min.y = 0.f;
	#ifdef HAVE_3D
		_pos_max.z = 0.f;
		_pos_min.z = 0.f;
		_vel_max.z = 0.f;
		_vel_min.z = 0.f;
	#endif
	S->addMessageF(1, "Bounds ready to work!\n");
}

Bounds::~Bounds()
{
	InputOutput::ScreenManager *S = InputOutput::ScreenManager::singleton();
	S->addMessageF(1, "Destroying maximum coordinates reduction processor...\n");
	if(_pos_max_reduction) delete _pos_max_reduction;
	S->addMessageF(1, "Destroying minimum coordinates reduction processor...\n");
	if(_pos_min_reduction)  delete _pos_min_reduction;
	S->addMessageF(1, "Destroying maximum velocity reduction processor...\n");
	if(_vel_max_reduction) delete _vel_max_reduction;
	S->addMessageF(1, "Destroying minimum velocity reduction processor...\n");
	if(_vel_min_reduction)  delete _vel_min_reduction;
	if(_device_mem)clReleaseMemObject(_device_mem); _device_mem=0;
	if(_pos_max_kernel)clReleaseKernel(_pos_max_kernel); _pos_max_kernel=0;
	if(_pos_min_kernel)clReleaseKernel(_pos_min_kernel); _pos_min_kernel=0;
	if(_vel_max_kernel)clReleaseKernel(_vel_max_kernel); _vel_max_kernel=0;
	if(_vel_min_kernel)clReleaseKernel(_vel_min_kernel); _vel_min_kernel=0;
	if(_program)clReleaseProgram(_program); _program=0;
	if(_path)delete[] _path; _path=0;
}

bool Bounds::execute()
{
    InputOutput::ScreenManager *S = InputOutput::ScreenManager::singleton();
	CalcServer *C = CalcServer::singleton();

    if(execute(&_pos_max, __BOUNDS_COORDS_MAX_OP__)){
        S->addMessage(0, "\tDuring the maximum coordinates computation.\n");
        return true;
    }
    if(execute(&_pos_min, __BOUNDS_COORDS_MIN_OP__)){
        S->addMessage(0, "\tDuring the minimum coordinates computation.\n");
        return true;
    }
    if(execute(&_vel_max, __BOUNDS_VEL_MAX_OP__)){
        S->addMessage(0, "\tDuring the maximum velocity computation.\n");
        return true;
    }
    if(execute(&_vel_min, __BOUNDS_VEL_MIN_OP__)){
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

	cl_kernel kernel = NULL;
	Reduction *reduction = NULL;
	if(op == __BOUNDS_COORDS_MAX_OP__){
        kernel = _pos_max_kernel;
        reduction = _pos_max_reduction;
	}
	if(op == __BOUNDS_COORDS_MIN_OP__){
        kernel = _pos_min_kernel;
        reduction = _pos_min_reduction;
	}
	if(op == __BOUNDS_VEL_MAX_OP__){
        kernel = _vel_max_kernel;
        reduction = _vel_max_reduction;
	}
	if(op == __BOUNDS_VEL_MIN_OP__){
        kernel = _vel_min_kernel;
        reduction = _vel_min_reduction;
	}

	#ifdef HAVE_GPUPROFILE
		cl_event event;
		cl_ulong end, start;
		profileTime(0.f);
		err_code = clEnqueueNDRangeKernel(C->command_queue,
                                          kernel,
                                          1,
                                          NULL,
                                          &_global_work_size,
                                          NULL,
                                          0,
                                          NULL,
                                          &event);
	#else
		err_code = clEnqueueNDRangeKernel(C->command_queue,
                                          kernel,
                                          1,
                                          NULL,
                                          &_global_work_size,
                                          NULL,
                                          0,
                                          NULL,
                                          NULL);
	#endif
	if(err_code != CL_SUCCESS) {
		S->addMessageF(3, "Imposible to execute bounds data preparation kernel.\n");
        S->printOpenCLError(err_code);
		return true;
	}
	#ifdef HAVE_GPUPROFILE
		err_code = clWaitForEvents(1, &event);
		if(err_code != CL_SUCCESS) {
			S->addMessageF(3, "I Cannot wait to the kernel end.\n");
            S->printOpenCLError(err_code);
			return true;
		}
		err_code = clGetEventProfilingInfo(event,
                                           CL_PROFILING_COMMAND_END,
                                           sizeof(cl_ulong),
                                           &end,
                                           0);
		if(err_code != CL_SUCCESS) {
			S->addMessageF(3, "I Cannot profile the kernel execution.\n");
            S->printOpenCLError(err_code);
			return true;
		}
		err_code = clGetEventProfilingInfo(event,
                                           CL_PROFILING_COMMAND_START,
                                           sizeof(cl_ulong),
                                           &start,
                                           0);
		if(err_code != CL_SUCCESS) {
			S->addMessageF(3, "I Cannot profile the kernel execution.\n");
            S->printOpenCLError(err_code);
			return true;
		}
		profileTime(profileTime() + (end - start)/1000.f);  // 10^-3 ms
	#endif

    cl_mem device_output = reduction->execute();
    if(!device_output)
        return true;
	if(C->getData((void *)output, device_output, sizeof(vec)))
		return true;
	return false;
}

bool Bounds::setupBounds()
{
	InputOutput::ScreenManager *S = InputOutput::ScreenManager::singleton();
	CalcServer *C = CalcServer::singleton();
	cl_int err_code = 0;
	char msg[1024];
	cl_device_id device;
	size_t local_work_size=0;

	_device_mem = C->allocMemory(C->n * sizeof( vec ));
	if(!_device_mem)
		return true;
	sprintf(msg,
            "\tAllocated memory = %u bytes\n",
            (unsigned int)C->allocated_mem);
	S->addMessage(0, msg);
	if(!loadKernelFromFile(&_pos_max_kernel,
                           &_program,
                           C->context,
                           C->device,
                           _path,
                           "MaximumCoordsFilter",
                           ""))
		return true;
	err_code  = sendArgument(_pos_max_kernel,
                             0,
                             sizeof(cl_mem),
                             (void*)&_device_mem);
	err_code |= sendArgument(_pos_max_kernel,
                             1,
                             sizeof(cl_mem),
                             (void*)&(C->imove));
	err_code |= sendArgument(_pos_max_kernel,
                             2,
                             sizeof(cl_mem),
                             (void*)&(C->pos));
	err_code |= sendArgument(_pos_max_kernel,
                             3,
                             sizeof(cl_uint),
                             (void*)&(C->n));
	if(err_code)
		return true;
	if(!loadKernelFromFile(&_pos_min_kernel,
                           &_program,
                           C->context,
                           C->device,
                           _path,
                           "MinimumCoordsFilter",
                           ""))
		return true;
	err_code  = sendArgument(_pos_min_kernel,
                             0,
                             sizeof(cl_mem),
                             (void*)&_device_mem);
	err_code |= sendArgument(_pos_min_kernel,
                             1,
                             sizeof(cl_mem),
                             (void*)&(C->imove));
	err_code |= sendArgument(_pos_min_kernel,
                             2,
                             sizeof(cl_mem),
                             (void*)&(C->pos));
	err_code |= sendArgument(_pos_min_kernel,
                             3,
                             sizeof(cl_uint),
                             (void*)&(C->n));
	if(err_code)
		return true;
	if(!loadKernelFromFile(&_vel_max_kernel,
                           &_program,
                           C->context,
                           C->device,
                           _path,
                           "MaximumVelFilter",
                           ""))
		return true;
	err_code  = sendArgument(_vel_max_kernel,
                             0,
                             sizeof(cl_mem),
                             (void*)&_device_mem);
	err_code |= sendArgument(_vel_max_kernel,
                             1,
                             sizeof(cl_mem),
                             (void*)&(C->imove));
	err_code |= sendArgument(_vel_max_kernel,
                             2,
                             sizeof(cl_mem),
                             (void*)&(C->v));
	err_code |= sendArgument(_vel_max_kernel,
                             3,
                             sizeof(cl_uint),
                             (void*)&(C->n));
	if(err_code)
		return true;
	if(!loadKernelFromFile(&_vel_min_kernel,
                           &_program,
                           C->context,
                           C->device,
                           _path,
                           "MinimumVelFilter",
                           ""))
		return true;
	err_code  = sendArgument(_vel_min_kernel,
                             0,
                             sizeof(cl_mem),
                             (void*)&_device_mem);
	err_code |= sendArgument(_vel_min_kernel,
                             1,
                             sizeof(cl_mem),
                             (void*)&(C->imove));
	err_code |= sendArgument(_vel_min_kernel,
                             2,
                             sizeof(cl_mem),
                             (void*)&(C->v));
	err_code |= sendArgument(_vel_min_kernel,
                             3,
                             sizeof(cl_uint),
                             (void*)&(C->n));
	if(err_code)
		return true;

	err_code = clGetCommandQueueInfo(C->command_queue,
                                     CL_QUEUE_DEVICE,
	                                 sizeof(cl_device_id),
                                     &device,
                                     NULL);
	if(err_code != CL_SUCCESS) {
		S->addMessageF(3, "I Cannot get the device from the command queue.\n");
	    S->printOpenCLError(err_code);
	    return true;
	}
	err_code = clGetKernelWorkGroupInfo(_pos_max_kernel,
                                        device,
                                        CL_KERNEL_WORK_GROUP_SIZE,
	                                    sizeof(size_t),
                                        &local_work_size,
                                        NULL);
	if(err_code != CL_SUCCESS) {
		S->addMessageF(3, "Failure retrieving the maximum local work size.\n");
	    S->printOpenCLError(err_code);
	    return true;
	}
	if(local_work_size < _local_work_size)
	    _local_work_size  = local_work_size;
	_global_work_size = globalWorkSize(_local_work_size);
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
        _pos_max_reduction = new Reduction(
            _device_mem,
            C->n,
            "vec",
            "(vec)(-INFINITY,-INFINITY,-INFINITY,0.f)",
            operation);
    #else
        _pos_max_reduction = new Reduction(
            _device_mem,
            C->n,
            "vec",
            "(vec)(-INFINITY,-INFINITY)",
            operation);
    #endif
    strcpy(operation, "");
    strcat(operation, "c.x = (a.x < b.x) ? a.x : b.x;\n");
    strcat(operation, "\tc.y = (a.y < b.y) ? a.y : b.y;\n");
    #ifdef HAVE_3D
        strcat(operation, "\tc.z = (a.z < b.z) ? a.z : b.z;\n");
        strcat(operation, "\tc.w = 0.f;\n");
        _pos_min_reduction = new Reduction(
            _device_mem,
            C->n,
            "vec",
            "(vec)(INFINITY,INFINITY,INFINITY,0.f)",
            operation);
    #else
        _pos_min_reduction = new Reduction(
            _device_mem,
            C->n,
            "vec",
            "(vec)(INFINITY,INFINITY)",
            operation);
    #endif
    _vel_max_reduction = new Reduction(
        _device_mem,
        C->n,
        "vec",
        "VEC_ZERO",
        "c = (dot(a,a) > dot(b,b)) ? a : b;\n");
    #ifdef HAVE_3D
        _vel_min_reduction = new Reduction(
            _device_mem,
            C->n,
            "vec",
            "(vec)(INFINITY,INFINITY,INFINITY,0.f)",
            "c = (dot(a,a) < dot(b,b)) ? a : b;\n");
    #else
        _vel_min_reduction = new Reduction(
            _device_mem,
            C->n,
            "vec",
            "(vec)(INFINITY,INFINITY)",
            "c = (dot(a,a) < dot(b,b)) ? a : b;\n");
    #endif

	return false;
}

}}  // namespace
