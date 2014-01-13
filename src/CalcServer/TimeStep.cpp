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
#include <CalcServer/Sensors.h>
#include <CalcServer.h>

namespace Aqua{ namespace CalcServer{

TimeStep::TimeStep()
	: Kernel("TimeStep")
	, _path(NULL)
	, _program(NULL)
	, _kernel(NULL)
	, _reduction(NULL)
	, _dt(0.f)
	, _is_dt_clamp(false)
{
	InputOutput::ScreenManager *S = InputOutput::ScreenManager::singleton();
	InputOutput::ProblemSetup *P  = InputOutput::ProblemSetup::singleton();
	int str_len;
	_dt = P->SPH_opts.h / P->SPH_opts.cs;
	_dt /= P->SPH_opts.dt_divisor;
	if(P->time_opts.dt_mode != __DT_VARIABLE__)
		return;

	str_len = strlen(P->OpenCL_kernels.time_step);
	if(str_len <= 0) {
		S->addMessageF(3, "The path of the kernel is empty.\n");
		exit(EXIT_FAILURE);
	}
	_path = new char[str_len+4];
	if(!_path) {
		S->addMessageF(3, "Memory cannot be allocated for the path.\n");
		exit(EXIT_FAILURE);
	}
	strcpy(_path, P->OpenCL_kernels.time_step);
	strcat(_path, ".cl");

	_local_work_size = 256;
	_global_work_size = globalWorkSize(_local_work_size);
	if(setupOpenCL()) {
		exit(EXIT_FAILURE);
	}
	S->addMessageF(1, "TimeStep ready to work!\n");
}

TimeStep::~TimeStep()
{
	InputOutput::ScreenManager *S = InputOutput::ScreenManager::singleton();
	InputOutput::ProblemSetup *P = InputOutput::ProblemSetup::singleton();
	if(_kernel)clReleaseKernel(_kernel); _kernel=NULL;
	if(_program)clReleaseProgram(_program); _program=NULL;
	if(_path)delete[] _path; _path=NULL;
	S->addMessageF(1, "Destroying time step reduction processor...\n");
	if(_reduction) delete _reduction; _reduction=NULL;
}

bool TimeStep::execute()
{
	InputOutput::ScreenManager *S = InputOutput::ScreenManager::singleton();
	InputOutput::ProblemSetup *P = InputOutput::ProblemSetup::singleton();
	CalcServer *C = CalcServer::singleton();
	unsigned int i;
	cl_int err_code = CL_SUCCESS;
	if(P->time_opts.dt_mode == __DT_FIX__){
		C->dt = P->time_opts.dt;
		return false;
	}
	else if(P->time_opts.dt_mode == __DT_FIXCALCULATED__){
		C->dt = _dt;
		return false;
	}
	err_code |= sendArgument(_kernel,  6, sizeof(cl_float), (void*)&(C->dt));
	err_code |= sendArgument(_kernel,  7, sizeof(cl_float), (void*)&(C->cs));
	if(err_code != CL_SUCCESS) {
		S->addMessageF(3, "I cannot send a variable to the kernel.\n");
		return true;
	}
	#ifdef HAVE_GPUPROFILE
		cl_event event;
		cl_ulong end, start;
		profileTime(0.f);
		err_code = clEnqueueNDRangeKernel(C->command_queue, _kernel, 1, NULL, &_global_work_size, NULL, 0, NULL, &event);
	#else
		err_code = clEnqueueNDRangeKernel(C->command_queue, _kernel, 1, NULL, &_global_work_size, NULL, 0, NULL, NULL);
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

    cl_mem reduced = _reduction->execute();
    if(!reduced)
        return true;
	if(C->getData((void *)&C->dt, reduced, sizeof(cl_float)))
		return true;
	C->dt /= C->dt_divisor;
	if(_dt > 10.f*C->dt)
	{
		char msg[256];
		sprintf(msg,
                "timestep has dramaticaly decreased! [%g -> %g]\n",
                _dt,
                C->dt);
		S->addMessageF(3, msg);
		_dt = C->dt;
	}
	if(C->dt < P->time_opts.dt_min){
		if(!_is_dt_clamp){
			char msg[256];
			sprintf(msg,
                    "timestep lower than minimum value [%g < %g], will be clamped therefore\n",
			        C->dt,
			        P->time_opts.dt_min);
			S->addMessageF(3, msg);
		}
		_is_dt_clamp = true;
		C->dt = P->time_opts.dt_min;
	}
	else
		_is_dt_clamp = false;
	return false;
}

bool TimeStep::setupOpenCL()
{
	CalcServer *C = CalcServer::singleton();
	int err_code;
	if(!loadKernelFromFile(&_kernel, &_program, C->context, C->device, _path, "TimeStep", ""))
		return true;
	err_code  = sendArgument(_kernel,
                             0,
                             sizeof(cl_mem),
                             (void*)&(C->dtconv));
	err_code |= sendArgument(_kernel,
                             1,
                             sizeof(cl_mem),
                             (void*)&(C->v));
	err_code |= sendArgument(_kernel,
                             2,
                             sizeof(cl_mem),
                             (void*)&(C->f));
	err_code |= sendArgument(_kernel,
                             3,
                             sizeof(cl_mem),
                             (void*)&(C->hp));
	err_code |= sendArgument(_kernel,
                             4,
                             sizeof(cl_mem),
                             (void*)&(C->sigma));
    err_code |= sendArgument(_kernel,
                             5,
                             sizeof(cl_uint),
                             (void*)&(C->n));

	if(err_code)
		return true;
    _reduction = new Reduction(C->dtconv,
                               C->n,
                               "float", "INFINITY", "c = (a < b) ? a : b;");
	return false;
}

}}  // namespace
