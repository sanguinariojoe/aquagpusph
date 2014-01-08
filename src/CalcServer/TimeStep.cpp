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
#include <CalcServer/Sensors.h>

// ----------------------------------------------------------------------------
// Include the calculation server
// ----------------------------------------------------------------------------
#include <CalcServer.h>

namespace Aqua{ namespace CalcServer{

TimeStep::TimeStep()
	: Kernel("TimeStep")
	, _path(NULL)
	, _program(NULL)
	, _kernel(NULL)
	, reduction(NULL)
	, MainDt(0.f)
	, dtClamp(0)
{
	InputOutput::ScreenManager *S = InputOutput::ScreenManager::singleton();
	InputOutput::ProblemSetup *P  = InputOutput::ProblemSetup::singleton();
	int str_len;
	MainDt = P->SPH_opts.h/P->SPH_opts.cs;
	MainDt /= P->SPH_opts.dt_divisor;
	if(P->time_opts.dt_mode != __DT_VARIABLE__)
		return;
	//! 1st.- Get data
	str_len = strlen(P->OpenCL_kernels.time_step);
	if(str_len <= 0) {
		S->addMessage(3, "(TimeStep::TimeStep): Path of TimeStep kernel is empty.\n");
		exit(EXIT_FAILURE);
	}
	_path = new char[str_len+4];
	if(!_path) {
		S->addMessage(3, "(TimeStep::TimeStep): Memory cannot be allocated for the path.\n");
		exit(EXIT_FAILURE);
	}
	strcpy(_path, P->OpenCL_kernels.time_step);
	strcat(_path, ".cl");
	//! 2nd.- Setup the kernel
	_local_work_size = 256;
	_global_work_size = globalWorkSize(_local_work_size);
	if(setupOpenCL()) {
		exit(EXIT_FAILURE);
	}
	S->addMessage(1, "(TimeStep::TimeStep): TimeStep ready to work!\n");
}

TimeStep::~TimeStep()
{
	InputOutput::ScreenManager *S = InputOutput::ScreenManager::singleton();
	InputOutput::ProblemSetup *P = InputOutput::ProblemSetup::singleton();
	if(_kernel)clReleaseKernel(_kernel); _kernel=NULL;
	if(_program)clReleaseProgram(_program); _program=NULL;
	if(_path)delete[] _path; _path=NULL;
	S->addMessage(1, "(TimeStep::~TimeStep): Destroying time step reduction processor...\n");
	if(reduction) delete reduction; reduction=NULL;
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
		C->dt = MainDt;
		return false;
	}
	err_code |= sendArgument(_kernel,  6, sizeof(cl_float), (void*)&(C->dt));
	err_code |= sendArgument(_kernel,  7, sizeof(cl_float), (void*)&(C->cs));
	if(err_code != CL_SUCCESS) {
		S->addMessage(3, "(TimeStep::Execute): I cannot send a variable to the kernel.\n");
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
		S->addMessage(3, "(TimeStep::Execute): I cannot execute the kernel.\n");
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
			S->addMessage(3, "(TimeStep::Execute): Impossible to wait for the kernels end.\n");
			return true;
		}
		err_code |= clGetEventProfilingInfo(event, CL_PROFILING_COMMAND_END, sizeof(cl_ulong), &end, 0);
		err_code |= clGetEventProfilingInfo(event, CL_PROFILING_COMMAND_START, sizeof(cl_ulong), &start, 0);
		if(err_code != CL_SUCCESS) {
			S->addMessage(3, "(TimeStep::Execute): I cannot profile the kernel execution.\n");
			return true;
		}
		profileTime(profileTime() + (end - start)/1000.f);  // 10^-3 ms
	#endif

    cl_mem reduced = reduction->execute();
    if(!reduced)
        return true;
	if(C->getData((void *)&C->dt, reduced, sizeof(cl_float)))
		return true;
	C->dt /= C->dt_divisor;
	if(MainDt > 10.f*C->dt)
	{
		char msg[256];
		sprintf(msg, "(TimeStep::Execute): timestep has dramaticaly decreased! [%g -> %g]\n",MainDt,C->dt);
		S->addMessage(3, msg);
		MainDt = C->dt;
	}
	if(C->dt < P->time_opts.dt_min){
		if(!dtClamp){
			char msg[256];
			sprintf(msg, "(TimeStep::Execute): timestep lower than minimum value [%g < %g], will be clamped therefore\n",
			        C->dt,P->time_opts.dt_min);
			S->addMessage(3, msg);
		}
		dtClamp = 1;
		C->dt = P->time_opts.dt_min;
	}
	else
		dtClamp = 0;
	return false;
}

bool TimeStep::setupOpenCL()
{
	CalcServer *C = CalcServer::singleton();
	int err_code;
	if(!loadKernelFromFile(&_kernel, &_program, C->context, C->device, _path, "TimeStep", ""))
		return true;
	err_code  = sendArgument(_kernel,  0, sizeof(cl_mem ), (void*)&(C->dtconv));
	err_code |= sendArgument(_kernel,  1, sizeof(cl_mem ), (void*)&(C->v));
	err_code |= sendArgument(_kernel,  2, sizeof(cl_mem ), (void*)&(C->f));
	err_code |= sendArgument(_kernel,  3, sizeof(cl_mem ), (void*)&(C->hp));
	err_code |= sendArgument(_kernel,  4, sizeof(cl_mem ), (void*)&(C->sigma));
    err_code |= sendArgument(_kernel,  5, sizeof(cl_uint), (void*)&(C->n));

	if(err_code)
		return true;
    reduction = new Reduction(C->dtconv, C->n, "float", "INFINITY", "c = (a < b) ? a : b;");
	return false;
}

}}  // namespace
