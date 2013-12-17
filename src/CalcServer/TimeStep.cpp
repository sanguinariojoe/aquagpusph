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
	, mPath(NULL)
	, clProgram(NULL)
	, clKernel(NULL)
	, reduction(NULL)
	, MainDt(0.f)
	, dtClamp(0)
{
	InputOutput::ScreenManager *S = InputOutput::ScreenManager::singleton();
	InputOutput::ProblemSetup *P  = InputOutput::ProblemSetup::singleton();
	int nChar;
	MainDt = P->SPHParameters.h/P->SPHParameters.cs;
	MainDt /= P->SPHParameters.DivDt;
	if(P->TimeParameters.dtMode != __DT_VARIABLE__)
		return;
	//! 1st.- Get data
	nChar = strlen(P->OpenCL_kernels.time_step);
	if(nChar <= 0) {
		S->addMessage(3, "(TimeStep::TimeStep): Path of TimeStep kernel is empty.\n");
		exit(EXIT_FAILURE);
	}
	mPath = new char[nChar+4];
	if(!mPath) {
		S->addMessage(3, "(TimeStep::TimeStep): Can't allocate memory for path.\n");
		exit(EXIT_FAILURE);
	}
	strcpy(mPath, P->OpenCL_kernels.time_step);
	strcat(mPath, ".cl");
	//! 2nd.- Setup the kernel
	clLocalWorkSize = 256;
	clGlobalWorkSize = globalWorkSize(clLocalWorkSize);
	if(setupOpenCL()) {
		exit(EXIT_FAILURE);
	}
	S->addMessage(1, "(TimeStep::TimeStep): TimeStep ready to work!\n");
}

TimeStep::~TimeStep()
{
	InputOutput::ScreenManager *S = InputOutput::ScreenManager::singleton();
	InputOutput::ProblemSetup *P = InputOutput::ProblemSetup::singleton();
	if(clKernel)clReleaseKernel(clKernel); clKernel=NULL;
	if(clProgram)clReleaseProgram(clProgram); clProgram=NULL;
	if(mPath)delete[] mPath; mPath=NULL;
	S->addMessage(1, "(TimeStep::~TimeStep): Destroying time step reduction processor...\n");
	if(reduction) delete reduction; reduction=NULL;
}

bool TimeStep::execute()
{
	InputOutput::ScreenManager *S = InputOutput::ScreenManager::singleton();
	InputOutput::ProblemSetup *P = InputOutput::ProblemSetup::singleton();
	CalcServer *C = CalcServer::singleton();
	unsigned int i;
	cl_int clFlag = CL_SUCCESS;
	if(P->TimeParameters.dtMode == __DT_FIX__){
		C->dt = P->TimeParameters.dt;
		return false;
	}
	else if(P->TimeParameters.dtMode == __DT_FIXCALCULATED__){
		C->dt = MainDt;
		return false;
	}
	clFlag |= sendArgument(clKernel,  6, sizeof(cl_float), (void*)&(C->dt));
	clFlag |= sendArgument(clKernel,  7, sizeof(cl_float), (void*)&(C->cs));
	if(clFlag != CL_SUCCESS) {
		S->addMessage(3, "(TimeStep::Execute): Can't send variable to kernel.\n");
		return true;
	}
	#ifdef HAVE_GPUPROFILE
		cl_event event;
		cl_ulong end, start;
		profileTime(0.f);
		clFlag = clEnqueueNDRangeKernel(C->clComQueue, clKernel, 1, NULL, &clGlobalWorkSize, NULL, 0, NULL, &event);
	#else
		clFlag = clEnqueueNDRangeKernel(C->clComQueue, clKernel, 1, NULL, &clGlobalWorkSize, NULL, 0, NULL, NULL);
	#endif
	if(clFlag != CL_SUCCESS) {
		S->addMessage(3, "(TimeStep::Execute): Can't execute the kernel.\n");
		if(clFlag == CL_INVALID_WORK_GROUP_SIZE)
			S->addMessage(0, "\tInvalid local work group size.\n");
		else if(clFlag == CL_OUT_OF_RESOURCES)
			S->addMessage(0, "\tDevice out of resources.\n");
		else if(clFlag == CL_MEM_OBJECT_ALLOCATION_FAILURE)
			S->addMessage(0, "\tAllocation error at device.\n");
		else if(clFlag == CL_OUT_OF_HOST_MEMORY)
			S->addMessage(0, "\tfailure to allocate resources required by the OpenCL implementation on the host.\n");
		return true;
	}
	#ifdef HAVE_GPUPROFILE
		clFlag = clWaitForEvents(1, &event);
		if(clFlag != CL_SUCCESS) {
			S->addMessage(3, "(TimeStep::Execute): Can't wait to kernels end.\n");
			return true;
		}
		clFlag |= clGetEventProfilingInfo(event, CL_PROFILING_COMMAND_END, sizeof(cl_ulong), &end, 0);
		clFlag |= clGetEventProfilingInfo(event, CL_PROFILING_COMMAND_START, sizeof(cl_ulong), &start, 0);
		if(clFlag != CL_SUCCESS) {
			S->addMessage(3, "(TimeStep::Execute): Can't profile kernel execution.\n");
			return true;
		}
		profileTime(profileTime() + (end - start)/1000.f);  // 10^-3 ms
	#endif

    cl_mem reduced = reduction->execute();
    if(!reduced)
        return true;
	if(C->getData((void *)&C->dt, reduced, sizeof(cl_float)))
		return true;
	C->dt /= C->DivDt;
	if(MainDt > 10.f*C->dt)
	{
		char Log[256];
		sprintf(Log, "(TimeStep::Execute): timestep has dramaticaly decreased! [%g -> %g]\n",MainDt,C->dt);
		S->addMessage(3, Log);
		MainDt = C->dt;
	}
	if(C->dt < P->TimeParameters.mindt){
		if(!dtClamp){
			char Log[256];
			sprintf(Log, "(TimeStep::Execute): timestep lower than minimum value [%g < %g], will be clamped therefore\n",
			        C->dt,P->TimeParameters.mindt);
			S->addMessage(3, Log);
		}
		dtClamp = 1;
		C->dt = P->TimeParameters.mindt;
	}
	else
		dtClamp = 0;
	return false;
}

bool TimeStep::setupOpenCL()
{
	CalcServer *C = CalcServer::singleton();
	int clFlag;
	if(!loadKernelFromFile(&clKernel, &clProgram, C->clContext, C->clDevice, mPath, "TimeStep", ""))
		return true;
	clFlag  = sendArgument(clKernel,  0, sizeof(cl_mem ), (void*)&(C->dtconv));
	clFlag |= sendArgument(clKernel,  1, sizeof(cl_mem ), (void*)&(C->v));
	clFlag |= sendArgument(clKernel,  2, sizeof(cl_mem ), (void*)&(C->f));
	clFlag |= sendArgument(clKernel,  3, sizeof(cl_mem ), (void*)&(C->hp));
	clFlag |= sendArgument(clKernel,  4, sizeof(cl_mem ), (void*)&(C->sigma));
    clFlag |= sendArgument(clKernel,  5, sizeof(cl_uint), (void*)&(C->n));

	if(clFlag)
		return true;
    reduction = new Reduction(C->dtconv, C->n, "float", "INFINITY", "c = (a < b) ? a : b;");
	return false;
}

}}  // namespace
