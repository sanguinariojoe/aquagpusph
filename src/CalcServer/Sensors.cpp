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
#include <TimeManager.h>

// ----------------------------------------------------------------------------
// Include the Problem setup manager header
// ----------------------------------------------------------------------------
#include <FileManager.h>

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

Sensors::Sensors()
	: Kernel("Sensors")
	, n(0)
	, _path(0)
	, Output(0)
	, OutputTime(0)
	, hPos(0)
	, hPress(0)
	, hDens(0)
	, hSumW(0)
	, _program(0)
	, _kernel(0)
	, _global_work_size(0)
	, _local_work_size(0)
{
	//! 1st.- Get data
	InputOutput::ScreenManager *S = InputOutput::ScreenManager::singleton();
	InputOutput::ProblemSetup *P  = InputOutput::ProblemSetup::singleton();
	n = P->SensorsParameters.pos.size();
	if(!n)          // No sensors
	    return;
	unsigned int str_len = strlen(P->SensorsParameters.script);
	if(str_len <= 0) {
	    S->addMessage(3, "(Sensors::Sensors): _path of rates kernel is empty.\n");
	    exit(EXIT_FAILURE);
	}
	_path = new char[str_len+4];
	if(!_path) {
	    S->addMessage(3, "(Sensors::Sensors): Memory cannot be allocated for the path.\n");
	    exit(EXIT_FAILURE);
	}
	strcpy(_path, P->SensorsParameters.script);
	strcat(_path, ".cl");
	initOutput();
	//! 2nd.- Allocate memory
	hPos = new vec[n];
	hPress = new float[n];
	hDens = new float[n];
	hSumW = new float[n];
	if((!hPos) || (!hPress) || (!hDens) || (!hSumW)){
	    S->addMessage(3, "(Sensors::Sensors): Can't allocate memory for output variables.\n");
	    exit(EXIT_FAILURE);
	}
	//! 2nd.- Setup the kernel
	_local_work_size  = localWorkSize();
	if(!_local_work_size){
	    S->addMessage(3, "(Sensors::Sensors): I cannot get a valid local work size for the required computation tool.\n");
	    exit(EXIT_FAILURE);
	}
	_global_work_size = globalWorkSize(_local_work_size);
	if(setupOpenCL()) {
	    exit(EXIT_FAILURE);
	}
	S->addMessage(1, "(Sensors::Sensors): Sensors ready to work!\n");
}

Sensors::~Sensors()
{
	if(!n)
	    return;
	if(Output)fclose(Output); Output=0;
	if(_path)delete[] _path; _path=0;
	if(hPress)delete[] hPress; hPress=0;
	if(hDens)delete[] hDens; hDens=0;
	if(hSumW)delete[] hSumW; hSumW=0;
	if(_kernel)clReleaseKernel(_kernel); _kernel=0;
	if(_program)clReleaseProgram(_program); _program=0;
}

bool Sensors::execute()
{
	InputOutput::ScreenManager *S = InputOutput::ScreenManager::singleton();
	CalcServer *C = CalcServer::singleton();
	cl_int err_code=0;
	if(!n)
	    return false;
	//! 1st.- Send data
	unsigned int i0 = C->N - n;
	err_code |= sendArgument(_kernel,  0, sizeof(cl_mem), (void*)&(C->f));
	err_code |= sendArgument(_kernel,  1, sizeof(cl_mem), (void*)&(C->drdt));
	err_code |= sendArgument(_kernel,  2, sizeof(cl_mem), (void*)&(C->press));
	err_code |= sendArgument(_kernel,  3, sizeof(cl_mem), (void*)&(C->pressin));
	err_code |= sendArgument(_kernel,  4, sizeof(cl_mem), (void*)&(C->dens));
	err_code |= sendArgument(_kernel,  5, sizeof(cl_mem), (void*)&(C->densin));
	err_code |= sendArgument(_kernel,  6, sizeof(cl_mem), (void*)&(C->refd));
	err_code |= sendArgument(_kernel,  7, sizeof(cl_mem), (void*)&(C->ifluid));
	err_code |= sendArgument(_kernel,  8, sizeof(cl_mem), (void*)&(C->gamma));
	err_code |= sendArgument(_kernel,  9, sizeof(cl_mem), (void*)&(C->shepard));
	err_code |= sendArgument(_kernel, 10, sizeof(cl_mem), (void*)&(C->sensor_mode));
	err_code |= sendArgument(_kernel, 11, sizeof(cl_float), (void*)&(C->cs));
	err_code |= sendArgument(_kernel, 12, sizeof(cl_uint), (void*)&(i0));
	err_code |= sendArgument(_kernel, 13, sizeof(cl_uint), (void*)&(n));
	if(err_code)
	    return true;
	//! 2nd.- Execute the kernel
	#ifdef HAVE_GPUPROFILE
	    cl_event event;
	    cl_ulong end, start;
	    profileTime(0.f);
	    err_code = clEnqueueNDRangeKernel(C->command_queue, _kernel, 1, NULL, &_global_work_size, NULL, 0, NULL, &event);
	#else
	    err_code = clEnqueueNDRangeKernel(C->command_queue, _kernel, 1, NULL, &_global_work_size, NULL, 0, NULL, NULL);
	#endif
	if(err_code != CL_SUCCESS) {
		S->addMessage(3, "(Sensors::Execute): I cannot execute the kernel.\n");
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
	// Profile the kernel execution
	#ifdef HAVE_GPUPROFILE
	    err_code = clWaitForEvents(1, &event);
	    if(err_code != CL_SUCCESS) {
	        S->addMessage(3, "(Rates::Execute): Impossible to wait for the sorting kernel ends.\n");
	        return true;
	    }
	    err_code |= clGetEventProfilingInfo(event, CL_PROFILING_COMMAND_END, sizeof(cl_ulong), &end, 0);
	    err_code |= clGetEventProfilingInfo(event, CL_PROFILING_COMMAND_START, sizeof(cl_ulong), &start, 0);
	    if(err_code != CL_SUCCESS) {
	        S->addMessage(3, "(Rates::Execute): I cannot profile the sorting kernel execution.\n");
	        return true;
	    }
	    profileTime(profileTime() + (end - start)/1000.f);  // 10^-3 ms
	#endif
	//! 3rd.- Print output.
	if(printOutput())
	    return true;
	return false;
}

bool Sensors::printOutput()
{
	InputOutput::ProblemSetup *P = InputOutput::ProblemSetup::singleton();
	InputOutput::TimeManager *T = InputOutput::TimeManager::singleton();
	CalcServer *C = CalcServer::singleton();
	int i;
	cl_int err_code=0;
	if(T->time() - OutputTime < 1.f / P->SensorsParameters.fps)    // Don't needed output
	    return false;
	OutputTime = T->time();
	//! 1st.- Retrieve data from server
	unsigned int i0 = C->N - n;
	err_code  = C->getData((void*)hPos, C->pos, n*sizeof(vec), i0*sizeof(vec));
	err_code |= C->getData((void*)hPress, C->press, n*sizeof(cl_float), i0*sizeof(cl_float));
	err_code |= C->getData((void*)hDens, C->dens, n*sizeof(cl_float), i0*sizeof(cl_float));
	err_code |= C->getData((void*)hSumW, C->shepard, n*sizeof(cl_float), i0*sizeof(cl_float));
	if(err_code){
	    return true;
	}
	//! 2nd.- Print data
	fprintf(Output, "%f", OutputTime);
	for(i=0;i<n;i++){
	    #ifndef HAVE_3D
	        fprintf(Output, "\t%g\t%g\t%g\t%g\t%g", hPos[i].x, hPos[i].y, hPress[i], hDens[i], hSumW[i]);
	    #else
	        fprintf(Output, "\t%g\t%g\t%g\t%g\t%g\t%g", hPos[i].x, hPos[i].y, hPos[i].z, hPress[i], hDens[i], hSumW[i]);
	    #endif
	}
	fprintf(Output, "\n");
	fflush(Output);
	return false;
}

bool Sensors::setupOpenCL()
{
	InputOutput::ScreenManager *S = InputOutput::ScreenManager::singleton();
	CalcServer *C = CalcServer::singleton();
	int err_code = 0;
	if(!loadKernelFromFile(&_kernel, &_program, C->context, C->device, _path, "Sensors", ""))
	    return true;
	//! Test for right work group size
	cl_device_id device;
	size_t local_work_size=0;
	err_code |= clGetCommandQueueInfo(C->command_queue,CL_QUEUE_DEVICE,
	                                sizeof(cl_device_id),&device, NULL);
	if(err_code != CL_SUCCESS) {
		S->addMessage(3, "(Sensors::setupOpenCL): I Cannot get the device from the command queue.\n");
	    return true;
	}
	err_code |= clGetKernelWorkGroupInfo(_kernel,device,CL_KERNEL_WORK_GROUP_SIZE,
	                                   sizeof(size_t), &local_work_size, NULL);
	if(err_code != CL_SUCCESS) {
		S->addMessage(3, "(Sensors::setupOpenCL): Failure retrieving the maximum local work size.\n");
	    return true;
	}
	if(local_work_size < _local_work_size)
	    _local_work_size  = local_work_size;
	_global_work_size = globalWorkSize(_local_work_size);
	return false;
}

bool Sensors::initOutput(){
	InputOutput::ScreenManager *S = InputOutput::ScreenManager::singleton();
	InputOutput::ProblemSetup *P  = InputOutput::ProblemSetup::singleton();
	InputOutput::FileManager *FileMgr = InputOutput::FileManager::singleton();
	int i, newFile=0;
	char fileName[256];
	if(strcmp(FileMgr->outputPrefix(), ""))
	    sprintf(fileName, "%sSensors.dat", FileMgr->outputPrefix());
	else
	    strcpy(fileName, "Sensors.dat");
	OutputTime = -P->SensorsParameters.fps;  // 1st frame print
	//! 1st.- Start reading the last file looking for the OutputTime
	Output = fopen(fileName, "r");
	if(Output){
	    newFile=1;
	    char Line[256];
	    while( fgets( Line, 256*sizeof(char), Output) )
	    {
	        int LineStartChar=0;
	        while( (Line[LineStartChar] == ' ') || (Line[LineStartChar] == '\t') )
	            LineStartChar++;
	        if(Line[LineStartChar] != '#'){
	            if(sscanf(Line, "%f",&OutputTime) != 1) {
	                S->addMessage(2, "(Sensors::InitOutput): Sensor file is bad formatted, and will be overwritten.\n");
	                S->addMessage(0, "\tConsider restart the simulation if you are starting from previous work.\n");
	                return 1;
	            }
	        }
	    }
	    fclose(Output);
	}
	//! 2nd.- Open the file
	if(Output)
	    Output = fopen(fileName, "wa");
	else
	    Output = fopen(fileName, "w");
	if(!Output){
	    S->addMessage(3, "(Sensors::InitOutput): Can't open the file.\n");
	    S->addMessage(0, fileName);
	}
	//! 3rd.- Header data
	if(newFile)
	    return false;
	fprintf(Output,"#########################################################\n");
	fprintf(Output,"#                                                       #\n");
	fprintf(Output,"#    #    ##   #  #   #                           #     #\n");
	fprintf(Output,"#   # #  #  #  #  #  # #                          #     #\n");
	fprintf(Output,"#  ##### #  #  #  # #####  ##  ###  #  #  ## ###  ###   #\n");
	fprintf(Output,"#  #   # #  #  #  # #   # #  # #  # #  # #   #  # #  #  #\n");
	fprintf(Output,"#  #   # #  #  #  # #   # #  # #  # #  #   # #  # #  #  #\n");
	fprintf(Output,"#  #   #  ## #  ##  #   #  ### ###   ### ##  ###  #  #  #\n");
	fprintf(Output,"#                            # #             #          #\n");
	fprintf(Output,"#                          ##  #             #          #\n");
	fprintf(Output,"#                                                       #\n");
	fprintf(Output,"#########################################################\n");
	fprintf(Output,"#\n");
	fprintf(Output,"# Another QUAlity GPU-SPH, by CEHINAV.\n");
	fprintf(Output,"#\thttp://canal.etsin.upm.es/\n");
	fprintf(Output,"#Authors:\n");
	fprintf(Output,"#\tCercós Pita, Jose Luis\n");
	fprintf(Output,"#\tMiguel Gonzalez, Leo\n");
	fprintf(Output,"#\tSouto Iglesias, Antonio\n");
	fprintf(Output,"#\n");
	fprintf(Output,"#########################################################\n");
	fprintf(Output,"#\n");
	fprintf(Output,"# Sensors output file.\n");
	fprintf(Output,"# Number of sensors: %d\n", n);
	fprintf(Output,"# Rows:\n");
	fprintf(Output,"# Time");
	for(i=0;i<n;i++)
	    #ifndef HAVE_3D
	        fprintf(Output,"\ts%d.X\ts%d.Y\ts%d.press\ts%d.dens\ts%d.sumW", i,i,i,i,i);
	    #else
	        fprintf(Output,"\ts%d.X\ts%d.Y\ts%d.Z\ts%d.press\ts%d.dens\ts%d.sumW", i,i,i,i,i,i);
	    #endif

	fprintf(Output,"\n#\n");
	fprintf(Output,"#  [s]");
	for(i=0;i<n;i++)
	    #ifndef HAVE_3D
	        fprintf(Output,"\t  [m]\t  [m]\t   [Pa/m]\t [kg/m3]\t      []");
	    #else
	        fprintf(Output,"\t  [m]\t  [m]\t  [m]\t   [Pa/m]\t [kg/m3]\t      []");
	    #endif
	fprintf(Output,"\n#\n");
	fprintf(Output,"#########################################################\n");
	fflush(Output);
	return false;
}

}}  // namespace
