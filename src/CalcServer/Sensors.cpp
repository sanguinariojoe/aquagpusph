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
#include <TimeManager.h>
#include <FileManager.h>
#include <ScreenManager.h>
#include <CalcServer/Sensors.h>
#include <CalcServer.h>

namespace Aqua{ namespace CalcServer{

Sensors::Sensors()
	: Kernel("Sensors")
	, _n(0)
	, _path(0)
	, _output(0)
	, _output_time(0)
	, _pos(0)
	, _press(0)
	, _dens(0)
	, _sum_W(0)
	, _program(0)
	, _kernel(0)
	, _global_work_size(0)
	, _local_work_size(0)
{
	InputOutput::ScreenManager *S = InputOutput::ScreenManager::singleton();
	InputOutput::ProblemSetup *P  = InputOutput::ProblemSetup::singleton();
	_n = P->SensorsParameters.pos.size();
	if(!_n)
	    return;
	unsigned int str_len = strlen(P->SensorsParameters.script);
	if(str_len <= 0) {
	    S->addMessageF(3, "The path of the kernel is empty.\n");
	    exit(EXIT_FAILURE);
	}
	_path = new char[str_len+4];
	if(!_path) {
	    S->addMessageF(3, "Memory cannot be allocated for the path.\n");
	    exit(EXIT_FAILURE);
	}
	strcpy(_path, P->SensorsParameters.script);
	strcat(_path, ".cl");
	initOutput();

	_pos = new vec[_n];
	_press = new float[_n];
	_dens = new float[_n];
	_sum_W = new float[_n];
	if((!_pos) || (!_press) || (!_dens) || (!_sum_W)){
	    S->addMessageF(3, "Failure allocating memory for the output variables.\n");
	    exit(EXIT_FAILURE);
	}
	//! 2nd.- Setup the kernel
	_local_work_size  = localWorkSize();
	if(!_local_work_size){
	    S->addMessageF(3, "I cannot get a valid local work size for the required computation tool.\n");
	    exit(EXIT_FAILURE);
	}
	_global_work_size = globalWorkSize(_local_work_size);
	if(setupOpenCL()) {
	    exit(EXIT_FAILURE);
	}
	S->addMessageF(1, "Sensors ready to work!\n");
}

Sensors::~Sensors()
{
	if(!_n)
	    return;
	if(_output)fclose(_output); _output=0;
	if(_path)delete[] _path; _path=0;
	if(_press)delete[] _press; _press=0;
	if(_dens)delete[] _dens; _dens=0;
	if(_sum_W)delete[] _sum_W; _sum_W=0;
	if(_kernel)clReleaseKernel(_kernel); _kernel=0;
	if(_program)clReleaseProgram(_program); _program=0;
}

bool Sensors::execute()
{
	InputOutput::ScreenManager *S = InputOutput::ScreenManager::singleton();
	CalcServer *C = CalcServer::singleton();
	cl_int err_code=0;
	if(!_n)
	    return false;

	unsigned int i0 = C->N - _n;
	err_code |= sendArgument(_kernel,
                             0,
                             sizeof(cl_mem),
                             (void*)&(C->f));
	err_code |= sendArgument(_kernel,
                             1,
                             sizeof(cl_mem),
                             (void*)&(C->drdt));
	err_code |= sendArgument(_kernel,
                             2,
                             sizeof(cl_mem),
                             (void*)&(C->press));
	err_code |= sendArgument(_kernel,
                             3,
                             sizeof(cl_mem),
                             (void*)&(C->pressin));
	err_code |= sendArgument(_kernel,
                             4,
                             sizeof(cl_mem),
                             (void*)&(C->dens));
	err_code |= sendArgument(_kernel,
                             5,
                             sizeof(cl_mem),
                             (void*)&(C->densin));
	err_code |= sendArgument(_kernel,
                             6,
                             sizeof(cl_mem),
                             (void*)&(C->refd));
	err_code |= sendArgument(_kernel,
                             7,
                             sizeof(cl_mem),
                             (void*)&(C->ifluid));
	err_code |= sendArgument(_kernel,
                             8,
                             sizeof(cl_mem),
                             (void*)&(C->gamma));
	err_code |= sendArgument(_kernel,
                             9,
                             sizeof(cl_mem),
                             (void*)&(C->shepard));
	err_code |= sendArgument(_kernel,
                             10,
                             sizeof(cl_mem),
                             (void*)&(C->sensor_mode));
	err_code |= sendArgument(_kernel,
                             11,
                             sizeof(cl_float),
                             (void*)&(C->cs));
	err_code |= sendArgument(_kernel,
                             12,
                             sizeof(cl_uint),
                             (void*)&(i0));
	err_code |= sendArgument(_kernel,
                             13,
                             sizeof(cl_uint),
                             (void*)&(_n));
	if(err_code)
	    return true;

	#ifdef HAVE_GPUPROFILE
	    cl_event event;
	    cl_ulong end, start;
	    profileTime(0.f);
	    err_code = clEnqueueNDRangeKernel(C->command_queue,
                                          _kernel,
                                          1,
                                          NULL,
                                          &_global_work_size,
                                          NULL,
                                          0,
                                          NULL,
                                          &event);
	#else
	    err_code = clEnqueueNDRangeKernel(C->command_queue,
                                          _kernel,
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

	if(printOutput())
	    return true;
	return false;
}

vec* Sensors::positions()
{
	CalcServer *C = CalcServer::singleton();

    if(!_n){
        return NULL;
    }
	cl_int err_code;
	unsigned int i0 = C->N - _n;
	err_code = C->getData((void*)_pos,
                          C->pos,
                          _n*sizeof(vec),
                          i0*sizeof(vec));
    return _pos;
}

bool Sensors::printOutput()
{
	InputOutput::ProblemSetup *P = InputOutput::ProblemSetup::singleton();
	InputOutput::TimeManager *T = InputOutput::TimeManager::singleton();
	CalcServer *C = CalcServer::singleton();
	int i;
	cl_int err_code=0;
	if(T->time() - _output_time < 1.f / P->SensorsParameters.fps)
	    return false;
	_output_time = T->time();

	unsigned int i0 = C->N - _n;
	err_code  = C->getData((void*)_pos,
                           C->pos,
                           _n*sizeof(vec),
                           i0*sizeof(vec));
	err_code |= C->getData((void*)_press,
                           C->press,
                           _n*sizeof(cl_float),
                           i0*sizeof(cl_float));
	err_code |= C->getData((void*)_dens,
                           C->dens,
                           _n*sizeof(cl_float),
                           i0*sizeof(cl_float));
	err_code |= C->getData((void*)_sum_W,
                           C->shepard,
                           _n*sizeof(cl_float),
                           i0*sizeof(cl_float));
	if(err_code){
	    return true;
	}
	//! 2nd.- Print data
	fprintf(_output, "%f", _output_time);
	for(i=0;i<_n;i++){
	    #ifndef HAVE_3D
	        fprintf(_output,
                    "\t%g\t%g\t%g\t%g\t%g",
                    _pos[i].x,
                    _pos[i].y,
                    _press[i],
                    _dens[i],
                    _sum_W[i]);
	    #else
	        fprintf(_output,
                    "\t%g\t%g\t%g\t%g\t%g\t%g",
                    _pos[i].x,
                    _pos[i].y,
                    _pos[i].z,
                    _press[i],
                    _dens[i],
                    _sum_W[i]);
	    #endif
	}
	fprintf(_output, "\n");
	fflush(_output);
	return false;
}

bool Sensors::setupOpenCL()
{
	InputOutput::ScreenManager *S = InputOutput::ScreenManager::singleton();
	CalcServer *C = CalcServer::singleton();
	int err_code = 0;
	if(!loadKernelFromFile(&_kernel,
                           &_program,
                           C->context,
                           C->device,
                           _path,
                           "Sensors",
                           ""))
	    return true;

	cl_device_id device;
	size_t local_work_size=0;
	err_code |= clGetCommandQueueInfo(C->command_queue,
                                      CL_QUEUE_DEVICE,
	                                  sizeof(cl_device_id),
                                      &device,
                                      NULL);
	if(err_code != CL_SUCCESS) {
		S->addMessageF(3, "I Cannot get the device from the command queue.\n");
        S->printOpenCLError(err_code);
	    return true;
	}
	err_code |= clGetKernelWorkGroupInfo(_kernel,
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

bool Sensors::initOutput(){
	InputOutput::ScreenManager *S = InputOutput::ScreenManager::singleton();
	InputOutput::ProblemSetup *P  = InputOutput::ProblemSetup::singleton();
	InputOutput::FileManager *FileMgr = InputOutput::FileManager::singleton();
	int i, new_file=0;
	char file_name[256];
    strcpy(file_name, "Sensors.dat");
	_output_time = -P->SensorsParameters.fps;  // 1st frame print

	// We must read the last file looking for the _output_time
	// =======================================================
	_output = fopen(file_name, "r");
	if(_output){
	    new_file=1;
	    char Line[256];
	    while( fgets( Line, 256*sizeof(char), _output) )
	    {
	        int LineStartChar=0;
	        while( (Line[LineStartChar] == ' ') || (Line[LineStartChar] == '\t') )
	            LineStartChar++;
	        if(Line[LineStartChar] != '#'){
	            if(sscanf(Line, "%f",&_output_time) != 1) {
	                S->addMessageF(2, "Sensor file is bad formatted, and will be overwritten.\n");
	                S->addMessage(0, "\tConsider restart the simulation if you are starting from previous work.\n");
	                return 1;
	            }
	        }
	    }
	    fclose(_output);
	}

	if(_output)
	    _output = fopen(file_name, "wa");
	else
	    _output = fopen(file_name, "w");
	if(!_output){
	    S->addMessage(3, "The output file cannot be opened.\n");
	    S->addMessage(0, file_name);
	}

	if(new_file)
	    return false;
	fprintf(_output,"#########################################################\n");
	fprintf(_output,"#                                                       #\n");
	fprintf(_output,"#    #    ##   #  #   #                           #     #\n");
	fprintf(_output,"#   # #  #  #  #  #  # #                          #     #\n");
	fprintf(_output,"#  ##### #  #  #  # #####  ##  ###  #  #  ## ###  ###   #\n");
	fprintf(_output,"#  #   # #  #  #  # #   # #  # #  # #  # #   #  # #  #  #\n");
	fprintf(_output,"#  #   # #  #  #  # #   # #  # #  # #  #   # #  # #  #  #\n");
	fprintf(_output,"#  #   #  ## #  ##  #   #  ### ###   ### ##  ###  #  #  #\n");
	fprintf(_output,"#                            # #             #          #\n");
	fprintf(_output,"#                          ##  #             #          #\n");
	fprintf(_output,"#                                                       #\n");
	fprintf(_output,"#########################################################\n");
	fprintf(_output,"#\n");
	fprintf(_output,"# Another QUAlity GPU-SPH, by CEHINAV.\n");
	fprintf(_output,"#\thttp://canal.etsin.upm.es/\n");
	fprintf(_output,"#Authors:\n");
	fprintf(_output,"#\tCercós Pita, Jose Luis\n");
	fprintf(_output,"#\tMiguel Gonzalez, Leo\n");
	fprintf(_output,"#\tSouto Iglesias, Antonio\n");
	fprintf(_output,"#\n");
	fprintf(_output,"#########################################################\n");
	fprintf(_output,"#\n");
	fprintf(_output,"# Sensors output file.\n");
	fprintf(_output,"# Number of sensors: %d\n", _n);
	fprintf(_output,"# Rows:\n");
	fprintf(_output,"# Time");
	for(i=0;i<_n;i++)
	    #ifndef HAVE_3D
	        fprintf(_output,"\ts%d.X\ts%d.Y\ts%d.press\ts%d.dens\ts%d.sumW", i,i,i,i,i);
	    #else
	        fprintf(_output,"\ts%d.X\ts%d.Y\ts%d.Z\ts%d.press\ts%d.dens\ts%d.sumW", i,i,i,i,i,i);
	    #endif

	fprintf(_output,"\n#\n");
	fprintf(_output,"#  [s]");
	for(i=0;i<_n;i++)
	    #ifndef HAVE_3D
	        fprintf(_output,"\t  [m]\t  [m]\t   [Pa/m]\t [kg/m3]\t      []");
	    #else
	        fprintf(_output,"\t  [m]\t  [m]\t  [m]\t   [Pa/m]\t [kg/m3]\t      []");
	    #endif
	fprintf(_output,"\n#\n");
	fprintf(_output,"#########################################################\n");
	fflush(_output);
	return false;
}

}}  // namespace
