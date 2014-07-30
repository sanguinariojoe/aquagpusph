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
#include <CalcServer/Boundary/GhostParticles.h>
#include <CalcServer.h>

namespace Aqua{ namespace CalcServer{ namespace Boundary{

GhostParticles::GhostParticles()
	: Kernel("GhostParticles")
	, _path(0)
	, _program(0)
	, _kernel(0)
	, _is_delta(false)
{
	InputOutput::ScreenManager *S = InputOutput::ScreenManager::singleton();
	InputOutput::ProblemSetup *P = InputOutput::ProblemSetup::singleton();
    unsigned int i;
	if(!P->ghost_particles.walls.size())
	    return;
	int str_len = strlen(P->OpenCL_kernels.ghost);
	if(str_len <= 0) {
	    S->addMessageF(3, "The path of the kernel is empty.\n");
	    exit(EXIT_FAILURE);
	}
	_path = new char[str_len+4];
	if(!_path) {
	    S->addMessageF(3, "Memory cannot be allocated for the path.\n");
	    exit(EXIT_FAILURE);
	}
	strcpy(_path, P->OpenCL_kernels.ghost);
	strcat(_path, ".cl");
    for(i=0;i<P->n_fluids;i++){
        if(P->fluids[i].delta > 0.f){
            _is_delta = true;
            break;
        }
    }

	_local_work_size  = localWorkSize();
	if(!_local_work_size){
	    S->addMessageF(3, "I cannot get a valid local work size for the required computation tool.\n");
	    exit(EXIT_FAILURE);
	}
	_global_work_size = globalWorkSize(_local_work_size);
	if(setupOpenCL()) {
	    exit(EXIT_FAILURE);
	}
	S->addMessageF(1, "GhostParticles ready to work!\n");
}

GhostParticles::~GhostParticles()
{
	CalcServer *C = CalcServer::singleton();
	if(_kernel)clReleaseKernel(_kernel); _kernel=0;
	if(_program)clReleaseProgram(_program); _program=0;
	if(_path) delete[] _path; _path=0;
}

bool GhostParticles::execute()
{
	InputOutput::ProblemSetup *P = InputOutput::ProblemSetup::singleton();
	if(!P->ghost_particles.walls.size())  // Have at least one wall
	    return false;
	InputOutput::ScreenManager *S = InputOutput::ScreenManager::singleton();
	CalcServer *C = CalcServer::singleton();
	cl_int err_code=0;
	// Loop over walls
	unsigned int i;
	for(i=0;i<P->ghost_particles.walls.size();i++){
	    // Send constant variables to the server
	    err_code |= sendArgument(_kernel, 0, sizeof(cl_mem),
                                 (void*)&(C->ifluid));
	    err_code |= sendArgument(_kernel, 1, sizeof(cl_mem),
                                 (void*)&(C->imove));
	    err_code |= sendArgument(_kernel, 2, sizeof(cl_mem),
                                 (void*)&(C->pos));
	    err_code |= sendArgument(_kernel, 3, sizeof(cl_mem),
                                 (void*)&(C->v));
	    err_code |= sendArgument(_kernel, 4, sizeof(cl_mem),
                                 (void*)&(C->dens));
	    err_code |= sendArgument(_kernel, 5, sizeof(cl_mem),
                                 (void*)&(C->mass));
	    err_code |= sendArgument(_kernel, 6, sizeof(cl_mem),
                                 (void*)&(C->press));
	    err_code |= sendArgument(_kernel, 7, sizeof(cl_mem),
                                 (void*)&(C->visc_dyn_corrected));
	    err_code |= sendArgument(_kernel, 8, sizeof(cl_mem),
                                 (void*)&(C->refd));
	    err_code |= sendArgument(_kernel, 9, sizeof(cl_mem),
                                 (void*)&(C->f));
	    err_code |= sendArgument(_kernel, 10, sizeof(cl_mem),
                                 (void*)&(C->drdt));
	    err_code |= sendArgument(_kernel, 11, sizeof(cl_mem),
                                 (void*)&(C->drdt_F));
	    err_code |= sendArgument(_kernel, 12, sizeof(cl_mem),
                                 (void*)&(C->shepard));
	    err_code |= sendArgument(_kernel, 13, sizeof(cl_mem),
                                 (void*)&(C->icell));
	    err_code |= sendArgument(_kernel, 14, sizeof(cl_mem),
                                 (void*)&(C->ihoc));
	    err_code |= sendArgument(_kernel, 15, sizeof(cl_uint),
                                 (void*)&(C->N));
	    err_code |= sendArgument(_kernel, 16, sizeof(uivec),
                                 (void*)&(C->num_cells_vec));
	    err_code |= sendArgument(_kernel, 17, sizeof(vec),
                                 (void*)&(C->g));
        #ifdef HAVE_3D
            #define WALL_ARGS 9
            err_code |= sendArgument(_kernel, 18, sizeof(vec),
                (void*)&(P->ghost_particles.walls.at(i)->p1));
            err_code |= sendArgument(_kernel, 19, sizeof(vec),
                (void*)&(P->ghost_particles.walls.at(i)->p2));
            err_code |= sendArgument(_kernel, 20, sizeof(vec),
                (void*)&(P->ghost_particles.walls.at(i)->p3));
            err_code |= sendArgument(_kernel, 21, sizeof(vec),
                (void*)&(P->ghost_particles.walls.at(i)->p4));
            err_code |= sendArgument(_kernel, 22, sizeof(vec),
                (void*)&(P->ghost_particles.walls.at(i)->n));
            err_code |= sendArgument(_kernel, 23, sizeof(vec),
                (void*)&(P->ghost_particles.walls.at(i)->v1));
            err_code |= sendArgument(_kernel, 24, sizeof(vec),
                (void*)&(P->ghost_particles.walls.at(i)->v2));
            err_code |= sendArgument(_kernel, 25, sizeof(vec),
                (void*)&(P->ghost_particles.walls.at(i)->v3));
            err_code |= sendArgument(_kernel, 26, sizeof(vec),
                (void*)&(P->ghost_particles.walls.at(i)->v4));
        #else
            #define WALL_ARGS 5
            err_code |= sendArgument(_kernel, 18, sizeof(vec),
                (void*)&(P->ghost_particles.walls.at(i)->p1));
            err_code |= sendArgument(_kernel, 19, sizeof(vec),
                (void*)&(P->ghost_particles.walls.at(i)->p2));
            err_code |= sendArgument(_kernel, 20, sizeof(vec),
                (void*)&(P->ghost_particles.walls.at(i)->n));
            err_code |= sendArgument(_kernel, 21, sizeof(vec),
                (void*)&(P->ghost_particles.walls.at(i)->v1));
            err_code |= sendArgument(_kernel, 22, sizeof(vec),
                (void*)&(P->ghost_particles.walls.at(i)->v2));
        #endif // HAVE_3D
        if(_is_delta) {
            err_code |= sendArgument(_kernel, 18 + WALL_ARGS, sizeof(cl_mem),
                                     (void*)&(C->delta));
            err_code |= sendArgument(_kernel, 19 + WALL_ARGS, sizeof(cl_float),
                                     (void*)&(C->dt));
            err_code |= sendArgument(_kernel, 20 + WALL_ARGS, sizeof(cl_float),
                                     (void*)&(C->cs));
        }
	    if(err_code != CL_SUCCESS) {
	        S->addMessageF(3, "Failure sending variables to GhostParticles effect computation kernel.\n");
	        return true;
	    }

	    size_t globalWorkSize = getGlobalWorkSize(C->n, _local_work_size);
	    #ifdef HAVE_GPUPROFILE
	        cl_event event;
	        cl_ulong end, start;
	        err_code = clEnqueueNDRangeKernel(C->command_queue,
                                              _kernel,
                                              1,
                                              NULL,
                                              &globalWorkSize,
                                              NULL,
                                              0,
                                              NULL,
                                              &event);
	    #else
	        err_code = clEnqueueNDRangeKernel(C->command_queue,
                                              _kernel,
                                              1,
                                              NULL,
                                              &globalWorkSize,
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
	}

	return false;
}

bool GhostParticles::setupOpenCL()
{
	InputOutput::ScreenManager *S = InputOutput::ScreenManager::singleton();
	InputOutput::ProblemSetup *P = InputOutput::ProblemSetup::singleton();
	char args[256];
	strcpy(args, "");
	sprintf(args,
            "-D__PRESS_MODEL__=%u -D__NORMAL_U_MODEL__=%u -D__TANGENT_U_MODEL__=%u",
	        P->ghost_particles.p_extension,
	        P->ghost_particles.vn_extension,
	        P->ghost_particles.vt_extension);
	if(_is_delta)
        strcat(args, " -D__DELTA_SPH__");
	CalcServer *C = CalcServer::singleton();
	char msg[1024];
	cl_int err_code;
	cl_device_id device;
	cl_ulong local_mem, required_local_mem;
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
	err_code = clGetDeviceInfo(device,
                               CL_DEVICE_LOCAL_MEM_SIZE,
                               sizeof(local_mem),
                               &local_mem,
                               NULL);
	if(err_code != CL_SUCCESS) {
		S->addMessageF(3, "Failure getting the local memory available on the device.\n");
        S->printOpenCLError(err_code);
	    return true;
	}
	if(!loadKernelFromFile(&_kernel,
                           &_program,
                           C->context,
                           C->device,
                           _path,
                           "Boundary",
                           args))
	    return true;
	if(_program)clReleaseProgram(_program); _program=0;

	err_code = clGetKernelWorkGroupInfo(_kernel,
                                        device,
                                        CL_KERNEL_LOCAL_MEM_SIZE,
	                                    sizeof(cl_ulong),
                                        &required_local_mem,
                                        NULL);
	if(err_code != CL_SUCCESS) {
		S->addMessageF(3, "Error retrieving the used local memory.\n");
        S->printOpenCLError(err_code);
	    return true;
	}
	if(local_mem < required_local_mem){
		S->addMessageF(3, "There are not enough local memory in the device.\n");
	    sprintf(msg, "\tNeeds %lu bytes, but only %lu bytes are available.\n",
	           required_local_mem, local_mem);
	    S->addMessage(0, msg);
	    return true;
	}
	// Test if the local work gorup size must be modified
	size_t local_work_size=0;
	err_code = clGetKernelWorkGroupInfo(_kernel,
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

	err_code = clGetKernelWorkGroupInfo(_kernel,
                                        device,
                                        CL_KERNEL_PREFERRED_WORK_GROUP_SIZE_MULTIPLE,
	                                    sizeof(size_t),
                                        &local_work_size,
                                        NULL);
	if(err_code != CL_SUCCESS) {
		S->addMessageF(3, "I cannot query the preferred local work size");
        S->printOpenCLError(err_code);
	    return true;
	}
	_local_work_size  = (_local_work_size/local_work_size) * local_work_size;
	_global_work_size = globalWorkSize(_local_work_size);

	required_local_mem += _local_work_size*(sizeof(vec     )
	                                + sizeof(cl_float)
	                                + sizeof(cl_float)
	                                + sizeof(cl_float));

	if(local_mem < required_local_mem){
		S->addMessageF(2, "Not enough local memory for boundary.\n");
	    sprintf(msg, "\tNeeds %lu bytes, but only %lu bytes are available.\n",
	           required_local_mem, local_mem);
	    S->addMessage(0, msg);
	    S->addMessage(0, "\tLocal memory usage will be avoided therefore.\n");
	}
	else{
        sprintf(args, "%s -DLOCAL_MEM_SIZE=%lu", args, _local_work_size);
	}
    if(_kernel)clReleaseKernel(_kernel); _kernel=0;
    if(!loadKernelFromFile(&_kernel,
                           &_program,
                           C->context,
                           C->device,
                           _path,
                           "Boundary",
                           args))
        return true;
    if(_program)clReleaseProgram(_program); _program=0;
	return false;
}

}}}  // namespaces
