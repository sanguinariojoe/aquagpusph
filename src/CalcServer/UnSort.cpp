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

/** @file
 * @brief UnSort Recover the original id of each particle.
 * (See Aqua::CalcServer::UnSort for details)
 * @note Hardcoded versions of the files CalcServer/UnSort.cl.in and
 * CalcServer/UnSort.hcl.in are internally included as a text array.
 */

#include <stdlib.h>
#include <math.h>
#include <vector>

#include <ProblemSetup.h>
#include <ScreenManager.h>
#include <CalcServer/UnSort.h>
#include <CalcServer.h>

namespace Aqua{ namespace CalcServer{

#ifndef DOXYGEN_SHOULD_SKIP_THIS
#include "CalcServer/UnSort.hcl"
#include "CalcServer/UnSort.cl"
#endif
const char* UNSORT_INC = (const char*)UnSort_hcl_in;
unsigned int UNSORT_INC_LEN = UnSort_hcl_in_len;
const char* UNSORT_SRC = (const char*)UnSort_cl_in;
unsigned int UNSORT_SRC_LEN = UnSort_cl_in_len;


UnSort::UnSort(const char *name, const char *var_name)
    : Tool(name)
    , _var_name(NULL)
    , _var(NULL)
    , _input(NULL)
    , _id_input(NULL)
    , _output(NULL)
    , _kernel(NULL)
    , _global_work_size(0)
    , _local_work_size(0)
    , _n(0)
{
    _var_name = new char[strlen(var_name) + 1];
    strcpy(_var_name, var_name);
}

UnSort::~UnSort()
{
    if(_var_name) delete[] _var_name; _var_name=NULL;
    if(_output) clReleaseMemObject(_output); _output=NULL;
    if(_kernel) clReleaseKernel(_kernel); _kernel=NULL;
}

bool UnSort::setup()
{
    if(variables()){
        return true;
    }
    if(setupMem()){
        return true;
    }

    _id_input = *(cl_mem*)_id_var->get();
    _input = *(cl_mem*)_var->get();
    _n = _id_var->size() / InputOutput::Variables::typeToBytes(_id_var->type());
    if(setupOpenCL())
        return true;
    return false;
}


bool UnSort::_execute()
{
    unsigned int i;
    cl_int err_code;
    char msg[1024];
    InputOutput::ScreenManager *S = InputOutput::ScreenManager::singleton();
    CalcServer *C = CalcServer::singleton();

    if(setVariables()){
        return true;
    }

    // Execute the kernels
    err_code = clEnqueueNDRangeKernel(C->command_queue(),
                                      _kernel,
                                      1,
                                      NULL,
                                      &_global_work_size,
                                      &_local_work_size,
                                      0,
                                      NULL,
                                      NULL);
    if(err_code != CL_SUCCESS) {
        sprintf(msg, "Failure executing the tool \"%s\".\n", name());
        S->addMessageF(3, msg);
        S->printOpenCLError(err_code);
        return true;
    }

    return false;
}

bool UnSort::variables()
{
    char msg[1024];
    InputOutput::ScreenManager *S = InputOutput::ScreenManager::singleton();
    CalcServer *C = CalcServer::singleton();
    InputOutput::Variables *vars = C->variables();
    if(!vars->get("id")){
        sprintf(msg,
                "The tool \"%s\" is using the undeclared variable \"%s\".\n",
                name(),
                "id");
        S->addMessageF(3, msg);
        return true;
    }
    if(strcmp(vars->get("id")->type(), "unsigned int*")){
        sprintf(msg,
                "The tool \"%s\" has found a wrong type for the variable \"%s\".\n",
                name(),
                "id");
        S->addMessageF(3, msg);
        sprintf(msg,
                "\t\"%s\" was expected, but \"%s\" was found.\n",
                "unsigned int*",
                vars->get("id")->type());
        S->addMessage(0, msg);
        return true;
    }
    _id_var = (InputOutput::ArrayVariable *)vars->get("id");

    if(!vars->get(_var_name)){
        sprintf(msg,
                "The tool \"%s\" is using the undeclared variable \"%s\".\n",
                name(),
                _var_name);
        S->addMessageF(3, msg);
        return true;
    }
    if(!strchr(vars->get(_var_name)->type(), '*')){
        sprintf(msg,
                "The tool \"%s\" has received the scalar variable \"%s\".\n",
                name(),
                _var_name);
        S->addMessageF(3, msg);
        return true;
    }
    _var = (InputOutput::ArrayVariable *)vars->get(_var_name);

    return false;
}

bool UnSort::setupMem()
{
    cl_int err_code;
    char msg[1024];
    size_t len_id, len_var;
    InputOutput::ScreenManager *S = InputOutput::ScreenManager::singleton();
    CalcServer *C = CalcServer::singleton();

    len_id = _id_var->size() / InputOutput::Variables::typeToBytes(_id_var->type());
    len_var = _var->size() / InputOutput::Variables::typeToBytes(_var->type());
    if(len_id > len_var){
        sprintf(msg,
                "Wrong variable length in the tool \"%s\".\n",
                name());
        S->addMessageF(3, msg);
        sprintf(msg,
                "\t\"%s\" has a length %lu.\n",
                "id",
                len_id);
        S->addMessage(0, msg);
        sprintf(msg,
                "\t\"%s\" has a length %lu.\n",
                _var_name,
                len_var);
        S->addMessage(0, msg);
        return true;
    }

    _output = clCreateBuffer(C->context(),
                             CL_MEM_WRITE_ONLY,
                             len_id * InputOutput::Variables::typeToBytes(_var->type()),
                             NULL,
                               &err_code);
    if(err_code != CL_SUCCESS){
        sprintf(msg,
                "Failure allocating output memory in the tool \"%s\".\n",
                name());
        S->addMessageF(3, msg);
        S->printOpenCLError(err_code);
    }
    allocatedMemory(len_id * InputOutput::Variables::typeToBytes(_var->type()));

    return false;
}

bool UnSort::setupOpenCL()
{
    cl_int err_code;
    char msg[1024];
    InputOutput::ScreenManager *S = InputOutput::ScreenManager::singleton();
    CalcServer *C = CalcServer::singleton();
    InputOutput::Variables *vars = C->variables();

    // Create a header for the source code where the operation will be placed
    char header[UNSORT_INC_LEN + 128];
    strcpy(header, "");
    strncat(header, UNSORT_INC, UNSORT_INC_LEN);
    strcat(header, "");

    // Setup the complete source code
    char source[strlen(header) + strlen(UNSORT_SRC) + 1];
    strcpy(source, header);
    strncat(source, UNSORT_SRC, UNSORT_SRC_LEN);
    strcat(source, "");

    // Starts a dummy kernel in order to study the local size that can be used
    _kernel = compile(source);
    if(!_kernel){
        return true;
    }
    err_code = clGetKernelWorkGroupInfo(_kernel,
                                        C->device(),
                                        CL_KERNEL_WORK_GROUP_SIZE,
                                        sizeof(size_t),
                                        &_local_work_size,
                                        NULL);
    if(err_code != CL_SUCCESS) {
        S->addMessageF(3, "Failure querying the work group size.\n");
        S->printOpenCLError(err_code);
        return true;
    }
    if(_local_work_size < __CL_MIN_LOCALSIZE__){
        S->addMessageF(3, "UnSort cannot be performed.\n");
        sprintf(msg,
                "\t%lu elements can be executed, but __CL_MIN_LOCALSIZE__=%lu\n",
                _local_work_size,
                __CL_MIN_LOCALSIZE__);
        S->addMessage(0, msg);
        return true;
    }

    _global_work_size = roundUp(_n, _local_work_size);
    err_code = clSetKernelArg(_kernel,
                              0,
                              _id_var->typesize(),
                              _id_var->get());
    if(err_code != CL_SUCCESS){
        S->addMessageF(3, "Failure sending the IDs argument\n");
        S->printOpenCLError(err_code);
        return true;
    }
    err_code = clSetKernelArg(_kernel,
                              1,
                              _var->typesize(),
                              _var->get());
    if(err_code != CL_SUCCESS){
        S->addMessageF(3, "Failure sending the input array argument\n");
        S->printOpenCLError(err_code);
        return true;
    }
    err_code = clSetKernelArg(_kernel,
                              2,
                              sizeof(cl_mem),
                              (void*)&_output);
    if(err_code != CL_SUCCESS){
        S->addMessageF(3, "Failure sending the output array argument\n");
        S->printOpenCLError(err_code);
        return true;
    }
    err_code = clSetKernelArg(_kernel,
                              3,
                              sizeof(unsigned int),
                              (void*)&_n);
    if(err_code != CL_SUCCESS){
        S->addMessageF(3, "Failure sending the array size argument\n");
        S->printOpenCLError(err_code);
        return true;
    }

    return false;
}

cl_kernel UnSort::compile(const char* source)
{
    cl_int err_code;
    cl_program program;
    cl_kernel kernel;
    char msg[1024];
    InputOutput::ScreenManager *S = InputOutput::ScreenManager::singleton();
    CalcServer *C = CalcServer::singleton();

    char flags[512];
    if(!strcmp(_var->type(), "unsigned int*")){
        sprintf(flags,
                "-DT=%s",
                "uint*");
    }
    else{
        sprintf(flags,
                "-DT=%s",
                _var->type());
    }
    strcpy(strchr(flags, '*'), "");
    #ifdef AQUA_DEBUG
        strcat(flags, " -DDEBUG ");
    #else
        strcat(flags, " -DNDEBUG ");
    #endif
    strcat(flags, " -cl-mad-enable -cl-fast-relaxed-math");
    #ifdef HAVE_3D
        strcat(flags, " -DHAVE_3D");
    #else
        strcat(flags, " -DHAVE_2D");
    #endif
    size_t source_length = strlen(source) + 1;
    program = clCreateProgramWithSource(C->context(),
                                        1,
                                        (const char **)&source,
                                        &source_length,
                                        &err_code);
    if(err_code != CL_SUCCESS) {
        S->addMessageF(3, "Failure creating the OpenCL program.\n");
        S->printOpenCLError(err_code);
        return NULL;
    }
    err_code = clBuildProgram(program, 0, NULL, flags, NULL, NULL);
    if(err_code != CL_SUCCESS) {
        S->addMessageF(3, "Error compiling the source code\n");
        S->printOpenCLError(err_code);
        S->addMessage(3, "--- Build log ---------------------------------\n");
        size_t log_size = 0;
        clGetProgramBuildInfo(program,
                              C->device(),
                              CL_PROGRAM_BUILD_LOG,
                              0,
                              NULL,
                              &log_size);
        char *log = (char*)malloc(log_size + sizeof(char));
        if(!log){
            sprintf(msg,
                    "Failure allocating %lu bytes for the building log\n",
                    log_size);
            S->addMessage(3, msg);
            S->addMessage(3, "--------------------------------- Build log ---\n");
            return NULL;
        }
        strcpy(log, "");
        clGetProgramBuildInfo(program,
                              C->device(),
                              CL_PROGRAM_BUILD_LOG,
                              log_size,
                              log,
                              NULL);
        strcat(log, "\n");
        S->addMessage(0, log);
        S->addMessage(3, "--------------------------------- Build log ---\n");
        free(log); log=NULL;
        clReleaseProgram(program);
        return NULL;
    }
    kernel = clCreateKernel(program, "unsort", &err_code);
    clReleaseProgram(program);
    if(err_code != CL_SUCCESS) {
        S->addMessageF(3, "Failure creating the kernel.\n");
        S->printOpenCLError(err_code);
        return NULL;
    }

    return kernel;
}

bool UnSort::setVariables()
{
    char msg[1024];
    cl_int err_code;
    InputOutput::ScreenManager *S = InputOutput::ScreenManager::singleton();

    if(_id_input != *(cl_mem*)_id_var->get()){
        err_code = clSetKernelArg(_kernel,
                                  0,
                                  _id_var->typesize(),
                                  _id_var->get());
        if(err_code != CL_SUCCESS) {
            sprintf(msg,
                    "Failure setting the input variable \"%s\" to the tool \"%s\".\n",
                    _id_var->name(),
                    name());
            S->addMessageF(3, msg);
            S->printOpenCLError(err_code);
            return true;
        }
        _id_input = *(cl_mem *)_id_var->get();
    }
    if(_input != *(cl_mem*)_var->get()){
        err_code = clSetKernelArg(_kernel,
                                  1,
                                  _var->typesize(),
                                  _var->get());
        if(err_code != CL_SUCCESS) {
            sprintf(msg,
                    "Failure setting the input variable \"%s\" to the tool \"%s\".\n",
                    _var->name(),
                    name());
            S->addMessageF(3, msg);
            S->printOpenCLError(err_code);
            return true;
        }
        _input = *(cl_mem *)_var->get();
    }

    return false;
}

}}  // namespaces
