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
 * @brief Set all the components of an array with the desired value.
 * (See Aqua::CalcServer::Set for details)
 * @note Hardcoded versions of the files CalcServer/Set.cl.in and
 * CalcServer/Set.hcl.in are internally included as a text array.
 */

#include <stdlib.h>
#include <math.h>
#include <vector>

#include <AuxiliarMethods.h>
#include <ProblemSetup.h>
#include <ScreenManager.h>
#include <CalcServer/Set.h>
#include <CalcServer.h>

namespace Aqua{ namespace CalcServer{

#ifndef DOXYGEN_SHOULD_SKIP_THIS
#include "CalcServer/Set.hcl"
#include "CalcServer/Set.cl"
#endif
const char* SET_INC = (const char*)Set_hcl_in;
unsigned int SET_INC_LEN = Set_hcl_in_len;
const char* SET_SRC = (const char*)Set_cl_in;
unsigned int SET_SRC_LEN = Set_cl_in_len;


Set::Set(const char *name, const char *var_name, const char *value)
    : Tool(name)
    , _var_name(NULL)
    , _value(NULL)
    , _var(NULL)
    , _input(NULL)
    , _kernel(NULL)
    , _global_work_size(0)
    , _local_work_size(0)
    , _n(0)
{
    _var_name = new char[strlen(var_name) + 1];
    strcpy(_var_name, var_name);
    _value = new char[strlen(value) + 1];
    strcpy(_value, value);
}

Set::~Set()
{
    if(_var_name) delete[] _var_name; _var_name=NULL;
    if(_value) delete[] _value; _value=NULL;
    if(_kernel) clReleaseKernel(_kernel); _kernel=NULL;
}

bool Set::setup()
{
    char msg[1024];
    InputOutput::ScreenManager *S = InputOutput::ScreenManager::singleton();

    sprintf(msg,
            "Loading the tool \"%s\"...\n",
            name());
    S->addMessageF(1, msg);

    if(variable()){
        return true;
    }

    _input = *(cl_mem*)_var->get();
    _n = _var->size() / InputOutput::Variables::typeToBytes(_var->type());
    if(setupOpenCL())
        return true;
    return false;
}


bool Set::_execute()
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

bool Set::variable()
{
    char msg[1024];
    InputOutput::ScreenManager *S = InputOutput::ScreenManager::singleton();
    CalcServer *C = CalcServer::singleton();
    InputOutput::Variables *vars = C->variables();
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

bool Set::setupOpenCL()
{
    cl_int err_code;
    cl_kernel kernel;
    char msg[1024];
    InputOutput::ScreenManager *S = InputOutput::ScreenManager::singleton();
    CalcServer *C = CalcServer::singleton();
    InputOutput::Variables *vars = C->variables();

    // Create a header for the source code where the operation will be placed
    char header[SET_INC_LEN + strlen(_value) + 128];
    strcpy(header, "");
    strncat(header, SET_INC, SET_INC_LEN);
    strcat(header, "");
    sprintf(header, "%s #define VALUE %s\n", header, _value);

    // Setup the complete source code
    char source[strlen(header) + strlen(SET_SRC) + 1];
    strcpy(source, header);
    strncat(source, SET_SRC, SET_SRC_LEN);
    strcat(source, "");

    // Starts a dummy kernel in order to study the local size that can be used
    kernel = compile(source);
    if(!kernel){
        return true;
    }
    err_code = clGetKernelWorkGroupInfo(kernel,
                                        C->device(),
                                        CL_KERNEL_WORK_GROUP_SIZE,
                                        sizeof(size_t),
                                        &_local_work_size,
                                        NULL);
    if(err_code != CL_SUCCESS) {
        S->addMessageF(3, "Failure querying the work group size.\n");
        S->printOpenCLError(err_code);
        clReleaseKernel(kernel);
        return true;
    }
    if(_local_work_size < __CL_MIN_LOCALSIZE__){
        S->addMessageF(3, "Set cannot be performed.\n");
        sprintf(msg,
                "\t%lu elements can be executed, but __CL_MIN_LOCALSIZE__=%lu\n",
                _local_work_size,
                __CL_MIN_LOCALSIZE__);
        S->addMessage(0, msg);
        return true;
    }

    _global_work_size = roundUp(_n, _local_work_size);
    _kernel = kernel;
    err_code = clSetKernelArg(kernel,
                              0,
                              _var->typesize(),
                              _var->get());
    if(err_code != CL_SUCCESS){
        S->addMessageF(3, "Failure sending the array argument\n");
        S->printOpenCLError(err_code);
        return true;
    }
    err_code = clSetKernelArg(kernel,
                              1,
                              sizeof(unsigned int),
                              (void*)&_n);
    if(err_code != CL_SUCCESS){
        S->addMessageF(3, "Failure sending the array size argument\n");
        S->printOpenCLError(err_code);
        return true;
    }

    return false;
}

cl_kernel Set::compile(const char* source)
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
        strcat(flags, " -g -DDEBUG ");
    #else
        strcat(flags, " -DNDEBUG ");
    #endif
    strcat(flags, " -cl-mad-enable -cl-no-signed-zeros -cl-finite-math-only -cl-fast-relaxed-math");
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
        S->addMessage(3, "Error compiling the source code\n");
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
    kernel = clCreateKernel(program, "set", &err_code);
    clReleaseProgram(program);
    if(err_code != CL_SUCCESS) {
        S->addMessageF(3, "Failure creating the kernel.\n");
        S->printOpenCLError(err_code);
        return NULL;
    }

    return kernel;
}

bool Set::setVariables()
{
    char msg[1024];
    cl_int err_code;
    InputOutput::ScreenManager *S = InputOutput::ScreenManager::singleton();

    if(_input == *(cl_mem*)_var->get()){
        return false;
    }

    err_code = clSetKernelArg(_kernel,
                              0,
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

    return false;
}

}}  // namespaces
