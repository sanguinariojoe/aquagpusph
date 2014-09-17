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

#include <stdlib.h>
#include <math.h>
#include <vector>

#include <ProblemSetup.h>
#include <ScreenManager.h>
#include <CalcServer/Copy.h>
#include <CalcServer.h>

namespace Aqua{ namespace CalcServer{

#include "CalcServer/Copy.hcl"
#include "CalcServer/Copy.cl"
const char* COPY_INC = (const char*)Copy_hcl_in;
unsigned int COPY_INC_LEN = Copy_hcl_in_len;
const char* COPY_SRC = (const char*)Copy_cl_in;
unsigned int COPY_SRC_LEN = Copy_cl_in_len;


Copy::Copy(const char *name, const char *input_name, const char *output_name)
    : Tool(name)
    , _input_name(NULL)
    , _output_name(NULL)
    , _input_var(NULL)
    , _output_var(NULL)
    , _input(NULL)
    , _output(NULL)
    , _kernel(NULL)
    , _global_work_size(0)
    , _local_work_size(0)
    , _n(0)
{
    _input_name = new char[strlen(input_name) + 1];
    strcpy(_input_name, input_name);
    _output_name = new char[strlen(output_name) + 1];
    strcpy(_output_name, output_name);
}

Copy::~Copy()
{
    if(_input_name) delete[] _input_name; _input_name=NULL;
    if(_output_name) delete[] _output_name; _output_name=NULL;
    if(_kernel) clReleaseKernel(_kernel); _kernel=NULL;
}

bool Copy::setup()
{
    char msg[1024];
    InputOutput::ScreenManager *S = InputOutput::ScreenManager::singleton();
    CalcServer *C = CalcServer::singleton();
    InputOutput::Variables *vars = C->variables();

    sprintf(msg,
            "Loading the tool \"%s\"...\n",
            name());
    S->addMessageF(1, msg);

    if(variables()){
        return true;
    }

    _input = (cl_mem*)_input_var->get();
    _output = (cl_mem*)_input_var->get();
    _n = _input_var->size() / vars->typeToBytes(_input_var->type());
    if(setupOpenCL())
        return true;
    return false;
}


bool Copy::execute()
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

bool Copy::variables()
{
    char msg[1024];
    InputOutput::ScreenManager *S = InputOutput::ScreenManager::singleton();
    CalcServer *C = CalcServer::singleton();
    InputOutput::Variables *vars = C->variables();
    if(!vars->get(_input_name)){
        sprintf(msg,
                "The tool \"%s\" has received undeclared variable \"%s\" as input.\n",
                name(),
                _input_name);
        S->addMessageF(3, msg);
        return true;
    }
    if(!strchr(vars->get(_input_name)->type(), '*')){
        sprintf(msg,
                "The tool \"%s\" has received the scalar variable \"%s\" as input.\n",
                name(),
                _input_name);
        S->addMessageF(3, msg);
        return true;
    }
    _input_var = (InputOutput::ArrayVariable *)vars->get(_input_name);
    size_t n_in = _input_var->size() / vars->typeToBytes(_input_var->type());
    if(!vars->get(_output_name)){
        sprintf(msg,
                "The tool \"%s\" has received undeclared variable \"%s\" as output.\n",
                name(),
                _output_name);
        S->addMessageF(3, msg);
        return true;
    }
    if(!strchr(vars->get(_output_name)->type(), '*')){
        sprintf(msg,
                "The tool \"%s\" has received the scalar variable \"%s\" as output.\n",
                name(),
                _output_name);
        S->addMessageF(3, msg);
        return true;
    }
    _output_var = (InputOutput::ArrayVariable *)vars->get(_output_name);
    size_t n_out = _input_var->size() / vars->typeToBytes(_input_var->type());
    if(!vars->isSameType(_input_var->type(), _output_var->type())){
        sprintf(msg,
                "The input and output types mismatch for the tool \"%s\".\n",
                name());
        S->addMessageF(3, msg);
        sprintf(msg,
                "\tInput variable \"%s\" is of type \"%s\".\n",
                _input_var->name(),
                _input_var->type());
        S->addMessageF(0, msg);
        sprintf(msg,
                "\tOutput variable \"%s\" is of type \"%s\".\n",
                _output_var->name(),
                _output_var->type());
        S->addMessageF(0, msg);
        return true;
    }
    if(n_in != n_out){
        sprintf(msg,
                "The input and output lengths mismatch for the tool \"%s\".\n",
                name());
        S->addMessageF(3, msg);
        sprintf(msg,
                "\tInput variable \"%s\" has a length n=%lu.\n",
                _input_var->name(),
                n_in);
        S->addMessageF(0, msg);
        sprintf(msg,
                "\tOutput variable \"%s\" has a length n=%lu.\n",
                _output_var->name(),
                n_out);
        S->addMessageF(0, msg);
        return true;
    }

    return false;
}

bool Copy::setupOpenCL()
{
    cl_int err_code;
    cl_kernel kernel;
    char msg[1024];
    InputOutput::ScreenManager *S = InputOutput::ScreenManager::singleton();
    CalcServer *C = CalcServer::singleton();
    InputOutput::Variables *vars = C->variables();

    // Create a header for the source code where the operation will be placed
    char header[COPY_INC_LEN + 128];
    strcpy(header, "");
    strncat(header, COPY_INC, COPY_INC_LEN);
    strcat(header, "");

    // Setup the complete source code
    char source[strlen(header) + strlen(COPY_SRC) + 1];
    strcpy(source, header);
    strncat(source, COPY_SRC, COPY_SRC_LEN);
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
        S->addMessageF(3, "Copy cannot be performed.\n");
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
                              _input_var->typesize(),
                              _input_var->get());
    if(err_code != CL_SUCCESS){
        S->addMessageF(3, "Failure sending the input array argument\n");
        S->printOpenCLError(err_code);
        return true;
    }
    err_code = clSetKernelArg(kernel,
                              1,
                              _output_var->typesize(),
                              _output_var->get());
    if(err_code != CL_SUCCESS){
        S->addMessageF(3, "Failure sending the output array argument\n");
        S->printOpenCLError(err_code);
        return true;
    }
    err_code = clSetKernelArg(kernel,
                              2,
                              sizeof(unsigned int),
                              (void*)&_n);
    if(err_code != CL_SUCCESS){
        S->addMessageF(3, "Failure sending the array size argument\n");
        S->printOpenCLError(err_code);
        return true;
    }

    return false;
}

cl_kernel Copy::compile(const char* source)
{
    cl_int err_code;
    cl_program program;
    cl_kernel kernel;
    char msg[1024];
    InputOutput::ScreenManager *S = InputOutput::ScreenManager::singleton();
    CalcServer *C = CalcServer::singleton();

    char flags[512];
    sprintf(flags,
            "-DT=%s",
            _input_var->type());
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
        size_t log_size;
        clGetProgramBuildInfo(program,
                              C->device(),
                              CL_PROGRAM_BUILD_LOG,
                              0,
                              NULL,
                              &log_size);
        char *log = (char*)malloc(log_size + sizeof(char));
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
    kernel = clCreateKernel(program, "main", &err_code);
    clReleaseProgram(program);
    if(err_code != CL_SUCCESS) {
        S->addMessageF(3, "Failure creating the kernel.\n");
        S->printOpenCLError(err_code);
        return NULL;
    }

    return kernel;
}

bool Copy::setVariables()
{
    char msg[1024];
    cl_int err_code;
    InputOutput::ScreenManager *S = InputOutput::ScreenManager::singleton();

    if((void*)_input != _input_var->get()){
        err_code = clSetKernelArg(_kernel,
                                  0,
                                  _input_var->typesize(),
                                  _input_var->get());
        if(err_code != CL_SUCCESS) {
            sprintf(msg,
                    "Failure setting the input variable \"%s\" to the tool \"%s\".\n",
                    _input_var->name(),
                    name());
            S->addMessageF(3, msg);
            S->printOpenCLError(err_code);
            return true;
        }

        _input = (cl_mem *)_input_var->get();
    }

    if((void*)_output != _output_var->get()){
        err_code = clSetKernelArg(_kernel,
                                  1,
                                  _output_var->typesize(),
                                  _output_var->get());
        if(err_code != CL_SUCCESS) {
            sprintf(msg,
                    "Failure setting the input variable \"%s\" to the tool \"%s\".\n",
                    _output_var->name(),
                    name());
            S->addMessageF(3, msg);
            S->printOpenCLError(err_code);
            return true;
        }

        _output = (cl_mem *)_output_var->get();
    }

    return false;
}

}}  // namespaces
