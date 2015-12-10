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
 * @brief Reductions, like scans, prefix sums, maximum or minimum, etc...
 * (See Aqua::CalcServer::Reduction for details)
 * @note Hardcoded versions of the files CalcServer/Reduction.cl.in and
 * CalcServer/Reduction.hcl.in are internally included as a text array.
 */

#include <stdlib.h>
#include <math.h>
#include <vector>

#include <ProblemSetup.h>
#include <ScreenManager.h>
#include <CalcServer/Reduction.h>
#include <CalcServer.h>

namespace Aqua{ namespace CalcServer{

#ifndef DOXYGEN_SHOULD_SKIP_THIS
#include "CalcServer/Reduction.hcl"
#include "CalcServer/Reduction.cl"
#endif
const char* REDUCTION_INC = (const char*)Reduction_hcl_in;
unsigned int REDUCTION_INC_LEN = Reduction_hcl_in_len;
const char* REDUCTION_SRC = (const char*)Reduction_cl_in;
unsigned int REDUCTION_SRC_LEN = Reduction_cl_in_len;


Reduction::Reduction(const char *name,
                     const char *input_name,
                     const char *output_name,
                     const char* operation,
                     const char* null_val)
    : Tool(name)
    , _input_name(NULL)
    , _output_name(NULL)
    , _operation(NULL)
    , _null_val(NULL)
    , _input_var(NULL)
    , _output_var(NULL)
    , _input(NULL)
{
    _input_name = new char[strlen(input_name) + 1];
    strcpy(_input_name, input_name);
    _output_name = new char[strlen(output_name) + 1];
    strcpy(_output_name, output_name);
    _operation = new char[strlen(operation) + 1];
    strcpy(_operation, operation);
    _null_val = new char[strlen(null_val) + 1];
    strcpy(_null_val, null_val);
}

Reduction::~Reduction()
{
    if(_input_name) delete[] _input_name; _input_name=NULL;
    if(_output_name) delete[] _output_name; _output_name=NULL;
    if(_operation) delete[] _operation; _operation=NULL;
    if(_null_val) delete[] _null_val; _null_val=NULL;

    unsigned int i;
    for(i=1;i<_mems.size();i++){
        if(_mems.at(i))
            clReleaseMemObject(_mems.at(i));
            _mems.at(i)=NULL;
    }
    for(i=0;i<_kernels.size();i++){
        if(_kernels.at(i))
            clReleaseKernel(_kernels.at(i));
        _kernels.at(i)=NULL;
    }
    _kernels.clear();
    _global_work_sizes.clear();
    _local_work_sizes.clear();
}

bool Reduction::setup()
{
    char msg[1024];
    InputOutput::ScreenManager *S = InputOutput::ScreenManager::singleton();

    sprintf(msg,
            "Loading the tool \"%s\"...\n",
            name());
    S->addMessageF(1, msg);

    if(variables()){
        return true;
    }

    _mems.push_back(*(cl_mem*)_input_var->get());
    _input = *(cl_mem*)_input_var->get();
    size_t n = _input_var->size() / InputOutput::Variables::typeToBytes(
        _input_var->type());
    _n.push_back(n);
    if(setupOpenCL())
        return true;
    return false;
}


bool Reduction::_execute()
{
    unsigned int i;
    cl_int err_code;
    char msg[1024];
    InputOutput::ScreenManager *S = InputOutput::ScreenManager::singleton();
    CalcServer *C = CalcServer::singleton();
    InputOutput::Variables *vars = C->variables();

    if(setVariables()){
        return true;
    }

    // Execute the kernels
    for(i = 0;i < _kernels.size(); i++){
        size_t _global_work_size = _global_work_sizes.at(i);
        size_t _local_work_size  = _local_work_sizes.at(i);
        err_code = clEnqueueNDRangeKernel(C->command_queue(),
                                          _kernels.at(i),
                                          1,
                                          NULL,
                                          &_global_work_size,
                                          &_local_work_size,
                                          0,
                                          NULL,
                                          NULL);
        if(err_code != CL_SUCCESS) {
            sprintf(msg,
                    "Failure executing the tool \"%s\" step %u.\n",
                    name(),
                    i);
            S->addMessageF(3, msg);
            S->printOpenCLError(err_code);
            return true;
        }
    }

    // Get back the result
    err_code = clEnqueueReadBuffer(C->command_queue(),
                                   _mems.at(_mems.size()-1),
                                   CL_TRUE,
                                   0,
                                   _output_var->typesize(),
                                   _output_var->get(),
                                   0,
                                   NULL,
                                   NULL);
    if(err_code != CL_SUCCESS) {
        sprintf(msg,
                "Failure in tool \"%s\" when reading back the reduced result.\n",
                name());
        S->addMessageF(3, msg);
        S->printOpenCLError(err_code);
        return true;
    }

    // Ensure that the variable is populated
    if(vars->populate(_output_var)){
        return true;
    }
    return false;
}

bool Reduction::variables()
{
    char msg[1024];
    InputOutput::ScreenManager *S = InputOutput::ScreenManager::singleton();
    CalcServer *C = CalcServer::singleton();
    InputOutput::Variables *vars = C->variables();
    if(!vars->get(_input_name)){
        sprintf(msg,
                "The tool \"%s\" has received the undeclared variable \"%s\" as input.\n",
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
    if(!vars->get(_output_name)){
        sprintf(msg,
                "The tool \"%s\" has received the undeclared variable \"%s\" as output.\n",
                name(),
                _output_name);
        S->addMessageF(3, msg);
        return true;
    }
    if(strchr(vars->get(_output_name)->type(), '*')){
        sprintf(msg,
                "The tool \"%s\" has received the array variable \"%s\" as output.\n",
                name(),
                _output_name);
        S->addMessageF(3, msg);
        return true;
    }
    _output_var = vars->get(_output_name);
    if(!vars->isSameType(_input_var->type(), _output_var->type())){
        sprintf(msg,
                "The input and output types mismatch for the tool \"%s\".\n",
                name());
        S->addMessageF(3, msg);
        sprintf(msg,
                "\tInput variable \"%s\" is of type \"%s\".\n",
                _input_var->name(),
                _input_var->type());
        S->addMessage(0, msg);
        sprintf(msg,
                "\tOutput variable \"%s\" is of type \"%s\".\n",
                _output_var->name(),
                _output_var->type());
        S->addMessage(0, msg);
        return true;
    }
    return false;
}

bool Reduction::setupOpenCL()
{
    size_t data_size, local_size, max_local_size;
    cl_int err_code;
    cl_kernel kernel;
    char msg[1024];
    InputOutput::ScreenManager *S = InputOutput::ScreenManager::singleton();
    CalcServer *C = CalcServer::singleton();
    InputOutput::Variables *vars = C->variables();

    // Get the elements data size to can allocate local memory later
    data_size = vars->typeToBytes(_input_var->type());

    // Create a header for the source code where the operation will be placed
    char header[REDUCTION_INC_LEN + strlen(_operation) + strlen(_null_val) + 128];
    strcpy(header, "");
    strncat(header, REDUCTION_INC, REDUCTION_INC_LEN);
    sprintf(header, "%s #define IDENTITY %s\n", header, _null_val);
    strcat(header, "T reduce(T a, T b) \n");
    strcat(header, "{ \n");
    strcat(header, "\tT c; \n");
    strcat(header, _operation);
    strcat(header, "\n");
    strcat(header, "\treturn c;\n");
    strcat(header, "} \n");

    // Setup the complete source code
    char source[strlen(header) + strlen(REDUCTION_SRC) + 1];
    strcpy(source, header);
    strncat(source, REDUCTION_SRC, REDUCTION_SRC_LEN);
    strcat(source, "");

    // Starts a dummy kernel in order to study the local size that can be used
    local_size = __CL_MAX_LOCALSIZE__;
    kernel = compile(source, local_size);
    if(!kernel){
        return true;
    }
    err_code = clGetKernelWorkGroupInfo(kernel,
                                        C->device(),
                                        CL_KERNEL_WORK_GROUP_SIZE,
                                        sizeof(size_t),
                                        &max_local_size,
                                        NULL);
    if(err_code != CL_SUCCESS) {
        S->addMessageF(3, "Failure querying the work group size.\n");
        S->printOpenCLError(err_code);
        clReleaseKernel(kernel);
        return true;
    }
    if(max_local_size < __CL_MIN_LOCALSIZE__){
        S->addMessageF(3, "Reduction cannot be performed (insufficient local memory)\n");
        sprintf(msg,
                "\t%lu elements can be executed, but __CL_MIN_LOCALSIZE__=%lu\n",
                max_local_size,
                __CL_MIN_LOCALSIZE__);
        S->addMessage(0, msg);
        return true;
    }
    local_size = max_local_size;
    if(!isPowerOf2(local_size)){
        local_size = nextPowerOf2(local_size) / 2;
    }

    if(kernel)clReleaseKernel(kernel); kernel=NULL;

    // Now we can start a loop while the amount of reduced data is greater than
    // one
    unsigned int n = _n.at(0);
    _n.clear();
    unsigned int i=0;
    while(n > 1){
        // Get work sizes
        _n.push_back(n);
        _local_work_sizes.push_back(local_size);
        _global_work_sizes.push_back(roundUp(n, local_size));
        _number_groups.push_back(
            _global_work_sizes.at(i) / _local_work_sizes.at(i)
        );
        // Build the output memory object
        cl_mem output = NULL;
        output = clCreateBuffer(C->context(),
                                CL_MEM_READ_WRITE,
                                _number_groups.at(i) * data_size,
                                NULL,
                                &err_code);
        if(err_code != CL_SUCCESS) {
            S->addMessageF(3, "Buffer memory allocation failure.\n");
            S->printOpenCLError(err_code);
            return true;
        }
        allocatedMemory(_number_groups.at(i) * data_size + allocatedMemory());
        _mems.push_back(output);
        // Build the kernel
        kernel = compile(source, local_size);
        if(!kernel){
            return true;
        }
        _kernels.push_back(kernel);

        err_code = clSetKernelArg(kernel,
                                  0,
                                  sizeof(cl_mem),
                                  (void*)&(_mems.at(i)));
        if(err_code != CL_SUCCESS){
            S->addMessageF(3, "Failure sending input argument\n");
            S->printOpenCLError(err_code);
            return true;
        }
        err_code = clSetKernelArg(kernel,
                                  1,
                                  sizeof(cl_mem),
                                  (void*)&(_mems.at(i+1)));
        if(err_code != CL_SUCCESS){
            S->addMessageF(3, "Failure sending output argument\n");
            S->printOpenCLError(err_code);
            return true;
        }
        err_code = clSetKernelArg(kernel,
                                  2,
                                  sizeof(cl_uint),
                                  (void*)&(n));
        if(err_code != CL_SUCCESS){
            S->addMessageF(3, "Failure sending number of threads argument\n");
            S->printOpenCLError(err_code);
            return true;
        }
        err_code = clSetKernelArg(kernel,
                                  3,
                                  local_size*data_size ,
                                  NULL);
        if(err_code != CL_SUCCESS){
            S->addMessageF(3, "Failure setting local memory\n");
            S->printOpenCLError(err_code);
            return true;
        }
        // Setup next step
        sprintf(msg,
                "\tStep %u, %u elements reduced to %u\n",
                i,
                n,
                _number_groups.at(i));
        S->addMessage(0, msg);
        n = _number_groups.at(i);
        i++;
    }
    return false;
}

cl_kernel Reduction::compile(const char* source, size_t local_work_size)
{
    cl_int err_code;
    cl_program program;
    cl_kernel kernel;
    char msg[1024];
    InputOutput::ScreenManager *S = InputOutput::ScreenManager::singleton();
    CalcServer *C = CalcServer::singleton();

    char flags[512];
    if(!strcmp(_output_var->type(), "unsigned int")){
        sprintf(flags,
                "-DT=%s -DLOCAL_WORK_SIZE=%luu",
                "uint",
                local_work_size);
    }
    else{
        sprintf(flags,
                "-DT=%s -DLOCAL_WORK_SIZE=%luu",
                _output_var->type(),
                local_work_size);
    }
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
    kernel = clCreateKernel(program, "reduction", &err_code);
    clReleaseProgram(program);
    if(err_code != CL_SUCCESS) {
        S->addMessageF(3, "Failure creating the kernel.\n");
        S->printOpenCLError(err_code);
        return NULL;
    }

    return kernel;
}

bool Reduction::setVariables()
{
    char msg[1024];
    cl_int err_code;
    InputOutput::ScreenManager *S = InputOutput::ScreenManager::singleton();

    if(_input == *(cl_mem*)_input_var->get()){
        return false;
    }

    err_code = clSetKernelArg(_kernels.at(0),
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

    _input = *(cl_mem *)_input_var->get();
    _mems.at(0) = _input;

    return false;
}


}}  // namespaces
