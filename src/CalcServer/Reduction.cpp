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
#include <InputOutput/Logger.h>
#include <CalcServer/Reduction.h>
#include <CalcServer.h>

namespace Aqua{ namespace CalcServer{

#ifndef DOXYGEN_SHOULD_SKIP_THIS
#include "CalcServer/Reduction.hcl"
#include "CalcServer/Reduction.cl"
#endif
std::string REDUCTION_INC = xxd2string(Reduction_hcl_in, Reduction_hcl_in_len);
std::string REDUCTION_SRC = xxd2string(Reduction_cl_in, Reduction_cl_in_len);


Reduction::Reduction(const std::string name,
                     const std::string input_name,
                     const std::string output_name,
                     const std::string operation,
                     const std::string null_val)
    : Tool(name)
    , _input_name(input_name)
    , _output_name(output_name)
    , _operation(operation)
    , _null_val(null_val)
    , _input_var(NULL)
    , _output_var(NULL)
    , _input(NULL)
{
}

Reduction::~Reduction()
{
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

void Reduction::setup()
{
    std::ostringstream msg;
    msg << "Loading the tool \"" << name() << "\"..." << std::endl;
    LOG(L_INFO, msg.str());

    variables();

    _mems.push_back(*(cl_mem*)_input_var->get());
    _input = *(cl_mem*)_input_var->get();
    size_t n = _input_var->size() / InputOutput::Variables::typeToBytes(
        _input_var->type());
    _n.push_back(n);
    setupOpenCL();
}

void Reduction::_execute()
{
    unsigned int i;
    cl_int err_code;
    CalcServer *C = CalcServer::singleton();
    InputOutput::Variables *vars = C->variables();

    setVariables();

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
            std::ostringstream msg;
            msg << "Failure executing the step " << i << " within the tool \""
                << name() << "\"." << std::endl;
            LOG(L_ERROR, msg.str());
            InputOutput::Logger::singleton()->printOpenCLError(err_code);
            throw std::runtime_error("OpenCL execution error");
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
            std::ostringstream msg;
            msg << "Failure reading back the result within the tool \""
                << name() << "\"." << std::endl;
            LOG(L_ERROR, msg.str());
            InputOutput::Logger::singleton()->printOpenCLError(err_code);
            throw std::runtime_error("OpenCL error");
    }

    // Ensure that the variable is populated
    vars->populate(_output_var);
}

void Reduction::variables()
{
    InputOutput::Variables *vars = CalcServer::singleton()->variables();
    if(!vars->get(_input_name)){
        std::stringstream msg;
        msg << "The tool \"" << name()
            << "\" is asking the undeclared input variable \""
            << _input_name << "\"." << std::endl;
        LOG(L_ERROR, msg.str());
        throw std::runtime_error("Invalid variable");
    }
    if(vars->get(_input_name)->type().find('*') == std::string::npos){
        std::stringstream msg;
        msg << "The tool \"" << name()
            << "\" is asking the input variable \"" << _input_name
            << "\", which is a scalar." << std::endl;
        LOG(L_ERROR, msg.str());
        throw std::runtime_error("Invalid variable type");
    }
    _input_var = (InputOutput::ArrayVariable *)vars->get(_input_name);
    if(!vars->get(_output_name)){
        std::stringstream msg;
        msg << "The tool \"" << name()
            << "\" is asking the undeclared output variable \""
            << _output_name << "\"." << std::endl;
        LOG(L_ERROR, msg.str());
        throw std::runtime_error("Invalid variable");
    }
    if(vars->get(_output_name)->type().find('*') != std::string::npos){
        std::stringstream msg;
        msg << "The tool \"" << name()
            << "\" is asking the output variable \"" << _output_name
            << "\", which is an array." << std::endl;
        LOG(L_ERROR, msg.str());
        throw std::runtime_error("Invalid variable type");
    }
    _output_var = vars->get(_output_name);
    if(!vars->isSameType(_input_var->type(), _output_var->type())){
        std::stringstream msg;
        msg << "Mismatching input and output types within the tool \"" << name()
            << "\"." << std::endl;
        LOG(L_ERROR, msg.str());
        msg.str("");
        msg << "\tInput variable \"" << _input_var->name()
            << "\" is of type \"" << _input_var->type()
            << "\"." << std::endl;
        LOG0(L_DEBUG, msg.str());
        msg << "\tOutput variable \"" << _output_var->name()
            << "\" is of type \"" << _output_var->type()
            << "\"." << std::endl;
        LOG0(L_DEBUG, msg.str());
        throw std::runtime_error("Invalid variable type");
    }
}

void Reduction::setupOpenCL()
{
    size_t data_size, local_size, max_local_size;
    cl_int err_code;
    cl_kernel kernel;
    CalcServer *C = CalcServer::singleton();
    InputOutput::Variables *vars = C->variables();

    // Get the elements data size to can allocate local memory later
    data_size = vars->typeToBytes(_input_var->type());

    std::ostringstream source;
    source << REDUCTION_INC << " #define IDENTITY " << _null_val << std::endl;
    source << "T reduce(T a, T b) " << std::endl;
    source << "{ " << std::endl;
    source << "    T c; " << std::endl;
    source << _operation << std::endl;
    source << "    return c; " << std::endl;
    source << "} " << std::endl;
    source << REDUCTION_SRC;

    // Starts a dummy kernel in order to study the local size that can be used
    local_size = __CL_MAX_LOCALSIZE__;
    kernel = compile(source.str(), local_size);
    err_code = clGetKernelWorkGroupInfo(kernel,
                                        C->device(),
                                        CL_KERNEL_WORK_GROUP_SIZE,
                                        sizeof(size_t),
                                        &max_local_size,
                                        NULL);
    if(err_code != CL_SUCCESS) {
        LOG(L_ERROR, "Failure querying the work group size.\n");
        InputOutput::Logger::singleton()->printOpenCLError(err_code);
        clReleaseKernel(kernel);
        throw std::runtime_error("OpenCL error");
    }
    if(max_local_size < __CL_MIN_LOCALSIZE__){
        LOG(L_ERROR, "insufficient local memory.\n");
        std::stringstream msg;
        msg << "\t" << max_local_size
            << " local work group size with __CL_MIN_LOCALSIZE__="
            << __CL_MIN_LOCALSIZE__ << std::endl;
        LOG0(L_DEBUG, msg.str());
        throw std::runtime_error("OpenCL error");
    }
    local_size = max_local_size;
    if(!isPowerOf2(local_size)){
        local_size = nextPowerOf2(local_size) / 2;
    }
    clReleaseKernel(kernel);

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
            std::stringstream msg;
            msg << "Failure allocating device memory in the tool \"" <<
                name() << "\"." << std::endl;
            LOG(L_ERROR, msg.str());
            InputOutput::Logger::singleton()->printOpenCLError(err_code);
            throw std::runtime_error("OpenCL allocation error");
        }
        allocatedMemory(_number_groups.at(i) * data_size + allocatedMemory());
        _mems.push_back(output);
        // Build the kernel
        kernel = compile(source.str(), local_size);
        _kernels.push_back(kernel);

        err_code = clSetKernelArg(kernel,
                                  0,
                                  sizeof(cl_mem),
                                  (void*)&(_mems.at(i)));
        if(err_code != CL_SUCCESS){
            LOG(L_ERROR, "Failure sending input argument\n");
            InputOutput::Logger::singleton()->printOpenCLError(err_code);
            throw std::runtime_error("OpenCL error");
        }
        err_code = clSetKernelArg(kernel,
                                  1,
                                  sizeof(cl_mem),
                                  (void*)&(_mems.at(i+1)));
        if(err_code != CL_SUCCESS){
            LOG(L_ERROR, "Failure sending output argument\n");
            InputOutput::Logger::singleton()->printOpenCLError(err_code);
            throw std::runtime_error("OpenCL error");
        }
        err_code = clSetKernelArg(kernel,
                                  2,
                                  sizeof(cl_uint),
                                  (void*)&(n));
        if(err_code != CL_SUCCESS){
            LOG(L_ERROR, "Failure sending number of threads argument\n");
            InputOutput::Logger::singleton()->printOpenCLError(err_code);
            throw std::runtime_error("OpenCL error");
        }
        err_code = clSetKernelArg(kernel,
                                  3,
                                  local_size*data_size ,
                                  NULL);
        if(err_code != CL_SUCCESS){
            LOG(L_ERROR, "Failure setting local memory\n");
            InputOutput::Logger::singleton()->printOpenCLError(err_code);
            throw std::runtime_error("OpenCL error");
        }
        // Setup next step
        std::stringstream msg;
        msg << "\tStep " << i << ", " << n << " elements reduced to "
            << _number_groups.at(i) << std::endl;
        LOG(L_DEBUG, msg.str());
        n = _number_groups.at(i);
        i++;
    }
}

cl_kernel Reduction::compile(const std::string source, size_t local_work_size)
{
    cl_int err_code;
    cl_program program;
    cl_kernel kernel;
    CalcServer *C = CalcServer::singleton();

    std::ostringstream flags;
    if(!_output_var->type().compare("unsigned int")){
        // Spaces are not a good business into definitions passed as args
        flags << "-DT=uint";
    }
    else{
        flags << "-DT=" << _output_var->type();
    }
    flags << " -DLOCAL_WORK_SIZE=" << local_work_size << "u";
    #ifdef AQUA_DEBUG
        flags << " -DDEBUG";
    #else
        flags << " -DNDEBUG";
    #endif
    flags << " -cl-mad-enable -cl-fast-relaxed-math";
    #ifdef HAVE_3D
        flags << " -DHAVE_3D";
    #else
        flags << " -DHAVE_2D";
    #endif

    size_t source_length = source.size();
    const char* source_cstr = source.c_str();
    program = clCreateProgramWithSource(C->context(),
                                        1,
                                        &source_cstr,
                                        &source_length,
                                        &err_code);
    if(err_code != CL_SUCCESS) {
        LOG(L_ERROR, "Failure creating the OpenCL program.\n");
        InputOutput::Logger::singleton()->printOpenCLError(err_code);
        throw std::runtime_error("OpenCL compilation error");
    }
    err_code = clBuildProgram(program, 0, NULL, flags.str().c_str(), NULL, NULL);
    if(err_code != CL_SUCCESS) {
        LOG0(L_ERROR, "Error compiling the source code\n");
        InputOutput::Logger::singleton()->printOpenCLError(err_code);
        LOG0(L_ERROR, "--- Build log ---------------------------------\n");
        size_t log_size = 0;
        clGetProgramBuildInfo(program,
                              C->device(),
                              CL_PROGRAM_BUILD_LOG,
                              0,
                              NULL,
                              &log_size);
        char *log = (char*)malloc(log_size + sizeof(char));
        if(!log){
            std::stringstream msg;
            msg << "Failure allocating " << log_size
                << " bytes for the building log" << std::endl;
            LOG0(L_ERROR, msg.str());
            LOG0(L_ERROR, "--------------------------------- Build log ---\n");
            throw std::bad_alloc();
        }
        strcpy(log, "");
        clGetProgramBuildInfo(program,
                              C->device(),
                              CL_PROGRAM_BUILD_LOG,
                              log_size,
                              log,
                              NULL);
        strcat(log, "\n");
        LOG0(L_DEBUG, log);
        LOG0(L_ERROR, "--------------------------------- Build log ---\n");
        free(log); log=NULL;
        clReleaseProgram(program);
        throw std::runtime_error("OpenCL compilation error");
    }
    kernel = clCreateKernel(program, "reduction", &err_code);
    clReleaseProgram(program);
    if(err_code != CL_SUCCESS) {
        LOG(L_ERROR, "Failure creating the kernel.\n");
        InputOutput::Logger::singleton()->printOpenCLError(err_code);
        throw std::runtime_error("OpenCL error");
    }

    return kernel;
}

void Reduction::setVariables()
{
    cl_int err_code;

    if(_input == *(cl_mem*)_input_var->get()){
        return;
    }

    err_code = clSetKernelArg(_kernels.at(0),
                              0,
                              _input_var->typesize(),
                              _input_var->get());
    if(err_code != CL_SUCCESS) {
        std::stringstream msg;
        msg << "Failure setting the input variable \"" << _input_var->name()
            << "\" to the tool \"" << name() << "\"." << std::endl;
        LOG(L_ERROR, msg.str());
        InputOutput::Logger::singleton()->printOpenCLError(err_code);
        throw std::runtime_error("OpenCL error");
    }

    _input = *(cl_mem *)_input_var->get();
    _mems.at(0) = _input;
}


}}  // namespaces
