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


Set::Set(const std::string name,
         const std::string var_name,
         const std::string value)
    : Tool(name)
    , _var_name(var_name)
    , _value(value)
    , _var(NULL)
    , _input(NULL)
    , _kernel(NULL)
    , _global_work_size(0)
    , _local_work_size(0)
    , _n(0)
{
}

Set::~Set()
{
    if(_kernel) clReleaseKernel(_kernel); _kernel=NULL;
}

void Set::setup()
{
    std::ostringstream msg;
    msg << "Loading the tool \"" << name() << "\"..." << std::endl;
    LOG(L_INFO, msg.str());

    variable();

    _input = *(cl_mem*)_var->get();
    _n = _var->size() / InputOutput::Variables::typeToBytes(_var->type());
    setupOpenCL();
}


void Set::_execute()
{
    unsigned int i;
    cl_int err_code;
    CalcServer *C = CalcServer::singleton();

    setVariables();

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
        std::stringstream msg;
        msg << "Failure executing the tool \"" <<
               name() << "\"." << std::endl;
        LOG(L_ERROR, msg.str());
        InputOutput::ScreenManager::singleton()->printOpenCLError(err_code);
        throw std::runtime_error("OpenCL execution error");
    }
}

void Set::variable()
{
    InputOutput::Variables vars = CalcServer::singleton()->variables();
    if(!vars.get(_var_name)){
        std::stringstream msg;
        msg << "The tool \"" << name()
            << "\" is asking the undeclared variable \""
            << _var_name << "\"." << std::endl;
        LOG(L_ERROR, msg.str());
        throw std::runtime_error("Invalid variable");
    }
    if(vars.get(_var_name)->type().find('*') == std::string::npos){
        std::stringstream msg;
        msg << "The tool \"" << name()
            << "\" is asking the variable \"" << _var_name
            << "\", which is a scalar." << std::endl;
        LOG(L_ERROR, msg.str());
        throw std::runtime_error("Invalid variable type");
    }
    _var = (InputOutput::ArrayVariable *)vars.get(_var_name);
}

void Set::setupOpenCL()
{
    cl_int err_code;
    cl_kernel kernel;
    CalcServer *C = CalcServer::singleton();

    // Create a header for the source code where the operation will be placed
    std::ostringstream source;
    source << SET_INC << " #define VALUE " << _value << SET_SRC;

    // Starts a dummy kernel in order to study the local size that can be used
    kernel = compile(source.str());
    err_code = clGetKernelWorkGroupInfo(kernel,
                                        C->device(),
                                        CL_KERNEL_WORK_GROUP_SIZE,
                                        sizeof(size_t),
                                        &_local_work_size,
                                        NULL);
    if(err_code != CL_SUCCESS) {
        LOG(L_ERROR, "Failure querying the work group size.\n");
        InputOutput::ScreenManager::singleton()->printOpenCLError(err_code);
        clReleaseKernel(kernel);
        throw std::runtime_error("OpenCL error");
    }
    if(_local_work_size < __CL_MIN_LOCALSIZE__){
        LOG(L_ERROR, "insufficient local memory.\n");
        std::stringstream msg;
        msg << "\t" << _local_work_size
            << " local work group size with __CL_MIN_LOCALSIZE__="
            << __CL_MIN_LOCALSIZE__ << std::endl;
        LOG0(L_DEBUG, msg.str());
        throw std::runtime_error("OpenCL error");
    }

    _global_work_size = roundUp(_n, _local_work_size);
    _kernel = kernel;
    err_code = clSetKernelArg(kernel,
                              0,
                              _var->typesize(),
                              _var->get());
    if(err_code != CL_SUCCESS){
        LOG(L_ERROR, "Failure sending the array argument\n");
        InputOutput::ScreenManager::singleton()->printOpenCLError(err_code);
        throw std::runtime_error("OpenCL error");
    }
    err_code = clSetKernelArg(kernel,
                              1,
                              sizeof(unsigned int),
                              (void*)&_n);
    if(err_code != CL_SUCCESS){
        LOG(L_ERROR, "Failure sending the array size argument\n");
        InputOutput::ScreenManager::singleton()->printOpenCLError(err_code);
        throw std::runtime_error("OpenCL error");
    }
}

cl_kernel Set::compile(const std::string source)
{
    cl_int err_code;
    cl_program program;
    cl_kernel kernel;
    CalcServer *C = CalcServer::singleton();

    std::ostringstream flags;
    if(!_var->type().compare("unsigned int*")){
        // Spaces are not a good business into definitions passed as args
        flags << "-DT=uint";
    }
    else{
        std::string t = trimCopy(_var->type());
        t.pop_back();  // Remove the asterisk
        flags << "-DT=" << t;
        flags << "-DT=" << _var->type();
    }
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
        InputOutput::ScreenManager::singleton()->printOpenCLError(err_code);
        throw std::runtime_error("OpenCL compilation error");
    }
    err_code = clBuildProgram(program, 0, NULL, flags.str().c_str(), NULL, NULL);
    if(err_code != CL_SUCCESS) {
        LOG0(L_ERROR, "Error compiling the source code\n");
        InputOutput::ScreenManager::singleton()->printOpenCLError(err_code);
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
    kernel = clCreateKernel(program, "set", &err_code);
    clReleaseProgram(program);
    if(err_code != CL_SUCCESS) {
        LOG(L_ERROR, "Failure creating the kernel.\n");
        InputOutput::ScreenManager::singleton()->printOpenCLError(err_code);
        throw std::runtime_error("OpenCL error");
    }

    return kernel;
}

void Set::setVariables()
{
    cl_int err_code;

    if(_input == *(cl_mem*)_var->get()){
        return;
    }

    err_code = clSetKernelArg(_kernel,
                              0,
                              _var->typesize(),
                              _var->get());
    if(err_code != CL_SUCCESS) {
        std::stringstream msg;
        msg << "Failure setting the variable \"" << _var->name()
            << "\" to the tool \"" << name() << "\"." << std::endl;
        LOG(L_ERROR, msg.str());
        InputOutput::ScreenManager::singleton()->printOpenCLError(err_code);
        throw std::runtime_error("OpenCL error");
    }

    _input = *(cl_mem *)_var->get();
}

}}  // namespaces
