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
 * @brief OpenCL kernel kernel based tool.
 * (see Aqua::CalcServer::Kernel for details)
 */

#include <clang-c/Index.h>
#include <clang-c/Platform.h>
#include <AuxiliarMethods.h>
#include <InputOutput/Logger.h>
#include <CalcServer.h>
#include <CalcServer/Kernel.h>

namespace Aqua{ namespace CalcServer{

Kernel::Kernel(const std::string tool_name,
               const std::string kernel_path,
               const std::string entry_point,
               const std::string n,
               bool once)
    : Tool(tool_name, once)
    , _path(kernel_path)
    , _entry_point(entry_point)
    , _n(n)
    , _kernel(NULL)
    , _work_group_size(0)
    , _global_work_size(0)
{
}

Kernel::~Kernel()
{
    if(_kernel) clReleaseKernel(_kernel); _kernel=NULL;

    for(auto it = _var_values.begin(); it < _var_values.end(); it++){
        free(*it);
    }
}

void Kernel::setup()
{
    std::ostringstream msg;
    msg << "Loading the tool \"" << name()
        << "\" from the file \"" << path() << "\"..." << std::endl;
    LOG(L_INFO, msg.str());

    Tool::setup();
    make(_entry_point);
    variables(_entry_point);
    setVariables();
    computeGlobalWorkSize();
}

cl_event Kernel::_execute(const std::vector<cl_event> events)
{
    cl_int err_code;
    cl_event event;
    CalcServer *C = CalcServer::singleton();

    setVariables();
    computeGlobalWorkSize();

    cl_uint num_events_in_wait_list = events.size();
    const cl_event *event_wait_list = events.size() ? events.data() : NULL;

    err_code = clEnqueueNDRangeKernel(C->command_queue(),
                                      _kernel,
                                      1,
                                      NULL,
                                      &_global_work_size,
                                      &_work_group_size,
                                      num_events_in_wait_list,
                                      event_wait_list,
                                      &event);
    if(err_code != CL_SUCCESS){
        std::stringstream msg;
        msg << "Failure executing the tool \"" <<
               name() << "\"." << std::endl;
        LOG(L_ERROR, msg.str());
        InputOutput::Logger::singleton()->printOpenCLError(err_code);
        throw std::runtime_error("OpenCL execution error");
    }

    return event;
}

void Kernel::make(const std::string entry_point,
                  const std::string add_flags,
                  const std::string header)
{
    unsigned int i;
    cl_program program;
    cl_kernel kernel;
    std::ostringstream source;
    std::ostringstream flags;
    size_t source_length = 0;
    cl_int err_code = CL_SUCCESS;
    size_t work_group_size = 0;
    CalcServer *C = CalcServer::singleton();

    // Read the script file
    std::ifstream script(path());
    if(!script) {
        std::stringstream msg;
        msg << "Failure reading the file \"" <<
               path() << "\"." << std::endl;
        LOG(L_ERROR, msg.str());
        throw std::ifstream::failure(msg.str());
    }
    source << header << script.rdbuf();

    // Setup the default flags
    flags << "-I" << getFolderFromFilePath(path()) << " ";
    if(C->base_path().compare("")){
        flags << "-I" << C->base_path() << " ";
    }
    // Setup the user registered flags
    for(auto def : C->definitions()) {
        flags << def << " ";
    }
    // Add the additionally specified flags
    flags << add_flags;

    // Try to compile without using local memory
    LOG(L_INFO, "Compiling without local memory... ");
    kernel = compile_kernel(source.str(),
                            entry_point,
                            flags.str());

    // Get the work group size
    err_code = clGetKernelWorkGroupInfo(kernel,
                                        C->device(),
                                        CL_KERNEL_WORK_GROUP_SIZE,
                                        sizeof(size_t),
                                        &work_group_size,
                                        NULL);
    if(err_code != CL_SUCCESS) {
        LOG(L_ERROR, "Failure querying the work group size.\n");
        InputOutput::Logger::singleton()->printOpenCLError(err_code);
        clReleaseKernel(kernel);
        throw std::runtime_error("OpenCL error");
    }
    LOG0(L_DEBUG, "OK\n");

    _kernel = kernel;
    _work_group_size = work_group_size;

    // Try to compile with local memory
    LOG(L_INFO, "Compiling with local memory... ");
    flags << " -DLOCAL_MEM_SIZE=" << work_group_size;
    try {
        kernel = compile_kernel(source.str(),
                                entry_point,
                                flags.str());
    } catch(std::runtime_error& e) {
        LOG(L_INFO, "Falling back to no local memory usage.\n");
        return;
    }
    cl_ulong used_local_mem;
    err_code = clGetKernelWorkGroupInfo(kernel,
                                        C->device(),
                                        CL_KERNEL_LOCAL_MEM_SIZE,
                                        sizeof(cl_ulong),
                                        &used_local_mem,
                                        NULL);
    if(err_code != CL_SUCCESS) {
        LOG(L_ERROR, "Failure querying the used local memory.\n");
        InputOutput::Logger::singleton()->printOpenCLError(err_code);
        clReleaseKernel(kernel);
        LOG(L_INFO, "Falling back to no local memory usage.\n");
        return;
    }
    cl_ulong available_local_mem;
    err_code = clGetDeviceInfo(C->device(),
                               CL_DEVICE_LOCAL_MEM_SIZE,
                               sizeof(cl_ulong),
                               &available_local_mem,
                               NULL);
    if(err_code != CL_SUCCESS) {
        LOG(L_ERROR, "Failure querying the device available local memory.\n");
        InputOutput::Logger::singleton()->printOpenCLError(err_code);
        clReleaseKernel(kernel);
        LOG(L_INFO, "Falling back to no local memory usage.\n");
        return;
    }

    if(available_local_mem < used_local_mem){
        LOG(L_ERROR, "Not enough available local memory.\n");
        InputOutput::Logger::singleton()->printOpenCLError(err_code);
        clReleaseKernel(kernel);
        LOG(L_INFO, "Falling back to no local memory usage.\n");
        return;
    }
    LOG0(L_DEBUG, "OK\n");

    // Swap kernels
    err_code = clReleaseKernel(_kernel);
    if(err_code != CL_SUCCESS) {
        LOG(L_WARNING, "Failure releasing non-local memory based kernel.\n");
        InputOutput::Logger::singleton()->printOpenCLError(err_code);
    }
    _kernel = kernel;
}

/** @brief Main traverse method, which will parse all tokens except functions
 * declarations.
 * @param cursor the cursor whose child may be visited. All kinds of cursors can
 * be visited, including invalid cursors (which, by definition, have no
 * children).
 * @param parent the visitor function that will be invoked for each child of
 * parent.
 * @param client_data pointer data supplied by the client, which will be passed
 * to the visitor each time it is invoked.
 * @return CXChildVisit_Continue if a function declaration is found,
 * CXChildVisit_Recurse otherwise.
 */
CXChildVisitResult cursorVisitor(CXCursor cursor,
                                 CXCursor parent,
                                 CXClientData client_data);

/** Method traverse method, which will parse the input arguments.
 * @param cursor the cursor whose child may be visited. All kinds of cursors can
 * be visited, including invalid cursors (which, by definition, have no
 * children).
 * @param parent the visitor function that will be invoked for each child of
 * parent.
 * @param client_data pointer data supplied by the client, which will be passed
 * to the visitor each time it is invoked.
 * @return CXChildVisit_Continue.
 */
CXChildVisitResult functionDeclVisitor(CXCursor cursor,
                                       CXCursor parent,
                                       CXClientData client_data);

/** @struct clientData
 * @brief Data structure to store the variables requested and a flag
 * to know if the entry point has been found.
 */
struct clientData{
    /// Entry point
    std::string entry_point;
    /// Number of instances of the entry point found.
    unsigned int entry_points;
    /// List of required variables
    std::vector<std::string> var_names;
};

void Kernel::variables(const std::string entry_point)
{
    InputOutput::Variables *vars = CalcServer::singleton()->variables();

    CXIndex index = clang_createIndex(0, 0);
    if(index == 0){
        LOG(L_ERROR, "Failure creating parser index.\n");
        throw std::runtime_error("clang initialization failure");
    }

    int argc = 2;
    const char* argv[2] = {"Kernel", _path.c_str()};
    CXTranslationUnit translation_unit = clang_parseTranslationUnit(
        index,
        0,
        argv,
        argc,
        NULL,
        0,
        CXTranslationUnit_None);
    if(translation_unit == 0){
        LOG(L_ERROR, "Failure parsing the source code.\n");
        throw std::runtime_error("clang parsing error");
    }

    CXCursor root_cursor = clang_getTranslationUnitCursor(translation_unit);
    struct clientData client_data;
    client_data.entry_point = entry_point;
    client_data.entry_points = 0;
    client_data.var_names = _var_names;
    clang_visitChildren(root_cursor, *cursorVisitor, &client_data);
    if(client_data.entry_points == 0){
        std::stringstream msg;
        msg << "The entry point \"" << entry_point
            << "\" cannot be found." << std::endl;
        LOG(L_ERROR, msg.str());
        throw std::runtime_error("Invalid entry point");
    }
    if(client_data.entry_points != 1){
        std::stringstream msg;
        msg << "The entry point \"" << entry_point
            << "\" has been found " << client_data.entry_points
            << "times." << std::endl;
        LOG(L_ERROR, msg.str());
        throw std::runtime_error("Invalid entry point");
    }
    _var_names = client_data.var_names;
    // Retain just the array variables as dependencies, provided that scalar
    // variables are synced when passed using clSetKernelArg()
    std::vector<InputOutput::Variable*> deps;
    for(auto var_name : _var_names) {
        InputOutput::Variable *var = vars->get(var_name);
        if(!var){
            std::stringstream msg;
            msg << "The tool \"" << name()
                << "\" is asking the undeclared variable \""
                << var_name << "\"." << std::endl;
            LOG(L_ERROR, msg.str());
            throw std::runtime_error("Invalid variable");
        }
        if(var->isArray())
            deps.push_back(var);
    }
    setDependencies(deps);
    
    for(unsigned int i = 0; i < _var_names.size(); i++){
        _var_values.push_back(NULL);
    }

    clang_disposeTranslationUnit(translation_unit);
    clang_disposeIndex(index);
}

CXChildVisitResult cursorVisitor(CXCursor cursor,
                                 CXCursor parent,
                                 CXClientData client_data)
{
    struct clientData *data = (struct clientData *)client_data;
    CXCursorKind kind = clang_getCursorKind(cursor);
    CXString name = clang_getCursorSpelling(cursor);
    if (kind == CXCursor_FunctionDecl ||
        kind == CXCursor_ObjCInstanceMethodDecl)
    {
        if(!data->entry_point.compare(clang_getCString(name))){
            data->entry_points++;
            clang_visitChildren(cursor, *functionDeclVisitor, client_data);
        }
        return CXChildVisit_Continue;
    }
    return CXChildVisit_Recurse;
}

CXChildVisitResult functionDeclVisitor(CXCursor cursor,
                                       CXCursor parent,
                                       CXClientData client_data)
{
    struct clientData *data = (struct clientData *)client_data;
    CXCursorKind kind = clang_getCursorKind(cursor);
    if (kind == CXCursor_ParmDecl){
        CXString name = clang_getCursorSpelling(cursor);
        data->var_names.push_back(clang_getCString(name));
    }
    return CXChildVisit_Continue;
}

void Kernel::setVariables()
{
    unsigned int i;
    cl_int err_code;
    InputOutput::Variables *vars = CalcServer::singleton()->variables();

    for(i = 0; i < _var_names.size(); i++){
        InputOutput::Variable *var = vars->get(_var_names.at(i));
        if(!var){
            std::stringstream msg;
            msg << "The tool \"" << name()
                << "\" is asking the undeclared variable \""
                << _var_names.at(i) << "\"." << std::endl;
            LOG(L_ERROR, msg.str());
            throw std::runtime_error("Invalid variable");
        }
        if(_var_values.at(i) == NULL){
            _var_values.at(i) = malloc(var->typesize());
        }
        else if(!memcmp(_var_values.at(i), var->get(), var->typesize())){
            // The variable still being valid
            continue;
        }

        // Update the variable
        err_code = clSetKernelArg(_kernel, i, var->typesize(), var->get());
        if(err_code != CL_SUCCESS) {
            std::stringstream msg;
            msg << "Failure setting the variable \"" << _var_names.at(i)
                << "\" (id=" << i
                << ") to the tool \"" << name() << "\"." << std::endl;
            LOG(L_ERROR, msg.str());
            InputOutput::Logger::singleton()->printOpenCLError(err_code);
            throw std::runtime_error("OpenCL error");
        }
        memcpy(_var_values.at(i), var->get(), var->typesize());
    }
}

void Kernel::computeGlobalWorkSize()
{
    unsigned int N;
    if(!_work_group_size){
        LOG(L_ERROR, "Work group size must be greater than 0.\n");
        throw std::runtime_error("Null work group size");
    }
    InputOutput::Variables *vars = CalcServer::singleton()->variables();
    try {
        vars->solve("unsigned int", _n, &N);
    } catch(...) {
        LOG(L_ERROR, "Failure evaluating the number of threads.\n");
        throw std::runtime_error("Invalid number of threads");
    }

    _global_work_size = (size_t)roundUp(N, (unsigned int)_work_group_size);
}

}}  // namespace
