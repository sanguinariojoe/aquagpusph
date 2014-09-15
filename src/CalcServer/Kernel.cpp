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
 * @brief Aqua::CalcServer tools base class.
 * (See Aqua::CalcServer::Kernel for details)
 */

#include <CalcServer/Kernel.h>
#include <CalcServer.h>
#include <ScreenManager.h>

namespace Aqua{ namespace CalcServer{

Kernel::Kernel(const char* tool_name, const char* kernel_path)
    : Tool(tool_name)
    , _path(NULL)
    , _kernel(NULL)
    , _work_group_size(0)
{
    path(kernel_path);
}

Kernel::~Kernel()
{
    unsigned int i;
    if(_path) delete[] _path; _path=NULL;
    if(_kernel) clReleaseKernel(_kernel); _kernel=NULL;
    for(i = 0; i < _var_names.size(); i++){
        delete[] _var_names.at(i);
    }
    _var_names.clear();
}

bool Kernel::setup()
{
    char msg[1024];
    InputOutput::ScreenManager *S = InputOutput::ScreenManager::singleton();

    sprintf(msg,
            "Loading the tool \"%s\" from the file \"%s\"...\n",
            name(),
            path());
    S->addMessageF(1, msg);

    if(compile()){
        return true;
    }

    if(variables()){
        return true;
    }

    if(setVariables()){
        return true;
    }

    return false;
}

bool Kernel::execute()
{
    return false;
}

void Kernel::path(const char* kernel_path)
{
    if(_path) delete[] _path; _path=NULL;
    _path = new char[strlen(kernel_path) + 1];
    strcpy(_path, kernel_path);
}

bool Kernel::compile(const char* entry_point,
                     const char* add_flags,
                     const char* header)
{
    cl_program program;
    cl_kernel kernel;
    char* source = NULL;
    char* flags = NULL;
    size_t source_length = 0;
    cl_int err_code = CL_SUCCESS;
    size_t work_group_size = 0;
    char msg[1024];
    InputOutput::ScreenManager *S = InputOutput::ScreenManager::singleton();
    CalcServer *C = CalcServer::singleton();

    // Allocate the required memory for the source code
    source_length = readFile(NULL, path());
    if(!source_length){
        return true;
    }
    source = new char[source_length + 1];
    if(!source){
        S->addMessageF(3, "Failure allocating memory for the source code.\n");
        return true;
    }
    // Read the file
    source_length = readFile(source, path());
    if(!source_length){
        delete[] source; source=NULL;
        return true;
    }
    // Append the header on top of the source code
    if(header){
        char *backup = source;
        source_length += strlen(header) * sizeof(char);
        source = new char[source_length + 1];
        if(!source) {
            S->addMessageF(3, "Failure allocate memory to append the header.\n");
            return true;
        }
        strcpy(source, header);
        strcat(source, backup);
        delete[] backup; backup=NULL;
    }

    // Setup the default flags
    flags = new char[1024];
    #ifdef AQUA_DEBUG
        strcpy(flags, "-g -DDEBUG ");
    #else
        strcpy(flags, "-DNDEBUG ");
    #endif
    strcat(flags, "-I");
    const char *folder = getFolderFromFilePath(path());
    strcat(flags, folder);

    strcat(flags, " -cl-mad-enable -cl-no-signed-zeros -cl-finite-math-only -cl-fast-relaxed-math ");
    #ifdef HAVE_3D
        strcat(flags, " -DHAVE_3D ");
    #else
        strcat(flags, " -DHAVE_2D ");
    #endif
    // Add the additionally specified flags
    if(add_flags)
        strcat(flags, add_flags);

    // Try to compile without using local memory
    S->addMessageF(1, "Compiling without local memory... ");
    program = clCreateProgramWithSource(C->context(),
                                        1,
                                        (const char **)&source,
                                        &source_length,
                                        &err_code);
    if(err_code != CL_SUCCESS) {
        S->addMessage(0, "FAIL\n");
        S->addMessageF(3, "Failure creating the OpenCL program.\n");
        S->printOpenCLError(err_code);
        delete[] flags; flags=NULL;
        delete[] source; source=NULL;
        return true;
    }
    err_code = clBuildProgram(program, 0, NULL, flags, NULL, NULL);
    if(err_code != CL_SUCCESS) {
        S->addMessage(0, "FAIL\n");
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
        delete[] flags; flags=NULL;
        delete[] source; source=NULL;
        clReleaseProgram(program);
        return true;
    }
    kernel = clCreateKernel(program, entry_point, &err_code);
    clReleaseProgram(program);
    if(err_code != CL_SUCCESS) {
        S->addMessage(0, "FAIL\n");
        S->addMessageF(3, "Failure creating the kernel.\n");
        S->printOpenCLError(err_code);
        delete[] flags; flags=NULL;
        return true;
    }

    // Get the work group size
    err_code = clGetKernelWorkGroupInfo(kernel,
                                        C->device(),
                                        CL_KERNEL_WORK_GROUP_SIZE,
                                        sizeof(size_t),
                                        &work_group_size,
                                        NULL);
    if(err_code != CL_SUCCESS) {
        S->addMessage(0, "FAIL\n");
        S->addMessageF(3, "Failure querying the work group size.\n");
        S->printOpenCLError(err_code);
        clReleaseKernel(kernel);
        delete[] flags; flags=NULL;
        return true;
    }
    S->addMessage(0, "OK\n");

    _kernel = kernel;
    _work_group_size = work_group_size;

    // Try to compile with local memory
    S->addMessageF(1, "Compiling with local memory... ");
    program = clCreateProgramWithSource(C->context(),
                                        1,
                                        (const char **)&source,
                                        &source_length,
                                        &err_code);
    delete[] source; source=NULL;
    if(err_code != CL_SUCCESS) {
        S->addMessage(0, "FAIL\n");
        S->addMessageF(3, "Failure creating the OpenCL program.\n");
        S->printOpenCLError(err_code);
        delete[] flags; flags=NULL;
        S->addMessageF(1, "Falling back to no local memory usage.\n");
        return false;
    }
    sprintf(flags, "%s -DLOCAL_MEM_SIZE=%lu", flags, work_group_size);
    err_code = clBuildProgram(program, 0, NULL, flags, NULL, NULL);
    delete[] flags; flags=NULL;
    if(err_code != CL_SUCCESS) {
        S->addMessage(0, "FAIL\n");
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
        S->addMessageF(1, "Falling back to no local memory usage.\n");
        return false;
    }
    kernel = clCreateKernel(program, entry_point, &err_code);
    clReleaseProgram(program);
    if(err_code != CL_SUCCESS) {
        S->addMessage(0, "FAIL\n");
        S->addMessageF(3, "Failure creating the kernel.\n");
        S->printOpenCLError(err_code);
        S->addMessageF(1, "Falling back to no local memory usage.\n");
        return false;
    }
    cl_ulong used_local_mem;
    err_code = clGetKernelWorkGroupInfo(kernel,
                                        C->device(),
                                        CL_KERNEL_LOCAL_MEM_SIZE,
                                        sizeof(cl_ulong),
                                        &used_local_mem,
                                        NULL);
    if(err_code != CL_SUCCESS) {
        S->addMessage(0, "FAIL\n");
        S->addMessageF(3, "Failure querying the used local memory.\n");
        S->printOpenCLError(err_code);
        S->addMessageF(1, "Falling back to no local memory usage.\n");
        return false;
    }
    cl_ulong available_local_mem;
    err_code = clGetDeviceInfo(C->device(),
                               CL_DEVICE_LOCAL_MEM_SIZE,
                               sizeof(cl_ulong),
                               &available_local_mem,
                               NULL);
    if(err_code != CL_SUCCESS) {
        S->addMessage(0, "FAIL\n");
        S->addMessageF(3, "Failure querying the available local memory.\n");
        S->printOpenCLError(err_code);
        S->addMessageF(1, "Falling back to no local memory usage.\n");
        return false;
    }

    if(available_local_mem < used_local_mem){
        S->addMessage(0, "FAIL\n");
        S->addMessageF(3, "Not enough available local memory.\n");
        S->printOpenCLError(err_code);
        S->addMessageF(1, "Falling back to no local memory usage.\n");
        return false;
    }
    S->addMessage(0, "OK\n");
    clReleaseKernel(_kernel);
    _kernel = kernel;

    return false;
}

/** @brief Main traverse method, which will parse all tokens except functions
 * declarations.
 * @param the cursor whose child may be visited. All kinds of cursors can be
 * visited, including invalid cursors (which, by definition, have no
 * children).
 * @param the visitor function that will be invoked for each child of
 * parent.
 * @param pointer data supplied by the client, which will be passed to the
 * visitor each time it is invoked.
 * @return CXChildVisit_Continue if a function declaration is found,
 * CXChildVisit_Recurse otherwise.
 */
CXChildVisitResult cursorVisitor(CXCursor cursor,
                                 CXCursor parent,
                                 CXClientData client_data);

/** @brief Method traverse method, which will parse the input arguments.
 * @param the cursor whose child may be visited. All kinds of cursors can be
 * visited, including invalid cursors (which, by definition, have no
 * children).
 * @param the visitor function that will be invoked for each child of
 * parent.
 * @param pointer data supplied by the client, which will be passed to the
 * visitor each time it is invoked.
 * @return CXChildVisit_Continue.
 */
CXChildVisitResult functionDeclVisitor(CXCursor cursor,
                                       CXCursor parent,
                                       CXClientData client_data);

/** @struct Data structure to store the variables requested and a flag
 * to know if the entry point has been found.
 */
struct clientData{
    /// Entry point
    const char* entry_point;
    /// Number of instances of the entry point found.
    unsigned int entry_points;
    /// List of required variables
    std::deque<char*> *var_names;
};

bool Kernel::variables(const char* entry_point)
{
	char msg[1024];
	InputOutput::ScreenManager *S = InputOutput::ScreenManager::singleton();

    CXIndex index = clang_createIndex(0, 0);
    if(index == 0){
        S->addMessageF(3, "Failure creating parser index.\n");
        return true;
    }

    int argc = 2;
    const char* argv[2] = {"Kernel", path()};
    CXTranslationUnit translation_unit = clang_parseTranslationUnit(
        index,
        0,
        argv,
        argc,
        0,
        0,
        CXTranslationUnit_None);
    if(translation_unit == 0){
        S->addMessageF(3, "Failure parsing the source code.\n");
        return true;
    }

    CXCursor root_cursor = clang_getTranslationUnitCursor(translation_unit);
    struct clientData client_data;
    client_data.entry_point = entry_point;
    client_data.entry_points = 0;
    client_data.var_names = &_var_names;
    clang_visitChildren(root_cursor, *cursorVisitor, &client_data);
    if(client_data.entry_points == 0){
        sprintf(msg, "The entry point \"%s\" cannot be found.\n", entry_point);
        S->addMessageF(3, msg);
        return true;
    }
    if(client_data.entry_points != 1){
        sprintf(msg,
                "Entry point \"%s\" found %u times.\n",
                entry_point,
                client_data.entry_points);
        S->addMessageF(3, msg);
        return true;
    }

    unsigned int i;
	for(i = 0; i < _var_names.size(); i++){
        _var_sizes.push_back(0);
        _var_values.push_back(NULL);
	}

    clang_disposeTranslationUnit(translation_unit);
    clang_disposeIndex(index);
    return false;
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
        if(!strcmp(clang_getCString(name), data->entry_point)){
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
        std::deque<char*> *var_names = data->var_names;
        char *var_name = new char[strlen(clang_getCString(name)) + 1];
        strcpy(var_name, clang_getCString(name));
        var_names->push_back(var_name);
    }
    return CXChildVisit_Continue;
}

bool Kernel::setVariables()
{
    unsigned int i;
    char msg[1024];
    cl_int err_code;
	InputOutput::ScreenManager *S = InputOutput::ScreenManager::singleton();
	CalcServer *C = CalcServer::singleton();
    InputOutput::Variables *vars = C->variables();

	for(i = 0; i < _var_names.size(); i++){
        if(!vars->get(_var_names.at(i))){
            sprintf(msg,
                    "The tool \"%s\" requires the undeclared variable \"%s\".\n",
                    name(),
                    _var_names.at(i));
            S->addMessageF(3, msg);
            return true;
        }
        InputOutput::Variable *var = vars->get(_var_names.at(i));
        if((_var_sizes.at(i) == var->typesize()) &&
           (_var_values.at(i) == var->get())){
            // The variable still being valid
            continue;
        }
        // Update the variable
        InputOutput::ScreenManager *S = InputOutput::ScreenManager::singleton();
        err_code = clSetKernelArg(_kernel, i, var->typesize(), var->get());
        if(err_code != CL_SUCCESS) {
            sprintf(msg,
                    "Failure setting the variable \"%s\" (id=%u) to the tool \"%s\".\n",
                    name(),
                    i,
                    _var_names.at(i));
            S->addMessageF(3, msg);
            S->printOpenCLError(err_code);
            return true;
        }
        _var_sizes.at(i) = var->typesize();
        _var_values.at(i) = var->get();
	}
	return false;
}

}}  // namespace
