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
 * @brief link-list based neighbours location algorithm.
 * (See Aqua::CalcServer::LinkList for details)
 * @note Hardcoded versions of the files CalcServer/LinkList.cl.in and
 * CalcServer/LinkList.hcl.in are internally included as a text array.
 */

#include <CalcServer/LinkList.h>
#include <CalcServer.h>
#include <ScreenManager.h>

namespace Aqua{ namespace CalcServer{

#ifndef DOXYGEN_SHOULD_SKIP_THIS
#include "CalcServer/LinkList.hcl"
#include "CalcServer/LinkList.cl"
#endif
const char* LINKLIST_INC = (const char*)LinkList_hcl_in;
unsigned int LINKLIST_INC_LEN = LinkList_hcl_in_len;
const char* LINKLIST_SRC = (const char*)LinkList_cl_in;
unsigned int LINKLIST_SRC_LEN = LinkList_cl_in_len;

LinkList::LinkList(const char* tool_name, const char* input)
    : Tool(tool_name)
    , _input_name(NULL)
    , _cell_length(0.f)
    , _min_pos(NULL)
    , _max_pos(NULL)
    , _ihoc(NULL)
    , _ihoc_lws(0)
    , _ihoc_gws(0)
    , _icell(NULL)
    , _icell_lws(0)
    , _icell_gws(0)
    , _ll(NULL)
    , _ll_lws(0)
    , _ll_gws(0)
{
    _input_name = new char[strlen(input) + 1];
    strcpy(_input_name, input);

    char min_pos_name[strlen(tool_name) + strlen("->Min. Pos.") + 1];
    strcpy(min_pos_name, tool_name);
    strcat(min_pos_name, "->Min. Pos.");
    const char *min_pos_op = "c.x = (a.x < b.x) ? a.x : b.x;\nc.y = (a.y < b.y) ? a.y : b.y;\n#ifdef HAVE_3D\nc.z = (a.z < b.z) ? a.z : b.z;\nc.w = 0.f;\n#endif\n";
    _min_pos = new Reduction(min_pos_name,
                             input,
                             "r_min",
                             min_pos_op,
                             "VEC_INFINITY");
    char max_pos_name[strlen(tool_name) + strlen("->Max. Pos.") + 1];
    strcpy(max_pos_name, tool_name);
    strcat(max_pos_name, "->Max. Pos.");
    const char *max_pos_op = "c.x = (a.x > b.x) ? a.x : b.x;\nc.y = (a.y > b.y) ? a.y : b.y;\n#ifdef HAVE_3D\nc.z = (a.z > b.z) ? a.z : b.z;\nc.w = 0.f;\n#endif\n";
    _max_pos = new Reduction(max_pos_name,
                             input,
                             "r_max",
                             max_pos_op,
                             "-VEC_INFINITY");
    char sort_name[strlen(tool_name) + strlen("->Radix-Sort") + 1];
    strcpy(sort_name, tool_name);
    strcat(sort_name, "->Radix-Sort");
    _sort = new RadixSort(sort_name);
}

LinkList::~LinkList()
{
    unsigned int i;
    if(_input_name) delete[] _input_name; _input_name=NULL;
    if(_min_pos) delete _min_pos; _min_pos=NULL;
    if(_max_pos) delete _max_pos; _max_pos=NULL;
    if(_sort) delete _sort; _sort=NULL;
    if(_ihoc) clReleaseKernel(_ihoc); _ihoc=NULL;
    if(_icell) clReleaseKernel(_icell); _icell=NULL;
    if(_ll) clReleaseKernel(_ll); _ll=NULL;
    for(i = 0; i < _ihoc_args.size(); i++){
        free(_ihoc_args.at(i));
    }
    _ihoc_args.clear();
    for(i = 0; i < _icell_args.size(); i++){
        free(_icell_args.at(i));
    }
    _icell_args.clear();
    for(i = 0; i < _ll_args.size(); i++){
        free(_ll_args.at(i));
    }
    _ll_args.clear();
}

bool LinkList::setup()
{
    char msg[1024];
    InputOutput::ScreenManager *S = InputOutput::ScreenManager::singleton();
    CalcServer *C = CalcServer::singleton();
    InputOutput::Variables *vars = C->variables();

    sprintf(msg,
            "Loading the tool \"%s\"...\n",
            name());
    S->addMessageF(L_INFO, msg);

    // Setup the reduction tools
    if(_min_pos->setup()){
        return true;
    }
    if(_max_pos->setup()){
        return true;
    }

    // Compute the cells length
    InputOutput::Variable *s = vars->get("support");
    InputOutput::Variable *h = vars->get("h");
    _cell_length = *(float*)s->get() * *(float*)h->get();

    // Setup the kernels
    if(setupOpenCL()){
        return true;
    }

    // Setup the radix-sort
    if(_sort->setup()){
        return true;
    }

    return false;
}

bool LinkList::_execute()
{
    cl_int err_code;
    char msg[1024];
    InputOutput::ScreenManager *S = InputOutput::ScreenManager::singleton();
    CalcServer *C = CalcServer::singleton();

    if(_min_pos->execute()){
        return true;
    }
    if(_max_pos->execute()){
        return true;
    }

    // Compute the number of cells, and allocate memory for ihoc
    if(nCells()){
        return true;
    }
    if(allocate()){
        return true;
    }

    // Check the validity of the variables
    if(setVariables()){
        return true;
    }

    // Compute the cell of each particle
    err_code = clEnqueueNDRangeKernel(C->command_queue(),
                                      _icell,
                                      1,
                                      NULL,
                                      &_icell_gws,
                                      &_icell_lws,
                                      0,
                                      NULL,
                                      NULL);
    if(err_code != CL_SUCCESS) {
        sprintf(msg,
                "Failure executing \"iCell\" from tool \"%s\".\n",
                name());
        S->addMessageF(L_ERROR, msg);
        S->printOpenCLError(err_code);
        return true;
    }

    // Sort the particles from the cells
    if(_sort->execute()){
        return true;
    }

    // Compute the head of cells
    err_code = clEnqueueNDRangeKernel(C->command_queue(),
                                      _ihoc,
                                      1,
                                      NULL,
                                      &_ihoc_gws,
                                      &_ihoc_lws,
                                      0,
                                      NULL,
                                      NULL);
    if(err_code != CL_SUCCESS) {
        sprintf(msg,
                "Failure executing \"iHoc\" from tool \"%s\".\n",
                name());
        S->addMessageF(L_ERROR, msg);
        S->printOpenCLError(err_code);
        return true;
    }

    err_code = clEnqueueNDRangeKernel(C->command_queue(),
                                      _ll,
                                      1,
                                      NULL,
                                      &_ll_gws,
                                      &_ll_lws,
                                      0,
                                      NULL,
                                      NULL);
    if(err_code != CL_SUCCESS) {
        sprintf(msg,
                "Failure executing \"linkList\" from tool \"%s\".\n",
                name());
        S->addMessageF(L_ERROR, msg);
        S->printOpenCLError(err_code);
        return true;
    }

    return false;
}

bool LinkList::setupOpenCL()
{
    unsigned int i;
    uivec4 n_cells;
    unsigned int n_radix, N;
    cl_int err_code;
    char msg[1024];
    InputOutput::ScreenManager *S = InputOutput::ScreenManager::singleton();
    CalcServer *C = CalcServer::singleton();
    InputOutput::Variables *vars = C->variables();

    // Create a header for the source code where the operation will be placed
    char header[LINKLIST_INC_LEN + 1];
    strcpy(header, "");
    strncat(header, LINKLIST_INC, LINKLIST_INC_LEN);
    strcat(header, "");

    // Setup the complete source code
    char source[strlen(header) + strlen(LINKLIST_SRC) + 1];
    strcpy(source, header);
    strncat(source, LINKLIST_SRC, LINKLIST_SRC_LEN);
    strcat(source, "");

    if(compile(source)){
        return true;
    }

    err_code = clGetKernelWorkGroupInfo(_ihoc,
                                        C->device(),
                                        CL_KERNEL_WORK_GROUP_SIZE,
                                        sizeof(size_t),
                                        &_ihoc_lws,
                                        NULL);
    if(err_code != CL_SUCCESS) {
        S->addMessageF(L_ERROR, "Failure querying the work group size (\"iHoc\").\n");
        S->printOpenCLError(err_code);
        return true;
    }
    if(_ihoc_lws < __CL_MIN_LOCALSIZE__){
        S->addMessageF(L_ERROR, "iHoc cannot be performed.\n");
        sprintf(msg,
                "\t%lu elements can be executed, but __CL_MIN_LOCALSIZE__=%lu\n",
                _ihoc_lws,
                __CL_MIN_LOCALSIZE__);
        S->addMessage(L_DEBUG, msg);
        return true;
    }
    n_cells = *(uivec4*)vars->get("n_cells")->get();
    _ihoc_gws = roundUp(n_cells.w, _ihoc_lws);
    const char *_ihoc_vars[3] = {"ihoc", "N", "n_cells"};
    for(i = 0; i < 3; i++){
        err_code = clSetKernelArg(_ihoc,
                                  i,
                                  vars->get(_ihoc_vars[i])->typesize(),
                                  vars->get(_ihoc_vars[i])->get());
        if(err_code != CL_SUCCESS){
            sprintf(msg,
                    "Failure sending \"%s\" argument to \"iHoc\".\n",
                    _ihoc_vars[i]);
            S->addMessageF(L_ERROR, msg);
            S->printOpenCLError(err_code);
            return true;
        }
        _ihoc_args.push_back(malloc(vars->get(_ihoc_vars[i])->typesize()));
        memcpy(_ihoc_args.at(i),
               vars->get(_ihoc_vars[i])->get(),
               vars->get(_ihoc_vars[i])->typesize());
    }

    err_code = clGetKernelWorkGroupInfo(_icell,
                                        C->device(),
                                        CL_KERNEL_WORK_GROUP_SIZE,
                                        sizeof(size_t),
                                        &_icell_lws,
                                        NULL);
    if(err_code != CL_SUCCESS) {
        S->addMessageF(L_ERROR, "Failure querying the work group size (\"iCell\").\n");
        S->printOpenCLError(err_code);
        return true;
    }
    if(_icell_lws < __CL_MIN_LOCALSIZE__){
        S->addMessageF(L_ERROR, "iCell cannot be performed.\n");
        sprintf(msg,
                "\t%lu elements can be executed, but __CL_MIN_LOCALSIZE__=%lu\n",
                _icell_lws,
                __CL_MIN_LOCALSIZE__);
        S->addMessage(L_DEBUG, msg);
        return true;
    }
    n_radix = *(unsigned int*)vars->get("n_radix")->get();
    _icell_gws = roundUp(n_radix, _icell_lws);
    const char *_icell_vars[8] = {"icell", _input_name, "N", "n_radix",
                                  "r_min", "support", "h", "n_cells"};
    for(i = 0; i < 8; i++){
        err_code = clSetKernelArg(_icell,
                                  i,
                                  vars->get(_icell_vars[i])->typesize(),
                                  vars->get(_icell_vars[i])->get());
        if(err_code != CL_SUCCESS){
            sprintf(msg,
                    "Failure sending \"%s\" argument to \"iCell\".\n",
                    _icell_vars[i]);
            S->addMessageF(L_ERROR, msg);
            S->printOpenCLError(err_code);
            return true;
        }
        _icell_args.push_back(malloc(vars->get(_icell_vars[i])->typesize()));
        memcpy(_icell_args.at(i),
               vars->get(_icell_vars[i])->get(),
               vars->get(_icell_vars[i])->typesize());
    }

    err_code = clGetKernelWorkGroupInfo(_ll,
                                        C->device(),
                                        CL_KERNEL_WORK_GROUP_SIZE,
                                        sizeof(size_t),
                                        &_ll_lws,
                                        NULL);
    if(err_code != CL_SUCCESS) {
        S->addMessageF(L_ERROR, "Failure querying the work group size (\"linkList\").\n");
        S->printOpenCLError(err_code);
        return true;
    }
    if(_ll_lws < __CL_MIN_LOCALSIZE__){
        S->addMessageF(L_ERROR, "linkList cannot be performed.\n");
        sprintf(msg,
                "\t%lu elements can be executed, but __CL_MIN_LOCALSIZE__=%lu\n",
                _ll_lws,
                __CL_MIN_LOCALSIZE__);
        S->addMessage(L_DEBUG, msg);
        return true;
    }
    N = *(unsigned int*)vars->get("N")->get();
    _ll_gws = roundUp(N, _ll_lws);
    const char *_ll_vars[3] = {"icell", "ihoc", "N"};
    for(i = 0; i < 3; i++){
        err_code = clSetKernelArg(_ll,
                                  i,
                                  vars->get(_ll_vars[i])->typesize(),
                                  vars->get(_ll_vars[i])->get());
        if(err_code != CL_SUCCESS){
            sprintf(msg,
                    "Failure sending \"%s\" argument to \"iCell\".\n",
                    _ll_vars[i]);
            S->addMessageF(L_ERROR, msg);
            S->printOpenCLError(err_code);
            return true;
        }
        _ll_args.push_back(malloc(vars->get(_ll_vars[i])->typesize()));
        memcpy(_ll_args.at(i),
               vars->get(_ll_vars[i])->get(),
               vars->get(_ll_vars[i])->typesize());
    }

    return false;
}

bool LinkList::compile(const char* source)
{
    cl_int err_code;
    cl_program program;
    char msg[1024];
    InputOutput::ScreenManager *S = InputOutput::ScreenManager::singleton();
    CalcServer *C = CalcServer::singleton();

    char flags[512];
    strcpy(flags, "");
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
        S->addMessageF(L_ERROR, "Failure creating the OpenCL program.\n");
        S->printOpenCLError(err_code);
        return true;
    }
    err_code = clBuildProgram(program, 0, NULL, flags, NULL, NULL);
    if(err_code != CL_SUCCESS) {
        S->addMessage(L_ERROR, "Error compiling the source code\n");
        S->printOpenCLError(err_code);
        S->addMessage(L_ERROR, "--- Build log ---------------------------------\n");
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
            S->addMessage(L_ERROR, msg);
            S->addMessage(L_ERROR, "--------------------------------- Build log ---\n");
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
        S->addMessage(L_DEBUG, log);
        S->addMessage(L_ERROR, "--------------------------------- Build log ---\n");
        free(log); log=NULL;
        clReleaseProgram(program);
        return true;
    }
    _ihoc = clCreateKernel(program, "iHoc", &err_code);
    if(err_code != CL_SUCCESS) {
        S->addMessageF(L_ERROR, "Failure creating the \"iHoc\" kernel.\n");
        S->printOpenCLError(err_code);
        clReleaseProgram(program);
        return true;
    }
    _icell = clCreateKernel(program, "iCell", &err_code);
    if(err_code != CL_SUCCESS) {
        S->addMessageF(L_ERROR, "Failure creating the \"iCell\" kernel.\n");
        S->printOpenCLError(err_code);
        clReleaseProgram(program);
        return true;
    }
    _ll = clCreateKernel(program, "linkList", &err_code);
    if(err_code != CL_SUCCESS) {
        S->addMessageF(L_ERROR, "Failure creating the \"linkList\" kernel.\n");
        S->printOpenCLError(err_code);
        clReleaseProgram(program);
        return true;
    }

    clReleaseProgram(program);
    return false;
}

bool LinkList::nCells()
{
    vec pos_min, pos_max;
    char msg[1024];
    InputOutput::ScreenManager *S = InputOutput::ScreenManager::singleton();
    CalcServer *C = CalcServer::singleton();
    InputOutput::Variables *vars = C->variables();

    if(!_cell_length){
        sprintf(msg,
                "Zero cell length detected in the tool \"%s\".\n",
                name());
        S->addMessageF(L_ERROR, msg);
        return true;
    }

    pos_min = *(vec*)vars->get("r_min")->get();
    pos_max = *(vec*)vars->get("r_max")->get();

    _n_cells.x = (unsigned int)((pos_max.x - pos_min.x) / _cell_length) + 6;
    _n_cells.y = (unsigned int)((pos_max.y - pos_min.y) / _cell_length) + 6;
    #ifdef HAVE_3D
        _n_cells.z = (unsigned int)((pos_max.z - pos_min.z) / _cell_length) + 6;
    #else
        _n_cells.z = 1;
    #endif
    _n_cells.w = _n_cells.x * _n_cells.y * _n_cells.z;

    return false;
}

bool LinkList::allocate()
{
    uivec4 n_cells;
    cl_int err_code;
    char msg[1024];
    InputOutput::ScreenManager *S = InputOutput::ScreenManager::singleton();
    CalcServer *C = CalcServer::singleton();
    InputOutput::Variables *vars = C->variables();

    if(vars->get("n_cells")->type().compare("uivec4")){
        sprintf(msg,
                "Wrong type found during the execution of the tool \"%s\".\n",
                name());
        S->addMessageF(L_ERROR, msg);
        sprintf(msg,
                "\tVariable \"%s\" type is \"%s\" (\"%s\" expected).\n",
                "n_cells",
                vars->get("n_cells")->type(),
                "uivec4");
        S->addMessage(L_DEBUG, msg);
        return true;
    }

    n_cells = *(uivec4*)vars->get("n_cells")->get();

    if(_n_cells.w <= n_cells.w){
        n_cells.x = _n_cells.x;
        n_cells.y = _n_cells.y;
        n_cells.z = _n_cells.z;
        vars->get("n_cells")->set(&n_cells);
        return false;
    }

    cl_mem mem = *(cl_mem*)vars->get("ihoc")->get();
    if(mem) clReleaseMemObject(mem); mem = NULL;

    mem = clCreateBuffer(C->context(),
                         CL_MEM_READ_WRITE,
                         _n_cells.w * sizeof(unsigned int),
                         NULL,
                         &err_code);
    if(err_code != CL_SUCCESS){
        sprintf(msg,
                "Buffer memory allocation failure during tool \"%s\" execution.\n",
                name());
        S->addMessageF(L_ERROR, msg);
        S->printOpenCLError(err_code);
        return true;
    }

    n_cells = _n_cells;
    vars->get("n_cells")->set(&n_cells);
    vars->get("ihoc")->set(&mem);
    _ihoc_gws = roundUp(n_cells.w, _ihoc_lws);

    return false;
}

bool LinkList::setVariables()
{
    unsigned int i;
    char msg[1024];
    cl_int err_code;
    InputOutput::ScreenManager *S = InputOutput::ScreenManager::singleton();
    CalcServer *C = CalcServer::singleton();
    InputOutput::Variables *vars = C->variables();

    const char *_ihoc_vars[3] = {"ihoc", "N", "n_cells"};
    for(i = 0; i < 3; i++){
        InputOutput::Variable *var = vars->get(_ihoc_vars[i]);
        if(!memcmp(var->get(), _ihoc_args.at(i), var->typesize())){
            continue;
        }
        err_code = clSetKernelArg(_ihoc,
                                  i,
                                  var->typesize(),
                                  var->get());
        if(err_code != CL_SUCCESS){
            sprintf(msg,
                    "Failure sending \"%s\" argument to \"iHoc\" in tool \"%s\".\n",
                    _ihoc_vars[i],
                    name());
            S->addMessageF(L_ERROR, msg);
            S->printOpenCLError(err_code);
            return true;
        }
        memcpy(_ihoc_args.at(i), var->get(), var->typesize());
    }

    const char *_icell_vars[8] = {"icell", _input_name, "N", "n_radix",
                                  "r_min", "support", "h", "n_cells"};
    for(i = 0; i < 8; i++){
        InputOutput::Variable *var = vars->get(_icell_vars[i]);
        if(!memcmp(var->get(), _icell_args.at(i), var->typesize())){
            continue;
        }
        err_code = clSetKernelArg(_icell,
                                  i,
                                  var->typesize(),
                                  var->get());
        if(err_code != CL_SUCCESS){
            sprintf(msg,
                    "Failure sending \"%s\" argument to \"iCell\" in tool \"%s\".\n",
                    _icell_vars[i],
                    name());
            S->addMessageF(L_ERROR, msg);
            S->printOpenCLError(err_code);
            return true;
        }
        memcpy(_icell_args.at(i), var->get(), var->typesize());
    }

    const char *_ll_vars[3] = {"icell", "ihoc", "N"};
    for(i = 0; i < 3; i++){
        InputOutput::Variable *var = vars->get(_ll_vars[i]);
        if(!memcmp(var->get(), _ll_args.at(i), var->typesize())){
            continue;
        }
        err_code = clSetKernelArg(_ll,
                                  i,
                                  var->typesize(),
                                  var->get());
        if(err_code != CL_SUCCESS){
            sprintf(msg,
                    "Failure sending \"%s\" argument to \"iCell\" in tool \"%s\".\n",
                    _ll_vars[i],
                    name());
            S->addMessageF(L_ERROR, msg);
            S->printOpenCLError(err_code);
            return true;
        }
        memcpy(_ll_args.at(i), var->get(), var->typesize());
    }

    return false;
}

}}  // namespace
