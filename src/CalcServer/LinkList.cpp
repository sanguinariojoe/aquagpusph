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
#include <InputOutput/Logger.h>

namespace Aqua{ namespace CalcServer{

#ifndef DOXYGEN_SHOULD_SKIP_THIS
#include "CalcServer/LinkList.hcl"
#include "CalcServer/LinkList.cl"
#endif
std::string LINKLIST_INC = xxd2string(LinkList_hcl_in, LinkList_hcl_in_len);
std::string LINKLIST_SRC = xxd2string(LinkList_cl_in, LinkList_cl_in_len);


LinkList::LinkList(const std::string tool_name, const std::string input)
    : Tool(tool_name)
    , _input_name(input)
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
    std::stringstream min_pos_name;
    min_pos_name << tool_name << "->Min. Pos.";
    std::string min_pos_op = "c.x = (a.x < b.x) ? a.x : b.x;\nc.y = (a.y < b.y) ? a.y : b.y;\n#ifdef HAVE_3D\nc.z = (a.z < b.z) ? a.z : b.z;\nc.w = 0.f;\n#endif\n";
    _min_pos = new Reduction(min_pos_name.str(),
                             input,
                             "r_min",
                             min_pos_op,
                             "VEC_INFINITY");
    std::stringstream max_pos_name;
    max_pos_name << tool_name << "->Max. Pos.";
    std::string max_pos_op = "c.x = (a.x > b.x) ? a.x : b.x;\nc.y = (a.y > b.y) ? a.y : b.y;\n#ifdef HAVE_3D\nc.z = (a.z > b.z) ? a.z : b.z;\nc.w = 0.f;\n#endif\n";
    _max_pos = new Reduction(max_pos_name.str(),
                             input,
                             "r_max",
                             max_pos_op,
                             "-VEC_INFINITY");
    std::stringstream sort_name;
    sort_name << tool_name << "->Radix-Sort";
    _sort = new RadixSort(sort_name.str());
}

LinkList::~LinkList()
{
    unsigned int i;
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

void LinkList::setup()
{
    InputOutput::Variables *vars = CalcServer::singleton()->variables();

    std::ostringstream msg;
    msg << "Loading the tool \"" << name() << "\"..." << std::endl;
    LOG(L_INFO, msg.str());

    // Setup the reduction tools
    _min_pos->setup();
    _max_pos->setup();

    // Compute the cells length
    InputOutput::Variable *s = vars->get("support");
    InputOutput::Variable *h = vars->get("h");
    _cell_length = *(float*)s->get() * *(float*)h->get();

    // Setup the kernels
    setupOpenCL();

    // Setup the radix-sort
    _sort->setup();
}

void LinkList::_execute()
{
    cl_int err_code;
    CalcServer *C = CalcServer::singleton();

    _min_pos->execute();
    _max_pos->execute();

    // Compute the number of cells, and allocate memory for ihoc
    nCells();
    allocate();

    // Check the validity of the variables
    setVariables();

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
        std::stringstream msg;
        msg << "Failure executing \"iCell\" from tool \"" <<
               name() << "\"." << std::endl;
        LOG(L_ERROR, msg.str());
        InputOutput::Logger::singleton()->printOpenCLError(err_code);
        throw std::runtime_error("OpenCL execution error");
    }

    // Sort the particles from the cells
    _sort->execute();

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
        std::stringstream msg;
        msg << "Failure executing \"iHoc\" from tool \"" <<
               name() << "\"." << std::endl;
        LOG(L_ERROR, msg.str());
        InputOutput::Logger::singleton()->printOpenCLError(err_code);
        throw std::runtime_error("OpenCL execution error");
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
        std::stringstream msg;
        msg << "Failure executing \"linkList\" from tool \"" <<
               name() << "\"." << std::endl;
        LOG(L_ERROR, msg.str());
        InputOutput::Logger::singleton()->printOpenCLError(err_code);
        throw std::runtime_error("OpenCL execution error");
    }
}

void LinkList::setupOpenCL()
{
    unsigned int i;
    uivec4 n_cells;
    unsigned int n_radix, N;
    cl_int err_code;
    CalcServer *C = CalcServer::singleton();
    InputOutput::Variables *vars = C->variables();

    // Create a header for the source code where the operation will be placed
    std::ostringstream source;
    source << LINKLIST_INC << LINKLIST_SRC;
    compile(source.str());

    err_code = clGetKernelWorkGroupInfo(_ihoc,
                                        C->device(),
                                        CL_KERNEL_WORK_GROUP_SIZE,
                                        sizeof(size_t),
                                        &_ihoc_lws,
                                        NULL);
    if(err_code != CL_SUCCESS) {
        LOG(L_ERROR, "Failure querying the work group size (\"iHoc\").\n");
        InputOutput::Logger::singleton()->printOpenCLError(err_code);
        throw std::runtime_error("OpenCL error");
    }
    if(_ihoc_lws < __CL_MIN_LOCALSIZE__){
        LOG(L_ERROR, "insufficient local memory for \"iHoc\".\n");
        std::stringstream msg;
        msg << "\t" << _ihoc_lws
            << " local work group size with __CL_MIN_LOCALSIZE__="
            << __CL_MIN_LOCALSIZE__ << std::endl;
        LOG0(L_DEBUG, msg.str());
        throw std::runtime_error("OpenCL error");
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
            std::stringstream msg;
            msg << "Failure sending \"" << _ihoc_vars[i]
                << "\" argument to \"iHoc\"." << std::endl;
            LOG(L_ERROR, msg.str());
            InputOutput::Logger::singleton()->printOpenCLError(err_code);
            throw std::runtime_error("OpenCL error");
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
        LOG(L_ERROR, "Failure querying the work group size (\"iCell\").\n");
        InputOutput::Logger::singleton()->printOpenCLError(err_code);
        throw std::runtime_error("OpenCL error");
    }
    if(_icell_lws < __CL_MIN_LOCALSIZE__){
        LOG(L_ERROR, "insufficient local memory for \"iCell\".\n");
        std::stringstream msg;
        msg << "\t" << _icell_lws
            << " local work group size with __CL_MIN_LOCALSIZE__="
            << __CL_MIN_LOCALSIZE__ << std::endl;
        LOG0(L_DEBUG, msg.str());
        throw std::runtime_error("OpenCL error");
    }
    n_radix = *(unsigned int*)vars->get("n_radix")->get();
    _icell_gws = roundUp(n_radix, _icell_lws);
    const char *_icell_vars[8] = {"icell", _input_name.c_str(), "N", "n_radix",
                                  "r_min", "support", "h", "n_cells"};
    for(i = 0; i < 8; i++){
        err_code = clSetKernelArg(_icell,
                                  i,
                                  vars->get(_icell_vars[i])->typesize(),
                                  vars->get(_icell_vars[i])->get());
        if(err_code != CL_SUCCESS){
            std::stringstream msg;
            msg << "Failure sending \"" << _icell_vars[i]
                << "\" argument to \"iCell\"." << std::endl;
            LOG(L_ERROR, msg.str());
            InputOutput::Logger::singleton()->printOpenCLError(err_code);
            throw std::runtime_error("OpenCL error");
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
        LOG(L_ERROR, "Failure querying the work group size (\"linkList\").\n");
        InputOutput::Logger::singleton()->printOpenCLError(err_code);
        throw std::runtime_error("OpenCL error");
    }
    if(_ll_lws < __CL_MIN_LOCALSIZE__){
        LOG(L_ERROR, "insufficient local memory for \"linkList\".\n");
        std::stringstream msg;
        msg << "\t" << _ll_lws
            << " local work group size with __CL_MIN_LOCALSIZE__="
            << __CL_MIN_LOCALSIZE__ << std::endl;
        LOG0(L_DEBUG, msg.str());
        throw std::runtime_error("OpenCL error");
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
            std::stringstream msg;
            msg << "Failure sending \"" << _ll_vars[i]
                << "\" argument to \"iCell\"." << std::endl;
            LOG(L_ERROR, msg.str());
            InputOutput::Logger::singleton()->printOpenCLError(err_code);
            throw std::runtime_error("OpenCL error");
        }
        _ll_args.push_back(malloc(vars->get(_ll_vars[i])->typesize()));
        memcpy(_ll_args.at(i),
               vars->get(_ll_vars[i])->get(),
               vars->get(_ll_vars[i])->typesize());
    }
}

void LinkList::compile(const std::string source)
{
    cl_int err_code;
    cl_program program;
    CalcServer *C = CalcServer::singleton();

    std::ostringstream flags;
    #ifdef AQUA_DEBUG
        flags << " -DDEBUG ";
    #else
        flags << " -DNDEBUG ";
    #endif
    flags << " -cl-mad-enable -cl-fast-relaxed-math";
    #ifdef HAVE_3D
        flags << " -DHAVE_3D ";
    #else
        flags << " -DHAVE_2D ";
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
        LOG(L_ERROR, "Error compiling the source code\n");
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
    _ihoc = clCreateKernel(program, "iHoc", &err_code);
    if(err_code != CL_SUCCESS) {
        LOG(L_ERROR, "Failure creating the \"iHoc\" kernel.\n");
        InputOutput::Logger::singleton()->printOpenCLError(err_code);
        clReleaseProgram(program);
        throw std::runtime_error("OpenCL error");
    }
    _icell = clCreateKernel(program, "iCell", &err_code);
    if(err_code != CL_SUCCESS) {
        LOG(L_ERROR, "Failure creating the \"iCell\" kernel.\n");
        InputOutput::Logger::singleton()->printOpenCLError(err_code);
        clReleaseProgram(program);
        throw std::runtime_error("OpenCL error");
    }
    _ll = clCreateKernel(program, "linkList", &err_code);
    if(err_code != CL_SUCCESS) {
        LOG(L_ERROR, "Failure creating the \"linkList\" kernel.\n");
        InputOutput::Logger::singleton()->printOpenCLError(err_code);
        clReleaseProgram(program);
        throw std::runtime_error("OpenCL error");
    }

    clReleaseProgram(program);
}

void LinkList::nCells()
{
    vec pos_min, pos_max;
    InputOutput::Variables *vars = CalcServer::singleton()->variables();

    if(!_cell_length){
        std::stringstream msg;
        msg << "Zero cell length detected in the tool \"" << name()
            << "\"." << std::endl;
        LOG(L_ERROR, msg.str());
        throw std::runtime_error("Invalid number of cells");
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
}

void LinkList::allocate()
{
    uivec4 n_cells;
    cl_int err_code;
    CalcServer *C = CalcServer::singleton();
    InputOutput::Variables *vars = C->variables();

    if(vars->get("n_cells")->type().compare("uivec4")){
        std::stringstream msg;
        msg << "\"n_cells\" has and invalid type for \"" << name()
            << "\"." << std::endl;
        LOG(L_ERROR, msg.str());
        msg.str("");
        msg << "\tVariable \"n_cells\" type is \"" << vars->get("n_cells")->type()
            << "\", while \"uivec4\" was expected" << std::endl;
        LOG0(L_DEBUG, msg.str());
        throw std::runtime_error("Invalid n_cells type");
    }

    n_cells = *(uivec4*)vars->get("n_cells")->get();

    if(_n_cells.w <= n_cells.w){
        n_cells.x = _n_cells.x;
        n_cells.y = _n_cells.y;
        n_cells.z = _n_cells.z;
        vars->get("n_cells")->set(&n_cells);
        return;
    }

    cl_mem mem = *(cl_mem*)vars->get("ihoc")->get();
    if(mem) clReleaseMemObject(mem); mem = NULL;

    mem = clCreateBuffer(C->context(),
                         CL_MEM_READ_WRITE,
                         _n_cells.w * sizeof(unsigned int),
                         NULL,
                         &err_code);
    if(err_code != CL_SUCCESS){
        std::stringstream msg;
        msg << "Failure allocating device memory in the tool \"" <<
               name() << "\"." << std::endl;
        LOG(L_ERROR, msg.str());
        InputOutput::Logger::singleton()->printOpenCLError(err_code);
        throw std::runtime_error("OpenCL allocation error");
    }

    n_cells = _n_cells;
    vars->get("n_cells")->set(&n_cells);
    vars->get("ihoc")->set(&mem);
    _ihoc_gws = roundUp(n_cells.w, _ihoc_lws);
}

void LinkList::setVariables()
{
    unsigned int i;
    cl_int err_code;
    InputOutput::Variables *vars = CalcServer::singleton()->variables();

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
            std::stringstream msg;
            msg << "Failure setting the variable \"" << _ihoc_vars[i]
                << "\" to the tool \"" << name()
                << "\" (\"iHoc\")." << std::endl;
            LOG(L_ERROR, msg.str());
            InputOutput::Logger::singleton()->printOpenCLError(err_code);
            throw std::runtime_error("OpenCL error");
        }
        memcpy(_ihoc_args.at(i), var->get(), var->typesize());
    }

    const char *_icell_vars[8] = {"icell", _input_name.c_str(), "N", "n_radix",
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
            std::stringstream msg;
            msg << "Failure setting the variable \"" << _ihoc_vars[i]
                << "\" to the tool \"" << name()
                << "\" (\"iCell\")." << std::endl;
            LOG(L_ERROR, msg.str());
            InputOutput::Logger::singleton()->printOpenCLError(err_code);
            throw std::runtime_error("OpenCL error");
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
            std::stringstream msg;
            msg << "Failure setting the variable \"" << _ihoc_vars[i]
                << "\" to the tool \"" << name()
                << "\" (\"linkList\")." << std::endl;
            LOG(L_ERROR, msg.str());
            InputOutput::Logger::singleton()->printOpenCLError(err_code);
            throw std::runtime_error("OpenCL error");
        }
        memcpy(_ll_args.at(i), var->get(), var->typesize());
    }
}

}}  // namespace
