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
 * @brief Methods to perform a radix sort using the GPU (or any device
 * supported by OpenCL).
 * (See Aqua::CalcServer::RadixSort for details)
 * @note Hardcoded versions of the files CalcServer/RadixSort.cl.in and
 * CalcServer/RadixSort.hcl.in are internally included as a text array.
 */

#include <AuxiliarMethods.h>
#include <InputOutput/Logger.h>
#include <CalcServer/RadixSort.h>

namespace Aqua{ namespace CalcServer{

#ifndef DOXYGEN_SHOULD_SKIP_THIS
#include "CalcServer/RadixSort.hcl"
#include "CalcServer/RadixSort.cl"
#endif
std::string RADIXSORT_INC = xxd2string(RadixSort_hcl_in, RadixSort_hcl_in_len);
std::string RADIXSORT_SRC = xxd2string(RadixSort_cl_in, RadixSort_cl_in_len);

RadixSort::RadixSort(const std::string tool_name,
                     const std::string variable,
                     const std::string permutations,
                     const std::string inv_permutations,
                     bool once)
    : Tool(tool_name, once)
    , _var_name(variable)
    , _perms_name(permutations)
    , _inv_perms_name(inv_permutations)
    , _var(NULL)
    , _perms(NULL)
    , _inv_perms(NULL)
    , _n(0)
    , _init_kernel(NULL)
    , _histograms_kernel(NULL)
    , _scan_kernel(NULL)
    , _paste_kernel(NULL)
    , _sort_kernel(NULL)
    , _inv_perms_kernel(NULL)
    , _in_keys(NULL)
    , _out_keys(NULL)
    , _in_permut(NULL)
    , _out_permut(NULL)
    , _histograms(NULL)
    , _global_sums(NULL)
    , _temp_mem(NULL)
    , _items(_ITEMS)
    , _groups(_GROUPS)
    , _bits(_STEPBITS)
    , _radix(_RADIX)
    , _histo_split(_HISTOSPLIT)
{
}

RadixSort::~RadixSort()
{
    if(_init_kernel) clReleaseKernel(_init_kernel); _init_kernel=NULL;
    if(_histograms_kernel) clReleaseKernel(_histograms_kernel); _histograms_kernel=NULL;
    if(_scan_kernel) clReleaseKernel(_scan_kernel); _scan_kernel=NULL;
    if(_paste_kernel) clReleaseKernel(_paste_kernel); _paste_kernel=NULL;
    if(_sort_kernel) clReleaseKernel(_sort_kernel); _sort_kernel=NULL;
    if(_inv_perms_kernel) clReleaseKernel(_inv_perms_kernel); _inv_perms_kernel=NULL;
    if(_in_keys) clReleaseMemObject(_in_keys); _in_keys=NULL;
    if(_out_keys) clReleaseMemObject(_out_keys); _out_keys=NULL;
    if(_in_permut) clReleaseMemObject(_in_permut); _in_permut=NULL;
    if(_out_permut) clReleaseMemObject(_out_permut); _out_permut=NULL;
    if(_histograms) clReleaseMemObject(_histograms); _histograms=NULL;
    if(_global_sums) clReleaseMemObject(_global_sums); _global_sums=NULL;
    if(_temp_mem) clReleaseMemObject(_temp_mem); _temp_mem=NULL;
}

void RadixSort::setup()
{
    std::ostringstream msg;
    msg << "Loading the tool \"" << name() << "\"..." << std::endl;
    LOG(L_INFO, msg.str());

    // Get the variables
    variables();

    // Setup the working tools
    setupOpenCL();
}

void RadixSort::_execute()
{
    cl_int err_code;
    unsigned int i, max_val;
    CalcServer *C = CalcServer::singleton();
    InputOutput::Variables *vars = C->variables();

    // Get maximum key bits, and needed pass
    max_val = UINT_MAX;
    if(!_var_name.compare("icell")){
        uivec4 n_cells = *(uivec4 *)vars->get("n_cells")->get();
        max_val = nextPowerOf2(n_cells.w);
    }
    else if(!isPowerOf2(max_val)){
        max_val = nextPowerOf2(max_val / 2);
    }
    for(i=0; (max_val&1) == 0; max_val >>= 1, i++);
    _key_bits = i;
    _key_bits = roundUp(_key_bits, _bits);
    if(_key_bits > __UINTBITS__){
        LOG(L_ERROR, "Resultant keys overflows unsigned int type.\n");
        throw std::runtime_error("Unsigned int overflow");
    }
    _n_pass = _key_bits / _bits;

    err_code = clEnqueueCopyBuffer(C->command_queue(),
                                   *(cl_mem *)_var->get(),
                                   _in_keys,
                                   0,
                                   0,
                                   _n * sizeof(cl_uint),
                                   0,
                                   NULL,
                                   NULL);
    if(err_code != CL_SUCCESS){
        std::ostringstream msg;
        msg << "Failure copying the keys to sort within the tool \"" << name()
            << "\"." << std::endl;
        LOG(L_ERROR, msg.str());
        InputOutput::Logger::singleton()->printOpenCLError(err_code);
        throw std::runtime_error("OpenCL error");
    }

    init();

    for(_pass = 0; _pass < _n_pass; _pass++){
        histograms();
        scan();
        reorder();
    }

    err_code = clEnqueueCopyBuffer(C->command_queue(),
                                   _in_keys,
                                   *(cl_mem *)_var->get(),
                                   0,
                                   0,
                                   _n * sizeof(cl_uint),
                                   0,
                                   NULL,
                                   NULL);
    if(err_code != CL_SUCCESS){
        std::ostringstream msg;
        msg << "Failure copying the sort keys within the tool \"" << name()
            << "\"." << std::endl;
        LOG(L_ERROR, msg.str());
        InputOutput::Logger::singleton()->printOpenCLError(err_code);
        throw std::runtime_error("OpenCL error");
    }
    err_code = clEnqueueCopyBuffer(C->command_queue(),
                                   _in_permut,
                                   *(cl_mem *)_perms->get(),
                                   0,
                                   0,
                                   _n * sizeof(cl_uint),
                                   0,
                                   NULL,
                                   NULL);
    if(err_code != CL_SUCCESS){
        std::ostringstream msg;
        msg << "Failure copying the permutations within the tool \"" << name()
            << "\"." << std::endl;
        LOG(L_ERROR, msg.str());
        InputOutput::Logger::singleton()->printOpenCLError(err_code);
        throw std::runtime_error("OpenCL error");
    }

    inversePermutations();
}

void RadixSort::init()
{
    cl_int err_code;
    CalcServer *C = CalcServer::singleton();

    err_code = clSetKernelArg(_init_kernel,
                              0,
                              sizeof(cl_mem),
                              (void*)&_in_permut);
    if(err_code != CL_SUCCESS){
        std::ostringstream msg;
        msg << "Failure sending argument 0 to \"init\" within the tool \""
            << name() << "\"." << std::endl;
        LOG(L_ERROR, msg.str());
        InputOutput::Logger::singleton()->printOpenCLError(err_code);
        throw std::runtime_error("OpenCL error");
    }

    err_code = clEnqueueNDRangeKernel(C->command_queue(),
                                      _init_kernel,
                                      1,
                                      NULL,
                                      &_global_work_size,
                                      NULL,
                                      0,
                                      NULL,
                                      NULL);
    if(err_code != CL_SUCCESS) {
        std::ostringstream msg;
        msg << "Failure executing \"init\" within the tool \""
            << name() << "\"." << std::endl;
        LOG(L_ERROR, msg.str());
        InputOutput::Logger::singleton()->printOpenCLError(err_code);
        throw std::runtime_error("OpenCL execution error");
    }
}

void RadixSort::histograms()
{
    cl_int err_code;
    CalcServer *C = CalcServer::singleton();
    size_t local_work_size = _items;
    size_t global_work_size = _groups * _items;

    err_code = clSetKernelArg(_histograms_kernel,
                              0,
                              sizeof(cl_mem),
                              (void*)&_in_keys);
    if(err_code != CL_SUCCESS){
        std::ostringstream msg;
        msg << "Failure sending argument 0 to \"histogram\" within the tool \""
            << name() << "\"." << std::endl;
        LOG(L_ERROR, msg.str());
        InputOutput::Logger::singleton()->printOpenCLError(err_code);
        throw std::runtime_error("OpenCL error");
    }
    err_code = clSetKernelArg(_histograms_kernel,
                              2,
                              sizeof(cl_uint),
                              (void*)&_pass);
    if(err_code != CL_SUCCESS){
        std::ostringstream msg;
        msg << "Failure sending argument 2 to \"histogram\" within the tool \""
            << name() << "\"." << std::endl;
        LOG(L_ERROR, msg.str());
        InputOutput::Logger::singleton()->printOpenCLError(err_code);
        throw std::runtime_error("OpenCL error");
    }

    err_code = clEnqueueNDRangeKernel(C->command_queue(),
                                      _histograms_kernel,
                                      1,
                                      NULL,
                                      &global_work_size,
                                      &local_work_size,
                                      0,
                                      NULL,
                                      NULL);
    if(err_code != CL_SUCCESS) {
        std::ostringstream msg;
        msg << "Failure executing \"histogram\" within the tool \""
            << name() << "\"." << std::endl;
        LOG(L_ERROR, msg.str());
        InputOutput::Logger::singleton()->printOpenCLError(err_code);
        throw std::runtime_error("OpenCL execution error");
    }
}

void RadixSort::scan()
{
    cl_int err_code;
    CalcServer *C = CalcServer::singleton();
    size_t global_work_size = _radix * _groups * _items / 2;
    size_t local_work_size = global_work_size / _histo_split;
    unsigned int maxmemcache=max(_histo_split,
                                 _items * _groups * _radix / _histo_split);

    // 1st scan
    // ========
    err_code = clSetKernelArg(_scan_kernel,
                              0,
                              sizeof(cl_mem),
                              (void*)&_histograms);
    if(err_code != CL_SUCCESS){
        std::ostringstream msg;
        msg << "Failure sending argument 0 to \"scan\" (1st call) "
            << "within the tool \"" << name() << "\"." << std::endl;
        LOG(L_ERROR, msg.str());
        InputOutput::Logger::singleton()->printOpenCLError(err_code);
        throw std::runtime_error("OpenCL error");
    }
    err_code = clSetKernelArg(_scan_kernel,
                              2,
                              sizeof(cl_mem),
                              (void*)&_global_sums);
    if(err_code != CL_SUCCESS){
        std::ostringstream msg;
        msg << "Failure sending argument 2 to \"scan\" (1st call) "
            << "within the tool \"" << name() << "\"." << std::endl;
        LOG(L_ERROR, msg.str());
        InputOutput::Logger::singleton()->printOpenCLError(err_code);
        throw std::runtime_error("OpenCL error");
    }

    err_code = clEnqueueNDRangeKernel(C->command_queue(),
                                      _scan_kernel,
                                      1,
                                      NULL,
                                      &global_work_size,
                                      &local_work_size,
                                      0,
                                      NULL,
                                      NULL);
    if(err_code != CL_SUCCESS) {
        std::ostringstream msg;
        msg << "Failure executing \"scan\" (1st call) within the tool \""
            << name() << "\"." << std::endl;
        LOG(L_ERROR, msg.str());
        InputOutput::Logger::singleton()->printOpenCLError(err_code);
        throw std::runtime_error("OpenCL execution error");
    }

    // 2nd scan
    // ========
    global_work_size = _histo_split / 2;
    local_work_size = global_work_size;
    err_code = clSetKernelArg(_scan_kernel,
                              0,
                              sizeof(cl_mem),
                              (void*)&_global_sums);
    if(err_code != CL_SUCCESS){
        std::ostringstream msg;
        msg << "Failure sending argument 0 to \"scan\" (2nd call) "
            << "within the tool \"" << name() << "\"." << std::endl;
        LOG(L_ERROR, msg.str());
        InputOutput::Logger::singleton()->printOpenCLError(err_code);
        throw std::runtime_error("OpenCL error");
    }
    err_code = clSetKernelArg(_scan_kernel,
                              2,
                              sizeof(cl_mem),
                              (void*)&_temp_mem);
    if(err_code != CL_SUCCESS){
        std::ostringstream msg;
        msg << "Failure sending argument 2 to \"scan\" (2nd call) "
            << "within the tool \"" << name() << "\"." << std::endl;
        LOG(L_ERROR, msg.str());
        InputOutput::Logger::singleton()->printOpenCLError(err_code);
        throw std::runtime_error("OpenCL error");
    }

    err_code = clEnqueueNDRangeKernel(C->command_queue(),
                                      _scan_kernel,
                                      1,
                                      NULL,
                                      &global_work_size,
                                      &local_work_size,
                                      0,
                                      NULL,
                                      NULL);
    if(err_code != CL_SUCCESS) {
        std::ostringstream msg;
        msg << "Failure executing \"scan\" (2nd call) within the tool \""
            << name() << "\"." << std::endl;
        LOG(L_ERROR, msg.str());
        InputOutput::Logger::singleton()->printOpenCLError(err_code);
        throw std::runtime_error("OpenCL execution error");
    }

    // Histograms paste
    // ================
    global_work_size = _radix * _groups * _items / 2;
    local_work_size = global_work_size / _histo_split;

    err_code = clEnqueueNDRangeKernel(C->command_queue(),
                                      _paste_kernel,
                                      1,
                                      NULL,
                                      &global_work_size,
                                      &local_work_size,
                                      0,
                                      NULL,
                                      NULL);
    if(err_code != CL_SUCCESS) {
        std::ostringstream msg;
        msg << "Failure executing \"paste\" within the tool \""
            << name() << "\"." << std::endl;
        LOG(L_ERROR, msg.str());
        InputOutput::Logger::singleton()->printOpenCLError(err_code);
        throw std::runtime_error("OpenCL execution error");
    }
}

void RadixSort::reorder()
{
    cl_int err_code;
    CalcServer *C = CalcServer::singleton();
    size_t local_work_size = _items;
    size_t global_work_size = _groups * _items;

    err_code = clSetKernelArg(_sort_kernel,
                              0,
                              sizeof(cl_mem),
                              (void*)&_in_keys);
    if(err_code != CL_SUCCESS){
        std::ostringstream msg;
        msg << "Failure sending argument 0 to \"sort\" within the tool \""
            << name() << "\"." << std::endl;
        LOG(L_ERROR, msg.str());
        InputOutput::Logger::singleton()->printOpenCLError(err_code);
        throw std::runtime_error("OpenCL error");
    }
    err_code = clSetKernelArg(_sort_kernel,
                              1,
                              sizeof(cl_mem),
                              (void*)&_out_keys);
    if(err_code != CL_SUCCESS){
        std::ostringstream msg;
        msg << "Failure sending argument 1 to \"sort\" within the tool \""
            << name() << "\"." << std::endl;
        LOG(L_ERROR, msg.str());
        InputOutput::Logger::singleton()->printOpenCLError(err_code);
        throw std::runtime_error("OpenCL error");
    }
    err_code = clSetKernelArg(_sort_kernel,
                              3,
                              sizeof(cl_uint),
                              (void*)&_pass);
    if(err_code != CL_SUCCESS){
        std::ostringstream msg;
        msg << "Failure sending argument 3 to \"sort\" within the tool \""
            << name() << "\"." << std::endl;
        LOG(L_ERROR, msg.str());
        InputOutput::Logger::singleton()->printOpenCLError(err_code);
        throw std::runtime_error("OpenCL error");
    }
    err_code = clSetKernelArg(_sort_kernel,
                              4,
                              sizeof(cl_mem),
                              (void*)&_in_permut);
    if(err_code != CL_SUCCESS){
        std::ostringstream msg;
        msg << "Failure sending argument 4 to \"sort\" within the tool \""
            << name() << "\"." << std::endl;
        LOG(L_ERROR, msg.str());
        InputOutput::Logger::singleton()->printOpenCLError(err_code);
        throw std::runtime_error("OpenCL error");
    }
    err_code = clSetKernelArg(_sort_kernel,
                              5,
                              sizeof(cl_mem),
                              (void*)&_out_permut);
    if(err_code != CL_SUCCESS){
        std::ostringstream msg;
        msg << "Failure sending argument 5 to \"sort\" within the tool \""
            << name() << "\"." << std::endl;
        LOG(L_ERROR, msg.str());
        InputOutput::Logger::singleton()->printOpenCLError(err_code);
        throw std::runtime_error("OpenCL error");
    }

    err_code = clEnqueueNDRangeKernel(C->command_queue(),
                                      _sort_kernel,
                                      1,
                                      NULL,
                                      &global_work_size,
                                      &local_work_size,
                                      0,
                                      NULL,
                                      NULL);
    if(err_code != CL_SUCCESS) {
        std::ostringstream msg;
        msg << "Failure executing \"sort\" within the tool \""
            << name() << "\"." << std::endl;
        LOG(L_ERROR, msg.str());
        InputOutput::Logger::singleton()->printOpenCLError(err_code);
        throw std::runtime_error("OpenCL execution error");
    }

    // Swap the memory identifiers for the next pass
    cl_mem d_temp;

    d_temp = _in_keys;
    _in_keys = _out_keys;
    _out_keys = d_temp;

    d_temp = _in_permut;
    _in_permut = _out_permut;
    _out_permut = d_temp;
}

void RadixSort::inversePermutations()
{
    cl_int err_code;
    CalcServer *C = CalcServer::singleton();

    err_code = clSetKernelArg(_inv_perms_kernel,
                              0,
                              sizeof(cl_mem),
                              (void*)&_in_permut);
    if(err_code != CL_SUCCESS){
        std::ostringstream msg;
        msg << "Failure sending argument 0 to \"inversePermutation\" within the tool \""
            << name() << "\"." << std::endl;
        LOG(L_ERROR, msg.str());
        InputOutput::Logger::singleton()->printOpenCLError(err_code);
        throw std::runtime_error("OpenCL error");
    }
    err_code = clSetKernelArg(_inv_perms_kernel,
                              1,
                              sizeof(cl_mem),
                              _inv_perms->get());
    if(err_code != CL_SUCCESS){
        std::ostringstream msg;
        msg << "Failure sending argument 1 to \"inversePermutation\" within the tool \""
            << name() << "\"." << std::endl;
        LOG(L_ERROR, msg.str());
        InputOutput::Logger::singleton()->printOpenCLError(err_code);
        throw std::runtime_error("OpenCL error");
    }

    err_code = clEnqueueNDRangeKernel(C->command_queue(),
                                      _inv_perms_kernel,
                                      1,
                                      NULL,
                                      &_global_work_size,
                                      NULL,
                                      0,
                                      NULL,
                                      NULL);
    if(err_code != CL_SUCCESS) {
        std::ostringstream msg;
        msg << "Failure executing \"inversePermutation\" within the tool \""
            << name() << "\"." << std::endl;
        LOG(L_ERROR, msg.str());
        InputOutput::Logger::singleton()->printOpenCLError(err_code);
        throw std::runtime_error("OpenCL execution error");
    }
}


void RadixSort::variables()
{
    size_t n;
    CalcServer *C = CalcServer::singleton();
    InputOutput::Variables *vars = C->variables();

    // Check and get the variables
    if(!vars->get(_var_name)){
        std::ostringstream msg;
        msg << "Tool \"" << name()
            << "\" is asking for the undeclared variable \"" << _var_name
            << "\"." << std::endl;
        LOG(L_ERROR, msg.str());
        throw std::runtime_error("Invalid variable");
    }
    if(vars->get(_var_name)->type().compare("unsigned int*")){
        std::ostringstream msg;
        msg << "Tool \"" << name()
            << "\" cannot process variable \"" << _var_name
            << "\"." << std::endl;
        LOG(L_ERROR, msg.str());
        msg.str("");
        msg << "\t\"unsigned int*\" type was expected, but \""
            << vars->get(_var_name)->type() << "\" has been received." << std::endl;
        LOG(L_DEBUG, msg.str());
        throw std::runtime_error("Invalid variable type");
    }
    _var = (InputOutput::ArrayVariable *)vars->get(_var_name);

    if(!vars->get(_perms_name)){
        std::ostringstream msg;
        msg << "Tool \"" << name()
            << "\" is asking for the undeclared permutations variable \""
            << _perms_name << "\"." << std::endl;
        LOG(L_ERROR, msg.str());
        throw std::runtime_error("Invalid variable");
    }
    if(vars->get(_perms_name)->type().compare("unsigned int*")){
        std::ostringstream msg;
        msg << "Tool \"" << name()
            << "\" cannot process permutations variable \"" << _perms_name
            << "\"." << std::endl;
        LOG(L_ERROR, msg.str());
        msg.str("");
        msg << "\t\"unsigned int*\" type was expected, but \""
            << vars->get(_perms_name)->type() << "\" has been received." << std::endl;
        LOG(L_DEBUG, msg.str());
        throw std::runtime_error("Invalid variable type");
    }
    _perms = (InputOutput::ArrayVariable *)vars->get(_perms_name);

    if(!vars->get(_inv_perms_name)){
        std::ostringstream msg;
        msg << "Tool \"" << name()
            << "\" is asking for the undeclared inverse permutations variable \""
            << _inv_perms_name << "\"." << std::endl;
        LOG(L_ERROR, msg.str());
        throw std::runtime_error("Invalid variable");
    }
    if(vars->get(_inv_perms_name)->type().compare("unsigned int*")){
        std::ostringstream msg;
        msg << "Tool \"" << name()
            << "\" cannot process inverse permutations variable \"" << _inv_perms_name
            << "\"." << std::endl;
        LOG(L_ERROR, msg.str());
        msg.str("");
        msg << "\t\"unsigned int*\" type was expected, but \""
            << vars->get(_inv_perms_name)->type() << "\" has been received." << std::endl;
        LOG(L_DEBUG, msg.str());
        throw std::runtime_error("Invalid variable type");
    }
    _inv_perms = (InputOutput::ArrayVariable *)vars->get(_inv_perms_name);

    // Check the lengths
    n = _var->size() / vars->typeToBytes(_var->type());
    if(!isPowerOf2(n)){
        std::ostringstream msg;
        msg << "Tool \"" << name()
            << "\" cannot process variable \"" << _var_name
            << "\"." << std::endl;
        LOG(L_ERROR, msg.str());
        msg.str("");
        msg << "\tThe variable has length, n=" << n
            << ", which is not power of 2." << std::endl;
        LOG(L_DEBUG, msg.str());
        throw std::runtime_error("Invalid variable length");
    }
    _n = n;
    n = _perms->size() / vars->typeToBytes(_perms->type());
    if(n != _n){
        std::ostringstream msg;
        msg << "Lengths mismatch in tool \"" << name()
            << "\"." << std::endl;
        LOG(L_ERROR, msg.str());
        msg.str("");
        msg << "\tVariable \"" << _var->name()
            << "\" has length, n=" << _n << std::endl;
        LOG(L_DEBUG, msg.str());
        msg.str("");
        msg << "\tVariable \"" << _perms->name()
            << "\" has length, n=" << n << std::endl;
        LOG(L_DEBUG, msg.str());
        throw std::runtime_error("Invalid variable length");
    }
    n = _inv_perms->size() / vars->typeToBytes(_inv_perms->type());
    if(n != _n){
        std::ostringstream msg;
        msg << "Lengths mismatch in tool \"" << name()
            << "\"." << std::endl;
        LOG(L_ERROR, msg.str());
        msg.str("");
        msg << "\tVariable \"" << _var->name()
            << "\" has length, n=" << _n << std::endl;
        LOG(L_DEBUG, msg.str());
        msg.str("");
        msg << "\tVariable \"" << _inv_perms->name()
            << "\" has length, n=" << n << std::endl;
        LOG(L_DEBUG, msg.str());
        throw std::runtime_error("Invalid variable length");
    }
}

void RadixSort::setupOpenCL()
{
    std::ostringstream source;
    source << RADIXSORT_INC << RADIXSORT_SRC;
    compile(source.str());

    // Check and correct _items, _groups and _histo_split
    setupDims();

    // Setup the memory objects
    setupMems();
    setupArgs();

    std::ostringstream msg;
    msg << "\titems: " << _items << std::endl;
    LOG0(L_DEBUG, msg.str());
    msg.str(""); msg << "\tgroups: " << _groups << std::endl;
    LOG0(L_DEBUG, msg.str());
    msg.str(""); msg << "\tsplits: " << _histo_split << std::endl;
    LOG0(L_DEBUG, msg.str());
}

void RadixSort::compile(const std::string source)
{
    cl_int err_code;
    cl_program program;
    CalcServer *C = CalcServer::singleton();

    std::ostringstream flags;
    flags << "-D_BITS=" << _bits << " -D_RADIX=" << _radix << " -DPERMUT";
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

    _init_kernel = clCreateKernel(program, "init", &err_code);
    if(err_code != CL_SUCCESS) {
        LOG(L_ERROR, "Failure creating the \"init\" kernel.\n");
        InputOutput::Logger::singleton()->printOpenCLError(err_code);
        clReleaseProgram(program);
        throw std::runtime_error("OpenCL error");
    }
    _histograms_kernel = clCreateKernel(program, "histogram", &err_code);
    if(err_code != CL_SUCCESS) {
        LOG(L_ERROR, "Failure creating the \"histogram\" kernel.\n");
        InputOutput::Logger::singleton()->printOpenCLError(err_code);
        clReleaseProgram(program);
        throw std::runtime_error("OpenCL error");
    }
    _scan_kernel = clCreateKernel(program, "scan", &err_code);
    if(err_code != CL_SUCCESS) {
        LOG(L_ERROR, "Failure creating the \"scan\" kernel.\n");
        InputOutput::Logger::singleton()->printOpenCLError(err_code);
        clReleaseProgram(program);
        throw std::runtime_error("OpenCL error");
    }
    _paste_kernel = clCreateKernel(program, "paste", &err_code);
    if(err_code != CL_SUCCESS) {
        LOG(L_ERROR, "Failure creating the \"paste\" kernel.\n");
        InputOutput::Logger::singleton()->printOpenCLError(err_code);
        clReleaseProgram(program);
        throw std::runtime_error("OpenCL error");
    }
    _sort_kernel = clCreateKernel(program, "sort", &err_code);
    if(err_code != CL_SUCCESS) {
        LOG(L_ERROR, "Failure creating the \"sort\" kernel.\n");
        InputOutput::Logger::singleton()->printOpenCLError(err_code);
        clReleaseProgram(program);
        throw std::runtime_error("OpenCL error");
    }
    _inv_perms_kernel = clCreateKernel(program, "inversePermutation", &err_code);
    if(err_code != CL_SUCCESS) {
        LOG(L_ERROR, "Failure creating the \"inversePermutation\" kernel.\n");
        InputOutput::Logger::singleton()->printOpenCLError(err_code);
        clReleaseProgram(program);
        throw std::runtime_error("OpenCL error");
    }
    clReleaseProgram(program);
}

void RadixSort::setupDims()
{
    cl_int err_code;
    size_t max_local_work_size, sort_local_work_size, scan_local_work_size;
    CalcServer *C = CalcServer::singleton();

    // For the _histograms_kernel and _sort_kernel _items can be used as the
    // upper bound
    err_code = clGetKernelWorkGroupInfo(_histograms_kernel,
                                        C->device(),
                                        CL_KERNEL_WORK_GROUP_SIZE,
                                        sizeof(size_t),
                                        &max_local_work_size,
                                        NULL);
    if(err_code != CL_SUCCESS) {
        LOG(L_ERROR, "Failure getting CL_KERNEL_WORK_GROUP_SIZE from \"histogram\".\n");
        InputOutput::Logger::singleton()->printOpenCLError(err_code);
        throw std::runtime_error("OpenCL error");
    }
    err_code = clGetKernelWorkGroupInfo(_sort_kernel,
                                        C->device(),
                                        CL_KERNEL_WORK_GROUP_SIZE,
                                        sizeof(size_t),
                                        &sort_local_work_size,
                                        NULL);
    if(err_code != CL_SUCCESS) {
        LOG(L_ERROR, "Failure getting CL_KERNEL_WORK_GROUP_SIZE from \"sort\".\n");
        InputOutput::Logger::singleton()->printOpenCLError(err_code);
        throw std::runtime_error("OpenCL error");
    }
    if(sort_local_work_size < max_local_work_size)
        max_local_work_size = sort_local_work_size;
    if(max_local_work_size < _items)
        _items = max_local_work_size;
    if(!isPowerOf2(_items))
        _items = nextPowerOf2(_items) / 2;

    // The _scan_kernel can be used to set an upper bound to the number of
    // histogram splits
    err_code = clGetKernelWorkGroupInfo(_scan_kernel,
                                        C->device(),
                                        CL_KERNEL_WORK_GROUP_SIZE,
                                        sizeof(size_t),
                                        &scan_local_work_size,
                                        NULL);
    if(err_code != CL_SUCCESS) {
        LOG(L_ERROR, "Failure getting CL_KERNEL_WORK_GROUP_SIZE from \"scan\".\n");
        InputOutput::Logger::singleton()->printOpenCLError(err_code);
        throw std::runtime_error("OpenCL error");
    }
    if(scan_local_work_size < _histo_split / 2)
        _histo_split = 2 * scan_local_work_size;
    if(!isPowerOf2(_histo_split))
        _histo_split = nextPowerOf2(_histo_split) / 2;

    // Finally using _scan_kernel and _paste_kernel we can adjust _groups,
    // _items and _radix
    err_code = clGetKernelWorkGroupInfo(_paste_kernel,
                                        C->device(),
                                        CL_KERNEL_WORK_GROUP_SIZE,
                                        sizeof(size_t),
                                        &max_local_work_size,
                                        NULL);
    if(err_code != CL_SUCCESS) {
        LOG(L_ERROR, "Failure getting CL_KERNEL_WORK_GROUP_SIZE from \"paste\".\n");
        InputOutput::Logger::singleton()->printOpenCLError(err_code);
        throw std::runtime_error("OpenCL error");
    }
    max_local_work_size = min(max_local_work_size, scan_local_work_size);
    while(max_local_work_size < _radix * _groups * _items / 2 / _histo_split){
        // We can't increase _histo_split, so we may start decreasing the number
        // of items
        _items /= 2;
        if(_items < __CL_MIN_LOCALSIZE__){
            _items = __CL_MIN_LOCALSIZE__;
            break;
        }
    }
    while(max_local_work_size < _radix * _groups * _items / 2 / _histo_split){
        // We have reached the minimum possible number of items, so we can
        // continue decreasing the number of groups
        _groups /= 2;
        if(!_groups){
            _groups = 1;
            break;
        }
    }
    if(max_local_work_size < _radix * _groups * _items / 2 / _histo_split){
        // We can try to reduce the radix, but it is a bad business
        LOG(L_ERROR, "Failure setting a number of items and groups compatible with \"scan\" and \"paste\".\n");
        LOG0(L_DEBUG, "\tYou can try to recompile the code decreasing __CL_MIN_LOCALSIZE__\n");
        throw std::runtime_error("OpenCL error");
    }

    _local_work_size = getLocalWorkSize(C->command_queue());
    _global_work_size = getGlobalWorkSize(_n, _local_work_size);
}

void RadixSort::setupMems()
{
    cl_int err_code;
    CalcServer *C = CalcServer::singleton();

    if(_in_keys) clReleaseMemObject(_in_keys); _in_keys=NULL;
    if(_out_keys) clReleaseMemObject(_out_keys); _out_keys=NULL;
    if(_in_permut) clReleaseMemObject(_in_permut); _in_permut=NULL;
    if(_out_permut) clReleaseMemObject(_out_permut); _out_permut=NULL;
    if(_histograms) clReleaseMemObject(_histograms); _histograms=NULL;
    if(_global_sums) clReleaseMemObject(_global_sums); _global_sums=NULL;
    if(_temp_mem) clReleaseMemObject(_temp_mem); _temp_mem=NULL;
    allocatedMemory(0);

    // Get the memory identifiers
    _in_keys = clCreateBuffer(C->context(),
                              CL_MEM_READ_WRITE,
                              _n * sizeof(unsigned int),
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
    _out_keys = clCreateBuffer(C->context(),
                              CL_MEM_READ_WRITE,
                              _n * sizeof(unsigned int),
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
    _in_permut = clCreateBuffer(C->context(),
                                CL_MEM_READ_WRITE,
                                _n * sizeof(unsigned int),
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
    _out_permut = clCreateBuffer(C->context(),
                                 CL_MEM_READ_WRITE,
                                 _n * sizeof(unsigned int),
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
    _histograms = clCreateBuffer(C->context(),
                                 CL_MEM_READ_WRITE,
                                 (_radix * _groups * _items) * sizeof(unsigned int),
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
    _global_sums = clCreateBuffer(C->context(),
                                  CL_MEM_READ_WRITE,
                                  _histo_split * sizeof(unsigned int),
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
    _temp_mem = clCreateBuffer(C->context(),
                               CL_MEM_READ_WRITE,
                               sizeof(unsigned int),
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

    allocatedMemory(4 * _n * sizeof(unsigned int) +
                    (_radix * _groups * _items) * sizeof(unsigned int) +
                    _histo_split * sizeof(unsigned int) +
                    sizeof(unsigned int));
}

void RadixSort::setupArgs()
{
    cl_int err_code;
    CalcServer *C = CalcServer::singleton();

    err_code = clSetKernelArg(_init_kernel,
                              1,
                              sizeof(cl_uint),
                              (void*)&_n);
    if(err_code != CL_SUCCESS){
        std::ostringstream msg;
        msg << "Failure sending argument 1 to \"init\" within the tool \""
            << name() << "\"." << std::endl;
        LOG(L_ERROR, msg.str());
        InputOutput::Logger::singleton()->printOpenCLError(err_code);
        throw std::runtime_error("OpenCL error");
    }

    err_code = clSetKernelArg(_histograms_kernel,
                              1,
                              sizeof(cl_mem),
                              (void*)&_histograms);
    if(err_code != CL_SUCCESS){
        std::ostringstream msg;
        msg << "Failure sending argument 1 to \"histogram\" within the tool \""
            << name() << "\"." << std::endl;
        LOG(L_ERROR, msg.str());
        InputOutput::Logger::singleton()->printOpenCLError(err_code);
        throw std::runtime_error("OpenCL error");
    }
    err_code = clSetKernelArg(_histograms_kernel,
                              3,
                              sizeof(cl_uint) * _radix * _items,
                              NULL);
    if(err_code != CL_SUCCESS){
        std::ostringstream msg;
        msg << "Failure sending argument 3 to \"histogram\" within the tool \""
            << name() << "\"." << std::endl;
        LOG(L_ERROR, msg.str());
        InputOutput::Logger::singleton()->printOpenCLError(err_code);
        throw std::runtime_error("OpenCL error");
    }
    err_code = clSetKernelArg(_histograms_kernel,
                              4,
                              sizeof(cl_uint),
                              (void*)&_n);
    if(err_code != CL_SUCCESS){
        std::ostringstream msg;
        msg << "Failure sending argument 4 to \"histogram\" within the tool \""
            << name() << "\"." << std::endl;
        LOG(L_ERROR, msg.str());
        InputOutput::Logger::singleton()->printOpenCLError(err_code);
        throw std::runtime_error("OpenCL error");
    }

    unsigned int maxmemcache = max(_histo_split,
                                   _items * _groups * _radix / _histo_split);
    err_code = clSetKernelArg(_scan_kernel,
                              1,
                              sizeof(cl_uint) * maxmemcache,
                              NULL);
    if(err_code != CL_SUCCESS){
        std::ostringstream msg;
        msg << "Failure sending argument 1 to \"scan\" within the tool \""
            << name() << "\"." << std::endl;
        LOG(L_ERROR, msg.str());
        InputOutput::Logger::singleton()->printOpenCLError(err_code);
        throw std::runtime_error("OpenCL error");
    }

    err_code = clSetKernelArg(_paste_kernel,
                              0,
                              sizeof(cl_mem),
                              (void*)&_histograms);
    if(err_code != CL_SUCCESS){
        std::ostringstream msg;
        msg << "Failure sending argument 0 to \"paste\" within the tool \""
            << name() << "\"." << std::endl;
        LOG(L_ERROR, msg.str());
        InputOutput::Logger::singleton()->printOpenCLError(err_code);
        throw std::runtime_error("OpenCL error");
    }
    err_code = clSetKernelArg(_paste_kernel,
                              1,
                              sizeof(cl_mem),
                              (void*)&_global_sums);
    if(err_code != CL_SUCCESS){
        std::ostringstream msg;
        msg << "Failure sending argument 1 to \"paste\" within the tool \""
            << name() << "\"." << std::endl;
        LOG(L_ERROR, msg.str());
        InputOutput::Logger::singleton()->printOpenCLError(err_code);
        throw std::runtime_error("OpenCL error");
    }

    err_code = clSetKernelArg(_sort_kernel,
                              2,
                              sizeof(cl_mem),
                              (void*)&_histograms);
    if(err_code != CL_SUCCESS){
        std::ostringstream msg;
        msg << "Failure sending argument 2 to \"sort\" within the tool \""
            << name() << "\"." << std::endl;
        LOG(L_ERROR, msg.str());
        InputOutput::Logger::singleton()->printOpenCLError(err_code);
        throw std::runtime_error("OpenCL error");
    }
    err_code = clSetKernelArg(_sort_kernel,
                              6,
                              sizeof(cl_uint) * _radix * _items,
                              NULL);
    if(err_code != CL_SUCCESS){
        std::ostringstream msg;
        msg << "Failure sending argument 6 to \"sort\" within the tool \""
            << name() << "\"." << std::endl;
        LOG(L_ERROR, msg.str());
        InputOutput::Logger::singleton()->printOpenCLError(err_code);
        throw std::runtime_error("OpenCL error");
    }
    err_code = clSetKernelArg(_sort_kernel,
                              7,
                              sizeof(cl_uint),
                              (void*)&_n);
    if(err_code != CL_SUCCESS){
        std::ostringstream msg;
        msg << "Failure sending argument 7 to \"sort\" within the tool \""
            << name() << "\"." << std::endl;
        LOG(L_ERROR, msg.str());
        InputOutput::Logger::singleton()->printOpenCLError(err_code);
        throw std::runtime_error("OpenCL error");
    }

    err_code = clSetKernelArg(_inv_perms_kernel,
                              1,
                              sizeof(cl_mem),
                              _inv_perms->get());
    if(err_code != CL_SUCCESS){
        LOG(L_ERROR, "Failure sending argument 1 to \"inversePermutation\"\n");
        InputOutput::Logger::singleton()->printOpenCLError(err_code);
        throw std::runtime_error("OpenCL error");
    }
    err_code = clSetKernelArg(_inv_perms_kernel,
                              2,
                              sizeof(cl_uint),
                              (void*)&_n);
    if(err_code != CL_SUCCESS){
        std::ostringstream msg;
        msg << "Failure sending argument 1 to \"inversePermutation\" within the tool \""
            << name() << "\"." << std::endl;
        LOG(L_ERROR, msg.str());
        InputOutput::Logger::singleton()->printOpenCLError(err_code);
        throw std::runtime_error("OpenCL error");
    }
}

}}  // namespace
