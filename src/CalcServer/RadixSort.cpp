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

#include <limits.h>
#include <CalcServer/RadixSort.h>
#include <ScreenManager.h>
#include <AuxiliarMethods.h>

namespace Aqua{ namespace CalcServer{

#ifndef DOXYGEN_SHOULD_SKIP_THIS
#include "CalcServer/RadixSort.hcl"
#include "CalcServer/RadixSort.cl"
#endif
const char* RADIXSORT_INC = (const char*)RadixSort_hcl_in;
unsigned int RADIXSORT_INC_LEN = RadixSort_hcl_in_len;
const char* RADIXSORT_SRC = (const char*)RadixSort_cl_in;
unsigned int RADIXSORT_SRC_LEN = RadixSort_cl_in_len;

RadixSort::RadixSort(const char* tool_name,
                     const char* variable,
                     const char* permutations,
                     const char* inv_permutations)
    : Tool(tool_name)
    , _var_name(NULL)
    , _perms_name(NULL)
    , _inv_perms_name(NULL)
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
    _var_name = new char[strlen(variable) + 1];
    strcpy(_var_name, variable);
    _perms_name = new char[strlen(permutations) + 1];
    strcpy(_perms_name, permutations);
    _inv_perms_name = new char[strlen(inv_permutations) + 1];
    strcpy(_inv_perms_name, inv_permutations);
}

RadixSort::~RadixSort()
{
    if(_var_name) delete[] _var_name; _var_name=NULL;
    if(_perms_name) delete[] _perms_name; _perms_name=NULL;
    if(_inv_perms_name) delete[] _inv_perms_name; _inv_perms_name=NULL;
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

bool RadixSort::setup()
{
    char msg[1024];
    InputOutput::ScreenManager *S = InputOutput::ScreenManager::singleton();
    CalcServer *C = CalcServer::singleton();
    InputOutput::Variables vars = C->variables();

    sprintf(msg,
            "Loading the tool \"%s\"...\n",
            name());
    S->addMessageF(L_INFO, msg);

    // Get the variables
    if(variables()){
        return true;
    }

    // Setup the working tools
    if(setupOpenCL()){
        return true;
    }

    return false;
}

bool RadixSort::_execute()
{
    cl_int err_code;
    unsigned int i, max_val;
    char msg[1024];
    InputOutput::ScreenManager *S = InputOutput::ScreenManager::singleton();
    CalcServer *C = CalcServer::singleton();
    InputOutput::Variables vars = C->variables();

    // Get maximum key bits, and needed pass
    max_val = UINT_MAX;
    if(!strcmp(_var_name, "icell")){
        uivec4 n_cells = *(uivec4 *)vars.get("n_cells")->get();
        max_val = nextPowerOf2(n_cells.w);
    }
    else if(!isPowerOf2(max_val)){
        max_val = nextPowerOf2(max_val / 2);
    }
    for(i=0; (max_val&1) == 0; max_val >>= 1, i++);
    _key_bits = i;
    _key_bits = roundUp(_key_bits, _bits);
    if(_key_bits > __UINTBITS__){
        S->addMessageF(L_ERROR, "Resultant keys overflows unsigned int type.\n");
        return true;
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
        sprintf(msg,
                "Failure copying the keys to sort with the tool \"%s\".\n",
                name());
        S->addMessageF(L_ERROR, msg);
        S->printOpenCLError(err_code);
        return true;
    }

    if(init())
        return true;

    for(_pass = 0; _pass < _n_pass; _pass++){
        if(histograms())
            return true;
        if(scan())
            return true;
        if(reorder())
            return true;
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
        sprintf(msg,
                "Failure copying the sorted keys with the tool \"%s\".\n",
                name());
        S->addMessageF(L_ERROR, msg);
        S->printOpenCLError(err_code);
        return true;
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
        sprintf(msg,
                "Failure copying the permutations with the tool \"%s\".\n",
                name());
        S->addMessageF(L_ERROR, msg);
        S->printOpenCLError(err_code);
        return true;
    }

    if(inversePermutations())
        return true;

    return false;
}

bool RadixSort::init()
{
    cl_int err_code;
    char msg[1024];
    InputOutput::ScreenManager *S = InputOutput::ScreenManager::singleton();
    CalcServer *C = CalcServer::singleton();

    err_code = clSetKernelArg(_init_kernel,
                              0,
                              sizeof(cl_mem),
                              (void*)&_in_permut);
    if(err_code != CL_SUCCESS){
        S->addMessageF(L_ERROR, "Failure sending argument 0 to \"init\"\n");
        S->printOpenCLError(err_code);
        return true;
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
        sprintf(msg,
                "Failure executing the tool \"%s\" kernel \"init\".\n",
                name());
        S->addMessageF(L_ERROR, msg);
        S->printOpenCLError(err_code);
        return true;
    }

    return false;
}

bool RadixSort::histograms()
{
    cl_int err_code;
    char msg[1024];
    InputOutput::ScreenManager *S = InputOutput::ScreenManager::singleton();
    CalcServer *C = CalcServer::singleton();
    size_t local_work_size = _items;
    size_t global_work_size = _groups * _items;

    err_code = clSetKernelArg(_histograms_kernel,
                              0,
                              sizeof(cl_mem),
                              (void*)&_in_keys);
    if(err_code != CL_SUCCESS){
        S->addMessageF(L_ERROR, "Failure sending argument 0 to \"histogram\"\n");
        S->printOpenCLError(err_code);
        return true;
    }
    err_code = clSetKernelArg(_histograms_kernel,
                              2,
                              sizeof(cl_uint),
                              (void*)&_pass);
    if(err_code != CL_SUCCESS){
        S->addMessageF(L_ERROR, "Failure sending argument 2 to \"histogram\"\n");
        S->printOpenCLError(err_code);
        return true;
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
        sprintf(msg,
                "Failure executing the tool \"%s\" kernel \"histogram\".\n",
                name());
        S->addMessageF(L_ERROR, msg);
        S->printOpenCLError(err_code);
        return true;
    }

    return false;
}

bool RadixSort::scan()
{
    cl_int err_code;
    char msg[1024];
    InputOutput::ScreenManager *S = InputOutput::ScreenManager::singleton();
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
        S->addMessageF(L_ERROR, "Failure sending argument 0 to \"scan\" (1st call)\n");
        S->printOpenCLError(err_code);
        return true;
    }
    err_code = clSetKernelArg(_scan_kernel,
                              2,
                              sizeof(cl_mem),
                              (void*)&_global_sums);
    if(err_code != CL_SUCCESS){
        S->addMessageF(L_ERROR, "Failure sending argument 2 to \"scan\" (1st call)\n");
        S->printOpenCLError(err_code);
        return true;
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
        sprintf(msg,
                "Failure executing the tool \"%s\" kernel \"scan\" (1st call).\n",
                name());
        S->addMessageF(L_ERROR, msg);
        S->printOpenCLError(err_code);
        return true;
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
        S->addMessageF(L_ERROR, "Failure sending argument 0 to \"scan\" (2nd call)\n");
        S->printOpenCLError(err_code);
        return true;
    }
    err_code = clSetKernelArg(_scan_kernel,
                              2,
                              sizeof(cl_mem),
                              (void*)&_temp_mem);
    if(err_code != CL_SUCCESS){
        S->addMessageF(L_ERROR, "Failure sending argument 2 to \"scan\" (2nd call)\n");
        S->printOpenCLError(err_code);
        return true;
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
        sprintf(msg,
                "Failure executing the tool \"%s\" kernel \"scan\" (2nd call).\n",
                name());
        S->addMessageF(L_ERROR, msg);
        S->printOpenCLError(err_code);
        return true;
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
        sprintf(msg,
                "Failure executing the tool \"%s\" kernel \"paste\".\n",
                name());
        S->addMessageF(L_ERROR, msg);
        S->printOpenCLError(err_code);
        return true;
    }

    return false;
}

bool RadixSort::reorder()
{
    cl_int err_code;
    char msg[1024];
    InputOutput::ScreenManager *S = InputOutput::ScreenManager::singleton();
    CalcServer *C = CalcServer::singleton();
    size_t local_work_size = _items;
    size_t global_work_size = _groups * _items;

    err_code = clSetKernelArg(_sort_kernel,
                              0,
                              sizeof(cl_mem),
                              (void*)&_in_keys);
    if(err_code != CL_SUCCESS){
        S->addMessageF(L_ERROR, "Failure sending argument 0 to \"sort\"\n");
        S->printOpenCLError(err_code);
        return true;
    }
    err_code = clSetKernelArg(_sort_kernel,
                              1,
                              sizeof(cl_mem),
                              (void*)&_out_keys);
    if(err_code != CL_SUCCESS){
        S->addMessageF(L_ERROR, "Failure sending argument 1 to \"sort\"\n");
        S->printOpenCLError(err_code);
        return true;
    }
    err_code = clSetKernelArg(_sort_kernel,
                              3,
                              sizeof(cl_uint),
                              (void*)&_pass);
    if(err_code != CL_SUCCESS){
        S->addMessageF(L_ERROR, "Failure sending argument 3 to \"sort\"\n");
        S->printOpenCLError(err_code);
        return true;
    }
    err_code = clSetKernelArg(_sort_kernel,
                              4,
                              sizeof(cl_mem),
                              (void*)&_in_permut);
    if(err_code != CL_SUCCESS){
        S->addMessageF(L_ERROR, "Failure sending argument 4 to \"sort\"\n");
        S->printOpenCLError(err_code);
        return true;
    }
    err_code = clSetKernelArg(_sort_kernel,
                              5,
                              sizeof(cl_mem),
                              (void*)&_out_permut);
    if(err_code != CL_SUCCESS){
        S->addMessageF(L_ERROR, "Failure sending argument 5 to \"sort\"\n");
        S->printOpenCLError(err_code);
        return true;
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
        sprintf(msg,
                "Failure executing the tool \"%s\" kernel \"paste\".\n",
                name());
        S->addMessageF(L_ERROR, msg);
        S->printOpenCLError(err_code);
        return true;
    }

    // Swap the memory identifiers for the next pass
    cl_mem d_temp;

    d_temp = _in_keys;
    _in_keys = _out_keys;
    _out_keys = d_temp;

    d_temp = _in_permut;
    _in_permut = _out_permut;
    _out_permut = d_temp;

    return false;
}

bool RadixSort::inversePermutations()
{
    cl_int err_code;
    char msg[1024];
    InputOutput::ScreenManager *S = InputOutput::ScreenManager::singleton();
    CalcServer *C = CalcServer::singleton();

    err_code = clSetKernelArg(_inv_perms_kernel,
                              0,
                              sizeof(cl_mem),
                              (void*)&_in_permut);
    if(err_code != CL_SUCCESS){
        S->addMessageF(L_ERROR, "Failure sending argument 0 to \"inversePermutation\"\n");
        S->printOpenCLError(err_code);
        return true;
    }
    err_code = clSetKernelArg(_inv_perms_kernel,
                              1,
                              sizeof(cl_mem),
                              _inv_perms->get());
    if(err_code != CL_SUCCESS){
        S->addMessageF(L_ERROR, "Failure sending argument 1 to \"inversePermutation\"\n");
        S->printOpenCLError(err_code);
        return true;
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
        sprintf(msg,
                "Failure executing the tool \"%s\" kernel \"inversePermutation\".\n",
                name());
        S->addMessageF(L_ERROR, msg);
        S->printOpenCLError(err_code);
        return true;
    }

    return false;
}


bool RadixSort::variables()
{
    size_t n;
    char msg[1024];
    InputOutput::ScreenManager *S = InputOutput::ScreenManager::singleton();
    CalcServer *C = CalcServer::singleton();
    InputOutput::Variables vars = C->variables();

    // Check and get the variables
    if(!vars.get(_var_name)){
        sprintf(msg,
                "Undeclared variable \"%s\" cannot be sorted by tool \"%s\".\n",
                _var_name,
                name());
        S->addMessageF(L_ERROR, msg);
        return true;
    }
    if(vars.get(_var_name)->type().compare("unsigned int*")){
        sprintf(msg,
                "Wrong type for the variable \"%s\" (tool: \"%s\").\n",
                _var_name,
                name());
        S->addMessageF(L_ERROR, msg);
        sprintf(msg,
                "\t\"%s\" was expected, but \"%s\" has been found.\n",
                "unsigned int*",
                vars.get(_var_name)->type());
        S->addMessage(L_DEBUG, msg);
        return true;
    }
    _var = (InputOutput::ArrayVariable *)vars.get(_var_name);

    if(!vars.get(_perms_name)){
        sprintf(msg,
                "Undeclared permutations variable \"%s\" (tool \"%s\").\n",
                _perms_name,
                name());
        S->addMessageF(L_ERROR, msg);
        return true;
    }
    if(vars.get(_perms_name)->type().compare("unsigned int*")){
        sprintf(msg,
                "Wrong type for the variable \"%s\" (tool: \"%s\").\n",
                _perms_name,
                name());
        S->addMessageF(L_ERROR, msg);
        sprintf(msg,
                "\t\"%s\" was expected, but \"%s\" has been found.\n",
                "unsigned int*",
                vars.get(_perms_name)->type());
        S->addMessage(L_DEBUG, msg);
        return true;
    }
    _perms = (InputOutput::ArrayVariable *)vars.get(_perms_name);

    if(!vars.get(_inv_perms_name)){
        sprintf(msg,
                "Undeclared permutations variable \"%s\" (tool \"%s\").\n",
                _inv_perms_name,
                name());
        S->addMessageF(L_ERROR, msg);
        return true;
    }
    if(vars.get(_inv_perms_name)->type().compare("unsigned int*")){
        sprintf(msg,
                "Wrong type for the variable \"%s\" (tool: \"%s\").\n",
                _inv_perms_name,
                name());
        S->addMessageF(L_ERROR, msg);
        sprintf(msg,
                "\t\"%s\" was expected, but \"%s\" has been found.\n",
                "unsigned int*",
                vars.get(_inv_perms_name)->type());
        S->addMessage(L_DEBUG, msg);
        return true;
    }
    _inv_perms = (InputOutput::ArrayVariable *)vars.get(_inv_perms_name);

    // Check the lengths
    n = _var->size() / vars.typeToBytes(_var->type());
    if(!isPowerOf2(n)){
        sprintf(msg,
                "Tool \"%s\" cannot sort the variable \"%s\".\n",
                name(),
                _var_name);
        S->addMessageF(L_ERROR, msg);
        sprintf(msg,
                "\tIt has a length %u, which is not power of 2.\n",
                n);
        S->addMessage(L_DEBUG, msg);
        return true;
    }
    _n = n;
    n = _perms->size() / vars.typeToBytes(_perms->type());
    if(n != _n){
        sprintf(msg,
                "Lengths mismatch in the tool \"%s\".\n",
                name());
        S->addMessageF(L_ERROR, msg);
        sprintf(msg,
                "\tVariable \"%s\" has a length n=%u.\n",
                _var->name(),
                _n);
        S->addMessage(L_DEBUG, msg);
        sprintf(msg,
                "\tVariable \"%s\" has a length n=%u.\n",
                _perms->name(),
                n);
        S->addMessage(L_DEBUG, msg);
        return true;
    }
    n = _inv_perms->size() / vars.typeToBytes(_inv_perms->type());
    if(n != _n){
        sprintf(msg,
                "Lengths mismatch in the tool \"%s\".\n",
                name());
        S->addMessageF(L_ERROR, msg);
        sprintf(msg,
                "\tVariable \"%s\" has a length n=%u.\n",
                _var->name(),
                _n);
        S->addMessage(L_DEBUG, msg);
        sprintf(msg,
                "\tVariable \"%s\" has a length n=%u.\n",
                _inv_perms->name(),
                n);
        S->addMessage(L_DEBUG, msg);
        return true;
    }

    return false;
}

bool RadixSort::setupOpenCL()
{
    cl_int err_code;
    cl_kernel kernel;
    char msg[1024];
    InputOutput::ScreenManager *S = InputOutput::ScreenManager::singleton();
    CalcServer *C = CalcServer::singleton();
    InputOutput::Variables vars = C->variables();

    // Create a header for the source code where the operation will be placed
    char header[RADIXSORT_INC_LEN + 128];
    strcpy(header, "");
    strncat(header, RADIXSORT_INC, RADIXSORT_INC_LEN);
    strcat(header, "");

    // Setup the complete source code
    char source[strlen(header) + strlen(RADIXSORT_SRC) + 1];
    strcpy(source, header);
    strncat(source, RADIXSORT_SRC, RADIXSORT_SRC_LEN);
    strcat(source, "");

    // Compile the kernels
    if(compile(source)){
        return true;
    }
    // Check and correct _items, _groups and _histo_split
    if(setupDims()){
        return true;
    }
    // Setup the memory objects
    if(setupMems()){
        return true;
    }
    if(setupArgs()){
        return true;
    }

    sprintf(msg, "\titems: %u\n", _items);
    S->addMessage(L_DEBUG, msg);
    sprintf(msg, "\tgroups: %u\n", _groups);
    S->addMessage(L_DEBUG, msg);
    sprintf(msg, "\tsplits: %u\n", _histo_split);
    S->addMessage(L_DEBUG, msg);

    return false;
}

bool RadixSort::compile(const char* source)
{
    cl_int err_code;
    cl_program program;
    char msg[1024];
    InputOutput::ScreenManager *S = InputOutput::ScreenManager::singleton();
    CalcServer *C = CalcServer::singleton();

    char flags[512];
    sprintf(flags,
            "-D_BITS=%u -D_RADIX=%u -DPERMUT",
            _bits,
            _radix);
    #ifdef AQUA_DEBUG
        strcat(flags, " -DDEBUG");
    #else
        strcat(flags, " -DNDEBUG");
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

    _init_kernel = clCreateKernel(program, "init", &err_code);
    if(err_code != CL_SUCCESS) {
        S->addMessageF(L_ERROR, "Failure creating the \"init\" kernel.\n");
        S->printOpenCLError(err_code);
        clReleaseProgram(program);
        return true;
    }
    _histograms_kernel = clCreateKernel(program, "histogram", &err_code);
    if(err_code != CL_SUCCESS) {
        S->addMessageF(L_ERROR, "Failure creating the \"histogram\" kernel.\n");
        S->printOpenCLError(err_code);
        clReleaseProgram(program);
        return true;
    }
    _scan_kernel = clCreateKernel(program, "scan", &err_code);
    if(err_code != CL_SUCCESS) {
        S->addMessageF(L_ERROR, "Failure creating the \"scan\" kernel.\n");
        S->printOpenCLError(err_code);
        clReleaseProgram(program);
        return true;
    }
    _paste_kernel = clCreateKernel(program, "paste", &err_code);
    if(err_code != CL_SUCCESS) {
        S->addMessageF(L_ERROR, "Failure creating the \"paste\" kernel.\n");
        S->printOpenCLError(err_code);
        clReleaseProgram(program);
        return true;
    }
    _sort_kernel = clCreateKernel(program, "sort", &err_code);
    if(err_code != CL_SUCCESS) {
        S->addMessageF(L_ERROR, "Failure creating the \"sort\" kernel.\n");
        S->printOpenCLError(err_code);
        clReleaseProgram(program);
        return true;
    }
    _inv_perms_kernel = clCreateKernel(program, "inversePermutation", &err_code);
    if(err_code != CL_SUCCESS) {
        S->addMessageF(L_ERROR, "Failure creating the \"inversePermutation\" kernel.\n");
        S->printOpenCLError(err_code);
        clReleaseProgram(program);
        return true;
    }
    clReleaseProgram(program);

    return false;
}

bool RadixSort::setupDims()
{
    cl_int err_code;
    size_t max_local_work_size, sort_local_work_size, scan_local_work_size;
    char msg[1024];
    InputOutput::ScreenManager *S = InputOutput::ScreenManager::singleton();
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
        S->addMessageF(L_ERROR, "Failure getting CL_KERNEL_WORK_GROUP_SIZE from \"histogram\".\n");
        S->printOpenCLError(err_code);
        return true;
    }
    err_code = clGetKernelWorkGroupInfo(_sort_kernel,
                                        C->device(),
                                        CL_KERNEL_WORK_GROUP_SIZE,
                                        sizeof(size_t),
                                        &sort_local_work_size,
                                        NULL);
    if(err_code != CL_SUCCESS) {
        S->addMessageF(L_ERROR, "Failure getting CL_KERNEL_WORK_GROUP_SIZE from \"sort\".\n");
        S->printOpenCLError(err_code);
        return true;
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
        S->addMessageF(L_ERROR, "Failure getting CL_KERNEL_WORK_GROUP_SIZE from \"scan\".\n");
        S->printOpenCLError(err_code);
        return true;
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
        S->addMessageF(L_ERROR, "Failure getting CL_KERNEL_WORK_GROUP_SIZE from \"paste\".\n");
        S->printOpenCLError(err_code);
        return true;
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
        S->addMessageF(L_ERROR, "Failure setting a number of items and groups compatible with \"scan\" and \"paste\".\n");
        S->addMessage(L_DEBUG, "\tYou can try to recompile the code decreasing __CL_MIN_LOCALSIZE__\n");
        return true;
    }

    _local_work_size = getLocalWorkSize(_n, C->command_queue());
    _global_work_size = getGlobalWorkSize(_n, _local_work_size);

    return false;
}

bool RadixSort::setupMems()
{
    cl_int err_code;
    char msg[1024];
    InputOutput::ScreenManager *S = InputOutput::ScreenManager::singleton();
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
        S->addMessageF(L_ERROR, "Buffer memory allocation failure.\n");
        S->printOpenCLError(err_code);
        return true;
    }
    _out_keys = clCreateBuffer(C->context(),
                              CL_MEM_READ_WRITE,
                              _n * sizeof(unsigned int),
                              NULL,
                              &err_code);
    if(err_code != CL_SUCCESS) {
        S->addMessageF(L_ERROR, "Buffer memory allocation failure.\n");
        S->printOpenCLError(err_code);
        return true;
    }
    _in_permut = clCreateBuffer(C->context(),
                                CL_MEM_READ_WRITE,
                                _n * sizeof(unsigned int),
                                NULL,
                                &err_code);
    if(err_code != CL_SUCCESS) {
        S->addMessageF(L_ERROR, "Buffer memory allocation failure.\n");
        S->printOpenCLError(err_code);
        return true;
    }
    _out_permut = clCreateBuffer(C->context(),
                                 CL_MEM_READ_WRITE,
                                 _n * sizeof(unsigned int),
                                 NULL,
                                 &err_code);
    if(err_code != CL_SUCCESS) {
        S->addMessageF(L_ERROR, "Buffer memory allocation failure.\n");
        S->printOpenCLError(err_code);
        return true;
    }
    _histograms = clCreateBuffer(C->context(),
                                 CL_MEM_READ_WRITE,
                                 (_radix * _groups * _items) * sizeof(unsigned int),
                                 NULL,
                                 &err_code);
    if(err_code != CL_SUCCESS) {
        S->addMessageF(L_ERROR, "Buffer memory allocation failure.\n");
        S->printOpenCLError(err_code);
        return true;
    }
    _global_sums = clCreateBuffer(C->context(),
                                  CL_MEM_READ_WRITE,
                                  _histo_split * sizeof(unsigned int),
                                  NULL,
                                  &err_code);
    if(err_code != CL_SUCCESS) {
        S->addMessageF(L_ERROR, "Buffer memory allocation failure.\n");
        S->printOpenCLError(err_code);
        return true;
    }
    _temp_mem = clCreateBuffer(C->context(),
                               CL_MEM_READ_WRITE,
                               sizeof(unsigned int),
                               NULL,
                               &err_code);
    if(err_code != CL_SUCCESS) {
        S->addMessageF(L_ERROR, "Buffer memory allocation failure.\n");
        S->printOpenCLError(err_code);
        return true;
    }

    allocatedMemory(4 * _n * sizeof(unsigned int) +
                    (_radix * _groups * _items) * sizeof(unsigned int) +
                    _histo_split * sizeof(unsigned int) +
                    sizeof(unsigned int));

    return false;
}

bool RadixSort::setupArgs()
{
    cl_int err_code;
    char msg[1024];
    InputOutput::ScreenManager *S = InputOutput::ScreenManager::singleton();
    CalcServer *C = CalcServer::singleton();

    err_code = clSetKernelArg(_init_kernel,
                              1,
                              sizeof(cl_uint),
                              (void*)&_n);
    if(err_code != CL_SUCCESS){
        S->addMessageF(L_ERROR, "Failure sending argument 1 to \"init\"\n");
        S->printOpenCLError(err_code);
        return true;
    }

    err_code = clSetKernelArg(_histograms_kernel,
                              1,
                              sizeof(cl_mem),
                              (void*)&_histograms);
    if(err_code != CL_SUCCESS){
        S->addMessageF(L_ERROR, "Failure sending argument 1 to \"histogram\"\n");
        S->printOpenCLError(err_code);
        return true;
    }
    err_code = clSetKernelArg(_histograms_kernel,
                              3,
                              sizeof(cl_uint) * _radix * _items,
                              NULL);
    if(err_code != CL_SUCCESS){
        S->addMessageF(L_ERROR, "Failure sending argument 3 to \"histogram\"\n");
        S->printOpenCLError(err_code);
        return true;
    }
    err_code = clSetKernelArg(_histograms_kernel,
                              4,
                              sizeof(cl_uint),
                              (void*)&_n);
    if(err_code != CL_SUCCESS){
        S->addMessageF(L_ERROR, "Failure sending argument 4 to \"histogram\"\n");
        S->printOpenCLError(err_code);
        return true;
    }

    unsigned int maxmemcache = max(_histo_split,
                                   _items * _groups * _radix / _histo_split);
    err_code = clSetKernelArg(_scan_kernel,
                              1,
                              sizeof(cl_uint) * maxmemcache,
                              NULL);
    if(err_code != CL_SUCCESS){
        S->addMessageF(L_ERROR, "Failure sending argument 1 to \"scan\"\n");
        S->printOpenCLError(err_code);
        return true;
    }

    err_code = clSetKernelArg(_paste_kernel,
                              0,
                              sizeof(cl_mem),
                              (void*)&_histograms);
    if(err_code != CL_SUCCESS){
        S->addMessageF(L_ERROR, "Failure sending argument 0 to \"paste\"\n");
        S->printOpenCLError(err_code);
        return true;
    }
    err_code = clSetKernelArg(_paste_kernel,
                              1,
                              sizeof(cl_mem),
                              (void*)&_global_sums);
    if(err_code != CL_SUCCESS){
        S->addMessageF(L_ERROR, "Failure sending argument 1 to \"paste\"\n");
        S->printOpenCLError(err_code);
        return true;
    }

    err_code = clSetKernelArg(_sort_kernel,
                              2,
                              sizeof(cl_mem),
                              (void*)&_histograms);
    if(err_code != CL_SUCCESS){
        S->addMessageF(L_ERROR, "Failure sending argument 2 to \"sort\"\n");
        S->printOpenCLError(err_code);
        return true;
    }
    err_code = clSetKernelArg(_sort_kernel,
                              6,
                              sizeof(cl_uint) * _radix * _items,
                              NULL);
    if(err_code != CL_SUCCESS){
        S->addMessageF(L_ERROR, "Failure sending argument 6 to \"sort\"\n");
        S->printOpenCLError(err_code);
        return true;
    }
    err_code = clSetKernelArg(_sort_kernel,
                              7,
                              sizeof(cl_uint),
                              (void*)&_n);
    if(err_code != CL_SUCCESS){
        S->addMessageF(L_ERROR, "Failure sending argument 7 to \"sort\"\n");
        S->printOpenCLError(err_code);
        return true;
    }

    err_code = clSetKernelArg(_inv_perms_kernel,
                              1,
                              sizeof(cl_mem),
                              _inv_perms->get());
    if(err_code != CL_SUCCESS){
        S->addMessageF(L_ERROR, "Failure sending argument 1 to \"inversePermutation\"\n");
        S->printOpenCLError(err_code);
        return true;
    }
    err_code = clSetKernelArg(_inv_perms_kernel,
                              2,
                              sizeof(cl_uint),
                              (void*)&_n);
    if(err_code != CL_SUCCESS){
        S->addMessageF(L_ERROR, "Failure sending argument 2 to \"inversePermutation\"\n");
        S->printOpenCLError(err_code);
        return true;
    }

    return false;
}

}}  // namespace
