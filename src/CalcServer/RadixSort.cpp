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
 * @brief Sorting permutations processor.
 * (See Aqua::CalcServer::RadixSort for details)
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

#include <AuxiliarMethods.h>

#include <ProblemSetup.h>
#include <ScreenManager.h>
#include <CalcServer/RadixSort.h>
#include <CalcServer.h>

namespace Aqua{ namespace CalcServer{

RadixSort::RadixSort()
    : _path(0)
    , _n(0)
    , _in_keys(0)
    , _out_keys(0)
    , _in_permut(0)
    , _out_permut(0)
    , _histograms(0)
    , _global_sums(0)
    , _temp_mem(0)
    , _init_kernel(0)
    , _transpose_kernel(0)
    , _histograms_kernel(0)
    , _histograms_scan_kernel(0)
    , _paste_histograms_kernel(0)
    , _sort_kernel(0)
    , _inv_permutations_kernel(0)
    , _program(0)
    #ifdef HAVE_GPUPROFILE
        , _time(0.f)
    #endif
    , _items(_ITEMS)
    , _groups(_GROUPS)
    , _bits(_STEPBITS)
    , _radix(_RADIX)
    , _histo_split(_HISTOSPLIT)
    , _tilesize_warn(false)
{
    InputOutput::ScreenManager *S = InputOutput::ScreenManager::singleton();
    S->addMessageF(1, "Initializating radix sort tool...\n");
    InputOutput::ProblemSetup *P = InputOutput::ProblemSetup::singleton();
    CalcServer *C = CalcServer::singleton();
    int str_len = 0;
    //! 1st.- Get data
    str_len = strlen(P->OpenCL_kernels.radix_sort);
    if(str_len <= 0){
        S->addMessageF(3, "Path of predictor kernel is empty.\n");
        exit(EXIT_FAILURE);
    }
    _path = new char[str_len+4];
    if(!_path) {
        S->addMessageF(3, "Memory cannot be allocated for the path.\n");
        exit(EXIT_FAILURE);
    }
    strcpy(_path, P->OpenCL_kernels.radix_sort);
    strcat(_path, ".cl");

    if(setN(C->num_icell))
        exit(EXIT_FAILURE);

    if(setupOpenCL())
        exit(EXIT_FAILURE);
    S->addMessageF(1, "Radix sort ready to work!\n");
}

RadixSort::~RadixSort()
{
    if(_path)delete[] _path;_path=0;
    _in_keys = NULL;                      // Let CalcServer destroy it
    _in_permut = NULL;                    // Let CalcServer destroy it
    if(_out_keys)
        clReleaseMemObject(_out_keys);
    _out_keys=0;
    if(_out_permut)
        clReleaseMemObject(_out_permut);
    _out_permut=0;
    if(_histograms)
        clReleaseMemObject(_histograms);
    _histograms=0;
    if(_global_sums)
        clReleaseMemObject(_global_sums);
    _global_sums=0;
    if(_temp_mem)
        clReleaseMemObject(_temp_mem);
    _temp_mem=0;
    if(_init_kernel)
        clReleaseKernel(_init_kernel);
    _init_kernel=0;
    if(_transpose_kernel)
        clReleaseKernel(_transpose_kernel);
    _transpose_kernel=0;
    if(_histograms_kernel)
        clReleaseKernel(_histograms_kernel);
    _histograms_kernel=0;
    if(_histograms_scan_kernel)
        clReleaseKernel(_histograms_scan_kernel);
    _histograms_scan_kernel=0;
    if(_paste_histograms_kernel)
        clReleaseKernel(_paste_histograms_kernel);
    _paste_histograms_kernel=0;
    if(_sort_kernel)
        clReleaseKernel(_sort_kernel);
    _sort_kernel=0;
    if(_inv_permutations_kernel)
        clReleaseKernel(_inv_permutations_kernel);
    _inv_permutations_kernel=0;
    if(_program)
        clReleaseProgram(_program);
    _program=0;
}

bool RadixSort::sort()
{
    InputOutput::ScreenManager *S = InputOutput::ScreenManager::singleton();
    CalcServer *C = CalcServer::singleton();
    #ifdef HAVE_GPUPROFILE
        _time = 0.f;
    #endif
    // Get maximum key bits, and needed pass
    unsigned int i, max_icell;
    max_icell = nextPowerOf2(C->num_cells);
    for(i=0; (max_icell&1) == 0; max_icell >>= 1, i++);
    _key_bits = i;
    _key_bits = roundUp(_key_bits, _bits);
    if(_key_bits > __UINTBITS__){
        S->addMessageF(3, "Resultant keys overflows unsigned int type.\n");
        return true;
    }
    _n_pass = _key_bits / _bits;
    unsigned int nbcol=_n/(_groups * _items);
    unsigned int nbrow=_groups * _items;
    if(init())
        return true;
    // if(transpose(nbrow,nbcol))
    //     return true;

    for(_pass=0;_pass<_n_pass;_pass++){
        // Create histograms
        if(histograms())
            return true;
        // Scan histograms and paste it
        if(scan())
            return true;
        // Reorder data
        if(reorder())
            return true;
    }

    // if(transpose(nbcol,nbrow))
    //     return true;

    C->icell = _in_keys;            // Take care, swaped into reorder
    C->permutation = _in_permut;    // Take care, swaped into reorder
    if(reversePermutations())
        return true;
    return false;
}

bool RadixSort::init()
{
    InputOutput::ScreenManager *S = InputOutput::ScreenManager::singleton();
    CalcServer *C = CalcServer::singleton();
    int err_code=0;

    err_code |= sendArgument(_init_kernel,
                             0,
                             sizeof(cl_mem),
                             (void*)&_in_permut);
    if(err_code != CL_SUCCESS) {
        S->addMessageF(3, "I cannot send a variable to the kernel.\n");
        return true;
    }
    //! 2nd.- Execute the kernel
    #ifdef HAVE_GPUPROFILE
        cl_event event;
        cl_ulong end, start;
        profileTime(0.f);
        err_code = clEnqueueNDRangeKernel(C->command_queue,
                                          _init_kernel,
                                          1,
                                          NULL,
                                          &_global_work_size,
                                          NULL,
                                          0,
                                          NULL,
                                          &event);
    #else
        err_code = clEnqueueNDRangeKernel(C->command_queue,
                                          _init_kernel,
                                          1,
                                          NULL,
                                          &_global_work_size,
                                          NULL,
                                          0,
                                          NULL,
                                          NULL);
    #endif
    if(err_code != CL_SUCCESS) {
        S->addMessageF(3, "I cannot execute the kernel.\n");
        S->printOpenCLError(err_code);
        return true;
    }
    #ifdef HAVE_GPUPROFILE
        err_code = clWaitForEvents(1, &event);
        if(err_code != CL_SUCCESS) {
            S->addMessage(3, "Impossible to wait for the kernels end.\n");
            S->printOpenCLError(err_code);
            return true;
        }
        err_code |= clGetEventProfilingInfo(event,
                                            CL_PROFILING_COMMAND_END,
                                            sizeof(cl_ulong),
                                            &end,
                                            0);
        if(err_code != CL_SUCCESS) {
            S->addMessage(3, "I cannot profile the kernel execution.\n");
            S->printOpenCLError(err_code);
            return true;
        }
        err_code |= clGetEventProfilingInfo(event,
                                            CL_PROFILING_COMMAND_START,
                                            sizeof(cl_ulong),
                                            &start,
                                            0);
        if(err_code != CL_SUCCESS) {
            S->addMessage(3, "I cannot profile the kernel execution.\n");
            S->printOpenCLError(err_code);
            return true;
        }
    #endif
    return false;
}

bool RadixSort::transpose(unsigned int nbrow, unsigned int nbcol)
{
    InputOutput::ScreenManager *S = InputOutput::ScreenManager::singleton();
    CalcServer *C = CalcServer::singleton();
    int err_code=0;
    unsigned int tilesize=__CL_MIN_LOCALSIZE__;
    if (nbrow%tilesize != 0) tilesize=1;
    if (nbcol%tilesize != 0) tilesize=1;
    if( _tilesize_warn ){
        tilesize=1;
    }
    else if (tilesize == 1) {
        S->addMessage(1, "List too small, local memory usage will be avoided.\n");
        _tilesize_warn=true;
    }

    err_code |= sendArgument(_transpose_kernel, 0, sizeof(cl_mem), (void*)&_in_keys);
    err_code |= sendArgument(_transpose_kernel, 1, sizeof(cl_mem), (void*)&_out_keys);
    err_code |= sendArgument(_transpose_kernel, 2, sizeof(cl_uint), (void*)&nbcol);
    err_code |= sendArgument(_transpose_kernel, 3, sizeof(cl_uint), (void*)&nbrow);
    err_code |= sendArgument(_transpose_kernel, 4, sizeof(cl_mem), (void*)&_in_permut);
    err_code |= sendArgument(_transpose_kernel, 5, sizeof(cl_mem), (void*)&_out_permut);
    err_code |= sendArgument(_transpose_kernel, 6, sizeof(cl_uint)*tilesize*tilesize, NULL);
    err_code |= sendArgument(_transpose_kernel, 7, sizeof(cl_uint)*tilesize*tilesize, NULL);
    err_code |= sendArgument(_transpose_kernel, 8, sizeof(cl_uint), (void*)&tilesize);
    if(err_code != CL_SUCCESS) {
        S->addMessageF(3, "I cannot send a variable to the kernel.\n");
        return 1;
    }

    size_t _global_work_size[2];
    size_t _local_work_size[2];
    _global_work_size[0]=nbrow/tilesize;
    _global_work_size[1]=nbcol;
    _local_work_size[0]=1;
    _local_work_size[1]=tilesize;
    if(_local_work_size[1] < __CL_MIN_LOCALSIZE__){
        _local_work_size[1]  *= __CL_MIN_LOCALSIZE__ / tilesize;
        _global_work_size[1] *= __CL_MIN_LOCALSIZE__ / tilesize;
    }
    #ifdef HAVE_GPUPROFILE
        cl_event event;
        cl_ulong end, start;
        err_code = clEnqueueNDRangeKernel(C->command_queue,
                                          _transpose_kernel,
                                          2,
                                          NULL,
                                          _global_work_size,
                                          _local_work_size,
                                          0,
                                          NULL,
                                          &event);
    #else
        err_code = clEnqueueNDRangeKernel(C->command_queue,
                                       _transpose_kernel,
                                       2,
                                       NULL,
                                       _global_work_size,
                                       _local_work_size,
                                       0,
                                       NULL,
                                       NULL);
    #endif
    if(err_code != CL_SUCCESS) {
        S->addMessage(3, "(RadixSort::_transpose): I cannot execute the kernel.\n");
        S->printOpenCLError(err_code);
        return true;
    }
    #ifdef HAVE_GPUPROFILE
        err_code = clWaitForEvents(1, &event);
        if(err_code != CL_SUCCESS) {
            S->addMessage(3, "Impossible to wait for the kernels end.\n");
            S->printOpenCLError(err_code);
            return true;
        }
        err_code |= clGetEventProfilingInfo(event,
                                            CL_PROFILING_COMMAND_END,
                                            sizeof(cl_ulong),
                                            &end,
                                            0);
        if(err_code != CL_SUCCESS) {
            S->addMessage(3, "I cannot profile the kernel execution.\n");
            S->printOpenCLError(err_code);
            return true;
        }
        err_code |= clGetEventProfilingInfo(event,
                                            CL_PROFILING_COMMAND_START,
                                            sizeof(cl_ulong),
                                            &start,
                                            0);
        if(err_code != CL_SUCCESS) {
            S->addMessage(3, "I cannot profile the kernel execution.\n");
            S->printOpenCLError(err_code);
            return true;
        }
        _time += (end - start)/1000.f;  // 10^-3 ms
    #endif
    //! 3rd.- Swap input/output arrays
    cl_mem d_temp;

    d_temp=_in_keys;
    _in_keys=_out_keys;
    _out_keys=d_temp;

    d_temp=_in_permut;
    _in_permut=_out_permut;
    _out_permut=d_temp;
    return false;
}

bool RadixSort::histograms()
{
    InputOutput::ScreenManager *S = InputOutput::ScreenManager::singleton();
    CalcServer *C = CalcServer::singleton();
    int err_code=0;
    size_t local_size = _items;
    size_t global_size = _groups*_items;
    //! 1st.- Send arguments
    err_code |= sendArgument(_histograms_kernel,
                             0,
                             sizeof(cl_mem),
                             (void*)&_in_keys);
    err_code |= sendArgument(_histograms_kernel,
                             2,
                             sizeof(cl_uint),
                             (void*)&_pass);
    err_code |= sendArgument(_histograms_kernel,
                             3,
                             sizeof(cl_uint)*_radix*_items,
                             NULL);
    if(err_code != CL_SUCCESS) {
        S->addMessageF(3, "I cannot send a variable to the kernel.\n");
        return true;
    }
    //! 2nd.- Execute the kernel
    #ifdef HAVE_GPUPROFILE
        cl_event event;
        cl_ulong end, start;
        err_code = clEnqueueNDRangeKernel(C->command_queue,
                                          _histograms_kernel,
                                          1,
                                          NULL,
                                          &global_size,
                                          &local_size,
                                          0,
                                          NULL,
                                          &event);
    #else
        err_code = clEnqueueNDRangeKernel(C->command_queue,
                                          _histograms_kernel,
                                          1,
                                          NULL,
                                          &global_size,
                                          &local_size,
                                          0,
                                          NULL,
                                          NULL);
    #endif
    if(err_code != CL_SUCCESS) {
        S->addMessageF(3, "I cannot execute the kernel.\n");
        S->printOpenCLError(err_code);
        return true;
    }
    #ifdef HAVE_GPUPROFILE
        err_code = clWaitForEvents(1, &event);
        if(err_code != CL_SUCCESS) {
            S->addMessage(3, "Impossible to wait for the kernels end.\n");
            S->printOpenCLError(err_code);
            return true;
        }
        err_code |= clGetEventProfilingInfo(event,
                                            CL_PROFILING_COMMAND_END,
                                            sizeof(cl_ulong),
                                            &end,
                                            0);
        if(err_code != CL_SUCCESS) {
            S->addMessage(3, "I cannot profile the kernel execution.\n");
            S->printOpenCLError(err_code);
            return true;
        }
        err_code |= clGetEventProfilingInfo(event,
                                            CL_PROFILING_COMMAND_START,
                                            sizeof(cl_ulong),
                                            &start,
                                            0);
        if(err_code != CL_SUCCESS) {
            S->addMessage(3, "I cannot profile the kernel execution.\n");
            S->printOpenCLError(err_code);
            return true;
        }
        _time += (end - start)/1000.f;  // 10^-3 ms
    #endif
    return false;
}

bool RadixSort::scan()
{
    InputOutput::ScreenManager *S = InputOutput::ScreenManager::singleton();
    CalcServer *C = CalcServer::singleton();
    int err_code=0;
    size_t global_size = _radix*_groups*_items/2;
    size_t local_size = global_size / _histo_split;
    unsigned int maxmemcache=max(_histo_split,
                                 _items * _groups * _radix / _histo_split);

    // 1st scan
    // ========
    err_code |= sendArgument(_histograms_scan_kernel,
                             0,
                             sizeof(cl_mem),
                             (void*)&_histograms);
    err_code |= sendArgument(_histograms_scan_kernel,
                             1,
                             sizeof(cl_uint) * maxmemcache,
                             NULL);
    err_code |= sendArgument(_histograms_scan_kernel,
                             2,
                             sizeof(cl_mem),
                             (void*)&_global_sums);
    if(err_code != CL_SUCCESS) {
        S->addMessageF(3, "Can't send one or more arguments to the first pass scan kernel.\n");
        return true;
    }

    #ifdef HAVE_GPUPROFILE
        cl_event event;
        cl_ulong end, start;
        err_code = clEnqueueNDRangeKernel(C->command_queue,
                                          _histograms_scan_kernel,
                                          1,
                                          NULL,
                                          &global_size,
                                          &local_size,
                                          0,
                                          NULL,
                                          &event);
    #else
        err_code = clEnqueueNDRangeKernel(C->command_queue,
                                          _histograms_scan_kernel,
                                          1,
                                          NULL,
                                          &global_size,
                                          &local_size,
                                          0,
                                          NULL,
                                          NULL);
    #endif
    if(err_code != CL_SUCCESS) {
        S->addMessageF(3, "I cannot execute the first pass kernel.\n");
        S->printOpenCLError(err_code);
        return true;
    }
    #ifdef HAVE_GPUPROFILE
        err_code = clWaitForEvents(1, &event);
        if(err_code != CL_SUCCESS) {
            S->addMessage(3, "Impossible to wait for the kernels end.\n");
            S->printOpenCLError(err_code);
            return true;
        }
        err_code |= clGetEventProfilingInfo(event,
                                            CL_PROFILING_COMMAND_END,
                                            sizeof(cl_ulong),
                                            &end,
                                            0);
        if(err_code != CL_SUCCESS) {
            S->addMessage(3, "I cannot profile the kernel execution.\n");
            S->printOpenCLError(err_code);
            return true;
        }
        err_code |= clGetEventProfilingInfo(event,
                                            CL_PROFILING_COMMAND_START,
                                            sizeof(cl_ulong),
                                            &start,
                                            0);
        if(err_code != CL_SUCCESS) {
            S->addMessage(3, "I cannot profile the kernel execution.\n");
            S->printOpenCLError(err_code);
            return true;
        }
        _time += (end - start)/1000.f;  // 10^-3 ms
    #endif

    // 2nd scan
    // ========
    global_size = _histo_split/2;
    local_size = global_size;
    err_code |= sendArgument(_histograms_scan_kernel,
                             0,
                             sizeof(cl_mem),
                             (void*)&_global_sums);
    err_code |= sendArgument(_histograms_scan_kernel,
                             2,
                             sizeof(cl_mem),
                             (void*)&_temp_mem);
    if(err_code != CL_SUCCESS) {
        S->addMessageF(3, "Can't send one or more arguments to the second pass scan kernel.\n");
        return true;
    }

    #ifdef HAVE_GPUPROFILE
        err_code = clEnqueueNDRangeKernel(C->command_queue,
                                          _histograms_scan_kernel,
                                          1,
                                          NULL,
                                          &global_size,
                                          &local_size,
                                          0,
                                          NULL,
                                          &event);
    #else
        err_code = clEnqueueNDRangeKernel(C->command_queue,
                                          _histograms_scan_kernel,
                                          1,
                                          NULL,
                                          &global_size,
                                          &local_size,
                                          0,
                                          NULL,
                                          NULL);
    #endif
    if(err_code != CL_SUCCESS) {
        S->addMessageF(3, "I cannot execute the second pass kernel.\n");
        S->printOpenCLError(err_code);
        return true;
    }
    #ifdef HAVE_GPUPROFILE
        err_code = clWaitForEvents(1, &event);
        if(err_code != CL_SUCCESS) {
            S->addMessage(3, "Impossible to wait for the kernels end.\n");
            S->printOpenCLError(err_code);
            return true;
        }
        err_code |= clGetEventProfilingInfo(event,
                                            CL_PROFILING_COMMAND_END,
                                            sizeof(cl_ulong),
                                            &end,
                                            0);
        if(err_code != CL_SUCCESS) {
            S->addMessage(3, "I cannot profile the kernel execution.\n");
            S->printOpenCLError(err_code);
            return true;
        }
        err_code |= clGetEventProfilingInfo(event,
                                            CL_PROFILING_COMMAND_START,
                                            sizeof(cl_ulong),
                                            &start,
                                            0);
        if(err_code != CL_SUCCESS) {
            S->addMessage(3, "I cannot profile the kernel execution.\n");
            S->printOpenCLError(err_code);
            return true;
        }
        _time += (end - start)/1000.f;  // 10^-3 ms
    #endif

    // Histograms paste
    // ================
    global_size = _radix*_groups*_items/2;
    local_size = global_size / _histo_split;

    #ifdef HAVE_GPUPROFILE
        err_code = clEnqueueNDRangeKernel(C->command_queue,
                                          _paste_histograms_kernel,
                                          1,
                                          NULL,
                                          &global_size,
                                          &local_size,
                                          0,
                                          NULL,
                                          &event);
    #else
        err_code = clEnqueueNDRangeKernel(C->command_queue,
                                          _paste_histograms_kernel,
                                          1,
                                          NULL,
                                          &global_size,
                                          &local_size,
                                          0,
                                          NULL,
                                          NULL);
    #endif
    if(err_code != CL_SUCCESS) {
        S->addMessageF(3, "I cannot execute pasting histograms kernel.\n");
        S->printOpenCLError(err_code);
        return true;
    }
    #ifdef HAVE_GPUPROFILE
        err_code = clWaitForEvents(1, &event);
        if(err_code != CL_SUCCESS) {
            S->addMessage(3, "Impossible to wait for the kernels end.\n");
            S->printOpenCLError(err_code);
            return true;
        }
        err_code |= clGetEventProfilingInfo(event,
                                            CL_PROFILING_COMMAND_END,
                                            sizeof(cl_ulong),
                                            &end,
                                            0);
        if(err_code != CL_SUCCESS) {
            S->addMessage(3, "I cannot profile the kernel execution.\n");
            S->printOpenCLError(err_code);
            return true;
        }
        err_code |= clGetEventProfilingInfo(event,
                                            CL_PROFILING_COMMAND_START,
                                            sizeof(cl_ulong),
                                            &start,
                                            0);
        if(err_code != CL_SUCCESS) {
            S->addMessage(3, "I cannot profile the kernel execution.\n");
            S->printOpenCLError(err_code);
            return true;
        }
        _time += (end - start)/1000.f;  // 10^-3 ms
    #endif
    return false;
}

bool RadixSort::reorder()
{
    InputOutput::ScreenManager *S = InputOutput::ScreenManager::singleton();
    CalcServer *C = CalcServer::singleton();
    int err_code=0;
    size_t local_size =  _items;
    size_t global_size = _groups*_items;

    err_code |= sendArgument(_sort_kernel,
                             0,
                             sizeof(cl_mem),
                             (void*)&_in_keys);
    err_code |= sendArgument(_sort_kernel,
                             1,
                             sizeof(cl_mem),
                             (void*)&_out_keys);
    err_code |= sendArgument(_sort_kernel,
                             3,
                             sizeof(cl_uint),
                             (void*)&_pass);
    err_code |= sendArgument(_sort_kernel,
                             4,
                             sizeof(cl_mem),
                             (void*)&_in_permut);
    err_code |= sendArgument(_sort_kernel,
                             5,
                             sizeof(cl_mem),
                             (void*)&_out_permut);
    err_code |= sendArgument(_sort_kernel,
                             6,
                             sizeof(cl_uint)*_radix*_items,
                             NULL);
    if(err_code != CL_SUCCESS) {
        S->addMessageF(3, "I cannot send a variable to the kernel.\n");
        return true;
    }

    #ifdef HAVE_GPUPROFILE
        cl_event event;
        cl_ulong end, start;
        err_code = clEnqueueNDRangeKernel(C->command_queue,
                                          _sort_kernel,
                                          1,
                                          NULL,
                                          &global_size,
                                          &local_size,
                                          0,
                                          NULL,
                                          &event);
    #else
        err_code = clEnqueueNDRangeKernel(C->command_queue,
                                          _sort_kernel,
                                          1,
                                          NULL,
                                          &global_size,
                                          &local_size,
                                          0,
                                          NULL,
                                          NULL);
    #endif
    if(err_code != CL_SUCCESS) {
        S->addMessageF(3, "I cannot execute the kernel.\n");
        S->printOpenCLError(err_code);
        return true;
    }
    #ifdef HAVE_GPUPROFILE
        err_code = clWaitForEvents(1, &event);
        if(err_code != CL_SUCCESS) {
            S->addMessage(3, "Impossible to wait for the kernels end.\n");
            S->printOpenCLError(err_code);
            return true;
        }
        err_code |= clGetEventProfilingInfo(event,
                                            CL_PROFILING_COMMAND_END,
                                            sizeof(cl_ulong),
                                            &end,
                                            0);
        if(err_code != CL_SUCCESS) {
            S->addMessage(3, "I cannot profile the kernel execution.\n");
            S->printOpenCLError(err_code);
            return true;
        }
        err_code |= clGetEventProfilingInfo(event,
                                            CL_PROFILING_COMMAND_START,
                                            sizeof(cl_ulong),
                                            &start,
                                            0);
        if(err_code != CL_SUCCESS) {
            S->addMessage(3, "I cannot profile the kernel execution.\n");
            S->printOpenCLError(err_code);
            return true;
        }
        _time += (end - start)/1000.f;  // 10^-3 ms
    #endif

    cl_mem d_temp;

    d_temp = _in_keys;
    _in_keys = _out_keys;
    _out_keys = d_temp;

    d_temp = _in_permut;
    _in_permut = _out_permut;
    _out_permut = d_temp;

    return false;
}

bool RadixSort::reversePermutations()
{
    InputOutput::ScreenManager *S = InputOutput::ScreenManager::singleton();
    CalcServer *C = CalcServer::singleton();
    int err_code=0;

    err_code = sendArgument(_inv_permutations_kernel,
                            0,
                            sizeof(cl_mem),
                            (void*)&_in_permut);
    if(err_code != CL_SUCCESS) {
        S->addMessageF(3, "I cannot send a variable to the kernel.\n");
        return true;
    }

    #ifdef HAVE_GPUPROFILE
        cl_event event;
        cl_ulong end, start;
        err_code = clEnqueueNDRangeKernel(C->command_queue,
                                          _inv_permutations_kernel,
                                          1,
                                          NULL,
                                          &_global_work_size,
                                          NULL,
                                          0,
                                          NULL,
                                          &event);
    #else
        err_code = clEnqueueNDRangeKernel(C->command_queue,
                                          _inv_permutations_kernel,
                                          1,
                                          NULL,
                                          &_global_work_size,
                                          NULL,
                                          0,
                                          NULL,
                                          NULL);
    #endif
    if(err_code != CL_SUCCESS) {
        S->addMessageF(3, "I cannot execute the kernel.\n");
        S->printOpenCLError(err_code);
        return true;
    }
    #ifdef HAVE_GPUPROFILE
        err_code = clWaitForEvents(1, &event);
        if(err_code != CL_SUCCESS) {
            S->addMessage(3, "Impossible to wait for the kernels end.\n");
            S->printOpenCLError(err_code);
            return true;
        }
        err_code |= clGetEventProfilingInfo(event,
                                            CL_PROFILING_COMMAND_END,
                                            sizeof(cl_ulong),
                                            &end,
                                            0);
        if(err_code != CL_SUCCESS) {
            S->addMessage(3, "I cannot profile the kernel execution.\n");
            S->printOpenCLError(err_code);
            return true;
        }
        err_code |= clGetEventProfilingInfo(event,
                                            CL_PROFILING_COMMAND_START,
                                            sizeof(cl_ulong),
                                            &start,
                                            0);
        if(err_code != CL_SUCCESS) {
            S->addMessage(3, "I cannot profile the kernel execution.\n");
            S->printOpenCLError(err_code);
            return true;
        }
        _time += (end - start)/1000.f;  // 10^-3 ms
    #endif
    return false;
}

bool RadixSort::setN(unsigned int N)
{
    InputOutput::ScreenManager *S = InputOutput::ScreenManager::singleton();
    CalcServer *C = CalcServer::singleton();
    char msg[1024];
    int err_code = 0;
    if(_n == N)
        return false;
    if(N <= 0){
        sprintf(msg, "No data to sort (n=%u)\n", N);
        S->addMessageF(3, msg);
        return true;
    }
    _local_work_size  = getLocalWorkSize(N, C->command_queue);
    if(!_local_work_size){
        sprintf(msg, "There are not a valid local work size for %u cells.\n", N);
        S->addMessageF(3, msg);
        return true;
    }
    _global_work_size = getGlobalWorkSize(N, _local_work_size);

    // Analise the amount of data
    unsigned int oldN = _n;
    _n = N;
    if(_n % (_items*_groups)){
        sprintf(msg, "n=%u is not divisible by _items*_groups (%u,%u).\n", _n, _items,_groups);
        S->addMessageF(3, msg);
        return true;
    }
    if(!isPowerOf2(_n)){
        sprintf(msg, "n=%u, which is not power of 2.\n", _n);
        S->addMessageF(3, msg);
        return true;
    }
    // Some checks
    if((_groups * _items * _radix) % _histo_split){
        S->addMessageF(3, "_groups * _items * _radix must be divisible by _histo_split\n");
        S->addMessage(0, "\tUnhandled exception\n");
        return true;
    }
    if(!isPowerOf2(_groups)){
        S->addMessageF(3, "_groups is not power of two\n");
        S->addMessage(0, "\tUnhandled exception\n");
        return true;
    }
    if(!isPowerOf2(_items)){
        S->addMessageF(3, "_items is not power of two\n");
        S->addMessage(0, "\tUnhandled exception\n");
        return true;
    }

    // Get the memory identifiers
    _in_keys = C->icell;
    _out_permut = C->permutation;

    // Keys and permutations
    if(_out_keys){
        C->allocated_mem -= oldN*sizeof( unsigned int );
        clReleaseMemObject(_out_keys);
        _out_keys=0;
    }
    _out_keys = C->allocMemory(_n*sizeof( unsigned int ));
    if(!_out_keys)
        return true;
    if(_in_permut){
        C->allocated_mem -= oldN*sizeof( unsigned int );
        clReleaseMemObject(_in_permut);
        _in_permut=0;
    }
    _in_permut = C->allocMemory(_n*sizeof( unsigned int ));
    if(!_in_permut)
        return true;
    C->permutation = _in_permut;
    // Auxiliar data
    if(_histograms){
        clReleaseMemObject(_histograms);
        _histograms=0;
        C->allocated_mem -= (_radix * _groups * _items)*sizeof( unsigned int );
    }
    _histograms = C->allocMemory((_radix * _groups * _items)*sizeof( unsigned int ));
    if(!_histograms)
        return true;
    if(_global_sums){
        clReleaseMemObject(_global_sums);
        _global_sums=0;
        C->allocated_mem -= (_histo_split)*sizeof( unsigned int );
    }
    _global_sums = C->allocMemory(_histo_split*sizeof( unsigned int ));
    if(!_global_sums)
        return true;
    if(_temp_mem){
        clReleaseMemObject(_temp_mem);
        _temp_mem=0;
        C->allocated_mem -= sizeof( unsigned int );
    }
    _temp_mem = C->allocMemory(sizeof( unsigned int ));
    if(!_temp_mem)
        return true;
    sprintf(msg, "Radix sort components = %u\n", _n);
    S->addMessageF(1, msg);
    sprintf(msg, "\tAllocated memory = %lu bytes\n", C->allocated_mem);
    S->addMessage(0, msg);
    return 0;
}

bool RadixSort::setupOpenCL()
{
    InputOutput::ScreenManager *S = InputOutput::ScreenManager::singleton();
    CalcServer *C = CalcServer::singleton();
    char CFlags[256];
    // sprintf(CFlags, "-D_BITS=%u -D_RADIX=%u -DTRANSPOSE -DPERMUT", _bits, _radix);
    sprintf(CFlags, "-D_BITS=%u -D_RADIX=%u -DPERMUT", _bits, _radix);
    char msg[1024];
    int err_code = 0;

    // Clean up all the previously generated kernels
    // =============================================
    if(_init_kernel)
        clReleaseKernel(_init_kernel);
    _init_kernel=0;
    if(_transpose_kernel)
        clReleaseKernel(_transpose_kernel);
    _transpose_kernel=0;
    if(_histograms_kernel)
        clReleaseKernel(_histograms_kernel);
    _histograms_kernel=0;
    if(_histograms_scan_kernel)
        clReleaseKernel(_histograms_scan_kernel);
    _histograms_scan_kernel=0;
    if(_paste_histograms_kernel)
        clReleaseKernel(_paste_histograms_kernel);
    _paste_histograms_kernel=0;
    if(_sort_kernel)
        clReleaseKernel(_sort_kernel);
    _sort_kernel=0;
    if(_inv_permutations_kernel)
        clReleaseKernel(_inv_permutations_kernel);
    _inv_permutations_kernel=0;
    if(_program)
        clReleaseProgram(_program);
    _program=0;

    // Load the new kernels
    if(!loadKernelFromFile(&_init_kernel,
                           &_program,
                           C->context,
                           C->device,
                           _path,
                           "init",
                           CFlags))
        return true;
    clReleaseProgram(_program);_program=0;
    if(!loadKernelFromFile(&_transpose_kernel,
                           &_program,
                           C->context,
                           C->device,
                           _path,
                           "transpose",
                           CFlags))
        return true;
    clReleaseProgram(_program);_program=0;
    if(!loadKernelFromFile(&_histograms_kernel,
                           &_program,
                           C->context,
                           C->device,
                           _path,
                           "histogram",
                           CFlags))
        return true;
    clReleaseProgram(_program);_program=0;
    if(!loadKernelFromFile(&_histograms_scan_kernel,
                           &_program,
                           C->context,
                           C->device,
                           _path,
                           "scanhistograms",
                           CFlags))
        return true;
    clReleaseProgram(_program);_program=0;
    if(!loadKernelFromFile(&_paste_histograms_kernel,
                           &_program,
                           C->context,
                           C->device,
                           _path,
                           "pastehistograms",
                           CFlags))
        return true;
    clReleaseProgram(_program);_program=0;
    if(!loadKernelFromFile(&_sort_kernel,
                           &_program,
                           C->context,
                           C->device,
                           _path,
                           "reorder",
                           CFlags))
        return true;
    clReleaseProgram(_program);_program=0;
    if(!loadKernelFromFile(&_inv_permutations_kernel,
                           &_program,
                           C->context,
                           C->device,
                           _path,
                           "reversePermutation",
                           CFlags))
        return true;
    // ------------------------------------------------------------------------
    // Send the fixed arguments
    // ------------------------------------------------------------------------
    err_code |= sendArgument(_init_kernel,
                             1,
                             sizeof(cl_uint),
                          (void*)&_n);
    if(err_code != CL_SUCCESS) {
        S->addMessageF(3, "Failure sending the fixed variables to the initialization kernel.\n");
        return true;
    }
    err_code |= sendArgument(_histograms_kernel,
                             1,
                             sizeof(cl_mem),
                             (void*)&_histograms);
    err_code |= sendArgument(_histograms_kernel,
                             4,
                             sizeof(cl_uint),
                             (void*)&_n);
    if(err_code != CL_SUCCESS) {
        S->addMessageF(3, "Failure sending the fixed variables to the histograms kernel.\n");
        return true;
    }
    err_code |= sendArgument(_paste_histograms_kernel,
                             0,
                             sizeof(cl_mem),
                             (void*)&_histograms);
    err_code |= sendArgument(_paste_histograms_kernel,
                             1,
                             sizeof(cl_mem),
                             (void*)&_global_sums);
    if(err_code != CL_SUCCESS) {
        S->addMessageF(3, "Failure sending the fixed variables to the histograms pasting kernel.\n");
        return true;
    }
    err_code |= sendArgument(_sort_kernel,
                             2,
                             sizeof(cl_mem),
                             (void*)&_histograms);
    err_code |= sendArgument(_sort_kernel,
                             7,
                             sizeof(cl_uint),
                             (void*)&_n);
    if(err_code != CL_SUCCESS) {
        S->addMessageF(3, "Failure sending the fixed variables to the sorting kernel.\n");
        return true;
    }
    err_code |= sendArgument(_inv_permutations_kernel,
                             1,
                             sizeof(cl_mem),
                             (void*)&C->permutation_inverse);
    err_code |= sendArgument(_inv_permutations_kernel,
                             2,
                             sizeof(cl_uint),
                             (void*)&_n);
    if(err_code != CL_SUCCESS) {
        S->addMessageF(3, "Failure sending the fixed variables to the reverse permutations kernel.\n");
        return true;
    }

    // Local work sizes
    // ================
    cl_device_id device;
    err_code = clGetCommandQueueInfo(C->command_queue,
                                      CL_QUEUE_DEVICE,
                                      sizeof(cl_device_id),
                                      &device,
                                      NULL);
    if(err_code != CL_SUCCESS) {
        S->addMessageF(3, "Can't get the device from the command queue.\n");
        S->printOpenCLError(err_code);
        return true;
    }
    cl_uint dims;
    err_code = clGetDeviceInfo(device,
                               CL_DEVICE_MAX_WORK_ITEM_DIMENSIONS,
                               sizeof(cl_uint),
                               &dims,
                               NULL);
    if(err_code != CL_SUCCESS){
        S->addMessageF(3, "Can't get the number of dimensions allowed in the device.\n");
        S->printOpenCLError(err_code);
        return true;
    }
    if(dims < 2){
        S->addMessageF(3, "The device can not execute 2D kernels.\n");
        return true;
    }

    // The transposition process requires __CL_MIN_LOCALSIZE__, otherwise the
    // device cannot be used
    size_t maxLocalSize = 0;
    err_code = clGetKernelWorkGroupInfo(_transpose_kernel,
                                         device,
                                         CL_KERNEL_WORK_GROUP_SIZE,
                                         sizeof(size_t),
                                         &maxLocalSize,
                                         NULL);
    if(err_code != CL_SUCCESS) {
        S->addMessageF(3, "Can't get maximum local work group size for the transposition kernel.\n");
        S->printOpenCLError(err_code);
        return true;
    }
    if(maxLocalSize < __CL_MIN_LOCALSIZE__) {
        S->addMessage(3, "This device is not able to execute the transposition kernel.\n");
        S->addMessage(0, "\tYou can try to recompile the code decreasing __CL_MIN_LOCALSIZE__\n");
        return true;
    }

    // With the _histograms_kernel and _sort_kernel _items can be used as
    // upper bound
    err_code = clGetKernelWorkGroupInfo(_sort_kernel,device,
                                         CL_KERNEL_WORK_GROUP_SIZE,
                                         sizeof(size_t),
                                         &maxLocalSize,
                                         NULL);
    if(err_code != CL_SUCCESS) {
        S->addMessageF(3, "Can't get maximum local work group size for the histograms kernel.\n");
        S->printOpenCLError(err_code);
        return true;
    }
    if(maxLocalSize < _items) _items = maxLocalSize;
    err_code = clGetKernelWorkGroupInfo(_histograms_kernel,
                                        device,
                                        CL_KERNEL_WORK_GROUP_SIZE,
                                        sizeof(size_t),
                                        &maxLocalSize,
                                        NULL);
    if(err_code != CL_SUCCESS) {
        S->addMessageF(3, "Can't get maximum local work group size for the histograms kernel.\n");
        S->printOpenCLError(err_code);
        return true;
    }
    if(maxLocalSize < _items) _items = maxLocalSize;

    // With the scan process we can easily set a bound with the number of
    // histograms splits
    err_code = clGetKernelWorkGroupInfo(_histograms_scan_kernel,
                                        device,
                                        CL_KERNEL_WORK_GROUP_SIZE,
                                        sizeof(size_t),
                                        &maxLocalSize,
                                        NULL);
    if(err_code != CL_SUCCESS) {
        S->addMessageF(3, "Can't get maximum local work group size for the histograms scan kernel.\n");
        S->printOpenCLError(err_code);
        return true;
    }
    if(maxLocalSize < _histo_split/2) _histo_split = 2*maxLocalSize;
    // With the scan histograms kernel and the paste one we must adjust it
    // with the GROUPS, ITEMS and RADIX
    size_t maxForScan = maxLocalSize;
    err_code = clGetKernelWorkGroupInfo(_sort_kernel,
                                        device,
                                        CL_KERNEL_WORK_GROUP_SIZE,
                                        sizeof(size_t),
                                        &maxLocalSize,
                                        NULL);
    if(err_code != CL_SUCCESS) {
        S->addMessageF(3, "Can't get maximum local work group size from the sorting kernel.\n");
        S->printOpenCLError(err_code);
        return true;
    }
    maxLocalSize = max(maxLocalSize, maxForScan);
    while(maxLocalSize < _radix*_groups*_items/2/_histo_split){
        // We can't increase _histo_split, so we start decreasing the number
        // of items
        _items /= 2;
        if(_items < __CL_MIN_LOCALSIZE__){
            _items = __CL_MIN_LOCALSIZE__;
            break;
        }
    }
    while(maxLocalSize < _radix*_groups*_items/2/_histo_split){
        // We have reached the minimum possible value for items, so we can
        // start decreasing the number of groups
        _groups /= 2;
        if(!_groups){
            _groups = 1;
            break;
        }
    }
    if(maxLocalSize < _radix*_groups*_items/2/_histo_split){
        // We can try to reduce the radix, but it is a bad bussiness
        S->addMessageF(3, "Can't be imposed a number of items and groups for this device.\n");
        S->addMessage(0, "\tYou can try to recompile the code decreasing __CL_MIN_LOCALSIZE__\n");
        return true;
    }

    cl_ulong usedMem = 0, availableMem = 0;
    err_code = clGetDeviceInfo(device, CL_DEVICE_LOCAL_MEM_SIZE, sizeof(cl_ulong), &availableMem, NULL);
    if(err_code != CL_SUCCESS) {
        S->addMessageF(3, "Can't get available local memory available on the device.\n");
        S->printOpenCLError(err_code);
        return true;
    }
    // The transpose kernel needs two matrices of tilesize*tilesize
    unsigned int tilesize = __CL_MIN_LOCALSIZE__;
    unsigned int nbcol    = _n/(_groups * _items);
    unsigned int nbrow    = _groups * _items;
    if (nbrow%tilesize != 0) tilesize = 1;
    if (nbcol%tilesize != 0) tilesize = 1;
    err_code = sendArgument(_transpose_kernel,
                             6,
                             sizeof(cl_uint)*tilesize*tilesize,
                             NULL);
    if(err_code != CL_SUCCESS) {
        S->addMessageF(3, "Failure setting the local memory for the transpose kernel.\n");
        S->printOpenCLError(err_code);
        return true;
    }
    err_code = sendArgument(_transpose_kernel,
                             7,
                             sizeof(cl_uint)*tilesize*tilesize,
                             NULL);
    if(err_code != CL_SUCCESS) {
        S->addMessageF(3, "Failure setting the local memory for the transpose kernel.\n");
        S->printOpenCLError(err_code);
        return true;
    }
    err_code = clGetKernelWorkGroupInfo(_transpose_kernel,
                                        device,
                                        CL_KERNEL_LOCAL_MEM_SIZE,
                                        sizeof(size_t),
                                        &usedMem,
                                        NULL);
    if(err_code != CL_SUCCESS) {
        S->addMessageF(3, "Impossible to get the local memory used by the transpose kernel.\n");
        S->printOpenCLError(err_code);
        return true;
    }
    if(availableMem < usedMem) {
        S->addMessageF(3, "This device is not able to execute the transposition kernel.\n");
        S->addMessage(0, "\tYou can try to recompile the code decreasing __CL_MIN_LOCALSIZE__\n");
        return true;
    }
    // The histograms and sorting kernels requires an array of radix*items
    err_code = sendArgument(_histograms_kernel,
                            3,
                            sizeof(cl_uint)*_radix*_items,
                            NULL);
    if(err_code != CL_SUCCESS) {
        S->addMessage(3, "Failure setting the local memory for the histograms kernel.\n");
        S->printOpenCLError(err_code);
        return true;
    }
    err_code = sendArgument(_sort_kernel,
                            6,
                            sizeof(cl_uint)*_radix*_items,
                            NULL);
    if(err_code != CL_SUCCESS) {
        S->addMessage(3, "Failure setting the local memory for the sorting kernel.\n");
        S->printOpenCLError(err_code);
        return true;
    }
    cl_ulong histoUsed = 0;
    err_code = clGetKernelWorkGroupInfo(_histograms_kernel,
                                        device,
                                        CL_KERNEL_LOCAL_MEM_SIZE,
                                        sizeof(size_t),
                                        &histoUsed,
                                        NULL);
    if(err_code != CL_SUCCESS) {
        S->addMessageF(3, "Impossible to get the local memory used by the histogram kernel.\n");
        S->printOpenCLError(err_code);
        return true;
    }
    err_code = clGetKernelWorkGroupInfo(_sort_kernel,
                                        device,
                                        CL_KERNEL_LOCAL_MEM_SIZE,
                                        sizeof(size_t),
                                        &usedMem,
                                        NULL);
    if(err_code != CL_SUCCESS) {
        S->addMessageF(3, "Impossible to get the local memory used by the sorting kernel.\n");
        S->printOpenCLError(err_code);
        return true;
    }
    usedMem = max(usedMem, histoUsed);
    while(availableMem < usedMem){
        // We can try to decrease the number of items
        _items /= 2;
        if(_items < __CL_MIN_LOCALSIZE__){
            _items = __CL_MIN_LOCALSIZE__;
            break;
        }
        err_code = sendArgument(_histograms_kernel,
                                3,
                                sizeof(cl_uint)*_radix*_items,
                                NULL);
        if(err_code != CL_SUCCESS) {
            S->addMessage(3, "Failure setting the local memory for the histogram kernel.\n");
            S->printOpenCLError(err_code);
            return true;
        }
        err_code = sendArgument(_sort_kernel,
                                6,
                                sizeof(cl_uint)*_radix*_items,
                                NULL);
        if(err_code != CL_SUCCESS) {
            S->addMessage(3, "Failure setting the local memory for the sorting kernel.\n");
            S->printOpenCLError(err_code);
            return true;
        }
        err_code = clGetKernelWorkGroupInfo(_histograms_kernel,
                                            device,
                                            CL_KERNEL_LOCAL_MEM_SIZE,
                                            sizeof(size_t),
                                            &histoUsed,
                                            NULL);
        if(err_code != CL_SUCCESS) {
            S->addMessageF(3, "Impossible to get the local memory used by the histogram kernel.\n");
            S->printOpenCLError(err_code);
            return true;
        }
        err_code = clGetKernelWorkGroupInfo(_sort_kernel,
                                            device,
                                            CL_KERNEL_LOCAL_MEM_SIZE,
                                            sizeof(size_t),
                                            &usedMem,
                                            NULL);
        if(err_code != CL_SUCCESS) {
            S->addMessageF(3, "Impossible to get the local memory used by the sorting kernel.\n");
            S->printOpenCLError(err_code);
            return true;
        }
        usedMem = max(usedMem, histoUsed);
    }
    if(availableMem < usedMem){
        // We can try to reduce the radix, but it is a bad bussiness
        S->addMessageF(3, "The device has not local memory enough for the histograms or sorting kernels.\n");
        S->addMessage(0, "\tYou can try to recompile the code decreasing __CL_MIN_LOCALSIZE__\n");
        return true;
    }
    // The scan steps requires an amount of memory of max(_histo_split,_items * _groups * _radix / _histo_split);
    unsigned int maxmemcache = max(_histo_split,
                                   _items*_groups*_radix / _histo_split);
    err_code = sendArgument(_histograms_scan_kernel,
                            1,
                            sizeof(cl_uint) * maxmemcache,
                            NULL);
    if(err_code != CL_SUCCESS) {
        S->addMessage(3, "Failure setting the local memory for the scanning kernel.\n");
        S->printOpenCLError(err_code);
        return true;
    }
    err_code = clGetKernelWorkGroupInfo(_histograms_scan_kernel,
                                        device,
                                        CL_KERNEL_LOCAL_MEM_SIZE,
                                        sizeof(size_t),
                                        &usedMem,
                                        NULL);
    if(err_code != CL_SUCCESS) {
        S->addMessageF(3, "Impossible to get the local memory used by the scanning kernel.\n");
        S->printOpenCLError(err_code);
        return true;
    }
    if(availableMem < usedMem){
        // We can try to decrease the number of splits
        S->addMessageF(3, "The device has not local memory enough for the scan kernel.\n");
        S->addMessage(0, "\tYou can try to recompile the code decreasing _histo_split\n");
        return true;
    }
    // It may happens that with the corrected data the local sizes are lower than the minimum local size
    if(    (_items < __CL_MIN_LOCALSIZE__)
        || (_histo_split/2 < __CL_MIN_LOCALSIZE__)
        || (_radix*_groups*_items/2/_histo_split < __CL_MIN_LOCALSIZE__) )
    {
        S->addMessageF(3, "I cannot find a valid set of values for this device.\n");
        S->addMessage(0, "\tYou can try to recompile the code decreasing __CL_MIN_LOCALSIZE__\n");
        return true;
    }
    // ------------------------------------------------------------------------
    // Report
    // ------------------------------------------------------------------------
    S->addMessageF(1, "OpenCL tools already built.\n");
    sprintf(msg,"\tITEMS = %u\n", _items);
    S->addMessage(0, msg);
    sprintf(msg,"\tGROUPS = %u\n", _groups);
    S->addMessage(0, msg);
    sprintf(msg,"\tBITS = %u\n", _bits);
    S->addMessage(0, msg);
    sprintf(msg,"\tRADIX = %u\n", _radix);
    S->addMessage(0, msg);
    sprintf(msg,"\tHISTOSPLIT = %u\n", _histo_split);
    S->addMessage(0, msg);
    return false;
}

}}  // namespace
