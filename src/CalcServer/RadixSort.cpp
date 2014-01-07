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

// ----------------------------------------------------------------------------
// Include the Problem setup manager header
// ----------------------------------------------------------------------------
#include <ProblemSetup.h>

// ----------------------------------------------------------------------------
// Include the Problem setup manager header
// ----------------------------------------------------------------------------
#include <ScreenManager.h>

// ----------------------------------------------------------------------------
// Include the main header
// ----------------------------------------------------------------------------
#include <CalcServer/RadixSort.h>

// ----------------------------------------------------------------------------
// Include the calculation server
// ----------------------------------------------------------------------------
#include <CalcServer.h>

namespace Aqua{ namespace CalcServer{

RadixSort::RadixSort()
	: Path(0)
	, n(0)
	, clInKeys(0)
	, clOutKeys(0)
	, clInPermut(0)
	, clOutPermut(0)
	, clHistograms(0)
	, clGlobalSums(0)
	, clTempMem(0)
	, ckInit(0)
	, ckTranspose(0)
	, ckHistogram(0)
	, ckScanHistogram(0)
	, ckPasteHistogram(0)
	, ckReorder(0)
	, ckReversePermutations(0)
	, _program(0)
	#ifdef HAVE_GPUPROFILE
	    , mTime(0.f)
	#endif
	, mItems(_ITEMS)
	, mGroups(_GROUPS)
	, mBits(_STEPBITS)
	, mRadix(_RADIX)
	, mHistoSplit(_HISTOSPLIT)
	, tilesizeWarning(false)
{
	InputOutput::ScreenManager *S = InputOutput::ScreenManager::singleton();
	S->addMessage(1, "INFO (RadixSort::Init): Initializating radix sort tool!\n");
	InputOutput::ProblemSetup *P = InputOutput::ProblemSetup::singleton();
	CalcServer *C = CalcServer::singleton();
	int str_len = 0;
	//! 1st.- Get data
	str_len = strlen(P->OpenCL_kernels.radix_sort);
	if(str_len <= 0) {
	    S->addMessage(3, "(RadixSort::Init): Path of predictor kernel is empty.\n");
	    exit(EXIT_FAILURE);
	}
	Path = new char[str_len+4];
	if(!Path) {
	    S->addMessage(3, "(RadixSort::Init): Memory cannot be allocated for the path.\n");
	    exit(EXIT_FAILURE);
	}
	strcpy(Path, P->OpenCL_kernels.radix_sort);
	strcat(Path, ".cl");
	//! 2nd.- Set the number of elements
	if(setN(C->num_icell))
	    exit(EXIT_FAILURE);
	//! 3rd.- Setup OpenCL
	if(setupOpenCL())
	    exit(EXIT_FAILURE);
	S->addMessage(1, "(RadixSort::Init): Radix sort ready to work!\n");
}

RadixSort::~RadixSort()
{
	if(Path)delete[] Path;Path=0;
	clInKeys = 0;                      // Let CalcServer destroy it
	clInPermut = 0;                    // Let CalcServer destroy it
	if(clOutKeys)clReleaseMemObject(clOutKeys);clOutKeys=0;
	if(clOutPermut)clReleaseMemObject(clOutPermut);clOutPermut=0;
	if(clHistograms)clReleaseMemObject(clHistograms);clHistograms=0;
	if(clGlobalSums)clReleaseMemObject(clGlobalSums);clGlobalSums=0;
	if(clTempMem)clReleaseMemObject(clTempMem);clTempMem=0;
	if(ckInit)clReleaseKernel(ckInit);ckInit=0;
	if(ckTranspose)clReleaseKernel(ckTranspose);ckTranspose=0;
	if(ckHistogram)clReleaseKernel(ckHistogram);ckHistogram=0;
	if(ckScanHistogram)clReleaseKernel(ckScanHistogram);ckScanHistogram=0;
	if(ckPasteHistogram)clReleaseKernel(ckPasteHistogram);ckPasteHistogram=0;
	if(ckReorder)clReleaseKernel(ckReorder);ckReorder=0;
	if(ckReversePermutations)clReleaseKernel(ckReversePermutations);ckReversePermutations=0;
	if(_program)clReleaseProgram(_program);_program=0;
}

bool RadixSort::sort()
{
	InputOutput::ScreenManager *S = InputOutput::ScreenManager::singleton();
	CalcServer *C = CalcServer::singleton();
	#ifdef HAVE_GPUPROFILE
	    mTime = 0.f;
	#endif
	// Get maximum key bits, and needed pass
	unsigned int i, maxLcell;
	maxLcell = nextPowerOf2(C->num_cells);                   // Any cell value cant be greather than this one (this value must be power of two)
	for(i=0; (maxLcell&1) == 0; maxLcell >>= 1, i++);  // Decompose value into power of two factors
	keyBits = i;
	keyBits = roundUp(keyBits, mBits);
	if(keyBits > __UINTBITS__){
	    S->addMessage(3, "(RadixSort::_scan): Resultant keys overflows unsigned int type.\n");
	    return true;
	}
	nPass = keyBits / mBits;
	unsigned int nbcol=n/(mGroups * mItems);
	unsigned int nbrow=mGroups * mItems;
	if(_init())
	    return true;
	// if(_transpose(nbrow,nbcol))
	//     return true;
	/// Perform radix sort
	for(pass=0;pass<nPass;pass++){
	    // Create histograms
	    if(_histograms())
	        return true;
	    // Scan histograms and paste it
	    if(_scan())
	        return true;
	    // Reorder data
	    if(_reorder())
	        return true;
	}
	/// Recover transposition
	// if(_transpose(nbcol,nbrow))
	//     return true;
	/// Send data to CalcServer
	C->icell = clInKeys;            // Take care, swaped into _reorder
	C->permutation = clInPermut;    // Take care, swaped into _reorder
	if(_reversePermutations())
	    return true;
	return false;
}

bool RadixSort::_init()
{
	InputOutput::ScreenManager *S = InputOutput::ScreenManager::singleton();
	CalcServer *C = CalcServer::singleton();
	int err_code=0;
	//! 1st.- Send arguments
	err_code |= sendArgument(ckInit, 0, sizeof(cl_mem), (void*)&clInPermut);
	if(err_code != CL_SUCCESS) {
		S->addMessage(3, "(RadixSort::_init): I cannot send a variable to the kernel.\n");
	    return true;
	}
	//! 2nd.- Execute the kernel
	#ifdef HAVE_GPUPROFILE
	    cl_event event;
	    cl_ulong end, start;
	    profileTime(0.f);
	    err_code = clEnqueueNDRangeKernel(C->command_queue, ckInit, 1, NULL, &_global_work_size, NULL, 0, NULL, &event);
	#else
	    err_code = clEnqueueNDRangeKernel(C->command_queue, ckInit, 1, NULL, &_global_work_size, NULL, 0, NULL, NULL);
	#endif
	if(err_code != CL_SUCCESS) {
		S->addMessage(3, "(RadixSort::_init): I cannot execute the kernel.\n");
	    if(err_code == CL_INVALID_WORK_GROUP_SIZE)
	        S->addMessage(0, "\tInvalid local work group size.\n");
	    else if(err_code == CL_OUT_OF_RESOURCES)
	        S->addMessage(0, "\tDevice out of resources.\n");
	    else if(err_code == CL_MEM_OBJECT_ALLOCATION_FAILURE)
	        S->addMessage(0, "\tAllocation error at device.\n");
	    else if(err_code == CL_OUT_OF_HOST_MEMORY)
	        S->addMessage(0, "\tfailure to allocate resources required by the OpenCL implementation on the host.\n");
	    return true;
	}
	#ifdef HAVE_GPUPROFILE
	    err_code = clWaitForEvents(1, &event);
	    if(err_code != CL_SUCCESS) {
	        S->addMessage(3, "(RadixSort::_init): Impossible to wait for the kernels end.\n");
	        return true;
	    }
	    err_code |= clGetEventProfilingInfo(event, CL_PROFILING_COMMAND_END, sizeof(cl_ulong), &end, 0);
	    err_code |= clGetEventProfilingInfo(event, CL_PROFILING_COMMAND_START, sizeof(cl_ulong), &start, 0);
	    if(err_code != CL_SUCCESS) {
	        S->addMessage(3, "(RadixSort::_init): I cannot profile the kernel execution.\n");
	        return true;
	    }
	    mTime += (end - start)/1000.f;  // 10^-3 ms
	#endif
	return false;
}

bool RadixSort::_transpose(unsigned int nbrow, unsigned int nbcol)
{
	InputOutput::ScreenManager *S = InputOutput::ScreenManager::singleton();
	CalcServer *C = CalcServer::singleton();
	int err_code=0;
	unsigned int tilesize=__CL_MIN_LOCALSIZE__;
	if (nbrow%tilesize != 0) tilesize=1;
	if (nbcol%tilesize != 0) tilesize=1;
	if( tilesizeWarning ){
	    tilesize=1;
	}
	else if (tilesize == 1) {
	    S->addMessage(1, "(RadixSort::_transpose): List too small, local memory usage will be avoided.\n");
	    tilesizeWarning=true;
	}
	//! 1st.- Send arguments
	err_code |= sendArgument(ckTranspose, 0, sizeof(cl_mem), (void*)&clInKeys);
	err_code |= sendArgument(ckTranspose, 1, sizeof(cl_mem), (void*)&clOutKeys);
	err_code |= sendArgument(ckTranspose, 2, sizeof(cl_uint), (void*)&nbcol);
	err_code |= sendArgument(ckTranspose, 3, sizeof(cl_uint), (void*)&nbrow);
	err_code |= sendArgument(ckTranspose, 4, sizeof(cl_mem), (void*)&clInPermut);
	err_code |= sendArgument(ckTranspose, 5, sizeof(cl_mem), (void*)&clOutPermut);
	err_code |= sendArgument(ckTranspose, 6, sizeof(cl_uint)*tilesize*tilesize, NULL);
	err_code |= sendArgument(ckTranspose, 7, sizeof(cl_uint)*tilesize*tilesize, NULL);
	err_code |= sendArgument(ckTranspose, 8, sizeof(cl_uint), (void*)&tilesize);
	if(err_code != CL_SUCCESS) {
		S->addMessage(3, "(RadixSort::_transpose): I cannot send a variable to the kernel.\n");
	    return 1;
	}
	//! 2nd.- Execute the kernel
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
	    err_code = clEnqueueNDRangeKernel(C->command_queue, ckTranspose, 2, NULL, _global_work_size, _local_work_size, 0, NULL, &event);
	#else
	    err_code = clEnqueueNDRangeKernel(C->command_queue, ckTranspose, 2, NULL, _global_work_size, _local_work_size, 0, NULL, NULL);
	#endif
	if(err_code != CL_SUCCESS) {
		S->addMessage(3, "(RadixSort::_transpose): I cannot execute the kernel.\n");
	    if(err_code == CL_INVALID_WORK_GROUP_SIZE)
	        S->addMessage(0, "\tInvalid local work group size.\n");
	    else if(err_code == CL_OUT_OF_RESOURCES)
	        S->addMessage(0, "\tDevice out of resources.\n");
	    else if(err_code == CL_MEM_OBJECT_ALLOCATION_FAILURE)
	        S->addMessage(0, "\tAllocation error at device.\n");
	    else if(err_code == CL_OUT_OF_HOST_MEMORY)
	        S->addMessage(0, "\tfailure to allocate resources required by the OpenCL implementation on the host.\n");
	    return true;
	}
	#ifdef HAVE_GPUPROFILE
	    err_code = clWaitForEvents(1, &event);
	    if(err_code != CL_SUCCESS) {
	        S->addMessage(3, "(RadixSort::_transpose): Impossible to wait for the kernels end.\n");
	        return true;
	    }
	    err_code |= clGetEventProfilingInfo(event, CL_PROFILING_COMMAND_END, sizeof(cl_ulong), &end, 0);
	    err_code |= clGetEventProfilingInfo(event, CL_PROFILING_COMMAND_START, sizeof(cl_ulong), &start, 0);
	    if(err_code != CL_SUCCESS) {
	        S->addMessage(3, "(RadixSort::_transpose): I cannot profile the kernel execution.\n");
	        return true;
	    }
	    mTime += (end - start)/1000.f;  // 10^-3 ms
	#endif
	//! 3rd.- Swap input/output arrays
	cl_mem d_temp;
	// keys
	d_temp=clInKeys;
	clInKeys=clOutKeys;
	clOutKeys=d_temp;
	// Permutations
	d_temp=clInPermut;
	clInPermut=clOutPermut;
	clOutPermut=d_temp;
	return false;
}

bool RadixSort::_histograms()
{
	InputOutput::ScreenManager *S = InputOutput::ScreenManager::singleton();
	CalcServer *C = CalcServer::singleton();
	int err_code=0;
	size_t localSize = mItems;
	size_t globalSize = mGroups*mItems;
	//! 1st.- Send arguments
	err_code |= sendArgument(ckHistogram, 0, sizeof(cl_mem), (void*)&clInKeys);
	err_code |= sendArgument(ckHistogram, 2, sizeof(cl_uint), (void*)&pass);
	err_code |= sendArgument(ckHistogram, 3, sizeof(cl_uint)*mRadix*mItems, NULL);
	if(err_code != CL_SUCCESS) {
		S->addMessage(3, "(RadixSort::_histograms): I cannot send a variable to the kernel.\n");
	    return true;
	}
	//! 2nd.- Execute the kernel
	#ifdef HAVE_GPUPROFILE
	    cl_event event;
	    cl_ulong end, start;
	    err_code = clEnqueueNDRangeKernel(C->command_queue, ckHistogram, 1, NULL, &globalSize, &localSize, 0, NULL, &event);
	#else
	    err_code = clEnqueueNDRangeKernel(C->command_queue, ckHistogram, 1, NULL, &globalSize, &localSize, 0, NULL, NULL);
	#endif
	if(err_code != CL_SUCCESS) {
		S->addMessage(3, "(RadixSort::_histograms): I cannot execute the kernel.\n");
	    if(err_code == CL_INVALID_WORK_GROUP_SIZE)
	        S->addMessage(0, "\tInvalid local work group size.\n");
	    else if(err_code == CL_OUT_OF_RESOURCES)
	        S->addMessage(0, "\tDevice out of resources.\n");
	    else if(err_code == CL_MEM_OBJECT_ALLOCATION_FAILURE)
	        S->addMessage(0, "\tAllocation error at device.\n");
	    else if(err_code == CL_OUT_OF_HOST_MEMORY)
	        S->addMessage(0, "\tfailure to allocate resources required by the OpenCL implementation on the host.\n");
	    return true;
	}
	#ifdef HAVE_GPUPROFILE
	    err_code = clWaitForEvents(1, &event);
	    if(err_code != CL_SUCCESS) {
	        S->addMessage(3, "(RadixSort::_histograms): Impossible to wait for the kernels end.\n");
	        return true;
	    }
	    err_code |= clGetEventProfilingInfo(event, CL_PROFILING_COMMAND_END, sizeof(cl_ulong), &end, 0);
	    err_code |= clGetEventProfilingInfo(event, CL_PROFILING_COMMAND_START, sizeof(cl_ulong), &start, 0);
	    if(err_code != CL_SUCCESS) {
	        S->addMessage(3, "(RadixSort::_histograms): I cannot profile the kernel execution.\n");
	        return true;
	    }
	    mTime += (end - start)/1000.f;  // 10^-3 ms
	#endif
	return false;
}

bool RadixSort::_scan()
{
	InputOutput::ScreenManager *S = InputOutput::ScreenManager::singleton();
	CalcServer *C = CalcServer::singleton();
	int err_code=0;
	size_t globalSize = mRadix*mGroups*mItems/2;
	size_t localSize = globalSize / mHistoSplit;
	unsigned int maxmemcache=max(mHistoSplit,mItems * mGroups * mRadix / mHistoSplit);
	//! 1st.- Send arguments to first scan pass
	err_code |= sendArgument(ckScanHistogram, 0, sizeof(cl_mem), (void*)&clHistograms);
	err_code |= sendArgument(ckScanHistogram, 1, sizeof(cl_uint)* maxmemcache, NULL);
	err_code |= sendArgument(ckScanHistogram, 2, sizeof(cl_mem), (void*)&clGlobalSums);
	if(err_code != CL_SUCCESS) {
		S->addMessage(3, "(RadixSort::_scan): Can't send arguments to first pass scan kernel.\n");
	    return true;
	}
	//! 2nd.- Execute the first scan pass kernel
	#ifdef HAVE_GPUPROFILE
	    cl_event event;
	    cl_ulong end, start;
	    err_code = clEnqueueNDRangeKernel(C->command_queue, ckScanHistogram, 1, NULL, &globalSize, &localSize, 0, NULL, &event);
	#else
	    err_code = clEnqueueNDRangeKernel(C->command_queue, ckScanHistogram, 1, NULL, &globalSize, &localSize, 0, NULL, NULL);
	#endif
	if(err_code != CL_SUCCESS) {
		S->addMessage(3, "(RadixSort::_scan): I cannot execute first pass kernel.\n");
	    if(err_code == CL_INVALID_WORK_GROUP_SIZE)
	        S->addMessage(0, "\tInvalid local work group size.\n");
	    else if(err_code == CL_OUT_OF_RESOURCES)
	        S->addMessage(0, "\tDevice out of resources.\n");
	    else if(err_code == CL_MEM_OBJECT_ALLOCATION_FAILURE)
	        S->addMessage(0, "\tAllocation error at device.\n");
	    else if(err_code == CL_OUT_OF_HOST_MEMORY)
	        S->addMessage(0, "\tfailure to allocate resources required by the OpenCL implementation on the host.\n");
	    return true;
	}
	#ifdef HAVE_GPUPROFILE
	    err_code = clWaitForEvents(1, &event);
	    if(err_code != CL_SUCCESS) {
	        S->addMessage(3, "(RadixSort::_scan): Impossible to wait for the first pass kernel end.\n");
	        return true;
	    }
	    err_code |= clGetEventProfilingInfo(event, CL_PROFILING_COMMAND_END, sizeof(cl_ulong), &end, 0);
	    err_code |= clGetEventProfilingInfo(event, CL_PROFILING_COMMAND_START, sizeof(cl_ulong), &start, 0);
	    if(err_code != CL_SUCCESS) {
	        S->addMessage(3, "(RadixSort::_scan): I cannot profile the first pass kernel execution.\n");
	        return true;
	    }
	    mTime += (end - start)/1000.f;  // 10^-3 ms
	#endif
	//! 3rd.- Send arguments to second scan pass
	globalSize = mHistoSplit/2;
	localSize = globalSize;
	err_code |= sendArgument(ckScanHistogram, 0, sizeof(cl_mem), (void*)&clGlobalSums);
	err_code |= sendArgument(ckScanHistogram, 2, sizeof(cl_mem), (void*)&clTempMem);
	if(err_code != CL_SUCCESS) {
		S->addMessage(3, "(RadixSort::_scan): Can't send arguments to second pass scan kernel.\n");
	    return true;
	}
	//! 4th.- Execute the second scan pass kernel
	#ifdef HAVE_GPUPROFILE
	    err_code = clEnqueueNDRangeKernel(C->command_queue, ckScanHistogram, 1, NULL, &globalSize, &localSize, 0, NULL, &event);
	#else
	    err_code = clEnqueueNDRangeKernel(C->command_queue, ckScanHistogram, 1, NULL, &globalSize, &localSize, 0, NULL, NULL);
	#endif
	if(err_code != CL_SUCCESS) {
		S->addMessage(3, "(RadixSort::_scan): I cannot execute second pass kernel.\n");
	    if(err_code == CL_INVALID_WORK_GROUP_SIZE)
	        S->addMessage(0, "\tInvalid local work group size.\n");
	    else if(err_code == CL_OUT_OF_RESOURCES)
	        S->addMessage(0, "\tDevice out of resources.\n");
	    else if(err_code == CL_MEM_OBJECT_ALLOCATION_FAILURE)
	        S->addMessage(0, "\tAllocation error at device.\n");
	    else if(err_code == CL_OUT_OF_HOST_MEMORY)
	        S->addMessage(0, "\tfailure to allocate resources required by the OpenCL implementation on the host.\n");
	    return true;
	}
	#ifdef HAVE_GPUPROFILE
	    err_code = clWaitForEvents(1, &event);
	    if(err_code != CL_SUCCESS) {
	        S->addMessage(3, "(RadixSort::_scan): Impossible to wait for the second pass kernel end.\n");
	        return true;
	    }
	    err_code |= clGetEventProfilingInfo(event, CL_PROFILING_COMMAND_END, sizeof(cl_ulong), &end, 0);
	    err_code |= clGetEventProfilingInfo(event, CL_PROFILING_COMMAND_START, sizeof(cl_ulong), &start, 0);
	    if(err_code != CL_SUCCESS) {
	        S->addMessage(3, "(RadixSort::_scan): I cannot profile the second pass kernel execution.\n");
	        return true;
	    }
	    mTime += (end - start)/1000.f;  // 10^-3 ms
	#endif
	//! 5th.- Send arguments to third scan pass
	globalSize = mRadix*mGroups*mItems/2;
	localSize = globalSize / mHistoSplit;
	//! 6th.- Execute the third scan pass kernel
	#ifdef HAVE_GPUPROFILE
	    err_code = clEnqueueNDRangeKernel(C->command_queue, ckPasteHistogram, 1, NULL, &globalSize, &localSize, 0, NULL, &event);
	#else
	    err_code = clEnqueueNDRangeKernel(C->command_queue, ckPasteHistogram, 1, NULL, &globalSize, &localSize, 0, NULL, NULL);
	#endif
	if(err_code != CL_SUCCESS) {
		S->addMessage(3, "(RadixSort::_scan): I cannot execute paste pass kernel.\n");
	    if(err_code == CL_INVALID_WORK_GROUP_SIZE)
	        S->addMessage(0, "\tInvalid local work group size.\n");
	    else if(err_code == CL_OUT_OF_RESOURCES)
	        S->addMessage(0, "\tDevice out of resources.\n");
	    else if(err_code == CL_MEM_OBJECT_ALLOCATION_FAILURE)
	        S->addMessage(0, "\tAllocation error at device.\n");
	    else if(err_code == CL_OUT_OF_HOST_MEMORY)
	        S->addMessage(0, "\tfailure to allocate resources required by the OpenCL implementation on the host.\n");
	    return true;
	}
	#ifdef HAVE_GPUPROFILE
	    err_code = clWaitForEvents(1, &event);
	    if(err_code != CL_SUCCESS) {
	        S->addMessage(3, "(RadixSort::_scan): Impossible to wait for the paste pass kernel end.\n");
	        return true;
	    }
	    err_code |= clGetEventProfilingInfo(event, CL_PROFILING_COMMAND_END, sizeof(cl_ulong), &end, 0);
	    err_code |= clGetEventProfilingInfo(event, CL_PROFILING_COMMAND_START, sizeof(cl_ulong), &start, 0);
	    if(err_code != CL_SUCCESS) {
	        S->addMessage(3, "(RadixSort::_scan): I cannot profile the paste pass kernel execution.\n");
	        return 4;
	    }
	    mTime += (end - start)/1000.f;  // 10^-3 ms
	#endif
	return false;
}

bool RadixSort::_reorder()
{
	InputOutput::ScreenManager *S = InputOutput::ScreenManager::singleton();
	CalcServer *C = CalcServer::singleton();
	int err_code=0;
	size_t localSize =  mItems;
	size_t globalSize = mGroups*mItems;
	//! 1st.- Send arguments
	err_code |= sendArgument(ckReorder, 0, sizeof(cl_mem), (void*)&clInKeys);
	err_code |= sendArgument(ckReorder, 1, sizeof(cl_mem), (void*)&clOutKeys);
	err_code |= sendArgument(ckReorder, 3, sizeof(cl_uint), (void*)&pass);
	err_code |= sendArgument(ckReorder, 4, sizeof(cl_mem), (void*)&clInPermut);
	err_code |= sendArgument(ckReorder, 5, sizeof(cl_mem), (void*)&clOutPermut);
	err_code |= sendArgument(ckReorder, 6, sizeof(cl_uint)*mRadix*mItems, NULL);
	if(err_code != CL_SUCCESS) {
		S->addMessage(3, (char*)"(RadixSort::_reorder): I cannot send a variable to the kernel.\n");
	    return true;
	}
	//! 2nd.- Execute the kernel
	#ifdef HAVE_GPUPROFILE
	    cl_event event;
	    cl_ulong end, start;
	    err_code = clEnqueueNDRangeKernel(C->command_queue, ckReorder, 1, NULL, &globalSize, &localSize, 0, NULL, &event);
	#else
	    err_code = clEnqueueNDRangeKernel(C->command_queue, ckReorder, 1, NULL, &globalSize, &localSize, 0, NULL, NULL);
	#endif
	if(err_code != CL_SUCCESS) {
		S->addMessage(3, "(RadixSort::_reorder): I cannot execute the kernel.\n");
	    if(err_code == CL_INVALID_PROGRAM_EXECUTABLE)
	        S->addMessage(0, "\tInvalid program (Compile errors maybe?).\n");
	    else if(err_code == CL_INVALID_COMMAND_QUEUE)
	        S->addMessage(0, "\tInvalid command queue.\n");
	    else if(err_code == CL_INVALID_KERNEL)
	        S->addMessage(0, "\tKernel is not a valid object.\n");
	    else if(err_code == CL_INVALID_CONTEXT)
	        S->addMessage(0, "\tContext associated to command queue don't match qith kernel context.\n");
	    else if(err_code == CL_INVALID_KERNEL_ARGS)
	        S->addMessage(0, "\tOne or more arguments are invalid (maybe don't specified).\n");
	    else if(err_code == CL_INVALID_WORK_DIMENSION)
	        S->addMessage(0, "\tDimension must be a value between 1 and 3.\n");
	    else if(err_code == CL_INVALID_WORK_GROUP_SIZE)
	        S->addMessage(0, "\tInvalid local work group size.\n");
	    else if(err_code == CL_INVALID_WORK_ITEM_SIZE)
	        S->addMessage(0, "\tLocal work group size is out of bounds.\n");
	    else if(err_code == CL_INVALID_GLOBAL_OFFSET)
	        S->addMessage(0, "\tGlobal offset must be NULL.\n");
	    else if(err_code == CL_OUT_OF_RESOURCES)
	        S->addMessage(0, "\tDevice out of resources.\n");
	    else if(err_code == CL_MEM_OBJECT_ALLOCATION_FAILURE)
	        S->addMessage(0, "\tAllocation error at device.\n");
	    else if(err_code == CL_INVALID_EVENT_WAIT_LIST)
	        S->addMessage(0, "\tInvalid event wait instruction.\n");
	    else if(err_code == CL_OUT_OF_HOST_MEMORY)
	        S->addMessage(0, "\tfailure to allocate resources required by the OpenCL implementation on the host.\n");
	    return true;
	}
	#ifdef HAVE_GPUPROFILE
	    err_code = clWaitForEvents(1, &event);
	    if(err_code != CL_SUCCESS) {
	        S->addMessage(3, (char*)"(RadixSort::_reorder): Impossible to wait for the kernels end.\n");
	        return true;
	    }
	    err_code |= clGetEventProfilingInfo(event, CL_PROFILING_COMMAND_END, sizeof(cl_ulong), &end, 0);
	    err_code |= clGetEventProfilingInfo(event, CL_PROFILING_COMMAND_START, sizeof(cl_ulong), &start, 0);
	    if(err_code != CL_SUCCESS) {
	        S->addMessage(3, (char*)"(RadixSort::_reorder): I cannot profile the kernel execution.\n");
	        return true;
	    }
	    mTime += (end - start)/1000.f;  // 10^-3 ms
	#endif
	//! 3rd.- Swap input/output arrays
	cl_mem d_temp;
	// keys
	d_temp=clInKeys;
	clInKeys=clOutKeys;
	clOutKeys=d_temp;
	// Permutations
	d_temp=clInPermut;
	clInPermut=clOutPermut;
	clOutPermut=d_temp;
	return false;
}

bool RadixSort::_reversePermutations()
{
	InputOutput::ScreenManager *S = InputOutput::ScreenManager::singleton();
	CalcServer *C = CalcServer::singleton();
	int err_code=0;
	//! 1st.- Send arguments
	err_code |= sendArgument(ckReversePermutations, 0, sizeof(cl_mem), (void*)&clInPermut);
	if(err_code != CL_SUCCESS) {
		S->addMessage(3, "(RadixSort::_reversePermutations): I cannot send a variable to the kernel.\n");
	    return true;
	}
	//! 2nd.- Execute the kernel
	#ifdef HAVE_GPUPROFILE
	    cl_event event;
	    cl_ulong end, start;
	    err_code = clEnqueueNDRangeKernel(C->command_queue, ckReversePermutations, 1, NULL, &_global_work_size, NULL, 0, NULL, &event);
	#else
	    err_code = clEnqueueNDRangeKernel(C->command_queue, ckReversePermutations, 1, NULL, &_global_work_size, NULL, 0, NULL, NULL);
	#endif
	if(err_code != CL_SUCCESS) {
		S->addMessage(3, "(RadixSort::_reversePermutations): I cannot execute the kernel.\n");
	    if(err_code == CL_INVALID_WORK_GROUP_SIZE)
	        S->addMessage(0, "\tInvalid local work group size.\n");
	    else if(err_code == CL_OUT_OF_RESOURCES)
	        S->addMessage(0, "\tDevice out of resources.\n");
	    else if(err_code == CL_MEM_OBJECT_ALLOCATION_FAILURE)
	        S->addMessage(0, "\tAllocation error at device.\n");
	    else if(err_code == CL_OUT_OF_HOST_MEMORY)
	        S->addMessage(0, "\tfailure to allocate resources required by the OpenCL implementation on the host.\n");
	    return true;
	}
	#ifdef HAVE_GPUPROFILE
	    err_code = clWaitForEvents(1, &event);
	    if(err_code != CL_SUCCESS) {
	        S->addMessage(3, "(RadixSort::_reversePermutations): Impossible to wait for the kernels end.\n");
	        return true;
	    }
	    err_code |= clGetEventProfilingInfo(event, CL_PROFILING_COMMAND_END, sizeof(cl_ulong), &end, 0);
	    err_code |= clGetEventProfilingInfo(event, CL_PROFILING_COMMAND_START, sizeof(cl_ulong), &start, 0);
	    if(err_code != CL_SUCCESS) {
	        S->addMessage(3, "(RadixSort::_reversePermutations): I cannot profile the kernel execution.\n");
	        return true;
	    }
	    mTime += (end - start)/1000.f;  // 10^-3 ms
	#endif
	return false;
}

bool RadixSort::setN(unsigned int N)
{
	InputOutput::ScreenManager *S = InputOutput::ScreenManager::singleton();
	CalcServer *C = CalcServer::singleton();
	char msg[1024];
	int err_code = 0;
	if(n == N)
	    return false;
	if(N <= 0){
	    sprintf(msg, "(RadixSort::setN): Any data to sort (n=%u)\n", N);
	    S->addMessage(3, msg);
	    return true;
	}
	_local_work_size  = getLocalWorkSize(N, C->command_queue);
	if(!_local_work_size){
	    sprintf(msg, "(RadixSort::setN): No valid local work size for %u cells.\n", N);
	    S->addMessage(3, msg);
	    return true;
	}
	_global_work_size = getGlobalWorkSize(N, _local_work_size);
	// Analise the amount of data
	unsigned int oldN = n;
	n = N;
	if(n % (mItems*mGroups)){
	    sprintf(msg, "(RadixSort::setN): n=%u is not divisible by mItems*mGroups (%u,%u).\n", n, mItems,mGroups);
	    S->addMessage(3, msg);
	    return true;
	}
    if(!isPowerOf2(n)){
	    sprintf(msg, "(RadixSort::setN): n=%u, that is not power of 2.\n", n);
	    S->addMessage(3, msg);
	    return true;
    }
	// Some checks
	if((mGroups * mItems * mRadix) % mHistoSplit){
	    S->addMessage(3, "(RadixSort::setN): mGroups * mItems * mRadix must be divisible by mHistoSplit\n");
	    S->addMessage(0, "\tUnhandled exception\n");
	    return true;
	}
	if(!isPowerOf2(mGroups)){
	    S->addMessage(3, "(RadixSort::setN): mGroups is not power of two\n");
	    S->addMessage(0, "\tUnhandled exception\n");
	    return true;
	}
	if(!isPowerOf2(mItems)){
	    S->addMessage(3, "(RadixSort::setN): mItems is not power of two\n");
	    S->addMessage(0, "\tUnhandled exception\n");
	    return true;
	}
	// Get the memory identifiers
	clInKeys = C->icell;
	clOutPermut = C->permutation;
	// Keys and permutations
	if(clOutKeys){
	    C->allocated_mem -= oldN*sizeof( unsigned int );
	    clReleaseMemObject(clOutKeys);
	    clOutKeys=0;
	}
	clOutKeys = C->allocMemory(n*sizeof( unsigned int ));
	if(!clOutKeys)
        return true;
	if(clInPermut){
	    C->allocated_mem -= oldN*sizeof( unsigned int );
	    clReleaseMemObject(clInPermut);
	    clInPermut=0;
	}
	clInPermut = C->allocMemory(n*sizeof( unsigned int ));
	if(!clInPermut)
        return true;
	C->permutation = clInPermut;
	// Auxiliar data
	if(clHistograms){
	    clReleaseMemObject(clHistograms);
	    clHistograms=0;
	    C->allocated_mem -= (mRadix * mGroups * mItems)*sizeof( unsigned int );
	}
	clHistograms = C->allocMemory((mRadix * mGroups * mItems)*sizeof( unsigned int ));
	if(!clHistograms)
        return true;
	if(clGlobalSums){
	    clReleaseMemObject(clGlobalSums);
	    clGlobalSums=0;
	    C->allocated_mem -= (mHistoSplit)*sizeof( unsigned int );
	}
	clGlobalSums = C->allocMemory(mHistoSplit*sizeof( unsigned int ));
	if(!clGlobalSums)
        return true;
	if(clTempMem){
	    clReleaseMemObject(clTempMem);
	    clTempMem=0;
	    C->allocated_mem -= sizeof( unsigned int );
	}
	clTempMem = C->allocMemory(sizeof( unsigned int ));
	if(!clTempMem)
        return true;
	sprintf(msg, "Radix sort components = %u\n", n);
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
	// sprintf(CFlags, "-D_BITS=%u -D_RADIX=%u -DTRANSPOSE -DPERMUT", mBits, mRadix);
	sprintf(CFlags, "-D_BITS=%u -D_RADIX=%u -DPERMUT", mBits, mRadix);
	char msg[1024];
	int err_code = 0;
	// ------------------------------------------------------------------------
	// Clean up all the previously generated kernels
	// ------------------------------------------------------------------------
	if(ckInit)clReleaseKernel(ckInit);ckInit=0;
	if(ckTranspose)clReleaseKernel(ckTranspose);ckTranspose=0;
	if(ckHistogram)clReleaseKernel(ckHistogram);ckHistogram=0;
	if(ckScanHistogram)clReleaseKernel(ckScanHistogram);ckScanHistogram=0;
	if(ckPasteHistogram)clReleaseKernel(ckPasteHistogram);ckPasteHistogram=0;
	if(ckReorder)clReleaseKernel(ckReorder);ckReorder=0;
	if(ckReversePermutations)clReleaseKernel(ckReversePermutations);ckReversePermutations=0;
	if(_program)clReleaseProgram(_program);_program=0;
    // Load the new kernels
	if(!loadKernelFromFile(&ckInit, &_program, C->context, C->device, Path, "init", CFlags))
	    return true;
	clReleaseProgram(_program);_program=0;
	if(!loadKernelFromFile(&ckTranspose, &_program, C->context, C->device, Path, "transpose", CFlags))
	    return true;
	clReleaseProgram(_program);_program=0;
	if(!loadKernelFromFile(&ckHistogram, &_program, C->context, C->device, Path, "histogram", CFlags))
	    return true;
	clReleaseProgram(_program);_program=0;
	if(!loadKernelFromFile(&ckScanHistogram, &_program, C->context, C->device, Path, "scanhistograms", CFlags))
	    return true;
	clReleaseProgram(_program);_program=0;
	if(!loadKernelFromFile(&ckPasteHistogram, &_program, C->context, C->device, Path, "pastehistograms", CFlags))
	    return true;
	clReleaseProgram(_program);_program=0;
	if(!loadKernelFromFile(&ckReorder, &_program, C->context, C->device, Path, "reorder", CFlags))
	    return true;
	clReleaseProgram(_program);_program=0;
	if(!loadKernelFromFile(&ckReversePermutations, &_program, C->context, C->device, Path, "reversePermutation", CFlags))
	    return true;
	// ------------------------------------------------------------------------
    // Send the fixed arguments
	// ------------------------------------------------------------------------
	err_code |= sendArgument(ckInit, 1, sizeof(cl_uint), (void*)&n);
	if(err_code != CL_SUCCESS) {
		S->addMessage(3, "(RadixSort::setupOpenCL): Fail sending fix variables to the initialization kernel.\n");
	    return true;
	}
	err_code |= sendArgument(ckHistogram, 1, sizeof(cl_mem), (void*)&clHistograms);
	err_code |= sendArgument(ckHistogram, 4, sizeof(cl_uint), (void*)&n);
	if(err_code != CL_SUCCESS) {
		S->addMessage(3, "(RadixSort::setupOpenCL): Fail sending fix variables to the histograms kernel.\n");
	    return true;
	}
	err_code |= sendArgument(ckPasteHistogram, 0, sizeof(cl_mem), (void*)&clHistograms);
	err_code |= sendArgument(ckPasteHistogram, 1, sizeof(cl_mem), (void*)&clGlobalSums);
	if(err_code != CL_SUCCESS) {
		S->addMessage(3, "(RadixSort::setupOpenCL): Fail sending fix variables to the histograms pasting kernel.\n");
	    return true;
	}
	err_code |= sendArgument(ckReorder, 2, sizeof(cl_mem), (void*)&clHistograms);
	err_code |= sendArgument(ckReorder, 7, sizeof(cl_uint), (void*)&n);
	if(err_code != CL_SUCCESS) {
		S->addMessage(3, "(RadixSort::setupOpenCL): Fail sending fix variables to the sorting kernel.\n");
	    return true;
	}
	err_code |= sendArgument(ckReversePermutations, 1, sizeof(cl_mem), (void*)&C->permutation_inverse);
	err_code |= sendArgument(ckReversePermutations, 2, sizeof(cl_uint), (void*)&n);
	if(err_code != CL_SUCCESS) {
		S->addMessage(3, "(RadixSort::setupOpenCL): Fail sending fix variables to the reverse permutations kernel.\n");
	    return true;
	}
	// ------------------------------------------------------------------------
	// Local work sizes
	// ------------------------------------------------------------------------
    cl_device_id device;
	err_code |= clGetCommandQueueInfo(C->command_queue,CL_QUEUE_DEVICE,
	                                sizeof(cl_device_id),&device, NULL);
	if(err_code != CL_SUCCESS) {
		S->addMessage(3, "(RadixSort::setupOpenCL): Can't get the device from the command queue.\n");
		S->addMessage(0, "\tUnhandled exception\n");
	    return true;
	}
	cl_uint dims;
	err_code |= clGetDeviceInfo(device,CL_DEVICE_MAX_WORK_ITEM_DIMENSIONS,sizeof(cl_uint),&dims, NULL);
	if(err_code != CL_SUCCESS){
		S->addMessage(3, "(RadixSort::setupOpenCL): Can't get the number of dimensions allowed in the device.\n");
		S->addMessage(0, "\tUnhandled exception\n");
	    return true;
	}
	if(dims < 2){
		S->addMessage(3, "(RadixSort::setupOpenCL): The device can not execute 2D kernels.\n");
	    return true;
	}
    // The transposition process requires __CL_MIN_LOCALSIZE__, otherwise the device can't be used
    size_t maxLocalSize = 0;
	err_code |= clGetKernelWorkGroupInfo(ckTranspose,device,CL_KERNEL_WORK_GROUP_SIZE,
	                                   sizeof(size_t), &maxLocalSize, NULL);
	if(err_code != CL_SUCCESS) {
		S->addMessage(3, "(LinkList::setupOpenCL): Can't get maximum local work group size from the transposition kernel.\n");
		S->addMessage(0, "\tUnhandled exception\n");
	    return true;
	}
	if(maxLocalSize < __CL_MIN_LOCALSIZE__) {
		S->addMessage(3, "(LinkList::setupOpenCL): This device is not able to execute the transposition kernel.\n");
		S->addMessage(0, "\tYou can try to recompile the code decreasing __CL_MIN_LOCALSIZE__\n");
	    return true;
    }
    // With the ckHistogram and ckReorder we can assign a maximum bound to mItems
	err_code |= clGetKernelWorkGroupInfo(ckReorder,device,CL_KERNEL_WORK_GROUP_SIZE,
	                                   sizeof(size_t), &maxLocalSize, NULL);
	if(err_code != CL_SUCCESS) {
		S->addMessage(3, "(LinkList::setupOpenCL): Can't get maximum local work group size from the histograms kernel.\n");
	    return true;
	}
	if(maxLocalSize < mItems) mItems = maxLocalSize;
	err_code |= clGetKernelWorkGroupInfo(ckHistogram,device,CL_KERNEL_WORK_GROUP_SIZE,
	                                   sizeof(size_t), &maxLocalSize, NULL);
	if(err_code != CL_SUCCESS) {
		S->addMessage(3, "(LinkList::setupOpenCL): Can't get maximum local work group size from the histograms kernel.\n");
		S->addMessage(0, "\tUnhandled exception\n");
	    return true;
	}
	if(maxLocalSize < mItems) mItems = maxLocalSize;
    // With the scan process we can easily set a bound for the number of histograms splits as well
	err_code |= clGetKernelWorkGroupInfo(ckScanHistogram,device,CL_KERNEL_WORK_GROUP_SIZE,
	                                   sizeof(size_t), &maxLocalSize, NULL);
	if(err_code != CL_SUCCESS) {
		S->addMessage(3, "(LinkList::setupOpenCL): Can't get maximum local work group size from the histograms scan kernel.\n");
		S->addMessage(0, "\tUnhandled exception\n");
	    return true;
	}
	if(maxLocalSize < mHistoSplit/2) mHistoSplit = 2*maxLocalSize;
    // With the scan histograms kernel and the paste them one we must adjust GROUPS, ITEMS and RADIX
    size_t maxForScan = maxLocalSize;
	err_code |= clGetKernelWorkGroupInfo(ckReorder,device,CL_KERNEL_WORK_GROUP_SIZE,
	                                   sizeof(size_t), &maxLocalSize, NULL);
	if(err_code != CL_SUCCESS) {
		S->addMessage(3, "(LinkList::setupOpenCL): Can't get maximum local work group size from the sorting kernel.\n");
		S->addMessage(0, "\tUnhandled exception\n");
	    return true;
	}
    maxLocalSize = max(maxLocalSize, maxForScan);
    while(maxLocalSize < mRadix*mGroups*mItems/2/mHistoSplit){
        // We can't increase mHistoSplit, so we start decreasing the number of items
        mItems /= 2;
        if(mItems < __CL_MIN_LOCALSIZE__){
            mItems = __CL_MIN_LOCALSIZE__;
            break;
        }
    }
    while(maxLocalSize < mRadix*mGroups*mItems/2/mHistoSplit){
        // We have reach the minimum possible value for items, so we can decrease the number of groups
        mGroups /= 2;
        if(!mGroups){
            mGroups = 1;
            break;
        }
    }
    if(maxLocalSize < mRadix*mGroups*mItems/2/mHistoSplit){
        // We can try to reduce the radix, but it is a bad bussiness
		S->addMessage(3, "(LinkList::setupOpenCL): Can't be imposed a number of items and groups for this device.\n");
		S->addMessage(0, "\tYou can try to recompile the code decreasing __CL_MIN_LOCALSIZE__\n");
	    return true;
    }
	// ------------------------------------------------------------------------
	// Local memory
	// ------------------------------------------------------------------------
    cl_ulong usedMem = 0, availableMem = 0;
	err_code |= clGetDeviceInfo(device, CL_DEVICE_LOCAL_MEM_SIZE, sizeof(cl_ulong), &availableMem, NULL);
	if(err_code != CL_SUCCESS) {
		S->addMessage(3, "(RadixSort::setupOpenCL): Can't get available local memory available on the device.\n");
		S->addMessage(0, "\tUnhandled exception\n");
	    return true;
	}
    // The transpose kernel needs two matrices of tilesize*tilesize
	unsigned int tilesize = __CL_MIN_LOCALSIZE__;
	unsigned int nbcol    = n/(mGroups * mItems);
	unsigned int nbrow    = mGroups * mItems;
	if (nbrow%tilesize != 0) tilesize = 1;
	if (nbcol%tilesize != 0) tilesize = 1;
	err_code |= sendArgument(ckTranspose, 6, sizeof(cl_uint)*tilesize*tilesize, NULL);
	err_code |= sendArgument(ckTranspose, 7, sizeof(cl_uint)*tilesize*tilesize, NULL);
	if(err_code != CL_SUCCESS) {
		S->addMessage(3, "(RadixSort::setupOpenCL): Can't set the local memory for the transpose kernel.\n");
		S->addMessage(0, "\tUnhandled exception\n");
	    return true;
	}
	err_code |= clGetKernelWorkGroupInfo(ckTranspose,device,CL_KERNEL_LOCAL_MEM_SIZE,
	                                   sizeof(size_t), &usedMem, NULL);
	if(err_code != CL_SUCCESS) {
		S->addMessage(3, "(RadixSort::setupOpenCL): Can't get the local memory used by the transpose kernel.\n");
		if(err_code == CL_INVALID_DEVICE)
            S->addMessage(0, "\tCL_INVALID_DEVICE\n");
		else if(err_code == CL_INVALID_VALUE)
            S->addMessage(0, "\tCL_INVALID_VALUE\n");
		else if(err_code == CL_INVALID_KERNEL)
            S->addMessage(0, "\tCL_INVALID_KERNEL\n");
        else{
            sprintf(msg, "\tUnhandled exception %i\n", err_code);
            S->addMessage(0, msg);
        }
	    return true;
	}
	if(availableMem < usedMem) {
		S->addMessage(3, "(LinkList::setupOpenCL): This device is not able to execute the transposition kernel.\n");
		S->addMessage(0, "\tYou can try to recompile the code decreasing __CL_MIN_LOCALSIZE__\n");
	    return true;
    }
    // The histograms and sorting kernels requires an array of radix*items
	err_code |= sendArgument(ckHistogram, 3, sizeof(cl_uint)*mRadix*mItems, NULL);
	err_code |= sendArgument(ckReorder, 6, sizeof(cl_uint)*mRadix*mItems, NULL);
	if(err_code != CL_SUCCESS) {
		S->addMessage(3, "(RadixSort::setupOpenCL): Can't set the local memory for histograms and sorting kernels.\n");
		S->addMessage(0, "\tUnhandled exception\n");
	    return true;
	}
	cl_ulong histoUsed = 0;
	err_code |= clGetKernelWorkGroupInfo(ckHistogram,device,CL_KERNEL_LOCAL_MEM_SIZE,
	                                   sizeof(size_t), &histoUsed, NULL);
	err_code |= clGetKernelWorkGroupInfo(ckReorder,device,CL_KERNEL_LOCAL_MEM_SIZE,
	                                   sizeof(size_t), &usedMem, NULL);
	if(err_code != CL_SUCCESS) {
		S->addMessage(3, "(RadixSort::setupOpenCL): Can't get the local memory used by the histograms and sorting kernels.\n");
		S->addMessage(0, "\tUnhandled exception\n");
	    return true;
	}
	usedMem = max(usedMem, histoUsed);
    while(availableMem < usedMem){
        // We can try to decrease the number of items
        mItems /= 2;
        if(mItems < __CL_MIN_LOCALSIZE__){
            mItems = __CL_MIN_LOCALSIZE__;
            break;
        }
        err_code |= sendArgument(ckHistogram, 3, sizeof(cl_uint)*mRadix*mItems, NULL);
        err_code |= sendArgument(ckReorder, 6, sizeof(cl_uint)*mRadix*mItems, NULL);
        if(err_code != CL_SUCCESS) {
            S->addMessage(3, "(RadixSort::setupOpenCL): Can't set the local memory for histograms and sorting kernels.\n");
            S->addMessage(0, "\tUnhandled exception\n");
            return true;
        }
        err_code |= clGetKernelWorkGroupInfo(ckHistogram,device,CL_KERNEL_LOCAL_MEM_SIZE,
                                           sizeof(size_t), &histoUsed, NULL);
        err_code |= clGetKernelWorkGroupInfo(ckReorder,device,CL_KERNEL_LOCAL_MEM_SIZE,
                                           sizeof(size_t), &usedMem, NULL);
        if(err_code != CL_SUCCESS) {
            S->addMessage(3, "(RadixSort::setupOpenCL): Can't get the local memory used by the histograms and sorting kernels.\n");
            S->addMessage(0, "\tUnhandled exception\n");
            return true;
        }
        usedMem = max(usedMem, histoUsed);
    }
    if(availableMem < usedMem){
        // We can try to reduce the radix, but it is a bad bussiness
		S->addMessage(3, "(LinkList::setupOpenCL): The device has not local memory enough for the histograms or sorting kernels.\n");
		S->addMessage(0, "\tYou can try to recompile the code decreasing __CL_MIN_LOCALSIZE__\n");
	    return true;
    }
    // The scan steps requires an amount of memory of max(mHistoSplit,mItems * mGroups * mRadix / mHistoSplit);
    unsigned int maxmemcache=max(mHistoSplit,mItems*mGroups*mRadix / mHistoSplit);
    err_code |= sendArgument(ckScanHistogram, 1, sizeof(cl_uint)* maxmemcache, NULL);
    if(err_code != CL_SUCCESS) {
        S->addMessage(3, "(RadixSort::setupOpenCL): Can't set the local memory for scanning kernel.\n");
        S->addMessage(0, "\tUnhandled exception\n");
        return true;
    }
    err_code |= clGetKernelWorkGroupInfo(ckScanHistogram,device,CL_KERNEL_LOCAL_MEM_SIZE,
                                       sizeof(size_t), &usedMem, NULL);
    if(err_code != CL_SUCCESS) {
        S->addMessage(3, "(RadixSort::setupOpenCL): Can't get the local memory used by the scanning kernel.\n");
        S->addMessage(0, "\tUnhandled exception\n");
        return true;
    }
    if(availableMem < usedMem){
        // We can try to decrease the number of splits
		S->addMessage(3, "(LinkList::setupOpenCL): The device has not local memory enough for the scan kernel.\n");
		S->addMessage(0, "\tYou can try to recompile the code decreasing mHistoSplit\n");
	    return true;
    }
    // It may happens that with the corrected data the local sizes are lower than the minimum local size
	if(    (mItems < __CL_MIN_LOCALSIZE__)
        || (mHistoSplit/2 < __CL_MIN_LOCALSIZE__)
        || (mRadix*mGroups*mItems/2/mHistoSplit < __CL_MIN_LOCALSIZE__) )
    {
		S->addMessage(3, "(LinkList::setupOpenCL): I can't find a valid set of values for this device.\n");
		S->addMessage(0, "\tYou can try to recompile the code decreasing __CL_MIN_LOCALSIZE__\n");
	    return true;
    }
	// ------------------------------------------------------------------------
	// Report
	// ------------------------------------------------------------------------
    S->addMessage(1, "(LinkList::setupOpenCL): OpenCL tools already built.\n");
    sprintf(msg,"\tITEMS = %u\n", mItems);
    S->addMessage(0, msg);
    sprintf(msg,"\tGROUPS = %u\n", mGroups);
    S->addMessage(0, msg);
    sprintf(msg,"\tBITS = %u\n", mBits);
    S->addMessage(0, msg);
    sprintf(msg,"\tRADIX = %u\n", mRadix);
    S->addMessage(0, msg);
    sprintf(msg,"\tHISTOSPLIT = %u\n", mHistoSplit);
    S->addMessage(0, msg);
	return false;
}

}}  // namespace
