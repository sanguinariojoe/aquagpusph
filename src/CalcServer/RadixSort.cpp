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
	, clProgram(0)
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
	int nChar = 0;
	//! 1st.- Get data
	nChar = strlen(P->OpenCL_kernels.radix_sort);
	if(nChar <= 0) {
	    S->addMessage(3, "(RadixSort::Init): Path of predictor kernel is empty.\n");
	    exit(EXIT_FAILURE);
	}
	Path = new char[nChar+4];
	if(!Path) {
	    S->addMessage(3, "(RadixSort::Init): Can't allocate memory for path.\n");
	    exit(EXIT_FAILURE);
	}
	strcpy(Path, P->OpenCL_kernels.radix_sort);
	strcat(Path, ".cl");
	//! 2nd.- Set the number of elements
	if(setN(C->nLcell))
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
	if(clProgram)clReleaseProgram(clProgram);clProgram=0;
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
	maxLcell = nextPowerOf2(C->lxy);                   // Any cell value cant be greather than this one (this value must be power of two)
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
	C->lcell = clInKeys;            // Take care, swaped into _reorder
	C->permutation = clInPermut;    // Take care, swaped into _reorder
	if(_reversePermutations())
	    return true;
	return false;
}

bool RadixSort::_init()
{
	InputOutput::ScreenManager *S = InputOutput::ScreenManager::singleton();
	CalcServer *C = CalcServer::singleton();
	int clFlag=0;
	//! 1st.- Send arguments
	clFlag |= sendArgument(ckInit, 0, sizeof(cl_mem), (void*)&clInPermut);
	if(clFlag != CL_SUCCESS) {
		S->addMessage(3, "(RadixSort::_init): Can't send variable to kernel.\n");
	    return true;
	}
	//! 2nd.- Execute the kernel
	#ifdef HAVE_GPUPROFILE
	    cl_event event;
	    cl_ulong end, start;
	    profileTime(0.f);
	    clFlag = clEnqueueNDRangeKernel(C->clComQueue, ckInit, 1, NULL, &clGlobalWorkSize, NULL, 0, NULL, &event);
	#else
	    clFlag = clEnqueueNDRangeKernel(C->clComQueue, ckInit, 1, NULL, &clGlobalWorkSize, NULL, 0, NULL, NULL);
	#endif
	if(clFlag != CL_SUCCESS) {
		S->addMessage(3, "(RadixSort::_init): Can't execute the kernel.\n");
	    if(clFlag == CL_INVALID_WORK_GROUP_SIZE)
	        S->addMessage(0, "\tInvalid local work group size.\n");
	    else if(clFlag == CL_OUT_OF_RESOURCES)
	        S->addMessage(0, "\tDevice out of resources.\n");
	    else if(clFlag == CL_MEM_OBJECT_ALLOCATION_FAILURE)
	        S->addMessage(0, "\tAllocation error at device.\n");
	    else if(clFlag == CL_OUT_OF_HOST_MEMORY)
	        S->addMessage(0, "\tfailure to allocate resources required by the OpenCL implementation on the host.\n");
	    return true;
	}
	#ifdef HAVE_GPUPROFILE
	    clFlag = clWaitForEvents(1, &event);
	    if(clFlag != CL_SUCCESS) {
	        S->addMessage(3, "(RadixSort::_init): Can't wait to kernels end.\n");
	        return true;
	    }
	    clFlag |= clGetEventProfilingInfo(event, CL_PROFILING_COMMAND_END, sizeof(cl_ulong), &end, 0);
	    clFlag |= clGetEventProfilingInfo(event, CL_PROFILING_COMMAND_START, sizeof(cl_ulong), &start, 0);
	    if(clFlag != CL_SUCCESS) {
	        S->addMessage(3, "(RadixSort::_init): Can't profile kernel execution.\n");
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
	int clFlag=0;
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
	clFlag |= sendArgument(ckTranspose, 0, sizeof(cl_mem), (void*)&clInKeys);
	clFlag |= sendArgument(ckTranspose, 1, sizeof(cl_mem), (void*)&clOutKeys);
	clFlag |= sendArgument(ckTranspose, 2, sizeof(cl_uint), (void*)&nbcol);
	clFlag |= sendArgument(ckTranspose, 3, sizeof(cl_uint), (void*)&nbrow);
	clFlag |= sendArgument(ckTranspose, 4, sizeof(cl_mem), (void*)&clInPermut);
	clFlag |= sendArgument(ckTranspose, 5, sizeof(cl_mem), (void*)&clOutPermut);
	clFlag |= sendArgument(ckTranspose, 6, sizeof(cl_uint)*tilesize*tilesize, NULL);
	clFlag |= sendArgument(ckTranspose, 7, sizeof(cl_uint)*tilesize*tilesize, NULL);
	clFlag |= sendArgument(ckTranspose, 8, sizeof(cl_uint), (void*)&tilesize);
	if(clFlag != CL_SUCCESS) {
		S->addMessage(3, "(RadixSort::_transpose): Can't send variable to kernel.\n");
	    return 1;
	}
	//! 2nd.- Execute the kernel
	size_t global_work_size[2];
	size_t local_work_size[2];
	global_work_size[0]=nbrow/tilesize;
	global_work_size[1]=nbcol;
	local_work_size[0]=1;
	local_work_size[1]=tilesize;
	if(local_work_size[1] < __CL_MIN_LOCALSIZE__){
	    local_work_size[1]  *= __CL_MIN_LOCALSIZE__ / tilesize;
	    global_work_size[1] *= __CL_MIN_LOCALSIZE__ / tilesize;
	}
	#ifdef HAVE_GPUPROFILE
	    cl_event event;
	    cl_ulong end, start;
	    clFlag = clEnqueueNDRangeKernel(C->clComQueue, ckTranspose, 2, NULL, global_work_size, local_work_size, 0, NULL, &event);
	#else
	    clFlag = clEnqueueNDRangeKernel(C->clComQueue, ckTranspose, 2, NULL, global_work_size, local_work_size, 0, NULL, NULL);
	#endif
	if(clFlag != CL_SUCCESS) {
		S->addMessage(3, "(RadixSort::_transpose): Can't execute the kernel.\n");
	    if(clFlag == CL_INVALID_WORK_GROUP_SIZE)
	        S->addMessage(0, "\tInvalid local work group size.\n");
	    else if(clFlag == CL_OUT_OF_RESOURCES)
	        S->addMessage(0, "\tDevice out of resources.\n");
	    else if(clFlag == CL_MEM_OBJECT_ALLOCATION_FAILURE)
	        S->addMessage(0, "\tAllocation error at device.\n");
	    else if(clFlag == CL_OUT_OF_HOST_MEMORY)
	        S->addMessage(0, "\tfailure to allocate resources required by the OpenCL implementation on the host.\n");
	    return true;
	}
	#ifdef HAVE_GPUPROFILE
	    clFlag = clWaitForEvents(1, &event);
	    if(clFlag != CL_SUCCESS) {
	        S->addMessage(3, "(RadixSort::_transpose): Can't wait to kernels end.\n");
	        return true;
	    }
	    clFlag |= clGetEventProfilingInfo(event, CL_PROFILING_COMMAND_END, sizeof(cl_ulong), &end, 0);
	    clFlag |= clGetEventProfilingInfo(event, CL_PROFILING_COMMAND_START, sizeof(cl_ulong), &start, 0);
	    if(clFlag != CL_SUCCESS) {
	        S->addMessage(3, "(RadixSort::_transpose): Can't profile kernel execution.\n");
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
	int clFlag=0;
	size_t localSize = mItems;
	size_t globalSize = mGroups*mItems;
	//! 1st.- Send arguments
	clFlag |= sendArgument(ckHistogram, 0, sizeof(cl_mem), (void*)&clInKeys);
	clFlag |= sendArgument(ckHistogram, 2, sizeof(cl_uint), (void*)&pass);
	clFlag |= sendArgument(ckHistogram, 3, sizeof(cl_uint)*mRadix*mItems, NULL);
	if(clFlag != CL_SUCCESS) {
		S->addMessage(3, "(RadixSort::_histograms): Can't send variable to kernel.\n");
	    return true;
	}
	//! 2nd.- Execute the kernel
	#ifdef HAVE_GPUPROFILE
	    cl_event event;
	    cl_ulong end, start;
	    clFlag = clEnqueueNDRangeKernel(C->clComQueue, ckHistogram, 1, NULL, &globalSize, &localSize, 0, NULL, &event);
	#else
	    clFlag = clEnqueueNDRangeKernel(C->clComQueue, ckHistogram, 1, NULL, &globalSize, &localSize, 0, NULL, NULL);
	#endif
	if(clFlag != CL_SUCCESS) {
		S->addMessage(3, "(RadixSort::_histograms): Can't execute the kernel.\n");
	    if(clFlag == CL_INVALID_WORK_GROUP_SIZE)
	        S->addMessage(0, "\tInvalid local work group size.\n");
	    else if(clFlag == CL_OUT_OF_RESOURCES)
	        S->addMessage(0, "\tDevice out of resources.\n");
	    else if(clFlag == CL_MEM_OBJECT_ALLOCATION_FAILURE)
	        S->addMessage(0, "\tAllocation error at device.\n");
	    else if(clFlag == CL_OUT_OF_HOST_MEMORY)
	        S->addMessage(0, "\tfailure to allocate resources required by the OpenCL implementation on the host.\n");
	    return true;
	}
	#ifdef HAVE_GPUPROFILE
	    clFlag = clWaitForEvents(1, &event);
	    if(clFlag != CL_SUCCESS) {
	        S->addMessage(3, "(RadixSort::_histograms): Can't wait to kernels end.\n");
	        return true;
	    }
	    clFlag |= clGetEventProfilingInfo(event, CL_PROFILING_COMMAND_END, sizeof(cl_ulong), &end, 0);
	    clFlag |= clGetEventProfilingInfo(event, CL_PROFILING_COMMAND_START, sizeof(cl_ulong), &start, 0);
	    if(clFlag != CL_SUCCESS) {
	        S->addMessage(3, "(RadixSort::_histograms): Can't profile kernel execution.\n");
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
	int clFlag=0;
	size_t globalSize = mRadix*mGroups*mItems/2;
	size_t localSize = globalSize / mHistoSplit;
	unsigned int maxmemcache=max(mHistoSplit,mItems * mGroups * mRadix / mHistoSplit);
	//! 1st.- Send arguments to first scan pass
	clFlag |= sendArgument(ckScanHistogram, 0, sizeof(cl_mem), (void*)&clHistograms);
	clFlag |= sendArgument(ckScanHistogram, 1, sizeof(cl_uint)* maxmemcache, NULL);
	clFlag |= sendArgument(ckScanHistogram, 2, sizeof(cl_mem), (void*)&clGlobalSums);
	if(clFlag != CL_SUCCESS) {
		S->addMessage(3, "(RadixSort::_scan): Can't send arguments to first pass scan kernel.\n");
	    return true;
	}
	//! 2nd.- Execute the first scan pass kernel
	#ifdef HAVE_GPUPROFILE
	    cl_event event;
	    cl_ulong end, start;
	    clFlag = clEnqueueNDRangeKernel(C->clComQueue, ckScanHistogram, 1, NULL, &globalSize, &localSize, 0, NULL, &event);
	#else
	    clFlag = clEnqueueNDRangeKernel(C->clComQueue, ckScanHistogram, 1, NULL, &globalSize, &localSize, 0, NULL, NULL);
	#endif
	if(clFlag != CL_SUCCESS) {
		S->addMessage(3, "(RadixSort::_scan): Can't execute first pass kernel.\n");
	    if(clFlag == CL_INVALID_WORK_GROUP_SIZE)
	        S->addMessage(0, "\tInvalid local work group size.\n");
	    else if(clFlag == CL_OUT_OF_RESOURCES)
	        S->addMessage(0, "\tDevice out of resources.\n");
	    else if(clFlag == CL_MEM_OBJECT_ALLOCATION_FAILURE)
	        S->addMessage(0, "\tAllocation error at device.\n");
	    else if(clFlag == CL_OUT_OF_HOST_MEMORY)
	        S->addMessage(0, "\tfailure to allocate resources required by the OpenCL implementation on the host.\n");
	    return true;
	}
	#ifdef HAVE_GPUPROFILE
	    clFlag = clWaitForEvents(1, &event);
	    if(clFlag != CL_SUCCESS) {
	        S->addMessage(3, "(RadixSort::_scan): Can't wait to first pass kernel end.\n");
	        return true;
	    }
	    clFlag |= clGetEventProfilingInfo(event, CL_PROFILING_COMMAND_END, sizeof(cl_ulong), &end, 0);
	    clFlag |= clGetEventProfilingInfo(event, CL_PROFILING_COMMAND_START, sizeof(cl_ulong), &start, 0);
	    if(clFlag != CL_SUCCESS) {
	        S->addMessage(3, "(RadixSort::_scan): Can't profile first pass kernel execution.\n");
	        return true;
	    }
	    mTime += (end - start)/1000.f;  // 10^-3 ms
	#endif
	//! 3rd.- Send arguments to second scan pass
	globalSize = mHistoSplit/2;
	localSize = globalSize;
	clFlag |= sendArgument(ckScanHistogram, 0, sizeof(cl_mem), (void*)&clGlobalSums);
	clFlag |= sendArgument(ckScanHistogram, 2, sizeof(cl_mem), (void*)&clTempMem);
	if(clFlag != CL_SUCCESS) {
		S->addMessage(3, "(RadixSort::_scan): Can't send arguments to second pass scan kernel.\n");
	    return true;
	}
	//! 4th.- Execute the second scan pass kernel
	#ifdef HAVE_GPUPROFILE
	    clFlag = clEnqueueNDRangeKernel(C->clComQueue, ckScanHistogram, 1, NULL, &globalSize, &localSize, 0, NULL, &event);
	#else
	    clFlag = clEnqueueNDRangeKernel(C->clComQueue, ckScanHistogram, 1, NULL, &globalSize, &localSize, 0, NULL, NULL);
	#endif
	if(clFlag != CL_SUCCESS) {
		S->addMessage(3, "(RadixSort::_scan): Can't execute second pass kernel.\n");
	    if(clFlag == CL_INVALID_WORK_GROUP_SIZE)
	        S->addMessage(0, "\tInvalid local work group size.\n");
	    else if(clFlag == CL_OUT_OF_RESOURCES)
	        S->addMessage(0, "\tDevice out of resources.\n");
	    else if(clFlag == CL_MEM_OBJECT_ALLOCATION_FAILURE)
	        S->addMessage(0, "\tAllocation error at device.\n");
	    else if(clFlag == CL_OUT_OF_HOST_MEMORY)
	        S->addMessage(0, "\tfailure to allocate resources required by the OpenCL implementation on the host.\n");
	    return true;
	}
	#ifdef HAVE_GPUPROFILE
	    clFlag = clWaitForEvents(1, &event);
	    if(clFlag != CL_SUCCESS) {
	        S->addMessage(3, "(RadixSort::_scan): Can't wait to second pass kernel end.\n");
	        return true;
	    }
	    clFlag |= clGetEventProfilingInfo(event, CL_PROFILING_COMMAND_END, sizeof(cl_ulong), &end, 0);
	    clFlag |= clGetEventProfilingInfo(event, CL_PROFILING_COMMAND_START, sizeof(cl_ulong), &start, 0);
	    if(clFlag != CL_SUCCESS) {
	        S->addMessage(3, "(RadixSort::_scan): Can't profile second pass kernel execution.\n");
	        return true;
	    }
	    mTime += (end - start)/1000.f;  // 10^-3 ms
	#endif
	//! 5th.- Send arguments to third scan pass
	globalSize = mRadix*mGroups*mItems/2;
	localSize = globalSize / mHistoSplit;
	//! 6th.- Execute the third scan pass kernel
	#ifdef HAVE_GPUPROFILE
	    clFlag = clEnqueueNDRangeKernel(C->clComQueue, ckPasteHistogram, 1, NULL, &globalSize, &localSize, 0, NULL, &event);
	#else
	    clFlag = clEnqueueNDRangeKernel(C->clComQueue, ckPasteHistogram, 1, NULL, &globalSize, &localSize, 0, NULL, NULL);
	#endif
	if(clFlag != CL_SUCCESS) {
		S->addMessage(3, "(RadixSort::_scan): Can't execute paste pass kernel.\n");
	    if(clFlag == CL_INVALID_WORK_GROUP_SIZE)
	        S->addMessage(0, "\tInvalid local work group size.\n");
	    else if(clFlag == CL_OUT_OF_RESOURCES)
	        S->addMessage(0, "\tDevice out of resources.\n");
	    else if(clFlag == CL_MEM_OBJECT_ALLOCATION_FAILURE)
	        S->addMessage(0, "\tAllocation error at device.\n");
	    else if(clFlag == CL_OUT_OF_HOST_MEMORY)
	        S->addMessage(0, "\tfailure to allocate resources required by the OpenCL implementation on the host.\n");
	    return true;
	}
	#ifdef HAVE_GPUPROFILE
	    clFlag = clWaitForEvents(1, &event);
	    if(clFlag != CL_SUCCESS) {
	        S->addMessage(3, "(RadixSort::_scan): Can't wait to paste pass kernel end.\n");
	        return true;
	    }
	    clFlag |= clGetEventProfilingInfo(event, CL_PROFILING_COMMAND_END, sizeof(cl_ulong), &end, 0);
	    clFlag |= clGetEventProfilingInfo(event, CL_PROFILING_COMMAND_START, sizeof(cl_ulong), &start, 0);
	    if(clFlag != CL_SUCCESS) {
	        S->addMessage(3, "(RadixSort::_scan): Can't profile paste pass kernel execution.\n");
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
	int clFlag=0;
	size_t localSize =  mItems;
	size_t globalSize = mGroups*mItems;
	//! 1st.- Send arguments
	clFlag |= sendArgument(ckReorder, 0, sizeof(cl_mem), (void*)&clInKeys);
	clFlag |= sendArgument(ckReorder, 1, sizeof(cl_mem), (void*)&clOutKeys);
	clFlag |= sendArgument(ckReorder, 3, sizeof(cl_uint), (void*)&pass);
	clFlag |= sendArgument(ckReorder, 4, sizeof(cl_mem), (void*)&clInPermut);
	clFlag |= sendArgument(ckReorder, 5, sizeof(cl_mem), (void*)&clOutPermut);
	clFlag |= sendArgument(ckReorder, 6, sizeof(cl_uint)*mRadix*mItems, NULL);
	if(clFlag != CL_SUCCESS) {
		S->addMessage(3, (char*)"(RadixSort::_reorder): Can't send variable to kernel.\n");
	    return true;
	}
	//! 2nd.- Execute the kernel
	#ifdef HAVE_GPUPROFILE
	    cl_event event;
	    cl_ulong end, start;
	    clFlag = clEnqueueNDRangeKernel(C->clComQueue, ckReorder, 1, NULL, &globalSize, &localSize, 0, NULL, &event);
	#else
	    clFlag = clEnqueueNDRangeKernel(C->clComQueue, ckReorder, 1, NULL, &globalSize, &localSize, 0, NULL, NULL);
	#endif
	if(clFlag != CL_SUCCESS) {
		S->addMessage(3, "(RadixSort::_reorder): Can't execute the kernel.\n");
	    if(clFlag == CL_INVALID_PROGRAM_EXECUTABLE)
	        S->addMessage(0, "\tInvalid program (Compile errors maybe?).\n");
	    else if(clFlag == CL_INVALID_COMMAND_QUEUE)
	        S->addMessage(0, "\tInvalid command queue.\n");
	    else if(clFlag == CL_INVALID_KERNEL)
	        S->addMessage(0, "\tKernel is not a valid object.\n");
	    else if(clFlag == CL_INVALID_CONTEXT)
	        S->addMessage(0, "\tContext associated to command queue don't match qith kernel context.\n");
	    else if(clFlag == CL_INVALID_KERNEL_ARGS)
	        S->addMessage(0, "\tOne or more arguments are invalid (maybe don't specified).\n");
	    else if(clFlag == CL_INVALID_WORK_DIMENSION)
	        S->addMessage(0, "\tDimension must be a value between 1 and 3.\n");
	    else if(clFlag == CL_INVALID_WORK_GROUP_SIZE)
	        S->addMessage(0, "\tInvalid local work group size.\n");
	    else if(clFlag == CL_INVALID_WORK_ITEM_SIZE)
	        S->addMessage(0, "\tLocal work group size is out of bounds.\n");
	    else if(clFlag == CL_INVALID_GLOBAL_OFFSET)
	        S->addMessage(0, "\tGlobal offset must be NULL.\n");
	    else if(clFlag == CL_OUT_OF_RESOURCES)
	        S->addMessage(0, "\tDevice out of resources.\n");
	    else if(clFlag == CL_MEM_OBJECT_ALLOCATION_FAILURE)
	        S->addMessage(0, "\tAllocation error at device.\n");
	    else if(clFlag == CL_INVALID_EVENT_WAIT_LIST)
	        S->addMessage(0, "\tInvalid event wait instruction.\n");
	    else if(clFlag == CL_OUT_OF_HOST_MEMORY)
	        S->addMessage(0, "\tfailure to allocate resources required by the OpenCL implementation on the host.\n");
	    return true;
	}
	#ifdef HAVE_GPUPROFILE
	    clFlag = clWaitForEvents(1, &event);
	    if(clFlag != CL_SUCCESS) {
	        S->addMessage(3, (char*)"(RadixSort::_reorder): Can't wait to kernels end.\n");
	        return true;
	    }
	    clFlag |= clGetEventProfilingInfo(event, CL_PROFILING_COMMAND_END, sizeof(cl_ulong), &end, 0);
	    clFlag |= clGetEventProfilingInfo(event, CL_PROFILING_COMMAND_START, sizeof(cl_ulong), &start, 0);
	    if(clFlag != CL_SUCCESS) {
	        S->addMessage(3, (char*)"(RadixSort::_reorder): Can't profile kernel execution.\n");
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
	int clFlag=0;
	//! 1st.- Send arguments
	clFlag |= sendArgument(ckReversePermutations, 0, sizeof(cl_mem), (void*)&clInPermut);
	if(clFlag != CL_SUCCESS) {
		S->addMessage(3, "(RadixSort::_reversePermutations): Can't send variable to kernel.\n");
	    return true;
	}
	//! 2nd.- Execute the kernel
	#ifdef HAVE_GPUPROFILE
	    cl_event event;
	    cl_ulong end, start;
	    clFlag = clEnqueueNDRangeKernel(C->clComQueue, ckReversePermutations, 1, NULL, &clGlobalWorkSize, NULL, 0, NULL, &event);
	#else
	    clFlag = clEnqueueNDRangeKernel(C->clComQueue, ckReversePermutations, 1, NULL, &clGlobalWorkSize, NULL, 0, NULL, NULL);
	#endif
	if(clFlag != CL_SUCCESS) {
		S->addMessage(3, "(RadixSort::_reversePermutations): Can't execute the kernel.\n");
	    if(clFlag == CL_INVALID_WORK_GROUP_SIZE)
	        S->addMessage(0, "\tInvalid local work group size.\n");
	    else if(clFlag == CL_OUT_OF_RESOURCES)
	        S->addMessage(0, "\tDevice out of resources.\n");
	    else if(clFlag == CL_MEM_OBJECT_ALLOCATION_FAILURE)
	        S->addMessage(0, "\tAllocation error at device.\n");
	    else if(clFlag == CL_OUT_OF_HOST_MEMORY)
	        S->addMessage(0, "\tfailure to allocate resources required by the OpenCL implementation on the host.\n");
	    return true;
	}
	#ifdef HAVE_GPUPROFILE
	    clFlag = clWaitForEvents(1, &event);
	    if(clFlag != CL_SUCCESS) {
	        S->addMessage(3, "(RadixSort::_reversePermutations): Can't wait to kernels end.\n");
	        return true;
	    }
	    clFlag |= clGetEventProfilingInfo(event, CL_PROFILING_COMMAND_END, sizeof(cl_ulong), &end, 0);
	    clFlag |= clGetEventProfilingInfo(event, CL_PROFILING_COMMAND_START, sizeof(cl_ulong), &start, 0);
	    if(clFlag != CL_SUCCESS) {
	        S->addMessage(3, "(RadixSort::_reversePermutations): Can't profile kernel execution.\n");
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
	int clFlag = 0;
	if(n == N)
	    return false;
	if(N <= 0){
	    sprintf(msg, "(RadixSort::setN): Any data to sort (n=%u)\n", N);
	    S->addMessage(3, msg);
	    return true;
	}
	clLocalWorkSize  = getLocalWorkSize(N, C->clComQueue);
	if(!clLocalWorkSize){
	    sprintf(msg, "(RadixSort::setN): No valid local work size for %u cells.\n", N);
	    S->addMessage(3, msg);
	    return true;
	}
	clGlobalWorkSize = getGlobalWorkSize(N, clLocalWorkSize);
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
	clInKeys = C->lcell;
	clOutPermut = C->permutation;
	// Keys and permutations
	if(clOutKeys){
	    C->AllocatedMem -= oldN*sizeof( unsigned int );
	    clReleaseMemObject(clOutKeys);
	    clOutKeys=0;
	}
	clOutKeys = C->allocMemory(n*sizeof( unsigned int ));
	if(!clOutKeys)
        return true;
	if(clInPermut){
	    C->AllocatedMem -= oldN*sizeof( unsigned int );
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
	    C->AllocatedMem -= (mRadix * mGroups * mItems)*sizeof( unsigned int );
	}
	clHistograms = C->allocMemory((mRadix * mGroups * mItems)*sizeof( unsigned int ));
	if(!clHistograms)
        return true;
	if(clGlobalSums){
	    clReleaseMemObject(clGlobalSums);
	    clGlobalSums=0;
	    C->AllocatedMem -= (mHistoSplit)*sizeof( unsigned int );
	}
	clGlobalSums = C->allocMemory(mHistoSplit*sizeof( unsigned int ));
	if(!clGlobalSums)
        return true;
	if(clTempMem){
	    clReleaseMemObject(clTempMem);
	    clTempMem=0;
	    C->AllocatedMem -= sizeof( unsigned int );
	}
	clTempMem = C->allocMemory(sizeof( unsigned int ));
	if(!clTempMem)
        return true;
	sprintf(msg, "Radix sort components = %u\n", n);
	S->addMessageF(1, msg);
	sprintf(msg, "\tAllocated memory = %lu bytes\n", C->AllocatedMem);
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
	int clFlag = 0;
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
	if(clProgram)clReleaseProgram(clProgram);clProgram=0;
    // Load the new kernels
	if(!loadKernelFromFile(&ckInit, &clProgram, C->clContext, C->clDevice, Path, "init", CFlags))
	    return true;
	clReleaseProgram(clProgram);clProgram=0;
	if(!loadKernelFromFile(&ckTranspose, &clProgram, C->clContext, C->clDevice, Path, "transpose", CFlags))
	    return true;
	clReleaseProgram(clProgram);clProgram=0;
	if(!loadKernelFromFile(&ckHistogram, &clProgram, C->clContext, C->clDevice, Path, "histogram", CFlags))
	    return true;
	clReleaseProgram(clProgram);clProgram=0;
	if(!loadKernelFromFile(&ckScanHistogram, &clProgram, C->clContext, C->clDevice, Path, "scanhistograms", CFlags))
	    return true;
	clReleaseProgram(clProgram);clProgram=0;
	if(!loadKernelFromFile(&ckPasteHistogram, &clProgram, C->clContext, C->clDevice, Path, "pastehistograms", CFlags))
	    return true;
	clReleaseProgram(clProgram);clProgram=0;
	if(!loadKernelFromFile(&ckReorder, &clProgram, C->clContext, C->clDevice, Path, "reorder", CFlags))
	    return true;
	clReleaseProgram(clProgram);clProgram=0;
	if(!loadKernelFromFile(&ckReversePermutations, &clProgram, C->clContext, C->clDevice, Path, "reversePermutation", CFlags))
	    return true;
	// ------------------------------------------------------------------------
    // Send the fixed arguments
	// ------------------------------------------------------------------------
	clFlag |= sendArgument(ckInit, 1, sizeof(cl_uint), (void*)&n);
	if(clFlag != CL_SUCCESS) {
		S->addMessage(3, "(RadixSort::setupOpenCL): Fail sending fix variables to the initialization kernel.\n");
	    return true;
	}
	clFlag |= sendArgument(ckHistogram, 1, sizeof(cl_mem), (void*)&clHistograms);
	clFlag |= sendArgument(ckHistogram, 4, sizeof(cl_uint), (void*)&n);
	if(clFlag != CL_SUCCESS) {
		S->addMessage(3, "(RadixSort::setupOpenCL): Fail sending fix variables to the histograms kernel.\n");
	    return true;
	}
	clFlag |= sendArgument(ckPasteHistogram, 0, sizeof(cl_mem), (void*)&clHistograms);
	clFlag |= sendArgument(ckPasteHistogram, 1, sizeof(cl_mem), (void*)&clGlobalSums);
	if(clFlag != CL_SUCCESS) {
		S->addMessage(3, "(RadixSort::setupOpenCL): Fail sending fix variables to the histograms pasting kernel.\n");
	    return true;
	}
	clFlag |= sendArgument(ckReorder, 2, sizeof(cl_mem), (void*)&clHistograms);
	clFlag |= sendArgument(ckReorder, 7, sizeof(cl_uint), (void*)&n);
	if(clFlag != CL_SUCCESS) {
		S->addMessage(3, "(RadixSort::setupOpenCL): Fail sending fix variables to the sorting kernel.\n");
	    return true;
	}
	clFlag |= sendArgument(ckReversePermutations, 1, sizeof(cl_mem), (void*)&C->reversePermutation);
	clFlag |= sendArgument(ckReversePermutations, 2, sizeof(cl_uint), (void*)&n);
	if(clFlag != CL_SUCCESS) {
		S->addMessage(3, "(RadixSort::setupOpenCL): Fail sending fix variables to the reverse permutations kernel.\n");
	    return true;
	}
	// ------------------------------------------------------------------------
	// Local work sizes
	// ------------------------------------------------------------------------
    cl_device_id device;
	clFlag |= clGetCommandQueueInfo(C->clComQueue,CL_QUEUE_DEVICE,
	                                sizeof(cl_device_id),&device, NULL);
	if(clFlag != CL_SUCCESS) {
		S->addMessage(3, "(RadixSort::setupOpenCL): Can't get the device from the command queue.\n");
		S->addMessage(0, "\tUnhandled exception\n");
	    return true;
	}
	cl_uint dims;
	clFlag |= clGetDeviceInfo(device,CL_DEVICE_MAX_WORK_ITEM_DIMENSIONS,sizeof(cl_uint),&dims, NULL);
	if(clFlag != CL_SUCCESS){
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
	clFlag |= clGetKernelWorkGroupInfo(ckTranspose,device,CL_KERNEL_WORK_GROUP_SIZE,
	                                   sizeof(size_t), &maxLocalSize, NULL);
	if(clFlag != CL_SUCCESS) {
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
	clFlag |= clGetKernelWorkGroupInfo(ckReorder,device,CL_KERNEL_WORK_GROUP_SIZE,
	                                   sizeof(size_t), &maxLocalSize, NULL);
	if(clFlag != CL_SUCCESS) {
		S->addMessage(3, "(LinkList::setupOpenCL): Can't get maximum local work group size from the histograms kernel.\n");
	    return true;
	}
	if(maxLocalSize < mItems) mItems = maxLocalSize;
	clFlag |= clGetKernelWorkGroupInfo(ckHistogram,device,CL_KERNEL_WORK_GROUP_SIZE,
	                                   sizeof(size_t), &maxLocalSize, NULL);
	if(clFlag != CL_SUCCESS) {
		S->addMessage(3, "(LinkList::setupOpenCL): Can't get maximum local work group size from the histograms kernel.\n");
		S->addMessage(0, "\tUnhandled exception\n");
	    return true;
	}
	if(maxLocalSize < mItems) mItems = maxLocalSize;
    // With the scan process we can easily set a bound for the number of histograms splits as well
	clFlag |= clGetKernelWorkGroupInfo(ckScanHistogram,device,CL_KERNEL_WORK_GROUP_SIZE,
	                                   sizeof(size_t), &maxLocalSize, NULL);
	if(clFlag != CL_SUCCESS) {
		S->addMessage(3, "(LinkList::setupOpenCL): Can't get maximum local work group size from the histograms scan kernel.\n");
		S->addMessage(0, "\tUnhandled exception\n");
	    return true;
	}
	if(maxLocalSize < mHistoSplit/2) mHistoSplit = 2*maxLocalSize;
    // With the scan histograms kernel and the paste them one we must adjust GROUPS, ITEMS and RADIX
    size_t maxForScan = maxLocalSize;
	clFlag |= clGetKernelWorkGroupInfo(ckReorder,device,CL_KERNEL_WORK_GROUP_SIZE,
	                                   sizeof(size_t), &maxLocalSize, NULL);
	if(clFlag != CL_SUCCESS) {
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
	clFlag |= clGetDeviceInfo(device, CL_DEVICE_LOCAL_MEM_SIZE, sizeof(cl_ulong), &availableMem, NULL);
	if(clFlag != CL_SUCCESS) {
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
	clFlag |= sendArgument(ckTranspose, 6, sizeof(cl_uint)*tilesize*tilesize, NULL);
	clFlag |= sendArgument(ckTranspose, 7, sizeof(cl_uint)*tilesize*tilesize, NULL);
	if(clFlag != CL_SUCCESS) {
		S->addMessage(3, "(RadixSort::setupOpenCL): Can't set the local memory for the transpose kernel.\n");
		S->addMessage(0, "\tUnhandled exception\n");
	    return true;
	}
	clFlag |= clGetKernelWorkGroupInfo(ckTranspose,device,CL_KERNEL_LOCAL_MEM_SIZE,
	                                   sizeof(size_t), &usedMem, NULL);
	if(clFlag != CL_SUCCESS) {
		S->addMessage(3, "(RadixSort::setupOpenCL): Can't get the local memory used by the transpose kernel.\n");
		if(clFlag == CL_INVALID_DEVICE)
            S->addMessage(0, "\tCL_INVALID_DEVICE\n");
		else if(clFlag == CL_INVALID_VALUE)
            S->addMessage(0, "\tCL_INVALID_VALUE\n");
		else if(clFlag == CL_INVALID_KERNEL)
            S->addMessage(0, "\tCL_INVALID_KERNEL\n");
        else{
            sprintf(msg, "\tUnhandled exception %i\n", clFlag);
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
	clFlag |= sendArgument(ckHistogram, 3, sizeof(cl_uint)*mRadix*mItems, NULL);
	clFlag |= sendArgument(ckReorder, 6, sizeof(cl_uint)*mRadix*mItems, NULL);
	if(clFlag != CL_SUCCESS) {
		S->addMessage(3, "(RadixSort::setupOpenCL): Can't set the local memory for histograms and sorting kernels.\n");
		S->addMessage(0, "\tUnhandled exception\n");
	    return true;
	}
	cl_ulong histoUsed = 0;
	clFlag |= clGetKernelWorkGroupInfo(ckHistogram,device,CL_KERNEL_LOCAL_MEM_SIZE,
	                                   sizeof(size_t), &histoUsed, NULL);
	clFlag |= clGetKernelWorkGroupInfo(ckReorder,device,CL_KERNEL_LOCAL_MEM_SIZE,
	                                   sizeof(size_t), &usedMem, NULL);
	if(clFlag != CL_SUCCESS) {
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
        clFlag |= sendArgument(ckHistogram, 3, sizeof(cl_uint)*mRadix*mItems, NULL);
        clFlag |= sendArgument(ckReorder, 6, sizeof(cl_uint)*mRadix*mItems, NULL);
        if(clFlag != CL_SUCCESS) {
            S->addMessage(3, "(RadixSort::setupOpenCL): Can't set the local memory for histograms and sorting kernels.\n");
            S->addMessage(0, "\tUnhandled exception\n");
            return true;
        }
        clFlag |= clGetKernelWorkGroupInfo(ckHistogram,device,CL_KERNEL_LOCAL_MEM_SIZE,
                                           sizeof(size_t), &histoUsed, NULL);
        clFlag |= clGetKernelWorkGroupInfo(ckReorder,device,CL_KERNEL_LOCAL_MEM_SIZE,
                                           sizeof(size_t), &usedMem, NULL);
        if(clFlag != CL_SUCCESS) {
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
    clFlag |= sendArgument(ckScanHistogram, 1, sizeof(cl_uint)* maxmemcache, NULL);
    if(clFlag != CL_SUCCESS) {
        S->addMessage(3, "(RadixSort::setupOpenCL): Can't set the local memory for scanning kernel.\n");
        S->addMessage(0, "\tUnhandled exception\n");
        return true;
    }
    clFlag |= clGetKernelWorkGroupInfo(ckScanHistogram,device,CL_KERNEL_LOCAL_MEM_SIZE,
                                       sizeof(size_t), &usedMem, NULL);
    if(clFlag != CL_SUCCESS) {
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
