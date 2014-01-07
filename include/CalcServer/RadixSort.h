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

#ifndef RADIXSORT_H_INCLUDED
#define RADIXSORT_H_INCLUDED

// ----------------------------------------------------------------------------
// Include Prerequisites
// ----------------------------------------------------------------------------
#include <sphPrerequisites.h>

// ----------------------------------------------------------------------------
// Include standar libraries
// ----------------------------------------------------------------------------
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

// ----------------------------------------------------------------------------
// Include OpenCL libraries
// ----------------------------------------------------------------------------
#include <CL/cl.h>

// ----------------------------------------------------------------------------
// Include auxiliar methods
// ----------------------------------------------------------------------------
#include <AuxiliarMethods.h>

/** @def _ITEMS Number of items in a group
 * @note Must be power of 2, and in some devices greather than 32.
 */
#ifndef _ITEMS
	#define _ITEMS  128
#endif
/** @def _GROUPS Number of groups (data must be divisible of _ITEMS*_GROUPS)
 * @note Must be power of 2
 */
#ifndef _GROUPS
	#define _GROUPS 32
#endif
/// @def __UINTBITS__ Bits of an unsigned integer variable
#ifndef __UINTBITS__
	#define __UINTBITS__ 32
#endif
/// @def _STEPBITS Bits that be sorted in each pass
#ifndef _STEPBITS
	#define _STEPBITS 4
#endif
/// @def _RADIX Bits that be sorted in each pass
#define _RADIX (1 << _STEPBITS)
/** @def _HISTOSPLIT Number of splits of the histogram
 * @remarks (_GROUPS * _ITEMS * _RADIX) % _HISTOSPLIT must be equal to zero
 * @remarks (2 * _HISTOSPLIT) % _RADIX must be equal to zero
 * @remarks Must be power of 2, and in some devices greather than 64.
 */
#ifndef _HISTOSPLIT
	#define _HISTOSPLIT 512
#endif

/* Modify data in order to impose that local size don't be lower than
* minimum allowed value.
*/
#if _ITEMS < __CL_MIN_LOCALSIZE__
	#undef _ITEMS
	#define _ITEMS __CL_MIN_LOCALSIZE__
#endif
#if _HISTOSPLIT/2 < __CL_MIN_LOCALSIZE__
	#undef _HISTOSPLIT
	#define _HISTOSPLIT 2*__CL_MIN_LOCALSIZE__
#endif
#if _RADIX*_GROUPS*_ITEMS/2/_HISTOSPLIT < __CL_MIN_LOCALSIZE__
	#undef _GROUPS
	#define _GROUPS 2*__CL_MIN_LOCALSIZE__*_HISTOSPLIT / (_RADIX*_ITEMS)
#endif

#ifndef NDEBUG
    #pragma message "_ITEMS: " XSTR(_ITEMS)
    #pragma message "_GROUPS: " XSTR(_GROUPS)
    #pragma message "_HISTOSPLIT: " XSTR(_HISTOSPLIT)
#endif

namespace Aqua{ namespace CalcServer{

/** @class RadixSort RadixSort.h CalcServer/RadixSort.h
 * @brief Methods to perform a radix sort using the GPU (or any device supported by OpenCL).
 * The code has 3 steps:
 * <ol><li> Create the histogram.</li>
 * <li> Scan the histogram to create accumulated histograms.</li>
 * <li> Permut the variables.</li></ul>
 * In order to improve the data access, the info is transposed. \n
 * Radix sort requires also several pass (VarBits / Radix). \n
 * To learn more about this code, please see also http://code.google.com/p/ocl-radix-sort/updates/list
 */
class RadixSort
{
public:
	/** Constructor
	 */
	RadixSort();

	/** Destructor
	 */
	~RadixSort();

	/** Sorts the icell components, returning permutations too.
	 * @return false if all gone right. \n true otherwise.
	 * @note This structure assumes number of particle as amount of
	 * data to sort, if not correct change it before calling this method.
	 * @warning After calling this method, CalcServer icell array, and permutation
	 * array may change (as memory direction), don't forgive resend it to the
	 * forwarded kernels.
	 */
	bool sort();

	/** Set the number of elements
	 * @param nElements Number of elements to sort.
	 * @return false if all gone right. \n true otherwise.
	 * @remarks After calling this method, probably you want to call
	 * SetupOpenCL too.
	 */
	bool setN(unsigned int nElements);

	#ifdef HAVE_GPUPROFILE
	    /** Set the kernel time consumed.
	     * @param t Kernel time consumed.
	     */
	    void profileTime(float t){_time = t;}
	    /** Get the kernel time consumed.
	     * @return Kernel time consumed.
	     */
	    float profileTime(){return _time;}
	#endif

private:
	/** Setup OpenCL kernel
	 * @return false if all gone right. \n true otherwise.
	 */
	bool setupOpenCL();

	/** Initialize permutations array.
	 * @return false if all gone right. \n true otherwise.
	 */
	bool _init();

	/** Transpose the permutations matrix to improve the data access.
	 * @param nbrow Number of rows.
	 * @param nbcol Number of columns.
	 * @return false if all gone right. \n true otherwise.
	 */
	bool _transpose(unsigned int nbrow, unsigned int nbcol);

	/** Perform histograms.
	 * @return false if all gone right. \n true otherwise.
	 */
	bool _histograms();

	/** Scan histograms.
	 * @return false if all gone right. \n true otherwise.
	 */
	bool _scan();

	/** Scan histograms.
	 * @return false if all gone right. \n true otherwise.
	 */
	bool _reorder();

	/** Build the reversed permutations vector.
	 * @return false if all gone right. \n true otherwise.
	 */
	bool _reversePermutations();

	/// Path of the kernels
	char* Path;
	/// Number of elements to sort (don't set manually)
	unsigned int n;
	/// Local work size (default value = 128)
	size_t _local_work_size;
	/// Global work size
	size_t _global_work_size;
	/// Key bits (maximum)
	unsigned int keyBits;
	/// Needed radix pass (keyBits / _STEPBITS)
	unsigned int nPass;
	/// Active pass (keyBits / _STEPBITS)
	unsigned int pass;

	/// Input keys (CalcServer icell array)
	cl_mem clInKeys;
	/// Output keys
	cl_mem clOutKeys;
	/// Input permutations
	cl_mem clInPermut;
	/// Output permutations (CalcServer permutations array)
	cl_mem clOutPermut;
	/// Histograms
	cl_mem clHistograms;
	/// Sums
	cl_mem clGlobalSums;
	/// Temporal memory
	cl_mem clTempMem;
	/// OpenCL initialization kernel
	cl_kernel ckInit;
	/// Transposition algorithm
	cl_kernel ckTranspose;
	/// OpenCL histogram kernel
	cl_kernel ckHistogram;
	/// OpenCL scan histogram kernel
	cl_kernel ckScanHistogram;
	/// OpenCL paste histogram kernel
	cl_kernel ckPasteHistogram;
	/// OpenCL permutations kernel
	cl_kernel ckReorder;
	/// OpenCL reverse permutations kernel
	cl_kernel ckReversePermutations;
	/// OpenCL program
	cl_program _program;

	#ifdef HAVE_GPUPROFILE
	    /// Kernel real time consumed
	    float _time;
	#endif

    /// Number of items in a group
    unsigned int mItems;
    /// Number of groups in a radix
    unsigned int mGroups;
    /// Bits of the Radix
    unsigned int mBits;
    /// Number of radices
    unsigned int mRadix;
    /// Splits of the histogram
    unsigned int mHistoSplit;

	/// Tile size warning shown
	bool tilesizeWarning;
};

}}  // namespace

#endif // RADIXSORT_H_INCLUDED
