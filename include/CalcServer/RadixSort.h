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

#include <sphPrerequisites.h>
#include <CL/cl.h>

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
 * @brief Methods to perform a radix sort using the GPU (or any device
 * supported by OpenCL).
 * The code has 3 steps:
 *   -# Create the histogram.
 *   -# Scan the histogram to create the accumulated one.
 *   -# Permut the variables.
 * To learn more about this code, please see also
 * http://code.google.com/p/ocl-radix-sort/updates/list
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

    /** Sorts the icell components, computing the permutations required too.
     * @return false if all gone right, true otherwise.
     * @warning After calling this method, CalcServer icell array, and
     * permutation array may change (including their memory direction), so
     ' don't forgive to resend it to the next kernels.
     */
    bool sort();

    /** Set the number of elements to sort
     * @param n_elements Number of elements to sort.
     * @return false if all gone right, true otherwise.
     */
    bool setN(unsigned int n_elements);

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
    /** Setup the OpenCL stuff
     * @return false if all gone right, true otherwise.
     */
    bool setupOpenCL();

    /** Initialize permutations array.
     * @return false if all gone right, true otherwise.
     */
    bool init();

    /** Transpose the permutations matrix to improve the data access.
     * @param nbrow Number of rows.
     * @param nbcol Number of columns.
     * @return false if all gone right, true otherwise.
     */
    bool transpose(unsigned int nbrow, unsigned int nbcol);

    /** Perform histograms.
     * @return false if all gone right, true otherwise.
     */
    bool histograms();

    /** Scan histograms.
     * @return false if all gone right, true otherwise.
     */
    bool scan();

    /** Scan histograms.
     * @return false if all gone right, true otherwise.
     */
    bool reorder();

    /** Build the reversed permutations vector.
     * @return false if all gone right, true otherwise.
     */
    bool reversePermutations();

    /// Path of the kernels
    char* _path;
    /// Number of elements to sort (don't set manually)
    unsigned int _n;
    /// Local work size (default value = 128)
    size_t _local_work_size;
    /// Global work size
    size_t _global_work_size;
    /// Key bits (maximum)
    unsigned int _key_bits;
    /// Needed radix pass (_key_bits / _STEPBITS)
    unsigned int _n_pass;
    /// Active pass (_key_bits / _STEPBITS)
    unsigned int _pass;

    /// Input keys (CalcServer icell array)
    cl_mem _in_keys;
    /// Output keys
    cl_mem _out_keys;
    /// Input permutations
    cl_mem _in_permut;
    /// Output permutations (CalcServer permutations array)
    cl_mem _out_permut;
    /// Histograms
    cl_mem _histograms;
    /// Sums
    cl_mem _global_sums;
    /// Temporal memory
    cl_mem _temp_mem;
    /// OpenCL initialization kernel
    cl_kernel _init_kernel;
    /// Transposition algorithm
    cl_kernel _transpose_kernel;
    /// OpenCL histogram kernel
    cl_kernel _histograms_kernel;
    /// OpenCL scan histogram kernel
    cl_kernel _histograms_scan_kernel;
    /// OpenCL paste histogram kernel
    cl_kernel _paste_histograms_kernel;
    /// OpenCL permutations kernel
    cl_kernel _sort_kernel;
    /// OpenCL reverse permutations kernel
    cl_kernel _inv_permutations_kernel;
    /// OpenCL program
    cl_program _program;

    #ifdef HAVE_GPUPROFILE
        /// Kernel real time consumed
        float _time;
    #endif

    /// Number of items in a group
    unsigned int _items;
    /// Number of groups in a radix
    unsigned int _groups;
    /// Bits of the Radix
    unsigned int _bits;
    /// Number of radices
    unsigned int _radix;
    /// Splits of the histogram
    unsigned int _histo_split;

    /// Tile size warning shown
    bool _tilesize_warn;
};

}}  // namespace

#endif // RADIXSORT_H_INCLUDED
