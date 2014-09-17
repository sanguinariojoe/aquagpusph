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
#include <CalcServer.h>
#include <CalcServer/Tool.h>

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
 * http://code.google.com/p/ocl-radix-sort/updates/list.
 */
class RadixSort : public Aqua::CalcServer::Tool
{
public:
    /** Constructor.
     * @param tool_name Tool name.
     * @param variable Variable to sort.
     * @param permutations Variable where the permutations will be stored.
     * @param inv_permutations Variable where the inverse permutations will be
     * stored.
     */
    RadixSort(const char* tool_name,
              const char* variable="icell",
              const char* permutations="id_unsorted",
              const char* inv_permutations="id_sorted");

    /** Destructor
     */
    ~RadixSort();

    /** Initialize the tool.
     * @return false if all gone right, true otherwise.
     */
    bool setup();

    /** Execute the tool.
     * @return false if all gone right, true otherwise.
     */
    bool execute();

private:
    /** Get the variables to compute.
     * @return false if all gone right, true otherwise.
     */
    bool variables();

    /** Setup the OpenCL stuff
     * @return false if all gone right, true otherwise.
     */
    bool setupOpenCL();

    /** Compile the source code and generate the corresponding kernels
     * @param source Source code to compile.
     * @return false if all gone right, true otherwise.
     */
    bool compile(const char* source);

    /** Setup the main computing dimensions _items, _groups and _histo_split
     * from the valid local work sizes per each kernel.
     * @return false if all gone right, true otherwise.
     */
    bool setupDims();

    /// Variable to sort name
    char *_var_name;

    /// Permutations array name
    char *_perms_name;

    /// Inverse permutations array name
    char *_inv_perms_name;

    /// Variable to sort
    InputOutput::ArrayVariable *_var;

    /// Permutations array
    InputOutput::ArrayVariable *_perms;

    /// Inverse permutations array
    InputOutput::ArrayVariable *_inv_perms;

    /// Number of keys to sort
    unsigned int _n;

    /// OpenCL initialization kernel
    cl_kernel _init_kernel;
    /// OpenCL histogram kernel
    cl_kernel _histograms_kernel;
    /// OpenCL scan histogram kernel
    cl_kernel _scan_kernel;
    /// OpenCL paste histogram kernel
    cl_kernel _paste_kernel;
    /// OpenCL permutations kernel
    cl_kernel _sort_kernel;
    /// OpenCL reverse permutations kernel
    cl_kernel _inv_perms_kernel;

    /// Number of items in a group
    unsigned int _items;
    /// Number of groups in a radix
    unsigned int _groups;
    /// Bits of the Radix
    unsigned int _bits;
    /// Number of Radices
    unsigned int _radix;
    /// Splits of the histogram
    unsigned int _histo_split;
};

}}  // namespace

#endif // RADIXSORT_H_INCLUDED
