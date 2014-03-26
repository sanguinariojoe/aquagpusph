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

/** @brief Methods to perform a radix sort using the GPU (or any device supported by OpenCL).
 * The code has 3 steps:
 * <ol><li> Create the histogram.</li>
 * <li> Scan the histogram to create accumulated histograms.</li>
 * <li> Permut the variables.</li></ul>
 * In order to improve the data access, the info is transposed. \n
 * Radix sort requires also several pass (VarBits / Radix). \n
 * To learn more about this code, please see also http://code.google.com/p/ocl-radix-sort/updates/list
 */

#ifndef HAVE_3D
    #include "types/2D.h"
#else
    #include "types/3D.h"
#endif

/// @def _BITS Bits of the Radix
#ifndef _BITS
	#define _BITS 8
#endif
/// @def _RADIX Resultant Radix
#ifndef _RADIX
	#define _RADIX (1 << _BITS)
#endif

#ifdef _g
	#error '_g' is already defined.
#endif
#define _g __global

#ifdef _l
	#error '_l' is already defined.
#endif
#define _l __local

#ifdef _c
	#error '_c' is already defined.
#endif
#define _c __constant

#pragma OPENCL EXTENSION cl_khr_byte_addressable_store : enable 

/** Initializes the permutations.
 * @param d_inPermut Initial permutations (null)
 * @param n Number of keys.
 */
__kernel void init(__global unsigned int* d_inPermut,
		      unsigned int n){
    // find position in global arrays
    unsigned int i = get_global_id(0);
	if(i >= n)
		return;

	d_inPermut[i] = i;
}

/** Transpose the input data to improve the data memory access in the next
 * steps. An array sorted such that the indexes are \n
 * [0, 1, ..., nbrow, 1+nbrow, ..., 2*nbrow, ..., n-1] \n
 * will be transposed as: \n
 * [0, nbrow, 2*nbrow, ..., 1, 1+nbrow, ..., n-1]
 * @param invect Input vector of keys
 * @param outvect Transposed vector of keys
 * @param nbcol nbrow / n
 * @param nbrow Number of rows
 * @param inperm Input permutations array
 * @param outperm Transposed permutations array
 * @param blockmat Shared memory keys matrix to accelerate the process
 * @param blockprem Shared memory permutations matrix to accelerate the process
 * @param tilesize Dimensions of the shared memory matrices
 * @remarks If the local size is not big enough, it will be artificially
 * increased due to some implementations will not support it otherwise
 * @note In the initialization nbrow = _ITEMS*_GROUPS, which of course will be
 * untransposed at the end of the radix sort process with
 * nbcol = _ITEMS*_GROUPS
 */
__kernel void transpose(const __global unsigned int* invect,
                        __global unsigned int* outvect,
                        const unsigned int nbcol,
                        const unsigned int nbrow,
                        const __global unsigned int* inperm,
                        __global unsigned int* outperm,
                        __local unsigned int* blockmat,
                        __local unsigned int* blockperm,
                        const unsigned int tilesize){
  
	unsigned int i0 = get_global_id(0)*tilesize;  // first row index
	unsigned int j  = get_global_id(1);           // column index

	unsigned int jloc = get_local_id(1);          // local column index

	// Discard out of tilesize threads, which has been artificially append
	if(jloc >= tilesize)
		return;
	// And therefore the j index must be corrected.
	j /= get_local_size(1) / tilesize;

	// fill the cache
	for(unsigned int iloc=0;iloc<tilesize;iloc++){
		unsigned int k=(i0+iloc)*nbcol+j;  // position in the matrix
		// Duplicate tilesize elements of original array into the columns of
		// local memory allocated matrix tilesize x tilesize (all colums are the same)
		blockmat[iloc*tilesize+jloc]=invect[k];
		#ifdef PERMUT 
			blockperm[iloc*tilesize+jloc]=inperm[k];
		#endif
	}

	barrier(CLK_LOCAL_MEM_FENCE);  

	// first row index in the transpose
	unsigned int j0=get_group_id(1)*tilesize;

	// put the cache at the good place
	for(unsigned int iloc=0;iloc<tilesize;iloc++){
		unsigned int kt=(j0+iloc)*nbrow+i0+jloc;  // position in the transpose
		outvect[kt]=blockmat[jloc*tilesize+iloc];
		#ifdef PERMUT 
			outperm[kt]=blockperm[jloc*tilesize+iloc];
		#endif
	}
}

/** Perform the local histograms. The histograms are the number of occurrences
 * of each radix. Since we are working in parallel, the number of ocurrences
 * of each radix will be splited in blocks of dimension _ITEMS*_GROUPS:
 *   | it(0)gr(0)ra(0) | it(1)gr(0)ra(0) | ... | it(items)gr(0)ra(0) |
 *   | it(0)gr(1)ra(0) | it(1)gr(1)ra(0) | ... | it(items)gr(1)ra(0) |
 *   | ... | | | |
 *   | it(0)gr(groups)ra(0) | it(1)gr(groups)ra(0) | ... | it(items)gr(groups)ra(0) |
 *   | it(0)gr(0)ra(1) | it(1)gr(0)ra(1) | ... | it(items)gr(0)ra(1) |
 *   | ... | | | |
 *   | it(0)gr(groups)ra(1) | it(1)gr(groups)ra(1) | ... | it(items)gr(groups)ra(1) |
 *   | ... | | | |
 *   | it(0)gr(groups)ra(radix) | it(1)gr(groups)ra(radix) | ... | it(items)gr(groups)ra(radix) |
 * where it is the thread, gr is the group, and ra is the radix.
 * @param d_Keys Input unsorted keys allocated into the device.
 * @param d_Histograms Input histograms.
 * @param pass Pass of the radix decomposition.
 * @param loc_histo Histogram local memory to speed up the process.
 * @param n Number of keys.
 */
__kernel void histogram(const __global unsigned int* d_Keys,
			__global unsigned int* d_Histograms,
			const unsigned int pass,
			__local unsigned int* loc_histo,
			const unsigned int n)
{
	unsigned int it = get_local_id(0);
	unsigned int ig = get_global_id(0);
	unsigned int gr = get_group_id(0);

	unsigned int groups = get_num_groups(0);
	unsigned int items  = get_local_size(0);

	// Initializate the histograms in each thread of this work group
	for(unsigned int ir=0;ir<_RADIX;ir++){
		loc_histo[ir * items + it] = 0;
	}

	barrier(CLK_LOCAL_MEM_FENCE);  

	// Set the keys analized by each thread
	unsigned int size = n/groups/items;
	#ifndef TRANSPOSE
	// If the data has not been transposed we must start reading from a
	// different place of ig
	unsigned int start= ig * size;
	#endif

	unsigned int key,radix,k;
	for(unsigned int j= 0; j< size;j++){
		// Get the key to count
		#ifdef TRANSPOSE
		k = groups * items * j + ig;
		#else
		k = j + start;
		#endif
		if(k >= n)
			return;
		key=d_Keys[k];   

		// Extract from the key the corresponding radix.
		// "key >> (pass * _BITS)" discards all the previously parsed data
		// and the comparation "& (_RADIX-1)" will return the radix in the
		// range (0 -> _RADIX-1)
		radix=(( key >> (pass * _BITS)) & (_RADIX-1));  

		// increment the local histogram of the radix
		loc_histo[radix * items + it ]++;
	}

	barrier(CLK_LOCAL_MEM_FENCE);  

	for(unsigned int ir=0;ir<_RADIX;ir++){
		d_Histograms[items * (ir * groups + gr) + it]=loc_histo[ir * items + it];
	}
  
	barrier(CLK_GLOBAL_MEM_FENCE);  
}

/** perform a parallel prefix sum (a scan) on the local histograms (see
 * Blelloch 1990), retrieving the accumulated histogram. Each workitem
 * worries about two memories.
 * See also http://http.developer.nvidia.com/GPUGems3/gpugems3_ch39.html
 * @param histo Local histograms, or number of ocurrences of each radix,
 * divided by blocks, as shown in histogram() method.
 * @param temp Local memory used to speed up the process.
 * @param globsum Total number of keys at each group (output).
 * @note In the first time that this kernel is called blocks of
 * _RADIX*_GROUPS*_ITEMS/_HISTOSPLIT accumulated histograms are generated
 * (_HISTOSPLIT number of them). Therefore _HISTOSPLIT global sums are
 * computed.
 * The second time that this kernel is called the accumulated histogram of
 * the _HISTOSPLIT global sums are computed, computing as well the total
 * sum of all the histograms (which is the number of keys), that will be
 * discarded.
 */
__kernel void scanhistograms( __global unsigned int* histo,__local unsigned int* temp,__global unsigned int* globsum)
{
	unsigned int it = get_local_id(0);
	unsigned int ig = get_global_id(0);
	unsigned int n  = get_local_size(0) * 2 ;
	unsigned int gr = get_group_id(0);
	unsigned int decale = 1; 

	temp[2*it]   = histo[2*ig];  
	temp[2*it+1] = histo[2*ig+1];  
 	
	// parallel prefix sum (algorithm of Blelloch 1990)
	// In each stage  (the first line is the input data):
	//   -# [h0, h1, h2, h3, ...]
	//   -# [h0, h0+h1, h2, h2+h3, ...]
	//   -# [h0, h0+h1, h2, h0+h1+h2+h3, ...]
	// Therefore at the end we have the global sum of the
	// local histograms in the last component.
	for (unsigned int d = n>>1; d > 0; d >>= 1){   
		barrier(CLK_LOCAL_MEM_FENCE);  
		if (it < d){  
			unsigned int ai = decale*(2*it+1)-1;  
			unsigned int bi = decale*(2*it+2)-1;  	
			temp[bi] += temp[ai];  
		}  
		decale *= 2; 
	}
  
	if (it == 0) {
		globsum[gr] = temp[n-1];
		temp[n - 1] = 0;
	}
                 
	// down sweep phase.
	// Now we perform the combination process, but in the
	// inverse order. For a example with 4 histograms, we
	// have the following process (the first line is the
	// input data):
	//   -# [h0, h0+h1, h2, 0]
	//   -# [h0, 0, h2, h0+h1]
	//   -# [0, h0, h0+h1, h0+h1+h2]
	// Which is the desired accumulated histogram.
	for (unsigned int d = 1; d < n; d *= 2){  
		decale >>= 1;  
		barrier(CLK_LOCAL_MEM_FENCE);

		if (it < d){  
			unsigned int ai = decale*(2*it+1)-1;  
			unsigned int bi = decale*(2*it+2)-1;  
			 
			unsigned int t = temp[ai];  
			temp[ai] = temp[bi];  
			temp[bi] += t;   
		}  
	}  
	barrier(CLK_LOCAL_MEM_FENCE);

	histo[2*ig] = temp[2*it];  
	histo[2*ig+1] = temp[2*it+1];  

	barrier(CLK_GLOBAL_MEM_FENCE);
}  

/** Use the _HISTOSPLIT accumulated global sums for updating the splited
 * accumulated histogram
 * @param histo Accumulated histogram. At the start the histogram is
 * reinitializated _HISTOSPLIT times.
 * @param globsum Accumulated global sums (with _HISTOSPLIT components).
 */
__kernel void pastehistograms( __global unsigned int* histo,__global unsigned int* globsum)
{
	unsigned int ig = get_global_id(0);
	unsigned int gr = get_group_id(0);
	unsigned int s;

	s=globsum[gr];

	histo[2*ig]   += s;  
	histo[2*ig+1] += s;  

	barrier(CLK_GLOBAL_MEM_FENCE);
}



/** Perform permutations using the accumulated histogram.
 * @param d_inKeys Input unsorted keys.
 * @param d_outKeys Output sorted keys.
 * @param d_Histograms Scanned histograms.
 * @param pass Pass of the radix decomposition.
 * @param d_inPermut Input permutations.
 * @param d_outPermut Output permutations.
 * @param loc_histo Histogram local memory to speed up the process.
 * @param n Number of keys.
 * @warning Radix sort needs several pass, so the output sorted keys of this
 * pass must be the input unsorted keys of the next pass. Don't forgive to swap
 * the OpenCL identifiers (for keys and permutations).
 * @note The output data from this kernel (keys and permutations) may need a
 * transposition.
 */
__kernel void reorder(const __global unsigned int* d_inKeys,
		      __global unsigned int* d_outKeys,
		      __global unsigned int* d_Histograms,
		      const unsigned int pass,
		      __global unsigned int* d_inPermut,
		      __global unsigned int* d_outPermut,
		      __local unsigned int* loc_histo,
		      const unsigned int n)
{
	unsigned int it = get_local_id(0);
	unsigned int ig = get_global_id(0);
	unsigned int gr = get_group_id(0);

	unsigned int groups = get_num_groups(0);
	unsigned int items  = get_local_size(0);

	// Set the keys analized by each thread
	unsigned int size = n/groups/items;
	#ifndef TRANSPOSE
	// If the data has not been transposed we must start reading from a
	// different place of ig
	unsigned int start= ig * (n/groups/items);
	#endif

	// take the accumulated histogram in the cache
	for(unsigned int ir=0;ir<_RADIX;ir++){
		loc_histo[ir*items + it] = d_Histograms[items*(ir*groups + gr) + it];
	}
	barrier(CLK_LOCAL_MEM_FENCE);  


	unsigned int newpos,key,radix,k,newpost;

	for(unsigned int j= 0; j< size;j++){
		// Get the key to sort
		#ifdef TRANSPOSE
		k = groups * items * j + ig;
		#else
		k = j + start;
		#endif
		if(k >= n)
			return;
		key   = d_inKeys[k];

		// Extract from the key the corresponding radix.
		// "key >> (pass * _BITS)" discards all the previously parsed data
		// and the comparation "& (_RADIX-1)" will return the radix in the
		// range (0 -> _RADIX-1)
		radix = ((key >> (pass * _BITS)) & (_RADIX-1)); 

		// Get the new position of the key from the histogram
		newpos = loc_histo[radix * items + it];

		// And transpose it (if proceed)
		#ifdef TRANSPOSE
		unsigned int ignew,jnew;
		ignew= newpos/(n/groups/items);
		jnew = newpos%(n/groups/items);
		newpost = jnew * (groups*items) + ignew;
		#else
		newpost=newpos;
		#endif

		// And set the new key
		d_outKeys[newpost]= key;  // killing line !!!

		#ifdef PERMUT 
		d_outPermut[newpost]=d_inPermut[k]; 
		#endif

		// The position is filled, modify the histogram to point to the
		// next available one
		newpos++;
		loc_histo[radix * items + it]=newpos;
	}
}

/** Compute the reversed permutations, which allows to know the original
 * position of a key from the sorted one.
 * @param d_directPermut Direct permutations (from the unsorted position to the
 * sorted one)
 * @param d_reversePermut Reverse permutations (from the sorted position to the
 * unsorted one)
 * @param n Number of keys.
 */
__kernel void reversePermutation(__global unsigned int* d_directPermut,
                                 __global unsigned int* d_reversePermut,
                                 unsigned int n)
{
	unsigned int i = get_global_id(0);
	if(i >= n)
		return;

	unsigned int permut = d_directPermut[i];
	d_reversePermut[permut] = i;
}

