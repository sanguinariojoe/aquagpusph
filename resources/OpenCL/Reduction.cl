/*
 *  This file is part of AQUA-gpusph, a free CFD program based on SPH.
 *  Copyright (C) 2012  Jose Luis Cercos Pita <jl.cercos@upm.es>
 *
 *  AQUA-gpusph is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  (at your option) any later version.
 *
 *  AQUA-gpusph is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with AQUA-gpusph.  If not, see <http://www.gnu.org/licenses/>.
 */

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

/** Reduction step. The objective of each step is obtain only one reduced value from each work group.
 * You can call this kernel recursively until only one work group will be computed.
 * @param input Input array where the reduced value is desired to be computed.
 * @param output Output array where reduced values will be stored.
 * @param N Number of input elements.
 * @param lmem local memory address array to store the output data while working.
 */
__kernel void Reduction( _g T *input, _g T *output,
                      unsigned int N,
                      _l T* lmem )
{
	unsigned int i;
	// Get the global index (to ensure not out of bounds reading operations)
	unsigned int gid = get_global_id(0);
	// Get id into the work group
	unsigned int tid = get_local_id(0);

	if(gid >= N)
		lmem[tid] = IDENTITY;	
	else
		lmem[tid] = input[gid];
	barrier(CLK_LOCAL_MEM_FENCE);

	// Start reducing the variables. The point is work in the first half of the
	// work group (that left to reduce), comparing it with the other half of data.
	// Since the number of threads into the work group can be obtained with the
	// get_local_size method, is better option use a defined variable because some
	// compilers can unroll the loop winning some performance.
	for(i=get_local_size(0)/2;i>0;i>>=1){
		// Ensure that we are not reading out of bounds
		if(tid < i)
			lmem[tid] = reduce(lmem[tid],lmem[tid+i]);
		barrier(CLK_LOCAL_MEM_FENCE);
	}

	if(tid == 0){
		output[get_group_id(0)] = lmem[0];
	}

	/*
	// Store the global data into the local memory address
	lmem[tid] = input[gid];
	barrier(CLK_LOCAL_MEM_FENCE);

	// Start reducing the variables. The point is work in the first half of the
	// work group (that left to reduce), comparing it with the other half of data.
	// Since the number of threads into the work group can be obtained with the
	// get_local_size method, is better option use a defined variable because some
	// compilers can unroll the loop winning some performance.
	for(i=get_local_size(0)/2;i>0;i>>=1){
		// Ensure that we are not reading out of bounds
		if(tid >= i)
			return;
		if(gid + i >= N)
			continue;
		lmem[tid] = reduce(lmem[tid],lmem[tid+i]);
		barrier(CLK_LOCAL_MEM_FENCE);
	}
	
	if(tid == 0) {
		output[get_group_id(0)] = lmem[0];
	}
	*/
}
