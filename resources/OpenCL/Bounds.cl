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

#ifndef HAVE_3D
    #include "types/2D.h"
#else
    #include "types/3D.h"
#endif

/** Filter the boundary particles/elements from the maximum position
 * computation. Filtering a particle imply to change its position by the
 * minimum known value such that any other position will be grater or equal
 * than its one.
 * @param output Output filtered positions.
 * @param imove Moving flag.
 * @param input Input positions to filter.
 * @param N Number of particles.
 */
__kernel void MaximumCoordsFilter(__global vec* output,
                                  __global int* imove,
                                  __global vec* input,
                                  unsigned int N)
{
	// find position in global arrays
	const unsigned int i = get_global_id(0);
	if(i >= N)
		return;

	// ---- | ------------------------ | ----
	// ---- V ---- Your code here ---- V ----

	if(imove[i] <= 0){
		#ifdef HAVE_3D
			output[i] = (vec)(-INFINITY,-INFINITY,-INFINITY,0.f);
		#else
			output[i] = (vec)(-INFINITY,-INFINITY);
		#endif
		return;
	}
	output[i] = input[i];

	// ---- A ---- Your code here ---- A ----
	// ---- | ------------------------ | ----
}

/** Filter the boundary particles/elements from the minimum position
 * computation. Filtering a particle imply to change its position by the
 * maximum known value such that any other position will be grater or equal
 * than its one.
 * @param output Output filtered positions.
 * @param imove Moving flag.
 * @param input Input positions to filter.
 * @param N Number of particles.
 */
__kernel void MinimumCoordsFilter(__global vec* output,
                                  __global int* imove,
                                  __global vec* input,
                                  unsigned int N)
{
	// find position in global arrays
	const unsigned int i = get_global_id(0);
	if(i >= N)
		return;

	// ---- | ------------------------ | ----
	// ---- V ---- Your code here ---- V ----

	if(imove[i] <= 0){
		#ifdef HAVE_3D
			output[i] = (vec)(INFINITY,INFINITY,INFINITY,0.f);
		#else
			output[i] = (vec)(INFINITY,INFINITY);
		#endif
		return;
	}
	output[i] = input[i];

	// ---- A ---- Your code here ---- A ----
	// ---- | ------------------------ | ----
}

/** This method discard the fixed particles for the maximum velocity computation.
 * @param output Output filtered data.
 * @param imove Particle moving flag.
 * @param input Input data to filter.
 * @param N Number of particles.
 */
__kernel void MaximumVelFilter(__global vec* output, __global int* imove, __global vec* input, unsigned int N)
{
	// find position in global arrays
	unsigned int i = get_global_id(0);
	if(i >= N)
		return;

	// ---- | ------------------------ | ----
	// ---- V ---- Your code here ---- V ----

	if(imove[i] <= 0){
		output[i] = VEC_ZERO;
		return;
	}
	output[i] = input[i];

	// ---- A ---- Your code here ---- A ----
	// ---- | ------------------------ | ----
}

/** This method discard the fixed particles for the minimum velocity computation.
 * @param output Output filtered data.
 * @param imove Particle moving flag.
 * @param input Input data to filter.
 * @param N Number of particles.
 */
__kernel void MinimumVelFilter(__global vec* output, __global int* imove, __global vec* input, unsigned int N)
{
	// find position in global arrays
	unsigned int i = get_global_id(0);
	if(i >= N)
		return;

	// ---- | ------------------------ | ----
	// ---- V ---- Your code here ---- V ----

	if(imove[i] <= 0){
		#ifdef HAVE_3D
			output[i] = (vec)(INFINITY,INFINITY,INFINITY,0.f);
		#else
			output[i] = (vec)(INFINITY,INFINITY);
		#endif
		return;
	}
	output[i] = input[i];

	// ---- A ---- Your code here ---- A ----
	// ---- | ------------------------ | ----
}
