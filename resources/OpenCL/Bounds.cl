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

#ifndef M_PI
	#define M_PI 3.14159265359
#endif

#ifdef _g
	#error '_g' is already defined.
#endif
#define _g __global

#ifdef _c
	#error '_c' is already defined.
#endif
#define _c __constant

/** This method discard the fixed particles for the maximum coordinates computation.
 * @param output Output filtered data.
 * @param imove Particle moving flag.
 * @param input Input data to filter.
 * @param N Number of particles.
 */
__kernel void MaximumCoordsFilter(_g vec* output, _g int* imove, _g vec* input, unsigned int N)
{
	// find position in global arrays
	unsigned int i = get_global_id(0);
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

/** This method discard the fixed particles for the minimum coordinates computation.
 * @param output Output filtered data.
 * @param imove Particle moving flag.
 * @param input Input data to filter.
 * @param N Number of particles.
 */
__kernel void MinimumCoordsFilter(_g vec* output, _g int* imove, _g vec* input, unsigned int N)
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

/** This method discard the fixed particles for the maximum velocity computation.
 * @param output Output filtered data.
 * @param imove Particle moving flag.
 * @param input Input data to filter.
 * @param N Number of particles.
 */
__kernel void MaximumVelFilter(_g vec* output, _g int* imove, _g vec* input, unsigned int N)
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
__kernel void MinimumVelFilter(_g vec* output, _g int* imove, _g vec* input, unsigned int N)
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
