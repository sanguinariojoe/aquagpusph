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
 * @brief OpenCL kernels to compute the fluid bounds box, as well as the
 * minimum and maximum velocities.
 * (See Aqua::CalcServer::Bounds for details)
 */

#ifndef HAVE_3D
    #include "types/2D.h"
#else
    #include "types/3D.h"
#endif

/** @brief Filter out the non fluid particles (boundaries or sensors) from the
 * maximum position computation.
 *
 * Filtering out a particle imply to change its position by the minimum known
 * value such that any other position will be grater or equal than its one.
 *
 * @param output Output filtered positions.
 * @param imove Moving flags.
 *   - imove > 0 for regular fluid particles.
 *   - imove = 0 for sensors.
 *   - imove < 0 for boundary elements/particles.
 * @param input Input positions to be filtered.
 * @param N Total number of particles.
 */
__kernel void MaximumCoordsFilter(__global vec* output,
                                  __global int* imove,
                                  __global vec* input,
                                  unsigned int N)
{
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

/** @brief Filter out the non fluid particles (boundaries or sensors) from the
 * minimum position computation.
 *
 * Filtering out a particle imply to change its position by the maximum known
 * value such that any other position will be lower or equal than its one.
 *
 * @param output Output filtered positions.
 * @param imove Moving flags.
 *   - imove > 0 for regular fluid particles.
 *   - imove = 0 for sensors.
 *   - imove < 0 for boundary elements/particles.
 * @param input Input positions to be filtered.
 * @param N Total number of particles.
 */
__kernel void MinimumCoordsFilter(__global vec* output,
                                  __global int* imove,
                                  __global vec* input,
                                  unsigned int N)
{
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

/** @brief Filter out the non fluid particles (boundaries or sensors) from the
 * maximum velocity computation.
 *
 * Filtering out a particle imply to change its velocity by the maximum known
 * value such that any other velocity will be lower or equal than its one.
 *
 * @param output Output filtered velocities.
 * @param imove Moving flags.
 *   - imove > 0 for regular fluid particles.
 *   - imove = 0 for sensors.
 *   - imove < 0 for boundary elements/particles.
 * @param input Input velocities to be filtered.
 * @param N Total number of particles.
 */
__kernel void MaximumVelFilter(__global vec* output,
                               __global int* imove,
                               __global vec* input,
                               unsigned int N)
{
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

/** @brief Filter out the non fluid particles (boundaries or sensors) from the
 * minimum velocity computation.
 *
 * Filtering out a particle imply to change its velocity by the minimum known
 * value such that any other velocity will be lower or equal than its one.
 *
 * @param output Output filtered velocities.
 * @param imove Moving flags.
 *   - imove > 0 for regular fluid particles.
 *   - imove = 0 for sensors.
 *   - imove < 0 for boundary elements/particles.
 * @param input Input velocities to be filtered.
 * @param N Total number of particles.
 */
__kernel void MinimumVelFilter(__global vec* output,
                               __global int* imove,
                               __global vec* input,
                               unsigned int N)
{
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
