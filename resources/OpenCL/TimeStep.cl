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

#ifndef M_PI
	#define M_PI 3,14159265359
#endif

#ifdef _g
	#error '_g' is already defined.
#endif
#define _g __global

#ifdef _c
	#error '_c' is already defined.
#endif
#define _c __constant

/** Compute the timestep requested by each particle.
 * @param dtconv Convection time step term
 * @param v Velocity of particles.
 * @param f Forces over particles.
 * @param hp Kernel height of particles.
 * @param N Number of particles.
 * @param dt Time step.
 * @param cs Sound speed
 */
__kernel void TimeStep(_g float* dtconv, _g vec* v, _g vec* f,
                      unsigned int N, float dt, float cs)
{
	// find position in global arrays
	unsigned int i = get_global_id(0);
	if(i >= N)
		return;

	// ---- | ------------------------ | ----
	// ---- V ---- Your code here ---- V ----

	const float vv = fast_length(v[i] + f[i] * dt);
	dtconv[i] = h / max(10.f * vv, cs);

	// ---- A ---- Your code here ---- A ----
	// ---- | ------------------------ | ----

}
