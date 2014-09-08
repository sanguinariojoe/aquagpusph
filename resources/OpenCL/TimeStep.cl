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
 * @brief Time step computation.
 * (See Aqua::CalcServer::TimeStep for details)
 */

#ifndef HAVE_3D
    #include "types/2D.h"
#else
    #include "types/3D.h"
#endif

/** @brief Compute the maximum valid time step for each particle.
 * @param dtconv Convective time step term.
 * @param v Velocity \f$ \mathbf{u} \f$.
 * @param dvdt Velocity rate of change \f$ \frac{d \mathbf{u}}{d t} \f$.
 * @param N Number of particles.
 * @param dt Time step \f$ \Delta t \f$.
 * @param cs Speed of sound \f$ c_s \f$.
 */
__kernel void TimeStep(__global float* dtconv,
                       __global vec* v,
                       __global vec* dvdt,
                       unsigned int N,
                       float dt,
                       float cs)
{
	// find position in global arrays
	unsigned int i = get_global_id(0);
	if(i >= N)
		return;

	// ---- | ------------------------ | ----
	// ---- V ---- Your code here ---- V ----

	const float vv = fast_length(v[i] + dvdt[i] * dt);
	dtconv[i] = h / max(10.f * vv, cs);

	// ---- A ---- Your code here ---- A ----
	// ---- | ------------------------ | ----

}
