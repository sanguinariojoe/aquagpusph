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
 * @brief Particles out of domain filter.
 * (See Aqua::CalcServer::Domain for details)
 */

#ifndef HAVE_3D
    #include "types/2D.h"
#else
    #include "types/3D.h"
#endif

/** @brief Check and destroy the particles out of the domain.
 * @param imove Moving flags.
 *   - imove > 0 for regular fluid particles.
 *   - imove = 0 for sensors.
 *   - imove < 0 for boundary elements/particles.
 * @param pos Position \f$ \mathbf{r} \f$.
 * @param v Velocity \f$ \mathbf{u} \f$.
 * @param dvdt Velocity rate of change \f$ \frac{d \mathbf{u}}{d t} \f$.
 * @param mass Mass \f$ m \f$.
 * @param N Number of particles.
 * @param min_bound Minimum position where a particle can be placed.
 * @param max_bound Maximum position where a particle can be placed.
 */
__kernel void Domain(__global int* imove,
                     __global vec* pos,
                     __global vec* v,
                     __global vec* dvdt,
                     __global float* mass,
                     unsigned int N,
                     vec min_bound,
                     vec max_bound)
{
	// find position in global arrays
	unsigned int i = get_global_id(0);
	if(i >= N)
		return;

	// ---- | ------------------------ | ----
	// ---- V ---- Your code here ---- V ----

	vec coords = pos[i];
	if(    (coords.x < min_bound.x)
	    || (coords.y < min_bound.y)
	    || (coords.x > max_bound.x)
	    || (coords.y > max_bound.y)
	    #ifdef HAVE_3D
	    || (coords.z < min_bound.z)
	    || (coords.z > max_bound.z)
	    #endif
	  )
	{
		// Set as fixed zero mass particle (sensor equivalent)
		imove[i] = 0;
		mass[i]  = 0.f;
		// Stop the particle
		v[i] = VEC_ZERO;
		dvdt[i] = VEC_ZERO;
		// Clamp the position
		pos[i].x = max(pos[i].x, min_bound.x);
		pos[i].x = min(pos[i].x, max_bound.x);
		pos[i].y = max(pos[i].y, min_bound.y);
		pos[i].y = min(pos[i].y, max_bound.y);
		#ifdef HAVE_3D
			pos[i].z = max(pos[i].z, min_bound.z);
			pos[i].z = min(pos[i].z, max_bound.z);
		#endif
	}

	// ---- A ---- Your code here ---- A ----
	// ---- | ------------------------ | ----
}
