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

/** Method called outside to do the predictor phase.
 * @param imove Fix particles flag.
 * @param pos Particle position.
 * @param v Particle velocity.
 * @param f Particle force.
 * @param mass Particle mass
 * @param N Number of particles.
 * @param minBound Minimum position that a particle can take.
 * @param maxBound Maximum position that a particle can take.
 */
__kernel void Domain(_g int* imove, _g vec* pos, _g vec* v, _g vec* f, _g float* mass,
                        unsigned int N, vec minBound, vec maxBound)
{
	// find position in global arrays
	unsigned int i = get_global_id(0);
	if(i >= N)
		return;

	// ---- | ------------------------ | ----
	// ---- V ---- Your code here ---- V ----

	vec coords = pos[i];
	if(    (coords.x < minBound.x)
	    || (coords.y < minBound.y)
	    || (coords.x > maxBound.x)
	    || (coords.y > maxBound.y)
	    #ifdef HAVE_3D
	    || (coords.z < minBound.z)
	    || (coords.z > maxBound.z)
	    #endif
	  )
	{
		// Set as fixed zero mass particle (sensor equivalent)
		imove[i] = 0;
		mass[i]  = 0.f;
		// Stop the particle
		v[i] = VEC_ZERO;
		f[i] = VEC_ZERO;
		// Clamp the position
		pos[i].x = max(pos[i].x, minBound.x);
		pos[i].x = min(pos[i].x, maxBound.x);
		pos[i].y = max(pos[i].y, minBound.y);
		pos[i].y = min(pos[i].y, maxBound.y);
		#ifdef HAVE_3D
			pos[i].z = max(pos[i].z, minBound.z);
			pos[i].z = min(pos[i].z, maxBound.z);
		#endif
	}

	// ---- A ---- Your code here ---- A ----
	// ---- | ------------------------ | ----
}
