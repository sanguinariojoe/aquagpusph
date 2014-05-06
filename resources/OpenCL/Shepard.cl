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

#ifndef M_PI
	#define M_PI 3.14159265359f
#endif
#ifndef iM_PI
	#define iM_PI 0.318309886f
#endif

#ifndef uint
	#define uint unsigned int
#endif

#ifdef _g
	#error '_g' is already defined.
#endif
#define _g __global

#ifdef _c
	#error '_c' is already defined.
#endif
#define _c __constant

#ifdef _l
	#error '_l' is already defined.
#endif
#define _l __local

/** Apply the shepard term.
 * @param iMove Movement flags.
 * @param f Forces over particles.
 * @param drdt Density evolution of particles.
 * @param normal Particle normal (represents ferrand effect direction).
 * @param shepard Shepard term (0th correction).
 * @param N Number of particles.
 */
__kernel void Shepard(_g int* iMove, _g vec* f, _g float* drdt,
                      _g float* drdt_F, _g float* shepard, uint N )
{
	uint i = get_global_id(0);
	if(i >= N)
		return;
	// Don't compute vertexes or sensors
	if(iMove[i]<=0)
		return;

	// ---- | ------------------------ | ----
	// ---- V ---- Your code here ---- V ----

	float iShepard = shepard[i];
	// Shepard bounds (alone particles, or clamped zones)
	if( (iShepard < 0.1f) || (iShepard > 1.f) )
		return;
	#ifdef __FORCE_CORRECTION__
		f[i] /= iShepard;
	#endif
	#ifdef __DENS_CORRECTION__
		drdt[i] /= iShepard;
		drdt_F[i] /= iShepard;
	#endif

	// ---- A ---- Your code here ---- A ----
	// ---- | ------------------------ | ----

}

