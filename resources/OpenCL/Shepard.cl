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
 * @brief Acceleration and density rate of change renormalization.
 * (See Aqua::CalcServer::Shepard for details)
 */

#ifndef HAVE_3D
    #include "types/2D.h"
#else
    #include "types/3D.h"
#endif

/** @brief Apply the Shepard renormalization term.
 * @param imove Moving flags.
 *   - imove > 0 for regular fluid particles.
 *   - imove = 0 for sensors.
 *   - imove < 0 for boundary elements/particles.
 * @param dvdt Velocity rate of change \f$ \frac{d \mathbf{u}}{d t} \f$.
 * @param drdt Density rate of change \f$ \frac{d \rho}{d t} \f$.
 * @param drdt_F Density rate of change restricted to the diffusive term
 * \f$ \left. \frac{d \rho}{d t} \right\vert_F \f$.
 * @param shepard Shepard term
 * \f$ \gamma(\mathbf{x}) = \int_{\Omega}
 *     W(\mathbf{y} - \mathbf{x}) \mathrm{d}\mathbf{x} \f$.
 * @param N Number of particles.
 * @see Aqua::CalcServer::Shepard
 * @see Aqua::CalcServer::Rates
 */
__kernel void Shepard(__global int* imove,
                      __global vec* dvdt,
                      __global float* drdt,
                      __global float* drdt_F,
                      __global float* shepard,
                      unsigned int N)
{
	unsigned int i = get_global_id(0);
	if(i >= N)
		return;
	// Don't compute vertexes or sensors
	if(imove[i]<=0)
		return;

	// ---- | ------------------------ | ----
	// ---- V ---- Your code here ---- V ----

	float iShepard = shepard[i];
	// Shepard bounds (alone particles, or clamped zones)
	if( (iShepard < 0.1f) || (iShepard > 1.f) )
		return;
	#ifdef __FORCE_CORRECTION__
		dvdt[i] /= iShepard;
	#endif
	#ifdef __DENS_CORRECTION__
		drdt[i] /= iShepard;
		drdt_F[i] /= iShepard;
	#endif

	// ---- A ---- Your code here ---- A ----
	// ---- | ------------------------ | ----

}

