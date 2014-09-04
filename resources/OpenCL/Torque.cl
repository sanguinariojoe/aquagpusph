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
 * @brief OpenCL kernels to compute the fluid global force and moment.
 * (See Aqua::CalcServer::Torque for details)
 */

#ifndef HAVE_3D
    #include "types/2D.h"
#else
    #include "types/3D.h"
#endif

/** @brief Tool to compute the global fluid force and moment.
 * @param moment Moment to be computed respect to @paramname{cor}
 * [N \f$ \cdot \f$ m].
 * @param force Force to be computed [N].
 * @param imove Moving flags.
 *   - imove > 0 for regular fluid particles.
 *   - imove = 0 for sensors.
 *   - imove < 0 for boundary elements/particles.
 * @param ifluid Fluid index.
 * @param pos Position \f$ \mathbf{r} \f$.
 * @param dvdt Velocity rate of change \f$ \frac{d \mathbf{u}}{d t} \f$.
 * @param mass Mass \f$ m \f$.
 * @param dens Density \f$ \rho \f$.
 * @param refd Density of reference of the fluid \f$ \rho_0 \f$.
 * @param N Number of particles.
 * @param grav Gravity acceleration \f$ \mathbf{g} \f$.
 * @param cor Center of rotation (point with respect the moment is computed).
 */
__kernel void Torque(__global vec* moment,
                     __global vec* force,
                     __global int* imove,
                     __global int* ifluid,
                     __global vec* pos,
                     __global vec* dvdt,
                     __global float* mass,
                     __global float* dens,
                     __global float* refd,
                     unsigned int N,
                     vec grav,
                     vec cor)
{
	// find position in global arrays
	unsigned int i = get_global_id(0);
	if(i >= N)
		return;
	if(imove[i] <= 0){
		moment[i] = VEC_ZERO;
		force[i] = VEC_ZERO;
		return;
	}

	// ---- | ------------------------ | ----
	// ---- V ---- Your code here ---- V ----

	vec arm  = pos[i] - cor;
	vec acc  = -dvdt[i];
	// float m  = mass[i]/dens[i]*refd[ifluid[i]]; // Volume based computation
	float m  = mass[i];                            // Mass based computation
	force[i] = m * acc;
	#ifndef HAVE_3D
		moment[i].x = m * (arm.y * acc.x - arm.x * acc.y);
		moment[i].y = m * (arm.y * acc.x - arm.x * acc.y);
	#else
		moment[i].x = m * (arm.y * acc.z - arm.z * acc.y);
		moment[i].y = m * (arm.z * acc.x - arm.x * acc.z);
		moment[i].z = m * (arm.x * acc.y - arm.y * acc.x);
		moment[i].w = 0.f;
	#endif

	// ---- A ---- Your code here ---- A ----
	// ---- | ------------------------ | ----
}
