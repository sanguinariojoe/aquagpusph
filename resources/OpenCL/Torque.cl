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

/** Method called outside to do the predictor phase.
 * @param torque Particle resultant moment.
 * @param force Particle resultant force.
 * @param imove Particle moving flag.
 * @param ifluid Particle fluid.
 * @param pos Position of particles.
 * @param f Particle acceleration.
 * @param mass Particle mass.
 * @param dens Particle density.
 * @param refd Fluids reference density.
 * @param N Number of particles.
 * @param cor Center of rotation.
 */
__kernel void Torque(_g vec* torque, _g vec* force, _g int* imove, _g int* ifluid, _g vec* pos, _g vec* f,
                     _g float* mass, _g float* dens, _g float* refd, unsigned int N, vec g, vec cor)
{
	// find position in global arrays
	unsigned int i = get_global_id(0);
	if(i >= N)
		return;
	if(imove[i] <= 0){
		torque[i] = VEC_ZERO;
		force[i] = VEC_ZERO;
		return;
	}

	// ---- | ------------------------ | ----
	// ---- V ---- Your code here ---- V ----

	vec arm  = pos[i] - cor;
	vec acc  = -f[i];
	// float m  = mass[i]/dens[i]*refd[ifluid[i]]; // Volume based computation
	float m  = mass[i];                            // Mass based computation
	force[i] = m*acc;
	#ifndef HAVE_3D
		torque[i].x = m*(arm.y * acc.x - arm.x * acc.y);
		torque[i].y = m*(arm.y * acc.x - arm.x * acc.y);
	#else
		torque[i].x = m*(arm.y * acc.z - arm.z * acc.y);
		torque[i].y = m*(arm.z * acc.x - arm.x * acc.z);
		torque[i].z = m*(arm.x * acc.y - arm.y * acc.x);
		torque[i].w = 0.f;
	#endif

	// ---- A ---- Your code here ---- A ----
	// ---- | ------------------------ | ----
}
