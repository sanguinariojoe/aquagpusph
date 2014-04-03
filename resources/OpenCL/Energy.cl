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
 * @param energy Particle resultant energy components:
 *     -# Mechanical energy: \f$ E = E_{pot} + E_{elas} + E_{kin} \f$
 *     -# Potential energy: \f$ E_{pot} = - m \mathbf{g} \cdot \mathbf{r} \f$
 *     -# Elastic energy: \f$ \frac{\mathrm{d} E_{elas}}{\mathrm{d} t} = 
 *                        m \frac{p - \rho \mathbf{g} \cdot \mathbf{r}}{\rho^2}
 *                        \frac{\mathrm{d} \rho}{\mathrm{d} t}\f$
 *     -# Kinetic energy: \f$ E_{kin} = \frac{1}{2} m \vert v \vert^2 \f$
 * @param imove Particle moving flag.
 * @param ifluid Particle fluid.
 * @param pos Position of particles.
 * @param v Particle velocity.
 * @param mass Particle mass.
 * @param dens Particle density.
 * @param press Particle pressure.
 * @param press Speed of sound.
 * @param drdt Density rate of change.
 * @param drdt_F Density rate of change (restricted to the diffusive term).
 * @param f Velocity rate of change.
 * @param refd Density of reference of the fluid.
 * @param gamma Eq. of state exponent.
 * @param cs Speed of sound.
 * @param grav Gravity acceleration.
 * @param N Number of particles.
 * @remarks Since the elastic energy should be integrated in time, it will
 * be excluded from the total energy computation for each particle.
 */
__kernel void Energy(_g vec4* energy, _g int* imove, _g int* ifluid,
                     _g vec* pos, _g vec* v, _g float* mass, _g float* dens,
                     _g float* press, _g float* drdt, _g float* drdt_F,
                     _g vec* f, _c float* refd, _c float* gamma,
                     float cs, vec grav, unsigned int N)
{
	// find position in global arrays
	unsigned int i = get_global_id(0);
	if(i >= N)
		return;
	if(imove[i] <= 0){
		energy[i] = (vec4)(0.f,0.f,0.f,0.f);
		return;
	}

	// ---- | ------------------------ | ----
	// ---- V ---- Your code here ---- V ----

	float m = mass[i];
	float d = dens[i];
	float p = press[i];
	float mp_dd = m * p/(d*d);
	// U = Internal energy (viscosity effect not implemented yet)
	energy[i].x = mp_dd * drdt[i];
	// H = Enthalpy
	energy[i].y = mp_dd * (drdt[i] + drdt_F[i]);
	// Epot = Potential energy
	energy[i].z = - m * dot(grav, pos[i]);
	// Ekin = Kinetic energy
	energy[i].w = 0.5f * m * dot(v[i], v[i]);

	// ---- A ---- Your code here ---- A ----
	// ---- | ------------------------ | ----
}
