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

/** @brief Tool to compute the fluid energy components.
 *
 * @param energy Particle resultant energy components:
 *   -# Internal energy: \f$ U = \int_0^t \sum_i \frac{p_i}{\rho_i^2}
     \left(
        \frac{\mathrm{d} \rho_i}{\mathrm{d} t}
        - \left. \frac{\mathrm{d} \rho_i}{\mathrm{d} t} \right\vert_F
     \right) m_i \mathrm{d}t \f$.
 *   -# Enthalpy: \f$ H = \int_0^t \sum_i \frac{p_i}{\rho_i^2}
     \frac{\mathrm{d} \rho_i}{\mathrm{d} t} m_i \mathrm{d}t \f$.
 *   -# Potential energy: \f$ E_{pot} = - \sum_i m_i
     \mathbf{g} \cdot \mathbf{r}_i \f$.
 *   -# Kinetic energy: \f$ E_{kin} = \sum_i \frac{1}{2} m_i
     \vert \mathbf{u}_i \vert^2 \f$.
 * @param imove Moving flags.
 *   - imove > 0 for regular fluid particles.
 *   - imove = 0 for sensors.
 *   - imove < 0 for boundary elements/particles.
 * @param ifluid Fluid index.
 * @param pos Position \f$ \mathbf{r} \f$.
 * @param v Velocity \f$ \mathbf{u} \f$.
 * @param mass Mass \f$ m \f$.
 * @param dens Density \f$ \rho \f$.
 * @param press Pressure \f$ p \f$.
 * @param drdt Density rate of change \f$ \frac{d \rho}{d t} \f$.
 * @param drdt_F Density rate of change restricted to the diffusive term
 * \f$ \left. \frac{d \rho}{d t} \right\vert_F \f$.
 * @param dvdt Velocity rate of change \f$ \frac{d \mathbf{u}}{d t} \f$.
 * @param refd Density of reference of the fluid \f$ \rho_0 \f$.
 * @param gamma Eq. of state exponent \f$ \gamma \f$.
 * @param cs Speed of sound \f$ c_s \f$.
 * @param grav Gravity acceleration \f$ \mathbf{g} \f$.
 * @param N Number of particles.
 */
__kernel void Energy(__global vec4* energy,
                     __global int* imove,
                     __global int* ifluid,
                     __global vec* pos,
                     __global vec* v,
                     __global float* mass,
                     __global float* dens,
                     __global float* press,
                     __global float* drdt,
                     __global float* drdt_F,
                     __global vec* dvdt,
                     __constant float* refd,
                     __constant float* gamma,
                     float cs,
                     vec grav,
                     unsigned int N)
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
