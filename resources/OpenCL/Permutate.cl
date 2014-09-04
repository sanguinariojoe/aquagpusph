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
 * @brief Particles sorting/unsorting tool.
 * (See Aqua::CalcServer::Permutate for details)
 */

#ifndef HAVE_3D
    #include "types/2D.h"
#else
    #include "types/3D.h"
#endif

/** Permutate the data arrays from the unsorted space to the sorted one, and
 * viceversa.
 * @param ifluid permuted fluid index.
 * @param ifluid_in Fluid index.
 * @param imove Permuted moving flags.
 *   - imove > 0 for regular fluid particles.
 *   - imove = 0 for sensors.
 *   - imove < 0 for boundary elements/particles.
 * @param imove_in Moving flags.
 *   - imove > 0 for regular fluid particles.
 *   - imove = 0 for sensors.
 *   - imove < 0 for boundary elements/particles.
 * @param pos Permuted position \f$ \mathbf{r} \f$.
 * @param pos_in Position \f$ \mathbf{r} \f$.
 * @param normal Permuted normal \f$ \mathbf{n} \f$.
 * @param normal_in Normal \f$ \mathbf{n} \f$.
 * @param v Permuted velocity \f$ \mathbf{u} \f$.
 * @param v_in Velocity \f$ \mathbf{u} \f$.
 * @param dvdt Permuted velocity rate of change \f$ \frac{d \mathbf{u}}{d t} \f$.
 * @param dvdt_in Velocity rate of change \f$ \frac{d \mathbf{u}}{d t} \f$.
 * @param dens Permuted density \f$ \rho \f$.
 * @param dens_in Density \f$ \rho \f$.
 * @param press Permuted pressure \f$ p \f$.
 * @param press_in Pressure \f$ p \f$.
 * @param mass Permuted mass \f$ m \f$.
 * @param mass_in Mass \f$ m \f$.
 * @param drdt Permuted density rate of change \f$ \frac{d \rho}{d t} \f$.
 * @param drdt_in Density rate of change \f$ \frac{d \rho}{d t} \f$.
 * @param drdt_F Permuted density rate of change restricted to the diffusive
 * term \f$ \left. \frac{d \rho}{d t} \right\vert_F \f$.
 * @param drdt_F_in Density rate of change restricted to the diffusive term
 * \f$ \left. \frac{d \rho}{d t} \right\vert_F \f$.
 * @param shepard Permutated shepard term
 * \f$ \gamma(\mathbf{x}) = \int_{\Omega}
 *     W(\mathbf{y} - \mathbf{x}) \mathrm{d}\mathbf{x} \f$.
 * @param shepard_in Shepard term
 * \f$ \gamma(\mathbf{x}) = \int_{\Omega}
 *     W(\mathbf{y} - \mathbf{x}) \mathrm{d}\mathbf{x} \f$.
 * @param permut Permutations array.
 * @param N Number of particles.
 */
__kernel void Permutate(__global int* ifluid, __global int* ifluid_in,
                        __global int* imove, __global int* imove_in,
                        __global vec* pos, __global vec* pos_in,
                        __global vec* normal, __global vec* normal_in,
                        __global vec* v, __global vec* v_in,
                        __global vec* dvdt, __global vec* dvdt_in,
                        __global float* dens, __global float* dens_in,
                        __global float* press, __global float* press_in,
                        __global float* mass, __global float* mass_in,
                        __global float* drdt, __global float* drdt_in,
                        __global float* drdt_F, __global float* drdt_F_in,
                        __global float* shepard, __global float* shepard_in,
                        __global unsigned int *permut,
                        unsigned int N)
{
	unsigned int i = get_global_id(0);
	if(i >= N)
		return;

	// ---- | ------------------------ | ----
	// ---- V ---- Your code here ---- V ----

	// We assume i in the unsorted space, and i_out in the sorted one 
	const unsigned int i_out = permut[i];

	ifluid[i_out] = ifluid_in[i];
	imove[i_out] = imove_in[i];
	pos[i_out] = pos_in[i];
	normal[i_out] = normal_in[i];
	v[i_out] = v_in[i];
	dvdt[i_out] = dvdt_in[i];
	dens[i_out] = dens_in[i];
	press[i_out] = press_in[i];
	mass[i_out] = mass_in[i];
	drdt[i_out] = drdt_in[i];
	drdt_F[i_out] = drdt_F_in[i];
	shepard[i_out] = shepard_in[i];

	// ---- A ---- Your code here ---- A ----
	// ---- | ------------------------ | ----
}
