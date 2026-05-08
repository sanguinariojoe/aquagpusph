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

/** @addtogroup cfd
 * @{
 */

/** @file
 * @brief Compute the fields on the inlet/outlet so they can transmit the waves
 * out the domain without reflecting them.
 *
 * This file is actually meant to be included by either Inlet.cl or Outlet.cl,
 * which are defining INWARD_NORMAL_SIGN before
 */

#ifndef INWARD_NORMAL_SIGN
#define INWARD_NORMAL_SIGN 1.f
#endif

#ifndef J_SHEPARD_LIMIT
#define J_SHEPARD_LIMIT FLT_EPSILON
#endif

#include "resources/Scripts/types/types.h"
#include "resources/Scripts/KernelFunctions/Kernel.h"

/** @brief Compute the characteristics of the outwards waves.
 * @param imove Moving flags.
 *   - imove > 0 for regular fluid particles.
 *   - imove = 0 for sensors.
 *   - imove < 0 for boundary elements/particles.
 * @param iset Set of particles index.
 * @param r Position \f$ \mathbf{r} \f$.
 * @param u Velocity \f$ \mathbf{u} \f$.
 * @param rho Density \f$ \rho \f$.
 * @param p Pressure \f$ p \f$.
 * @param j1 First characteristic \f$ J_1 \f$.
 * @param j2 Second characteristic \f$ J_2 \f$.
 * @param j3 Third characteristic \f$ J_3 \f$.
 * @param refd Density of reference of the fluid \f$ \rho_0 \f$.
 * @param N Number of particles.
 * @param dt Time step \f$ \Delta t \f$.
 * @param cs Speed of sound \f$ c_s \f$.
 * @param p0 Background pressure \f$ p_0 \f$.
 * @param g Gravity acceleration \f$ \mathbf{g} \f$.
 * @param io_r Lower corner of the inlet/outlet square.
 * @param io_n = Velocity direction.
 * @param io_U = Constant velocity magnitude.
 * @param io_rFS The point where the pressure is the reference one
 * (\f$ p_0 \f$).
 */
__kernel void
characteristics(const __global int* restrict imove,
                const __global unsigned int* restrict iset,
                const __global vec* restrict r,
                const __global vec* restrict u,
                const __global float* restrict rho,
                const __global float* restrict p,
                __global float* restrict j1,
                __global float* restrict j2,
                __global float* restrict j3,
                const __constant float* restrict refd,
                usize N,
                usize nbuffer,
                float dt,
                float cs,
                float p0,
                vec g,
                vec io_r,
                vec io_n,
                float io_U,
                vec io_rFS)
{
	const usize i = get_global_id(0);
	if (i >= N)
		return;
	if (imove[i] != 1)
		return;

	// Discard the particles at the inlet/outlet
	if (dot(r[i] - io_r, INWARD_NORMAL_SIGN * io_n) < 0.f)
		return;

	const float cs2 = cs * cs;
	const float un = dot(u[i], io_n);

	// Get the reference values
	const float uref = io_U;
	const float pref = refd[iset[i]] * dot(g, r[i] - io_rFS) + p0;
	const float rhoref = refd[iset[i]] + (p[i] - p0) / cs2;

	j1[i] = -cs2 * (rho[i] - rhoref) + p[i] - pref;
	j2[i] = rho[i] * cs * (un - uref) + p[i] - pref;
	j3[i] = -rho[i] * cs * (un - uref) + p[i] - pref;
}

/** @brief Extrapolate the characteristics to the particles at the
 * inlet/outlet.
 * @param imove Moving flags.
 *   - imove > 0 for regular fluid particles.
 *   - imove = 0 for sensors.
 *   - imove < 0 for boundary elements/particles.
 * @param r Position \f$ \mathbf{r} \f$.
 * @param rho Density \f$ \rho \f$.
 * @param m Mass \f$ m \f$.
 * @param jhoc Head and tail of chains for each cell.
 * @param j1 First characteristic \f$ J_1 \f$.
 * @param j2 Second characteristic \f$ J_2 \f$.
 * @param j3 Third characteristic \f$ J_3 \f$.
 * @param j1_tmp First characteristic (smoothed) \f$ J_1 \f$.
 * @param j2_tmp Second characteristic (smoothed) \f$ J_2 \f$.
 * @param j3_tmp Third characteristic (smoothed) \f$ J_3 \f$.
 * @param shepard Shepard renormalization factor \f$ \gamma \f$.
 * @param io_r Lower corner of the inlet/outlet square.
 * @param io_n = Velocity direction.
 * @param N Number of particles.
 */
__kernel void
extrapolate(const __global int* restrict imove,
            const __global vec* restrict r,
            const __global float* restrict rho,
            const __global float* restrict m,
            const __global svec2* restrict jhoc,
            const __global float* restrict j1,
            const __global float* restrict j2,
            const __global float* restrict j3,
            __global float* restrict j1_tmp,
            __global float* restrict j2_tmp,
            __global float* restrict j3_tmp,
            __global float* restrict shepard,
            usize N,
            vec io_r,
            vec io_n)
{
	const usize i = get_global_id(0);
	const usize it = get_local_id(0);
	if (i >= N)
		return;
	if (imove[i] != 1)
		return;

	const vec_xyz r_i = r[i].XYZ;

	// Discard the particles that already passed through the inlet/outlet
	if (dot(r_i - io_r.XYZ, INWARD_NORMAL_SIGN * io_n.XYZ) > 0.f)
		return;

	__private float __j1 = 0.f;
	__private float __j2 = 0.f;
	__private float __j3 = 0.f;
	__private float __shepard = 0.f;

	FOR_NEIGHS(N, jhoc)
	{
		if (imove[j] != 1)
			continue;
		if (dot(r[j] - io_r, INWARD_NORMAL_SIGN * io_n) <= 0.f) {
			// Do not use other inlet/outlet particles to interpolate
			continue;
		}

		const vec_xyz r_ij = r[j].XYZ - r_i;
		const float q = length(r_ij) / H;
		if (q >= SUPPORT) {
			continue;
		}
		{
			const float w_ij = kernelW(q) * CONW * m[j] / rho[j];

			__shepard += w_ij;
			__j1 += j1[j] * w_ij;
			__j2 += j2[j] * w_ij;
			__j3 += j3[j] * w_ij;
		}
	}
	END_FOR_NEIGHS()

	const float div = __shepard > J_SHEPARD_LIMIT ? 1.f / __shepard : 1.f;
	j1_tmp[i] = __j1 * div;
	j2_tmp[i] = __j2 * div;
	j3_tmp[i] = __j3 * div;
	shepard[i] = __shepard;
}

/** @brief Set the field values at the inlet
 * @param imove Moving flags.
 *   - imove > 0 for regular fluid particles.
 *   - imove = 0 for sensors.
 *   - imove < 0 for boundary elements/particles.
 * @param iset Set of particles index.
 * @param r Position \f$ \mathbf{r} \f$.
 * @param u Velocity \f$ \mathbf{u} \f$.
 * @param rho Density \f$ \rho \f$.
 * @param p Pressure \f$ p \f$.
 * @param j1 First characteristic \f$ J_1 \f$.
 * @param j2 Second characteristic \f$ J_2 \f$.
 * @param j3 Third characteristic \f$ J_3 \f$.
 * @param refd Density of reference of the fluid \f$ \rho_0 \f$.
 * @param N Number of particles.
 * @param dt Time step \f$ \Delta t \f$.
 * @param cs Speed of sound \f$ c_s \f$.
 * @param p0 Background pressure \f$ p_0 \f$.
 * @param g Gravity acceleration \f$ \mathbf{g} \f$.
 * @param io_r Lower corner of the inlet/outlet square.
 * @param io_n = Velocity direction.
 * @param io_U = Constant velocity magnitude.
 * @param io_rFS The point where the pressure is the reference one
 * (\f$ p_0 \f$).
 */
__kernel void
values(const __global int* restrict imove,
       const __global unsigned int* restrict iset,
       const __global vec* restrict r,
       __global vec* restrict u,
       __global float* restrict rho,
       __global float* restrict p,
       const __global float* restrict j1,
       const __global float* restrict j2,
       const __global float* restrict j3,
       const __constant float* restrict refd,
       usize N,
       float dt,
       float cs,
       float p0,
       vec g,
       vec io_r,
       vec io_n,
       float io_U,
       vec io_rFS)
{
	const usize i = get_global_id(0);
	if (i >= N)
		return;
	if (imove[i] != 1)
		return;

	// Discard the particles that already passed through the inlet/outlet
	if (dot(r[i] - io_r, INWARD_NORMAL_SIGN * io_n) > 0.f)
		return;

	const float cs2 = cs * cs;

	// Get the reference values
	const float uref = io_U;
	const float pref = refd[iset[i]] * dot(g, r[i] - io_rFS) + p0;
	const float rhoref = refd[iset[i]] + (p[i] - p0) / cs2;

	rho[i] = rhoref + 1.f / cs2 * (-j1[i] + 0.5f * j2[i] + 0.5f * j3[i]);
	u[i] = (uref + 1.f / (2.f * rho[i] * cs) * (j2[i] - j3[i])) * io_n;
	p[i] = pref + 0.5f * (j2[i] + j3[i]);
}

/*
 * @}
 */
