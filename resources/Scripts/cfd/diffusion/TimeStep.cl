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
 * @brief Variable time step computation.
 */

#include "resources/Scripts/types/types.h"
#include "resources/Scripts/KernelFunctions/Kernel.h"
#include "resources/Scripts/cfd/ideal_gas/sound_speed.hcl"

/** @brief Compute the maximum time step for each particle.
 *
 * In SPH the time step is selected to enforce the particles may not move more
 * than \f$ 0.1 h \f$, where the Courant factor is not taken into account yet.
 *
 * Along this line, the distance moved by a particle can be written as follows:
 *
 * \f$ \vert \mathbf{r}_{n+1} - \mathbf{r}_{n} \vert =
 *     \vert \mathbf{u} \vert \Delta t +
 *     \frac{1}{2} \left\vert
 *                     \frac{\mathrm{d} \mathbf{u}}{\mathrm{d} t}
 *                 \right\vert {\Delta t}^2 +
       \mathcal{O}({\Delta t}^3) \f$
 *
 * Such that, taking maximums, and rearraging the equation:
 *
 * \f$ \Delta t = \frac{1}{20} \min \left(
 *     \frac{h}{\vert \mathbf{u} \vert},
 *     \sqrt{\frac{2 h}{\left\vert
 *                          \frac{\mathrm{d} \mathbf{u}}{\mathrm{d} t}
 *                      \right\vert}}
 * \right) \f$
 *
 * @param dt_var Variable time step \f$ \mathrm{min} \left(
 * C_f \frac{h}{c_s}, C_f \frac{h}{10 \vert \mathbf{u} \vert}\right)\f$.
 * @param u Velocity \f$ \mathbf{u}_{n+1/2} \f$.
 * @param dudt Velocity rate of change \f$ \frac{d \mathbf{u}}{d t} \f$.
 * @param N Number of particles.
 * @param dt Fixed time step \f$ \Delta t = C_f \frac{h}{c_s} \f$.
 * @param dt_min Minimum time step \f$ \Delta t_{\mathrm{min}} \f$.
 * @param courant Courant factor \f$ C_f \f$.
 * @param h Kernel characteristic length \f$ h \f$.
 * @param div_u divergence of u * rho
 * @param grad_p grad of p / rho
 * @param gamma politropic coeficient
 * @param lambda thermal conductivity
 * @param C_p heat t constat pressure
 * @param D_xx diffusion coeficcient
 */
/*float
givemax(const local vec16 ddss)
{
	return max(max(max(max(max(max(max(max(max(max(max(max(max(max(max(ddss.s0,
	                                                                   ddss.s1),
	                                                               ddss.s2),
	                                                           ddss.s3),
	                                                       ddss.s4),
	                                                   ddss.s5),
	                                               ddss.s6),
	                                           ddss.s7),
	                                       ddss.s8),
	                                   ddss.s9),
	                               ddss.sA),
	                           ddss.sB),
	                       ddss.sC),
	                   ddss.sD),
	               ddss.sE),
	           ddss.sF);
}
*/

float
givemax(const local vec16 ddss)
{

const vec8 tmp1 = max(dt_u7.s01234567, dt_u7.s89ABCDEF);
const vec4 tmp2 = max(tmp1.s0123, tmp1.s4567);
const vec2 tmp3 = max(tmp2.s01, tmp1.s23);
const float tmp4 = max(tmp2.s0, tmp1.s1);

return tmp4;
}

/*
float
givemax(const local vec16 ddss)
{

	float all[16] = { 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f,
		              0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f };

	all[0] = ddss.s0;
	all[1] = ddss.s1;
	all[2] = ddss.s2;
	all[3] = ddss.s3;
	all[4] = ddss.s4;
	all[5] = ddss.s5;
	all[6] = ddss.s6;
	all[7] = ddss.s7;
	all[8] = ddss.s8;
	all[9] = ddss.s9;
	all[10] = ddss.sA;
	all[11] = ddss.sB;
	all[12] = ddss.sC;
	all[13] = ddss.sD;
	all[14] = ddss.sE;
	all[15] = ddss.sF;

	return max(all);
}
*/

__kernel void
entry(__global float* dt_var,
      const __global int* imove,
      const __global unsigned int* iset,
      const __global vec* u,
      const __global vec* dudt,
      const __global float* rho,
      const __global float* p,
      const __global float* m,
      const usize N,
      const float dt,
      const float dt_min,
      const float courant,
      const float h,
      const __global float* div_u,
      const __global vec* grad_p,
      __constant float* gamma,
      const __global float* lambda,
      const __global float* cp,
      const __global float* D_H2,
      const __global float* D_O2,
      const __global float* D_N2,
      const __global float* D_H2O)
{
	const usize i = get_global_id(0);
	if (i >= N)
		return;
	if (imove[i] <= 0) {
		dt_var[i] = dt;
		return;
	}

	float dxx = H;

	float s_i = sound_speed_perfect_gas(gamma[iset[i]], p[i], rho[i]);

	const float dt_u1 = courant * 0.4f * dxx /
	                    sqrt((4.0f * dxx * div_u[i] / rho[i]) *
	                             (4.0f * dxx * div_u[i] / rho[i]) +
	                         s_i * s_i);
	const float dt_u2 = courant * sqrt(dxx / (length(grad_p[i]) + 1.0e-12f));
	const float dt_u3 =
	    courant * 0.4f * dxx / sqrt(length(u[i]) * length(u[i]) + s_i * s_i);

	const float dt_u4 = courant * dxx / (length(u[i]) + 1.0e-12f);
	// const float dt_u5 = 0.1f / zeta_dot[i];
	// float xi_loc = lambda[i] / (rho[i] * cp[i]);
	const float dt_u5 = courant * dxx * dxx / (lambda[i] / (rho[i] * cp[i]));
	const float dt_u6 =
	    courant * dxx / (s_i + dxx * sqrt(div_u[i] * div_u[i]) / rho[i]);

	/*const float dt_u7 = courant * dxx * dxx / (D_H2[i] / rho[i]);
	const float dt_u8 = courant * dxx * dxx / (D_O2[i] / rho[i]);
	const float dt_u9 = courant * dxx * dxx / (D_N2[i] / rho[i]);
	const float dt_u10 = courant * dxx * dxx / (D_H2O[i] / rho[i]);*/

	const vec16 dt_many = courant * dxx * dxx / (Ds[i] / rho[i]);
	const float dt_u7 = givemax(dt_many);

	/*	const float dt_u = min(
	        min(min(min(min(min(min(min(min(dt_u1, dt_u2), dt_u3), dt_u4),
	   dt_u5), dt_u6), dt_u7), dt_u8), dt_u9), dt_u10);
	        */

	const float dt_u =
	    min(min(min(min(min(min(dt_u1, dt_u2), dt_u3), dt_u4), dt_u5), dt_u6),
	        dt_u7);

	dt_var[i] = max(min(dt, dt_u), dt_min);
}
