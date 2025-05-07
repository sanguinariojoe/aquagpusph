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
 * @brief Fluid particles interactions computation.
 */

#if defined(LOCAL_MEM_SIZE) && defined(NO_LOCAL_MEM)
#error NO_LOCAL_MEM has been set.
#endif

#include "resources/Scripts/types/types.h"
#include "resources/Scripts/KernelFunctions/Kernel.h"

/** @brief Fluid particles interactions computation.
 *
 * computation of diffusion of the different speceis
 *
 * @param imove Moving flags.
 *   - imove > 0 for regular fluid particles.
 *   - imove = 0 for sensors.
 *   - imove < 0 for boundary elements/particles.
 * @param r Position \f$ \mathbf{r} \f$.
 * @param u Velocity \f$ \mathbf{u} \f$.
 * @param rho Density \f$ \rho \f$.
 * @param m Mass \f$ m \f$.
 * @param N Number of particles.
 * @param icell Cell where each particle is located.
 * @param ihoc Head of chain for each cell (first particle found).
 * @param n_cells Number of cells in each direction
 * @param y_xx mass fraction
 * @param D_xx diffusion coeficient
 *
 */
__kernel void
entry(const __global int* imove,
      const __global vec* r,
      // const __global vec* u,
      const __global float* rho,
      const __global float* m,
      // const __global float* p,
      __global float* rhs_yh2,
      __global float* rhs_yo2,
      __global float* rhs_yn2,
      __global float* rhs_yh2o,
      const __global float* y_H2,
      const __global float* y_O2,
      const __global float* y_N2,
      const __global float* y_H2O,
      const __global float* D_H2,
      const __global float* D_O2,
      const __global float* D_N2,
      const __global float* D_H2O,
      usize N,
      LINKLIST_LOCAL_PARAMS)
{
	const usize i = get_global_id(0);
	const usize it = get_local_id(0);
	if (i >= N)
		return;
	if (imove[i] != 1) {
		return;
	}

	const vec_xyz r_i = r[i].XYZ;

	const float yh2_i = y_H2[i];
	const float yh2o_i = y_H2O[i];
	const float yn2_i = y_N2[i];
	const float yo2_i = y_O2[i];

	const float Dh2_i = D_H2[i];
	const float Dh2o_i = D_H2O[i];
	const float Dn2_i = D_N2[i];
	const float Do2_i = D_O2[i];

	const float rho_i = rho[i];

// printf("D in i is %f\n", D_H2[i]);

// Initialize the output
#ifndef LOCAL_MEM_SIZE

#define _RHS_YH2_ rhs_yh2[i]
#define _RHS_YO2_ rhs_yo2[i]
#define _RHS_YN2_ rhs_yn2[i]
#define _RHS_YH2O_ rhs_yh2o[i]

#else

#define _RHS_YH2_ rhs_yh2_l[it]
#define _RHS_YO2_ rhs_yo2_l[it]
#define _RHS_YN2_ rhs_yn2_l[it]
#define _RHS_YH2O_ rhs_yh2o_l[it]

	__local float rhs_yh2_l[LOCAL_MEM_SIZE];
	__local float rhs_yo2_l[LOCAL_MEM_SIZE];
	__local float rhs_yn2_l[LOCAL_MEM_SIZE];
	__local float rhs_yh2o_l[LOCAL_MEM_SIZE];

	_RHS_YH2_ = 0.f;
	_RHS_YO2_ = 0.f;
	_RHS_YN2_ = 0.f;
	_RHS_YH2O_ = 0.f;

#endif

	const usize c_i = icell[i];
	BEGIN_NEIGHS(c_i, N, n_cells, icell, ihoc)
	{
		if (i == j) {
			j++;
			continue;
		}
		if (imove[j] != 1) {
			j++;
			continue;
		}
		const vec_xyz r_ij = r[j].XYZ - r_i;
		const float q = length(r_ij) / H;
		if (q >= SUPPORT) {
			j++;
			continue;
		}
		{
			const float rho_j = rho[j];
			const float f_ij = kernelF(q) * CONF * m[j];

			const float yh2_j = y_H2[j];
			const float yh2o_j = y_H2O[j];
			const float yn2_j = y_N2[j];
			const float yo2_j = y_O2[j];

			const float Dh2_j = D_H2[j];
			const float Dh2o_j = D_H2O[j];
			const float Dn2_j = D_N2[j];
			const float Do2_j = D_O2[j];

			//_RHS_YH2_ += 4.0f * Dh2_i * Dh2_j / (rho_i * Dh2_i + rho_j *
			//Dh2_j)*(yh2_i-yh2_j)*f_ij; _RHS_YO2_ += 4.0f * Do2_i * Do2_j /
			//(rho_i * Do2_i + rho_j * Do2_j)*(yo2_i-yo2_j)*f_ij; _RHS_YN2_
			//+= 4.0f * Dn2_i * Dn2_j / (rho_i * Dn2_i + rho_j *
			//Dn2_j)*(yn2_i-yn2_j)*f_ij; _RHS_YH2O_ += 4.0f * Dh2o_i * Dh2o_j /
			//(rho_i * Dh2o_i + rho_j * Dh2o_j) * (yh2o_i - yh2o_j)*f_ij;

			_RHS_YH2_ += -4.0f * Dh2_i * Dh2_j /
			             ((rho_i * rho_j) * (Dh2_i + Dh2_j)) * (yh2_i - yh2_j) *
			             f_ij;
			_RHS_YO2_ += -4.0f * Do2_i * Do2_j /
			             ((rho_i * rho_j) * (Do2_i + Do2_j)) * (yo2_i - yo2_j) *
			             f_ij;
			_RHS_YN2_ += -4.0f * Dn2_i * Dn2_j /
			             ((rho_i * rho_j) * (Dn2_i + Dn2_j)) * (yn2_i - yn2_j) *
			             f_ij;
			_RHS_YH2O_ += -4.0f * Dh2o_i * Dh2o_j /
			              ((rho_i * rho_j) * (Dh2o_i + Dh2o_j)) *
			              (yh2o_i - yh2o_j) * f_ij;
		}
	}
	END_NEIGHS()

#ifdef LOCAL_MEM_SIZE

	rhs_yh2[i] = _RHS_YH2_;
	rhs_yo2[i] = _RHS_YO2_;
	rhs_yn2[i] = _RHS_YN2_;
	rhs_yh2o[i] = _RHS_YH2O_;

#endif
}
