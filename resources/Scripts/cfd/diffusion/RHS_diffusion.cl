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
 * @param ys mass fraction
 * @param Ds diffusion coeficient
 * @param rhsy right hand sides of mass fractions equations
 */
__kernel void
entry(const __global int* imove,
      const __global vec* r,
      const __global float* rho,
      const __global float* m,
      __global vec16* rhsy,
      const __global vec16* ys,
      const __global vec16* Ds,
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

	const vec16 ys_i = ys[i];
	const vec16 Ds_i = Ds[i];

	const float rho_i = rho[i];

#ifndef LOCAL_MEM_SIZE

#define _RHS_YS_ rhsy[i]

#else
#define _RHS_YS_ rhsyl[it]

	__local vec16 rhsyl[LOCAL_MEM_SIZE];

	_RHS_YS_ = VEC16_ZERO;

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

			const vec16 ys_j = ys[j];

			const vec16 Ds_j = Ds[j];

			_RHS_YS_ += -4.0f * Ds_i * Ds_j /
			            ((rho_i * rho_j) * (Ds_i + Ds_j)) * (ys_i - ys_j) *
			            f_ij;
		}
	}
	END_NEIGHS()

#ifdef LOCAL_MEM_SIZE

	rhsy[i] = _RHS_YS_;

#endif
}
