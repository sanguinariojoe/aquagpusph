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

/**
 * \addtogroup ideal_gas
 * @{
 */

/** @file
 * @brief Sort the internal energy by the cell indexes
 *
 * This is an extension of resources/Scripts/basic/Sort.cl
 */

#include "resources/Scripts/types/types.h"

/** @brief Sort the internal energy.
 *
 * @param D_xx unsorted diffusion coeficient
 * @param D_xx_in sorted diffusion coeficient
 * @param id_sorted Permutations list from the unsorted space to the sorted
 * one.
 * @param N Number of particles.
 */
__kernel void
entry(const __global float* D_H2_in,
      const __global float* D_O2_in,
      const __global float* D_N2_in,
      const __global float* D_H2O_in,
      __global float* D_H2,
      __global float* D_O2,
      __global float* D_N2,
      __global float* D_H2O,
      const __global usize* id_sorted,
      usize N)
{
	usize i = get_global_id(0);
	if (i >= N)
		return;

	const usize i_out = id_sorted[i];

	D_H2[i_out] = D_H2_in[i];
	D_O2[i_out] = D_O2_in[i];
	D_N2[i_out] = D_N2_in[i];
	D_H2O[i_out] = D_H2O_in[i];

	/*    rhs_yh2[i_out] = rhs_yh2_in[i];
	    rhs_yo2[i_out] = rhs_yo2_in[i];
	    rhs_yn2[i_out] = rhs_yn2_in[i];
	    rhs_yh2o[i_out] = rhs_yh2o_in[i];*/
}

/*
 * @}
 */
