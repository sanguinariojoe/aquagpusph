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

#include "resources/Scripts/types/types.h"

/** @brief Sort the internal energy.
 *

 * @param N Number of particles.
 */

/// @param xi unordered thermal difusivity
/// @param xi_in ordered thermal difusivity
/// @param nu unordered kinematic viscosity
/// @param nu_in ordered kinematic viscosity

__kernel void
entry(const __global usize* id_sorted,
      __global float* xi,
      const __global float* xi_in,
      usize N)
{
	usize i = get_global_id(0);
	if (i >= N)
		return;

	const usize i_out = id_sorted[i];

	xi[i_out] = xi_in[i];

}

/*
 * @}
 */
