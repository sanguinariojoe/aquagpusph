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
 * @brief Mirroring process for the symmetry boundary condition.
 */

#include "resources/Scripts/types/types.h"


/// @brief 
/// @param z_in first component... Tracer unreacting
/// @param ys_in ordered mass fraction 
/// @param dysdt_in rate of change ordered mass fraction
 
__kernel void
feed(const __global usize* mirror_src,
     usize N,
     __global float* z_in,
	__global vec16* ys_in,
	__global vec16* dysdt_in
	 )
{

	const usize ii = get_global_id(0);
	if (ii >= N)
		return;
	const usize i = mirror_src[ii];
	if (i >= N)
		return;

	z_in[ii] = z_in[i];
	ys_in[ii] = ys_in[i];
	dysdt_in[ii] = dysdt_in[i];

}