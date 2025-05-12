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
/// @param y_xx_in ordered mass fraction 
/// @param dy_xxdt_in rate of change ordered mass fraction
 
__kernel void
feed(const __global usize* mirror_src,
     usize N,
     __global float* z_in,
     __global float* y_H2_in,
     __global float* dy_H2dt_in,
     __global float* y_O2_in,
     __global float* dy_O2dt_in,
     __global float* y_N2_in,
     __global float* dy_N2dt_in,
     __global float* y_H2O_in,
     __global float* dy_H2Odt_in)
{

	const usize ii = get_global_id(0);
	if (ii >= N)
		return;
	const usize i = mirror_src[ii];
	if (i >= N)
		return;

	z_in[ii] = z_in[i];

	y_H2_in[ii] = y_H2_in[i];
	dy_H2dt_in[ii] = dy_H2dt_in[i];

	y_O2_in[ii] = y_O2_in[i];
	dy_O2dt_in[ii] = dy_O2dt_in[i];

	y_N2_in[ii] = y_N2_in[i];
	dy_N2dt_in[ii] = dy_N2dt_in[i];

	y_H2O_in[ii] = y_H2O_in[i];
	dy_H2Odt_in[ii] = dy_H2Odt_in[i];
}