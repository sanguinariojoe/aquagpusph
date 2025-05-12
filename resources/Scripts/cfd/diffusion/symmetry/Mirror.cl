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


/**
 * @brief reflection routine for the species and diffusion
 * mirroring
 *
 */
/// @param D_H2_in Diffusion coefficient H2
/// @param D_O2_in Diffusion coefficient O2
/// @param D_N2_in Diffusion coefficient N2
/// @param D_H2O_in Diffusion coefficient H2O
__kernel void
feed(const __global usize* mirror_src,
     usize N,
     __global float* D_H2_in,
     __global float* D_O2_in,
     __global float* D_N2_in,
     __global float* D_H2O_in)
{
	const usize ii = get_global_id(0);
	if (ii >= N)
		return;
	const usize i = mirror_src[ii];
	if (i >= N)
		return;
		
	D_H2_in[ii] = D_H2_in[i];
	D_O2_in[ii] = D_O2_in[i];
	D_N2_in[ii] = D_N2_in[i];
	D_H2O_in[ii] = D_H2O_in[i];
}
