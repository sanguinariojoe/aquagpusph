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
 * @brief 1st order Euler integration scheme for the internal energy.
 *
 * This is an extension of resources/Scripts/basic/time_scheme/euler.cl
 */

/** @brief 1st order Euler time integration scheme predictor stage
 * @param z first gas component
 * @param ys_in ordered mass fraction
 * @param ys unordered mass fraction
 * @param dysdt_in ordered mass fraction rate of change
 * @param dysdt unordered mass fraction rate of change
 * @param N Number of particles.
 */

#include "resources/Scripts/types/types.h"

__kernel void
predictor(const __global float* z,
          const __global species_t* ys,
          const __global species_t* dysdt,
          __global float* z_in,
          __global species_t* ys_in,
          __global species_t* dysdt_in,
          const usize N)
{
	const usize i = get_global_id(0);
	if (i >= N)
		return;
	z_in[i] = z[i];
	ys_in[i] = ys[i];
	dysdt_in[i] = dysdt[i];
}

/** @brief 1st order Euler time integration scheme corrector stage
 * @param imove Moving flags.
 *   - imove > 0 for regular fluid particles.
 *   - imove = 0 for sensors.
 *   - imove < 0 for boundary elements/particles.

 * @param N Number of particles.
 * @param dt Time step \f$ \Delta t \f$.
 * @param y_xx unordered mass fraction
 * @param dy_xxdt unordered mass fraction rate of change
 * @param z first gas component
 * @param dz_dt first gas component rate of change
 */
__kernel void
corrector(const __global int* imove,
          __global species_t* ys,
          const __global species_t* dysdt,
          __global float* z,
          const __global float* dz_dt,
          const unsigned int N,
          const float dt)
{
	usize i = get_global_id(0);
	if (i >= N)
		return;

	if (imove[i] > 0) {
		z[i] += dt * dz_dt[i];
		ys[i] += dt * dysdt[i];
	}
}

/*
 * @}
 */
