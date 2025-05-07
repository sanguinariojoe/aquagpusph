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


/** @brief 1st order Euler time integration scheme corrector stage
 * @param imove Moving flags.
 *   - imove > 0 for regular fluid particles.
 *   - imove = 0 for sensors.
 *   - imove < 0 for boundary elements/particles.
 * @param rhs_yxx increments in the gases calculated in previous routine
 * (diffusion)
 * @param dy_xxdt add to the increments of the composition variation
 * \f$ \left. \frac{d e}{d t} \right\vert_{n+1/2} \f$.
 * @param N Number of particles.
 * @param dt Time step \f$ \Delta t \f$.
 */
__kernel void
add(const __global int* imove,
    __global float* dz_dt,
    __global float* dy_H2dt,
    __global float* dy_O2dt,
    __global float* dy_N2dt,
    __global float* dy_H2Odt,
    const __global float* rhs_yh2,
    const __global float* rhs_yo2,
    const __global float* rhs_yn2,
    const __global float* rhs_yh2o,
    const unsigned int N,
    const float dt)
{
	usize i = get_global_id(0);
	if (i >= N)
		return;

	if (imove[i] > 0) {

		dz_dt[i] += rhs_yh2[i];
		dy_H2dt[i] += rhs_yh2[i];
		dy_O2dt[i] += rhs_yo2[i];
		dy_N2dt[i] += rhs_yn2[i];
		dy_H2Odt[i] += rhs_yh2o[i];
	}
}

/*
 * @}
 */
