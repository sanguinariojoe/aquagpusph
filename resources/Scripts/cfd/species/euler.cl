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
 * @param y_xx_in ordered mass fraction
 * @param y_xx unordered mass fraction
 * @param dy_xxdt_in ordered mass fraction rate of change
 * @param dy_xxdt unordered mass fraction rate of change
 * @param N Number of particles.
 */
__kernel void
predictor(const __global float* z,
          const __global vec16* ys,
          const __global vec16* dysdt,
          /*const __global float* y_H2,
          const __global float* dy_H2dt,
          const __global float* y_O2,
          const __global float* dy_O2dt,
          const __global float* y_N2,
          const __global float* dy_N2dt,
          const __global float* y_H2O,
          const __global float* dy_H2Odt,*/
          __global float* z_in,
          __global vec16* ys_in,
          __global vec16* dysdt_in,
          /*__global float* y_H2_in,
          __global float* dy_H2dt_in,
          __global float* y_O2_in,
          __global float* dy_O2dt_in,
          __global float* y_N2_in,
          __global float* dy_N2dt_in,
          __global float* y_H2O_in,
          __global float* dy_H2Odt_in,*/
          const usize N)
{
	const usize i = get_global_id(0);
	if (i >= N)
		return;
	z_in[i] = z[i];
    ysin[i]= ys[i];
    dysdt_in[i] = dysdt[i];

	/*ys_in[i].H2 = y_H2[i];
	dysdt_in[i].H2 = dy_H2dt[i];

	ys_O2_in[i] = y_O2[i];
	dys_O2dt_in[i] = dy_O2dt[i];

	ys_N2_in[i] = y_N2[i];
	dys_N2dt_in[i] = dy_N2dt[i];

	ys_H2O_in[i] = y_H2O[i];
	dys_H2Odt_in[i] = dy_H2Odt[i];*/
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
          __global vec16* ys,
          const __global vec16* dysdt,/*
          __global float* y_H2,
          const __global float* dy_H2dt,
          __global float* y_O2,
          const __global float* dy_O2dt,
          __global float* y_N2,
          const __global float* dy_N2dt,
          __global float* y_H2O,
          const __global float* dy_H2Odt,*/
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

        /*
		y_H2[i] += dt * dy_H2dt[i];
		y_O2[i] += dt * dy_O2dt[i];
		y_N2[i] += dt * dy_N2dt[i];
		y_H2O[i] += dt * dy_H2Odt[i];
        */
	}
}

/*
 * @}
 */
