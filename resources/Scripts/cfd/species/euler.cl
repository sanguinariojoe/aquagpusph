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
 * @param eint Internal energy \f$ \e_{n+1/2} \f$.
 * @param deintdt Internal energy rate of change
 * \f$ \left. \frac{d e}{d t} \right\vert_{n+1/2} \f$.
 * @param eint_in Internal energy \f$ e_{n+1} \f$.
 * @param deintdt_in Internal energy rate of change
 * \f$ \left. \frac{d e}{d t} \right\vert_{n+1} \f$.
 * @param N Number of particles.
 */
__kernel void predictor(const __global float* z,
                        const __global float* y_H2,
                        const __global float* dy_H2dt,
                        const __global float* y_O2,
                        const __global float* dy_O2dt,
                        const __global float* y_N2,
                        const __global float* dy_N2dt,
                        const __global float* y_H2O,
                        const __global float* dy_H2Odt,
                        __global float* z_in,
                        __global float* y_H2_in,
                        __global float* dy_H2dt_in,
                        __global float* y_O2_in,
                        __global float* dy_O2dt_in,
                        __global float* y_N2_in,
                        __global float* dy_N2dt_in,
                        __global float* y_H2O_in,
                        __global float* dy_H2Odt_in,
                        const usize N)
{
    const usize i = get_global_id(0);
    if(i >= N)
        return;
    z_in[i] = z[i];

    y_H2_in[i] = y_H2[i];
    dy_H2dt_in[i] = dy_H2dt[i];

    y_O2_in[i] = y_O2[i];
    dy_O2dt_in[i] = dy_O2dt[i];
    
    y_N2_in[i] = y_N2[i];
    dy_N2dt_in[i] = dy_N2dt[i];
    
    y_H2O_in[i] = y_H2O[i];
    dy_H2Odt_in[i] = dy_H2Odt[i];

}

/** @brief 1st order Euler time integration scheme corrector stage
 * @param imove Moving flags.
 *   - imove > 0 for regular fluid particles.
 *   - imove = 0 for sensors.
 *   - imove < 0 for boundary elements/particles.
 * @param eint Internal energy \f$ \e_{n+1/2} \f$.
 * @param deintdt Internal energy rate of change
 * \f$ \left. \frac{d e}{d t} \right\vert_{n+1/2} \f$.
 * @param N Number of particles.
 * @param dt Time step \f$ \Delta t \f$.
 */
__kernel void corrector(const __global int* imove,
                        __global float* y_H2,
                        const __global float* dy_H2dt,
                        __global float* y_O2,
                        const __global float* dy_O2dt,
                        __global float* y_N2,
                        const __global float* dy_N2dt,
                        __global float* y_H2O,
                        const __global float* dy_H2Odt,
                        __global float* z,
                        const __global float* dz_dt,                     
                        const unsigned int N,
                        const float dt)
{
    usize i = get_global_id(0);
    if(i >= N)
        return;

    if(imove[i] > 0) {
        z[i] += dt * dz_dt[i];
        y_H2[i] += dt * dy_H2dt[i];
        y_O2[i] += dt * dy_O2dt[i];
        y_N2[i] += dt * dy_N2dt[i];
        y_H2O[i] += dt * dy_H2Odt[i];
    }
}

/*
 * @}
 */
