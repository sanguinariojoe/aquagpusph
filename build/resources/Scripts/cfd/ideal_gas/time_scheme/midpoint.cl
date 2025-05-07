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
 *  @brief Semi-implicit Midpoint Euler time integration scheme for the
 * internal energy.
 *
 * This is an extension of resources/Scripts/basic/time_scheme/midpoint.cl
 */

#include "resources/Scripts/types/types.h"

/** @brief Semi-implicit Midpoint Euler time integration scheme predictor stage
 * @param imove Moving flags.
 *   - imove > 0 for regular fluid/solid particles.
 *   - imove = 0 for sensors.
 *   - imove < 0 for boundary elements/particles.
 * @param eint Internal energy \f$ \e_{n} \f$.
 * @param deintdt Internal energy rate of change
 * \f$ \left. \frac{d e}{d t} \right\vert_{n} \f$.
 * @param eint_in Internal energy \f$ e_{n} \f$.
 * @param deintdt_in Internal energy rate of change
 * \f$ \left. \frac{d e}{d t} \right\vert_{n+1/2} \f$.
 * @param N Number of particles.
 */
__kernel void predictor(const __global float* eint,
                        const __global float* deintdt,
                        __global float* eint_in,
                        __global float* deintdt_in,
                        usize N)
{
    const usize i = get_global_id(0);
    if(i >= N)
        return;

    deintdt_in[i] = deintdt[i];
    eint_in[i] = eint[i];
}

/** @brief Advance to the time step midpoint
 * @param imove Moving flags.
 *   - imove > 0 for regular fluid particles.
 *   - imove = 0 for sensors.
 *   - imove < 0 for boundary elements/particles.
 * @param eint_in Internal energy \f$ e_{n+1/2} \f$.
 * @param deintdt Internal energy rate of change
 * \f$ \left. \frac{d e}{d t} \right\vert_{n+1/2} \f$.
 * @param eint Internal energy \f$ \e_{n} \f$.
 * @param N Number of particles.
 * @param dt Time step \f$ \Delta t \f$.
 */
__kernel void midpoint(const __global int* imove,
                       const __global float* eint_in,
                       const __global float* deintdt,
                       __global float* eint,
                       usize N,
                       float dt)
{
    const usize i = get_global_id(0);
    if(i >= N)
        return;
    if(imove[i] <= 0)
        return;

    eint[i] = eint_in[i] + 0.5f * dt * deintdt[i];
}

/** @brief Relax the obtained output to avoid diverging inner iterator
 * @param imove Moving flags.
 *   - imove > 0 for regular fluid particles.
 *   - imove = 0 for sensors.
 *   - imove < 0 for boundary elements/particles.
 * @param deintdt_in Internal energy rate of change
 * \f$ \left. \frac{d e}{d t} \right\vert_{n+1/2} \f$.
 * @param deintdt Internal energy rate of change
 * \f$ \left. \frac{d e}{d t} \right\vert_{n+1/2} \f$.
 * @param N Number of particles.
 * @param dt Time step \f$ \Delta t \f$.
 */
__kernel void relax(const __global int* imove,
                    const __global float* deintdt_in,
                    __global float* deintdt,
                    usize N,
                    float relax_midpoint)
{
    const usize i = get_global_id(0);
    if(i >= N)
        return;
    if(imove[i] <= 0)
        return;

    deintdt[i] = relax_midpoint * deintdt_in[i] +
        (1.f - relax_midpoint) * deintdt[i];
}

/** @brief 1st order Semi-implicit Euler time integration scheme corrector
 * stage
 * @param imove Moving flags.
 *   - imove > 0 for regular fluid particles.
 *   - imove = 0 for sensors.
 *   - imove < 0 for boundary elements/particles.
 * @param eint_in Internal energy \f$ e_{n+1/2} \f$.
 * @param deintdt Internal energy rate of change
 * \f$ \left. \frac{d e}{d t} \right\vert_{n+1/2} \f$.
 * @param eint Internal energy \f$ \e_{n+1} \f$.
 * @param N Number of particles.
 * @param dt Time step \f$ \Delta t \f$.
 */
__kernel void corrector(const __global int* imove,
                        const __global float* eint_in,
                        const __global float* deintdt,
                        __global float* eint,
                        usize N,
                        float dt)
{
    const usize i = get_global_id(0);
    if(i >= N)
        return;
    if(imove[i] <= 0)
        return;

    eint[i] = eint_in[i] + dt * deintdt[i];
}

/*
 * @}
 */
 
