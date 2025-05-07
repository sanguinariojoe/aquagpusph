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
 *  @brief Improved Euler integration scheme for the internal energy.
 *
 * This is an extension of resources/Scripts/basic/time_scheme/improved_euler.cl
 */

#include "resources/Scripts/types/types.h"

/** @brief Improved Euler time integration scheme predictor stage
 * @param imove Moving flags.
 *   - imove > 0 for regular fluid/solid particles.
 *   - imove = 0 for sensors.
 *   - imove < 0 for boundary elements/particles.
 * @param eint Internal energy \f$ \e_{n+1/2} \f$.
 * @param deintdt Internal energy rate of change
 * \f$ \left. \frac{d e}{d t} \right\vert_{n+1/2} \f$.
 * @param eint_in Internal energy \f$ e_{n+1} \f$.
 * @param deintdt_in Internal energy rate of change
 * \f$ \left. \frac{d e}{d t} \right\vert_{n+1} \f$.
 * @param N Number of particles.
 * @param dt Time step \f$ \Delta t \f$.
 */
__kernel void predictor(__global int* imove,
                        const __global float* eint,
                        const __global float* deintdt,
                        __global float* eint_in,
                        __global float* deintdt_in,
                        usize N,
                        float dt)
{
    usize i = get_global_id(0);
    if(i >= N)
        return;

    float DT = dt;
    if(imove[i] <= 0)
        DT = 0.f;

    deintdt_in[i] = deintdt[i];
    eint_in[i] = eint[i] + DT * deintdt[i];
}

/** @brief Improved Euler time integration scheme corrector stage
 * @param imove Moving flags.
 *   - imove > 0 for regular fluid particles.
 *   - imove = 0 for sensors.
 *   - imove < 0 for boundary elements/particles.
 * @param deintdt Internal energy rate of change
 * \f$ \left. \frac{d e}{d t} \right\vert_{n+1/2} \f$.
 * @param deintdt_in Internal energy rate of change
 * \f$ \left. \frac{d e}{d t} \right\vert_{n} \f$.
 * @param eint Internal energy \f$ \e_{n+1/2} \f$.
 * @param N Number of particles.
 * @param dt Time step \f$ \Delta t \f$.
 */
__kernel void corrector(const __global int* imove,
                        const __global float* deintdt,
                        const __global float* deintdt_in,
                        __global float* eint,
                        usize N,
                        float dt)
{
    usize i = get_global_id(0);
    if(i >= N)
        return;

    if(imove[i] > 0) {
        const float DT = 0.5f * dt;
        eint[i] += DT * (deintdt[i] - deintdt_in[i]);
    }
}

/*
 * @}
 */
 
