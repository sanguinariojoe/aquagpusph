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
/*__kernel void predictor(const __global float* xi,
                        __global float* xi_in,
                        const usize N)
{
    const usize i = get_global_id(0);
    if(i >= N)
        return;

    //rhs_qdot_in[i] = rhs_qdot[i];
    xi_in[i] = xi[i];

}
*/

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
#include "resources/Scripts/types/types.h"

__kernel void add(const __global int* imove,
                        //__global float* eint,
                        __global vec* dudt,
                        __global float* deintdt,
                        const __global vec* rhs_dissipation_impulse,
                        const __global float* rhs_dissipation_energy,
                        const unsigned int N)
{
    usize i = get_global_id(0);
    if(i >= N)
        return;

    if(imove[i] > 0) {
        //eint[i] += dt * rhs_qdot[i];
        dudt[i] += rhs_dissipation_impulse[i];
        deintdt[i] += rhs_dissipation_energy[i];
    }
}

/*
 * @}
 */
