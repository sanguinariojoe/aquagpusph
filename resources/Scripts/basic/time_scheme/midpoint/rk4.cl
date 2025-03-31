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

/** @defgroup basic Basic preset
 *
 * @brief Basic preset of tools to build more complex sets of tools later
 * 
 * @{
 */

/** @file
 *  @brief RK4 implicit midpoint triggerer
 */

#include "resources/Scripts/types/types.h"

/** @brief Compute the final derivative value for the RK4 triggering scheme
 * @param dudt0 Velocity rate of change
 * \f$ \left. \frac{d \mathbf{u}}{d t} \right\vert_{0} \f$.
 * @param drhodt0 Density rate of change
 * \f$ \left. \frac{d \rho}{d t} \right\vert_{0} \f$.
 * @param dudt1 Velocity rate of change
 * \f$ \left. \frac{d \mathbf{u}}{d t} \right\vert_{1} \f$.
 * @param drhodt1 Density rate of change
 * \f$ \left. \frac{d \rho}{d t} \right\vert_{1} \f$.
 * @param dudt2 Velocity rate of change
 * \f$ \left. \frac{d \mathbf{u}}{d t} \right\vert_{2} \f$.
 * @param drhodt2 Density rate of change
 * \f$ \left. \frac{d \rho}{d t} \right\vert_{2} \f$.
 * @param dudt Velocity rate of change
 * \f$ \left. \frac{d \mathbf{u}}{d t} \right\vert_{3} \f$.
 * @param drhodt Density rate of change
 * \f$ \left. \frac{d \rho}{d t} \right\vert_{3} \f$.
 * @param N Number of particles.
 */
__kernel void entry(const __global vec* dudt0,
                    const __global float* drhodt0,
                    const __global vec* dudt1,
                    const __global float* drhodt1,
                    const __global vec* dudt2,
                    const __global float* drhodt2,
                    __global vec* dudt,
                    __global float* drhodt,
                    usize N)
{
    const usize i = get_global_id(0);
    if(i >= N)
        return;

    dudt[i] = (dudt0[i] + dudt[i]) / 6.f + (dudt1[i] + dudt2[i]) / 3.f;
    drhodt[i] = (drhodt0[i] + drhodt[i]) / 6.f + (drhodt1[i] + drhodt2[i]) / 3.f;
}

/*
 * @}
 */
 
