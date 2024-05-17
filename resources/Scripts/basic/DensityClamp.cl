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

/** @addtogroup basic
 * @{
 */

/** @file
 * @brief Particles out of domain filter.
 */

#include "resources/Scripts/types/types.h"

/** @brief Clamp the density
 *
 * In some situations it can be convenient to define a minimum or maximum
 * density values in order to avoid that a crazy particle may cause the
 * simulation blow up.
 * Use this tool only if you know exactly what you are doing
 *
 * @param rho_in Density \f$ \rho_{n+1/2} \f$.
 * @param N Number of particles.
 * @param rho_min Minimum tolerated density value \f$ \rho_{min} \f$.
 * @param rho_max Maximum tolerated density value \f$ \rho_{max} \f$.
 */
__kernel void entry(__global float* rho_in,
                    usize N,
                    float rho_min,
                    float rho_max)
{
    usize i = get_global_id(0);
    if(i >= N)
        return;

    if(rho_in[i] < rho_min) rho_in[i] = rho_min;
    if(rho_in[i] > rho_max) rho_in[i] = rho_max;
}

/*
 * @}
 */
