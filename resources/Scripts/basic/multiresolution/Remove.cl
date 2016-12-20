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
 * @brief Set the mass of the particles according to the iterations since they
 * were split/coalesced
 */

#include "resources/Scripts/types/types.h"

#ifndef M_ITERS
    #define M_ITERS 10
#endif


/** @brief Set the mass of the particles according to their gamma value.
 *
 * When a particle is split in a set of daughter particles, in order to avoid
 * shocks, its effect is smoothly transfered to the daughters, using for that a
 * \f$ \gamma_m \f$ mass multiplier.
 *
 * @param imove Moving flags.
 *   - imove > 0 for regular fluid/solid particles.
 *   - imove = 0 for sensors.
 *   - imove < 0 for boundary elements/particles.
 * @param miter Mass transfer iteration (Positive for growing particles,
 * negative for shrinking particles).
 * @param r Position \f$ \mathbf{r} \f$.
 * @param domain_max Maximum point of the computational domain.
 * @param N Number of particles.
 */
__kernel void entry(__global int* imove,
                    __global const int* miter,
                    __global vec* r,
                    vec domain_max,
                    unsigned int N)
{
    unsigned int i = get_global_id(0);
    if(i >= N)
        return;

    // Neglect boundary elements/particles
    if(imove[i] <= 0){
        return;
    }

    if(miter[i] < -M_ITERS){
        // Outdated partner particle
        r[i] = domain_max + VEC_ONE;
        imove[i] = -255;
    }
}

/*
 * @}
 */