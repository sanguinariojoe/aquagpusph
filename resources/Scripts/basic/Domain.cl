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

/** @brief Check and destroy the particles out of the domain.
 *
 * Usually is a good methodology to impose a computational domain, such that a
 * single lost particle will not cause the simulation blow up.
 * Since the particles out of the domain will take 0 mass, you can control how
 * many particles are lost controlling the total mass 
 *
 * @param imove Moving flags.
 *   - imove > 0 for regular fluid/solid particles.
 *   - imove = 0 for sensors.
 *   - imove < 0 for boundary elements/particles.
 * @param r Position \f$ \mathbf{r} \f$.
 * @param u Velocity \f$ \mathbf{u} \f$.
 * @param dudt Velocity rate of change \f$ \frac{d \mathbf{u}}{d t} \f$.
 * @param m Mass \f$ m \f$.
 * @param N Number of particles.
 * @param domain_min Minimum point of the domain.
 * @param domain_max Maximum point of the domain.
 */
__kernel void entry(__global int* imove,
                    __global vec* r,
                    __global vec* u,
                    __global vec* dudt,
                    __global float* m,
                    uint N,
                    vec domain_min,
                    vec domain_max)
{
    unsigned int i = get_global_id(0);
    if(i >= N)
        return;
    if(imove[i] <= -255)
        return;

    const vec coords = r[i];
    if(    (coords.x < domain_min.x)
        || (coords.y < domain_min.y)
        || (coords.x > domain_max.x)
        || (coords.y > domain_max.y)
        #ifdef HAVE_3D
        || (coords.z < domain_min.z)
        || (coords.z > domain_max.z)
        #endif
      )
    {
        // Set as buffer particle (non-interacting one)
        imove[i] = -256;
        m[i] = 0.f;
        // Stop the particle
        u[i] = VEC_ZERO;
        dudt[i] = VEC_ZERO;
        // Move the particle to a more convenient position (in order to can use
        // it as buffer)
        r[i] = domain_max;
    }
}

/*
 * @}
 */
