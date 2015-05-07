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

/** @file
 * @brief Particles out of domain filter.
 */

#ifndef HAVE_3D
    #include "types/2D.h"
#else
    #include "types/3D.h"
#endif

/** @brief Check and destroy the particles out of the domain.
 * @param imove Moving flags.
 *   - imove > 0 for regular fluid particles.
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
__kernel void main(__global int* imove,
                   __global vec* r,
                   __global vec* u,
                   __global vec* dudt,
                   __global float* m,
                   uint N,
                   vec domain_min,
                   vec domain_max)
{
    // find position in global arrays
    unsigned int i = get_global_id(0);
    if(i >= N)
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
        // Set as fixed zero mass particle (sensor equivalent)
        imove[i] = 0;
        m[i] = 0.f;
        // Stop the particle
        u[i] = VEC_ZERO;
        dudt[i] = VEC_ZERO;
        // Move the particle to a more convenient position (in order to can use
        // it as buffer)
        r[i] = domain_max;
    }
}
