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
 *  @brief Splitting particles methods
 */

#ifndef HAVE_3D
    #include "../../types/2D.h"
#else
    #include "../../types/3D.h"
#endif

/** @brief Remove the daughter particles which should coalesce.
 *
 * The particles should coalesce just when the refinement target level is lower
 * than the level of the particle.
 * This tool is decreasing the ilevel value of the coalescing mother particles.
 *
 * @param imove Moving flags.
 *   - imove > 0 for regular fluid/solid particles.
 *   - imove = 0 for sensors.
 *   - imove < 0 for boundary elements/particles.
 * @param ilevel0 Level of refinement of the particle, by construction.
 * @param ilevel Current refinement level of the particle.
 * @param level Target refinement level of the particle.
 * @param r Position \f$ \mathbf{r} \f$.
 * @param u Velocity \f$ \mathbf{u} \f$.
 * @param dudt Velocity rate of change \f$ \frac{d \mathbf{u}}{d t} \f$.
 * @param m Mass \f$ m \f$.
 * @param N Number of particles.
 * @see P.N. Sun, A. Colagrossi, S. Marrone, A.M. Zhang. Multi-resolution
 * delta-SPH model. 2016
 */
__kernel void entry(__global int* imove,
                    __global const unsigned int* ilevel0,
                    __global unsigned int* ilevel,
                    __global const unsigned int* level,
                    __global vec* r,
                    __global vec* u,
                    __global vec* dudt,
                    __global float* m0,
                    unsigned int N,
                    vec domain_max)
{
    unsigned int i = get_global_id(0);
    if(i >= N)
        return;

    if((imove[i] <= 0) ||           // Neglect boundary elements/particles
       (ilevel[i] <= level[i])) {   // Not asked to coalesce
        return;
    }

    ilevel[i]--;
    if(ilevel[i] > ilevel0[i]) {
        // It is just a mother particle asked to coalesce
        return;
    }
    
    // A daughter particles which should be removed
    imove[i] = -255;
    m0[i] = 0.f;
    // Stop the particle
    u[i] = VEC_ZERO;
    dudt[i] = VEC_ZERO;
    // Move the particle to a more convenient position (in order to can use
    // it as buffer)
    r[i] = domain_max;
}

/*
 * @}
 */