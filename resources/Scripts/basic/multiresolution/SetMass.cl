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
 *  @brief Set the mass of the particles according to their gamma value
 */

#ifndef HAVE_3D
    #include "../../types/2D.h"
#else
    #include "../../types/3D.h"
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
 * @param gamma_m Mass multiplier \f$ \gamma_m \f$.
 * @param m Mass \f$ m \f$.
 * @param N Number of particles.
 * @see P.N. Sun, A. Colagrossi, S. Marrone, A.M. Zhang. Multi-resolution
 * delta-SPH model. 2016
 */
__kernel void set_mass(__global const int* imove,
                       __global const float* gamma_m,
                       __global float* m,
                       unsigned int N)
{
    unsigned int i = get_global_id(0);
    if(i >= N)
        return;

    // Neglect boundary elements/particles
    if(imove[i] <= 0)
        return;

    m[i] *= gamma_m[i];
}

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
 * @param ilevel0 Level of refinement of the particle, by construction.
 * @param ilevel Current refinement level of the particle.
 * @param gamma_m Mass multiplier \f$ \gamma_m \f$.
 * @param m Mass \f$ m \f$.
 * @param N Number of particles.
 * @see P.N. Sun, A. Colagrossi, S. Marrone, A.M. Zhang. Multi-resolution
 * delta-SPH model. 2016
 */
__kernel void init_gamma(__global const int* imove,
                         __global const unsigned int* ilevel0,
                         __global const unsigned int* ilevel,
                         __global float* gamma_m,
                         unsigned int N)
{
    unsigned int i = get_global_id(0);
    if(i >= N)
        return;

    // Neglect boundary elements/particles
    if(imove[i] <= 0)
        return;

    if((ilevel0[i] <= ilevel[i]) && (ilevel0[i] != 0)){
        // non-mother particles
        gamma_m[i] = 0.f;
        return;
    }
    gamma_m[i] = 1.f;
}

/*
 * @}
 */