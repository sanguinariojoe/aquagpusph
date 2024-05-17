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
 * @brief Variable time step computation.
 */

#include "resources/Scripts/types/types.h"

/** @brief Compute the maximum time step for each particle.
 *
 * In SPH the time step is selected to enforce the particles may not move more
 * than \f$ 0.1 h \f$, where the Courant factor is not taken into account yet.
 *
 * Along this line, the distance moved by a particle can be written as follows:
 *
 * \f$ \vert \mathbf{r}_{n+1} - \mathbf{r}_{n} \vert = 
 *     \vert \mathbf{u} \vert \Delta t +
 *     \frac{1}{2} \left\vert
 *                     \frac{\mathrm{d} \mathbf{u}}{\mathrm{d} t}
 *                 \right\vert {\Delta t}^2 +
       \mathcal{O}({\Delta t}^3) \f$
 *
 * Such that, taking maximums, and rearraging the equation:
 *
 * \f$ \Delta t = \frac{1}{20} \min \left(
 *     \frac{h}{\vert \mathbf{u} \vert},
 *     \sqrt{\frac{2 h}{\left\vert
 *                          \frac{\mathrm{d} \mathbf{u}}{\mathrm{d} t}
 *                      \right\vert}}
 * \right) \f$
 *
 * @param dt_var Variable time step \f$ \mathrm{min} \left(
 * C_f \frac{h}{c_s}, C_f \frac{h}{10 \vert \mathbf{u} \vert}\right)\f$.
 * @param u Velocity \f$ \mathbf{u}_{n+1/2} \f$.
 * @param dudt Velocity rate of change \f$ \frac{d \mathbf{u}}{d t} \f$.
 * @param N Number of particles.
 * @param dt Fixed time step \f$ \Delta t = C_f \frac{h}{c_s} \f$.
 * @param dt_min Minimum time step \f$ \Delta t_{\mathrm{min}} \f$.
 * @param courant Courant factor \f$ C_f \f$.
 * @param h Kernel characteristic length \f$ h \f$.
 */
__kernel void entry(__global float* dt_var,
                    __global vec* u,
                    __global vec* dudt,
                    usize N,
                    float dt,
                    float dt_min,
                    float courant,
                    float h)
{
    const usize i = get_global_id(0);
    if(i >= N)
        return;

    const float dr_max = 0.1f * h;
    const float dt_u = courant * min(dr_max / length(u[i]),
                                     sqrt(dr_max / (0.5f * length(dudt[i]))));
    dt_var[i] = max(min(dt, dt_u), dt_min);
}
