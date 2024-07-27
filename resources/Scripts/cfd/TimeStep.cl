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
 * than \f$ Ma_{\Delta t} h \f$, where the Courant factor is not taken into
 * account yet.
 *
 * Along this line, the distance moved by a particle can be written as follows:
 *
 * \f$ \vert \mathbf{r}_{n+1} - \mathbf{r}_{n} \vert = 
 *     \vert \mathbf{u} \vert \Delta t +
 *     \mathcal{O}({\Delta t}^2) \f$
 *
 * Such that rearraging the equation:
 *
 * \f$ \Delta t = Ma \frac{h}{\vert \mathbf{u} \vert} \f$
 *
 * @param imove Moving flags.
 *   - imove > 0 for regular fluid particles.
 *   - imove = 0 for sensors.
 *   - imove < 0 for boundary elements/particles.
 * @param u Velocity \f$ \mathbf{u}_{n+1/2} \f$.
 * @param dt_var Variable time step \f$ \mathrm{min} \left(
 * C_f \frac{h}{c_s}, C_f \frac{h}{10 \vert \mathbf{u} \vert}\right)\f$.
 * @param dt_Ma The maximum Mach number, \f$ Ma_{\Delta t} \f$
 * @param dt Fixed time step \f$ \Delta t = C_f \frac{h}{c_s} \f$.
 * @param dt_min Minimum time step \f$ \Delta t_{\mathrm{min}} \f$.
 * @param courant Courant factor \f$ C_{\Delta t} \f$.
 * @param h Kernel characteristic length \f$ h \f$.
 * @param Ma Mach number \f$ Ma \f$.
 * @param N Number of particles.
 */
__kernel void entry(const __global int* imove,
                    const __global vec* u,
                    __global float* dt_var,
                    const usize N,
                    const float dt,
                    const float dt_min,
                    const float courant,
                    const float dt_Ma,
                    const float h)
{
    const usize i = get_global_id(0);
    if(i >= N)
        return;
    if(imove[i] <= 0) {
        dt_var[i] = dt;
        return;
    }

    const float dr_max = dt_Ma * h;
    const float dt_u = courant * dr_max / length(u[i]);
    dt_var[i] = max(min(dt, dt_u), dt_min);
}
