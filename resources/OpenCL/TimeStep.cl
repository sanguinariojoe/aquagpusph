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

#ifndef HAVE_3D
    #include "types/2D.h"
#else
    #include "types/3D.h"
#endif

/** @brief Compute the maximum time step for each particle.
 * @param dt_var Variable time step \f$ \mathrm{min} \left(
 * C_f \frac{h}{c_s}, C_f \frac{h}{10 \vert \mathbf{u} \vert}\right)\f$.
 * @param v Velocity \f$ \mathbf{u}_{n+1/2} \f$.
 * @param dvdt Velocity rate of change \f$ \frac{d \mathbf{u}}{d t} \f$.
 * @param N Number of particles.
 * @param dt Fixed time step \f$ \Delta t = C_f \frac{h}{c_s} \f$.
 * @param dt_min Minimum time step \f$ \Delta t_{\mathrm{min}} \f$.
 * @param courant Courant factor \f$ C_f \f$.
 * @param h Kernel characteristic length \f$ h \f$.
 */
__kernel void main(__global float* dt_var,
                   __global vec* v,
                   __global vec* dvdt,
                   unsigned int N,
                   float dt,
                   float dt_min,
                   float courant,
                   float h)
{
    unsigned int i = get_global_id(0);
    if(i >= N)
        return;

    const float vv = 10.f * fast_length(v[i] + dvdt[i] * dt);
    dt_var[i] = max(min(dt, h / vv), dt_min);
}
