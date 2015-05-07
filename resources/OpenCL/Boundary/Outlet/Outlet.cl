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
 * @brief Vanish the density varaition rate of the particles crossing the outlet
 * infinite plane, destroying (moving out the domain) the particles far away.
 */

#ifndef HAVE_3D
    #include "../../types/2D.h"
#else
    #include "../../types/3D.h"
#endif

/** @brief Outlet boundary condition, which is consisting into vanishing the
 * density variation rate of the particles tresspassing the outlet plane, as
 * well as moving out of the computational domain the particles far away, such
 * that they can be used as a buffer.
 *
 * @param imove Moving flags.
 *   - imove > 0 for regular fluid particles.
 *   - imove = 0 for sensors.
 *   - imove < 0 for boundary elements/particles.
 * @param r Position \f$ \mathbf{r} \f$.
 * @param dudt Velocity rate of change \f$ \frac{d \mathbf{u}}{d t} \f$.
 * @param drhodt Density rate of change \f$ \frac{d \rho}{d t} \f$.
 * @param N Number of particles.
 * @param domain_max Maximum point of the computational domain.
 * @param inlet_r Lower corner of the inlet square.
 * @param inlet_n = Velocity direction of the generated particles.
 */
__kernel void main(__global int* imove,
                   __global vec* r,
                   __global vec* dudt,
                   __global float* drhodt,
                   unsigned int N,
                   vec domain_max,
                   vec outlet_r,
                   vec outlet_n)
{
    // find position in global arrays
    unsigned int i = get_global_id(0);
    if(i >= N)
        return;
    if(imove[i] <= 0)
        return;

    // Compute the distance to the outlet plane
    const float dist = dot(r[i] - outlet_r, outlet_n);
    if(dist < 0.f)
        return;

    drhodt[i] = 0.f;
    dudt[i] = VEC_ZERO;
    // Destroy the particles far away from the outlet plane
    if(dist > SUPPORT * H)
        r[i] = domain_max + VEC_ONE;
}
