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
 * @brief Vanish the velocity and desnity rates of variation of the velocity
 * and density for the dummy particles of the inlet.
 */

#ifndef HAVE_3D
    #include "../../types/2D.h"
#else
    #include "../../types/3D.h"
#endif

/** @brief Enforce the paticles
 *
 * @param imove Moving flags.
 *   - imove > 0 for regular fluid particles.
 *   - imove = 0 for sensors.
 *   - imove < 0 for boundary elements/particles.
 * @param r Position \f$ \mathbf{r} \f$.
 * @param dudt Velocity rate of change \f$ \frac{d \mathbf{u}}{d t} \f$.
 * @param drhodt Density rate of change \f$ \frac{d \rho}{d t} \f$.
 * @param N Number of particles.
 * @param inlet_r Lower corner of the inlet square.
 * @param inlet_n = Velocity direction of the generated particles.
 */
__kernel void main(__global int* imove,
                   __global vec* r,
                   __global vec* u,
                   __global vec* dudt,
                   __global float* drhodt,
                   unsigned int N,
                   vec inlet_r,
                   float inlet_U,
                   vec inlet_n)
{
    // find position in global arrays
    unsigned int i = get_global_id(0);
    if(i >= N)
        return;
    if(imove[i] <= 0)
        return;

    // Discard the particles already passed through the inlet
    if(dot(r[i] - inlet_r, inlet_n) > 0.f)
        return;

    u[i] = inlet_U * inlet_n;
    dudt[i] = VEC_ZERO;
    drhodt[i] = 0.f;
}
