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
 * @brief OpenCL tool to compute the fluid global force and moment.
 */

#include "resources/Scripts/types/types.h"

/** @brief Tool to compute the pressure force and moment of each particle.
 *
 * In this approach, since the global force and motion is looked for, the force
 * is computed as:
 * \f$$ f_a = m_a \cdot \frac{\nabla}{\rho_a} \f$$
 * And the moment as:
 * \f$$ m_a = m_a \cdot \mathbf{r}_0 \times \frac{\nabla}{\rho_a} \f$$
 *
 * @param pressureForces_f Force of each particle to be computed [N].
 * @param pressureForces_m Moment of each particle to be computed
 * [N \f$ \cdot \f$ m].
 * @param imove Moving flags.
 *   - imove > 0 for regular fluid particles.
 *   - imove = 0 for sensors.
 *   - imove < 0 for boundary elements/particles.
 * @param r Position \f$ \mathbf{r} \f$.
 * @param grad_p Pressure gradient \f$ \frac{\nabla p}{rho} \f$.
 * @param m Mass \f$ m \f$.
 * @param N Number of particles.
 * @param pressureForces_r Point with respect the moments are computed
 * \f$ \mathbf{r}_0 \f$.
 */
__kernel void entry(__global vec* pressureForces_f,
                    __global vec4* pressureForces_m,
                    const __global int* imove,
                    const __global vec* r,
                    const __global vec* grad_p,
                    const __global float* m,
                    unsigned int N,
                    vec pressureForces_r)
{
    // find position in global arrays
    unsigned int i = get_global_id(0);
    if(i >= N)
        return;
    if(imove[i] != 1){
        pressureForces_f[i] = VEC_ZERO;
        pressureForces_m[i] = (vec4)(0.f, 0.f, 0.f, 0.f);
        return;
    }

    const vec arm = r[i] - pressureForces_r;
    const vec acc = grad_p[i];  // Already divided by rho_i
    const float mass = m[i];
    pressureForces_f[i] = mass * acc;
    pressureForces_m[i].z = mass * (arm.x * acc.y - arm.y * acc.x);
    pressureForces_m[i].w = 0.f;
    #ifdef HAVE_3D
        pressureForces_m[i].x = mass * (arm.y * acc.z - arm.z * acc.y);
        pressureForces_m[i].y = mass * (arm.z * acc.x - arm.x * acc.z);
    #else
        pressureForces_m[i].x = 0.f;
        pressureForces_m[i].y = 0.f;
    #endif
}
