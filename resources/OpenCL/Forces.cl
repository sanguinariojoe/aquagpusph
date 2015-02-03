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

#ifndef HAVE_3D
    #include "types/2D.h"
#else
    #include "types/3D.h"
#endif

/** @brief Tool to compute the force and moment of each particle.
 *
 * In this approach, since the global force and motion is looked for, the force
 * is computed as:
 * \f$$ f_a = - m_a \cdot \left(
 *   \frac{d \mathbf{u}}{d t} - \mathbf{g} \right) \f$$
 * And the moment as:
 * \f$$ f_a = - m_a \cdot \mathbf{r}_0 \times \left(
 *   \frac{d \mathbf{u}}{d t} - \mathbf{g} \right) \f$$
 *
 * @param forces_f Force of each particle to be computed [N].
 * @param forces_m Moment of each particle to be computed
 * [N \f$ \cdot \f$ m].
 * @param imove Moving flags.
 *   - imove > 0 for regular fluid particles.
 *   - imove = 0 for sensors.
 *   - imove < 0 for boundary elements/particles.
 * @param r Position \f$ \mathbf{r} \f$.
 * @param dudt Velocity rate of change
 * \f$ \frac{d \mathbf{u}}{d t} \f$.
 * @param m Mass \f$ m \f$.
 * @param g Gravity acceleration \f$ \mathbf{g} \f$.
 * @param forces_r Point with respect the moments are computed
 * \f$ \mathbf{r}_0 \f$.
 */
__kernel void main(__global vec* forces_f,
                   __global vec4* forces_m,
                   __global int* imove,
                   __global vec* r,
                   __global vec* dudt,
                   __global float* m,
                   unsigned int N,
                   vec g,
                   vec forces_r)
{
    // find position in global arrays
    unsigned int i = get_global_id(0);
    if(i >= N)
        return;
    if(imove[i] <= 0){
        forces_f[i] = VEC_ZERO;
        forces_m[i] = (vec4)(0.f, 0.f, 0.f, 0.f);
        return;
    }

    const vec arm = r[i] - forces_r;
    const vec acc = g - dudt[i];
    const float mass = m[i];
    forces_f[i] = mass * acc;
    forces_m[i].z = mass * (arm.x * acc.y - arm.y * acc.x);
    forces_m[i].w = 0.f;
    #ifdef HAVE_3D
        forces_m[i].x = mass * (arm.y * acc.z - arm.z * acc.y);
        forces_m[i].y = mass * (arm.z * acc.x - arm.x * acc.z);
    #else
        forces_m[i].x = 0.f;
        forces_m[i].y = 0.f;
    #endif
}
