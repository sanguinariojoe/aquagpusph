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
 * @brief The simplest boundary technique to assert the non-tresspasable
 * boundary condition.
 */

#include "resources/Scripts/types/types.h"
#include "resources/Scripts/KernelFunctions/Kernel.h"

#ifndef __DR_FACTOR__
    /** @def __DR_FACTOR__
     * @brief The boundary elements effect is restricted to a quadrangular area
     * of \f$ R \times R \f$, where \f$ R = DR_FACTOR \cdot \Delta r \f$.
     */
    #define __DR_FACTOR__ 0.5f
#endif

#ifndef __MIN_BOUND_DIST__
    /** @def __MIN_BOUND_DIST__
     * @brief The elastic bounce is not tolerating that a particle becomes
     * closer than this distance (multiplied by \f$ \Delta r \f$).
     */
    #define __MIN_BOUND_DIST__ 0.0f
#endif

/** @brief Performs the boundary effect on the fluid particles.
 * @param imove Moving flags.
 *   - imove > 0 for regular fluid particles.
 *   - imove = 0 for sensors.
 *   - imove < 0 for boundary elements/particles.
 * @param r Position \f$ \mathbf{r} \f$.
 * @param normal Normal \f$ \mathbf{n} \f$.
 * @param m Area \f$ s \f$.
 * @param u_in Velocity \f$ \mathbf{u} \f$.
 * @param jhoc Head and tail of chains for each cell.
 * @param dudt Velocity rate of change
 * \f$ \left. \frac{d \mathbf{u}}{d t} \right\vert_{n+1} \f$.
 * @param N Number of particles.
 * @param dt Time step \f$ \Delta t \f$.
 */
__kernel void entry(const __global int* restrict imove,
                    const __global vec* restrict r_in,
                    const __global vec* restrict normal,
                    const __global float* restrict m,
                    const __global vec* restrict u_in,
                    const __global svec2* restrict jhoc,
                    __global vec* restrict dudt,
                    usize N,
                    float dt)
{
    const usize i = get_global_id(0);
    const usize it = get_local_id(0);
    if(i >= N)
        return;
    if(imove[i] != 1)
        return;
    if(!dt)
        return;

    const vec_xyz r_i = r_in[i].XYZ;
    __private vec_xyz __dudt = dudt[i].XYZ;
    __private vec_xyz u_i = u_in[i].XYZ + 0.5f * dt * __dudt;

    FOR_NEIGHS(N, jhoc){
        if (imove[j] != -3)
            continue;
        const vec_xyz r_ij = r_in[j].XYZ - r_i;
        const vec_xyz n_j = normal[j].XYZ;
        const float rn = dot(r_ij, n_j);
        if(rn < 0.f){
            // The particle is on the "wrong" side of the wall.
            continue;
        }
#ifdef HAVE_3D
        const float dr = sqrt(m[j]);
#else
        const float dr = m[j];
#endif
        const float R = __DR_FACTOR__ * dr;
        const vec_xyz rt = r_ij - rn * n_j;
        if(dot(rt, rt) >= R * R){
            // The particle is passing too far from the boundary element
            continue;
        }

        {
            const float drn = dt * dot(u_i, n_j);
            if(drn < 0.f){
                // The particle is already running away from the boundary
                continue;
            }

            // The particle should be corrected if:
            //   - It is already placed in the effect zone.
            //   - It is entering inside the effect zone.
            if(rn - drn <= __MIN_BOUND_DIST__ * dr){
                // Reflect the particle velocity, so its module remains
                // constant
                const vec_xyz u = u_in[i].XYZ + dt * __dudt;
                const vec_xyz u_r = u - 2.f * dot(u, n_j) * n_j;
                // Modify the values for the next wall tests.
                __dudt = (u_r - u_in[i].XYZ) / dt;
                u_i = u_in[i].XYZ + 0.5f * dt * __dudt;
            }
        }
    }END_FOR_NEIGHS()

    dudt[i].XYZ = __dudt;
}

/** @brief Compute the force of each fluid particle on the boundary due to the
 * elastic bounce.
 *
 * @param imove Moving flags.
 *   - imove > 0 for regular fluid particles.
 *   - imove = 0 for sensors.
 *   - imove < 0 for boundary elements/particles.
 * @param m Mass \f$ m \f$.
 * @param dudt_preelastic Velocity rate of change before the elastic bounce
 * \f$ \left. \frac{d \mathbf{u}}{d t} \right\vert_{n+1} \f$.
 * @param dudt_elastic Velocity rate of change after the elatic bounce
 * \f$ \left. \frac{d \mathbf{u}}{d t} \right\vert_{n+1} \f$.
 * @param N Number of particles.
 */
__kernel void force_bound(const __global int* restrict imove,
                          const __global float* restrict m,
                          const __global vec* restrict dudt_preelastic,
                          const __global vec* restrict dudt_elastic,
                          __global vec* restrict force_elastic,
                          usize N)
{
    const usize i = get_global_id(0);
    if(i >= N)
        return;
    if(imove[i] != 1) {
        force_elastic[i] = VEC_ZERO;
        return;
    }

    force_elastic[i] = -m[i] * (dudt_elastic[i] - dudt_preelastic[i]);
}
