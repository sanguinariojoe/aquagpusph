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

#if defined(LOCAL_MEM_SIZE) && defined(NO_LOCAL_MEM)
    #error NO_LOCAL_MEM has been set.
#endif

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
    #define __MIN_BOUND_DIST__ 0.1f
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
 * @param dudt Velocity rate of change
 * \f$ \left. \frac{d \mathbf{u}}{d t} \right\vert_{n+1} \f$.
 * @param N Number of particles.
 * @param dt Time step \f$ \Delta t \f$.
 * @param icell Cell where each particle is located.
 * @param ihoc Head of chain for each cell (first particle found).
 * @param n_cells Number of cells in each direction
 */
__kernel void entry(const __global int* imove,
                    const __global vec* r_in,
                    const __global vec* normal,
                    const __global float* m,
                    const __global vec* u_in,
                    __global vec* dudt,
                    usize N,
                    float dt,
                    LINKLIST_LOCAL_PARAMS)
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

    // Initialize the output
    #ifndef LOCAL_MEM_SIZE
        #define _U_ u_i
        #define _DUDT_ dudt[i].XYZ
        vec_xyz u_i;
    #else
        #define _U_ u_i[it]
        #define _DUDT_ dudt_l[it]
        __local vec_xyz u_i[LOCAL_MEM_SIZE];
        __local vec_xyz dudt_l[LOCAL_MEM_SIZE];
        _DUDT_ = dudt[i].XYZ;
    #endif
    _U_ = u_in[i].XYZ + 0.5f * dt * dudt[i].XYZ;

    const usize c_i = icell[i];
    BEGIN_NEIGHS(c_i, N, n_cells, icell, ihoc){
        if (imove[j] != -3) {
            j++;
            continue;
        }
        const vec_xyz r_ij = r_in[j].XYZ - r_i;
        const vec_xyz n_j = normal[j].XYZ;
        const float rn = dot(r_ij, n_j);
        if(rn < 0.f){
            // The particle is on the "wrong" side of the wall.
            j++;
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
            j++;
            continue;
        }

        {
            const float drn = dt * dot(_U_, n_j);
            if(drn < 0.f){
                // The particle is already running away from the boundary
                j++;
                continue;
            }

            // The particle should be corrected if:
            //   - It is already placed in the effect zone.
            //   - It is entering inside the effect zone.
            if(rn - drn <= __MIN_BOUND_DIST__ * dr){
                // Reflect the particle velocity, so its module remains
                // constant
                const vec_xyz u_1 = u_in[i].XYZ + dt * dudt[i].XYZ; // full step velocity
                const vec_xyz u_r = u_in[i].XYZ - 2.f * dot(u_1, n_j) * n_j;  // reversed-normal initial velocity
                // Modify the values for the next wall tests.
                _DUDT_ = (u_r - u_1) / dt;

//                _U_ = u_r;  // I think this would apply the elastic force twice!

            }
        }
    }END_NEIGHS()

    #ifdef LOCAL_MEM_SIZE
        dudt[i].XYZ = _DUDT_;
    #endif
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
__kernel void force_bound(const __global int* imove,
                          const __global float* m,
                          const __global vec* dudt_preelastic,
                          const __global vec* dudt_elastic,
                          __global vec* force_elastic,
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
