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
    #define __DR_FACTOR__ 0.6f
#endif

#ifndef __MIN_BOUND_DIST__
    /** @def __MIN_BOUND_DIST__
     * @brief The elastic bounce is not tolerating that a particle becomes
     * closer than this distance (multiplied by \f$ \Delta r \f$).
     */
    #define __MIN_BOUND_DIST__ 0.f
#endif

#ifndef __ELASTIC_FACTOR__
    /** @def __ELASTIC_FACTOR__
     * @brief The amount of kinetic energy conserved in the interaction.
     *
     * A factor of 1 imply that the velocity of the particle will be preserved
     * (except for the direction), while a factor of 0 imply that the particle
     * will loss all its normal to the boundary velocity.
     *
     * The tangential velocity is not affected.
     */
    #define __ELASTIC_FACTOR__ 1.0f
#endif

/** @brief Performs the boundary effect on the fluid particles.
 * @param imove Moving flags.
 *   - imove > 0 for regular fluid particles.
 *   - imove = 0 for sensors.
 *   - imove < 0 for boundary elements/particles.
 * @param r Position \f$ \mathbf{r} \f$.
 * @param normal Normal \f$ \mathbf{n} \f$.
 * @param u_in Velocity \f$ \mathbf{u} \f$.
 * @param dudt Velocity rate of change
 * \f$ \left. \frac{d \mathbf{u}}{d t} \right\vert_{n+1} \f$.
 * @param N Number of particles.
 * @param dr Distance between particles \f$ \Delta r \f$.
 * @param dt Time step \f$ \Delta t \f$.
 * @param icell Cell where each particle is located.
 * @param ihoc Head of chain for each cell (first particle found).
 * @param n_cells Number of cells in each direction
 */
__kernel void entry(const __global int* imove,
                    const __global vec* r,
                    const __global vec* normal,
                    const __global vec* u_in,
                    __global vec* dudt,
                    uint N,
                    float dr,
                    float dt,
                    LINKLIST_LOCAL_PARAMS)
{
    const uint i = get_global_id(0);
    const uint it = get_local_id(0);
    if(i >= N)
        return;
    if(imove[i] != 1)
        return;

    const float R = __DR_FACTOR__ * dr;
    const vec_xyz r_i = r[i].XYZ;
    vec_xyz dudt_i = dudt[i].XYZ;
    const vec_xyz u_i = u_in[i].XYZ + 0.5f * dt * dudt_i;

    const unsigned int c_i = icell[i];
    BEGIN_NEIGHS(c_i, N, n_cells, icell, ihoc){
        if (imove[j] != -3) {
            j++;
            continue;
        }
        const vec_xyz r_ij = r[j].XYZ - r_i;
        const vec_xyz n_j = normal[j].XYZ;
        const float r0 = dot(r_ij, n_j);
        if(r0 < 0.f){
            // The particle is on the "wrong" side of the wall.
            j++;
            continue;
        }
        const vec_xyz rt = r_ij - r0 * n_j;
        if(dot(rt, rt) >= R * R){
            // The particle is passing too far from the boundary element
            j++;
            continue;
        }

        {
            if(dot(dudt_i, n_j) < 0.f){
                // The particle is accelerating away from the boundary, we
                // cannot do anything more
                j++;
                continue;
            }

            const float dist = dt * dot(u_i - u_in[j].XYZ, n_j);
            // ------------------------------------------------------------------
            // The particle should be corrected if:
            //   - It is already placed in the effect zone.
            //   - It is entering inside the effect zone.
            // ------------------------------------------------------------------
            if(r0 - dist <= __MIN_BOUND_DIST__ * dr){
                // ------------------------------------------------------------------
                // Reflect particle acceleration (using elastic factor)
                // ------------------------------------------------------------------
                dudt[i].XYZ = dudt_i -
                    (1.f + __ELASTIC_FACTOR__) * dot(dudt_i, n_j) * n_j;

                // Modify the value for the next walls test.
                dudt_i = dudt[i].XYZ;
            }
        }
    }END_NEIGHS()
}
