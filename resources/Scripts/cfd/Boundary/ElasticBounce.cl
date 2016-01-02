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

#ifndef HAVE_3D
    #include "../../types/2D.h"
    #include "../../KernelFunctions/Wendland2D.hcl"
#else
    #include "../../types/3D.h"
    #include "../../KernelFunctions/Wendland3D.hcl"
#endif

#ifndef __DR_FACTOR__
    /** @def __DR_FACTOR__
     * @brief The boundary elements effect is restricted to a quadrangular area
     * of \f$ R \times R \f$, where \f$ R = DR_FACTOR \cdot \Delta r \f$.
     */
    #define __DR_FACTOR__ 1.5f
#endif

#ifndef __MIN_BOUND_DIST__
    /** @def __MIN_BOUND_DIST__
     * @brief The elastic bounce is not tolerating that a particle becomes
     * closer than this distance (multiplied by \f$ \Delta r \f$).
     */
    #define __MIN_BOUND_DIST__ 0.3f
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
    #define __ELASTIC_FACTOR__ 0.0f
#endif

/** @brief Performs the boundary effect on the fluid particles.
 * @param imove Moving flags.
 *   - imove > 0 for regular fluid particles.
 *   - imove = 0 for sensors.
 *   - imove < 0 for boundary elements/particles.
 * @param r Position \f$ \mathbf{r} \f$.
 * @param normal Normal \f$ \mathbf{n} \f$.
 * @param u Velocity \f$ \mathbf{u} \f$.
 * @param dudt Velocity rate of change
 * \f$ \left. \frac{d \mathbf{u}}{d t} \right\vert_{n+1} \f$.
 * @param icell Cell where each particle is located.
 * @param ihoc Head of chain for each cell (first particle found).
 * @param N Number of particles.
 * @param n_cells Number of cells in each direction
 * @param dr Distance between particle \f$ \Delta r \f$.
 * @param dt Time step \f$ \Delta t \f$.
 */
__kernel void entry(const __global int* imove,
                    const __global vec* r,
                    const __global vec* normal,
                    __global vec* u,
                    __global vec* dudt,
                    // Link-list data
                    __global uint *icell,
                    __global uint *ihoc,
                    // Simulation data
                    uint N,
                    uivec4 n_cells,
                    float dr,
                    float dt)
{
    const uint i = get_global_id(0);
    const uint it = get_local_id(0);
    if(i >= N)
        return;
    if(imove[i] != 1)
        return;

    const float R = __DR_FACTOR__ * dr;
    const vec_xyz r_i = r[i].XYZ;
    vec_xyz u_i = u[i].XYZ;
    vec_xyz dudt_i = dudt[i].XYZ;

    BEGIN_LOOP_OVER_NEIGHS(){
        if((imove[j] != -2) && (imove[j] != -3)){
            j++;
            continue;
        }
        const vec_xyz r_ij = r[j].XYZ - r_i;
        const vec_xyz n_j = normal[j].XYZ;
        const float r0 = dot(r_ij, n_j);
        if(r0 < 0.f){
            // The particle is in the "wrong" side of the wall.
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
            const float u_n = dot(u_i - u[j].XYZ, n_j);
            const float dudt_n = dot(dudt_i - dudt[j].XYZ, n_j);
            const float dist = dt * u_n + 0.5f * dt * dt * dudt_n;
            if(dist < 0.f){
                // The particle is already running away from the boundary
                j++;
                continue;
            }

            // ------------------------------------------------------------------
            // The particle should be corrected if:
            //   - It is already placed in the effect zone.
            //   - It is entering inside the effect zone.
            // ------------------------------------------------------------------
            if(r0 - dist <= __MIN_BOUND_DIST__ * dr){
                // ------------------------------------------------------------------
                // Reflect particle velocity (using elastic factor)
                // ------------------------------------------------------------------
                dudt[i].XYZ = dudt_i - (1.f + __ELASTIC_FACTOR__) * dudt_n * n_j;
                u[i].XYZ = u_i - (1.f + __ELASTIC_FACTOR__) * u_n * n_j;

                // Modify the value for the next walls test.
                u_i = u[i].XYZ;
                dudt_i = dudt[i].XYZ;
            }
        }
    }END_LOOP_OVER_NEIGHS()
}