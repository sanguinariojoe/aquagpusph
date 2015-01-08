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

#ifndef HAVE_3D
    #include "../types/2D.h"
    #include "../KernelFunctions/Wendland2D.hcl"
#else
    #include "../types/3D.h"
    #include "../KernelFunctions/Wendland3D.hcl"
#endif

#ifndef HAVE_3D
    #ifndef NEIGH_CELLS
        /** @def NEIGH_CELLS
         * @brief Number of neigh cells.
         *
         * In 2D cases 9 cells must be computed, while in 3D simulations 27
         * cells must be computed.
         */
        #define NEIGH_CELLS 9
    #endif
#else
    #ifndef NEIGH_CELLS
        /** @def NEIGH_CELLS
         * @brief Number of neigh cells.
         *
         * In 2D cases 9 cells must be computed, while in 3D simulations 27
         * cells must be computed.
         */
        #define NEIGH_CELLS 27
    #endif
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
 * @param pos Position \f$ \mathbf{r} \f$.
 * @param normal Normal \f$ \mathbf{n} \f$.
 * @param v Velocity \f$ \mathbf{u} \f$.
 * @param dvdt Velocity rate of change
 * \f$ \left. \frac{d \mathbf{u}}{d t} \right\vert_{n+1} \f$.
 * @param icell Cell where each particle is located.
 * @param ihoc Head of chain for each cell (first particle found).
 * @param N Number of particles.
 * @param n_cells Number of cells in each direction
 * @param dr Distance between particle \f$ \Delta r \f$.
 * @param dt Time step \f$ \Delta t \f$.
 */
__kernel void main(const __global int* imove,
                   const __global vec* pos,
                   const __global vec* normal,
                   __global vec* v,
                   __global vec* dvdt,
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
    if(imove[i] <= 0)
        return;

    const float R = __DR_FACTOR__ * dr;
    const uint c_i = icell[i];
    const vec_xyz pos_i = pos[i].XYZ;
    vec_xyz v_i = v[i].XYZ;
    vec_xyz dvdt_i = dvdt[i].XYZ;

    // Loop over neighbour particles
    // =============================
    for(int ci = -1; ci <= 1; ci++) {
        for(int cj = -1; cj <= 1; cj++) {
            #ifdef HAVE_3D
            for(int ck = -1; ck <= 1; ck++) {
            #else
            const int ck = 0; {
            #endif
                const uint c_j = c_i +
                                ci +
                                cj * n_cells.x +
                                ck * n_cells.x * n_cells.y;
                uint j = ihoc[c_j];
                while((j < N) && (icell[j] == c_j)) {
                    if(imove[j] >= 0){
                        j++;
                        continue;
                    }
                    const vec_xyz r = pos[j].XYZ - pos_i;
                    const vec_xyz n_j = normal[j].XYZ;
                    const float r0 = dot(r, n_j);
                    if(r0 < 0.f){
                        // The boundary element is not well oriented
                        j++;
                        continue;
                    }
                    const vec_xyz rt = r - r0 * n_j;
                    if(dot(rt, rt) >= R * R){
                        // The particle is passing too far from the boundary element
                        j++;
                        continue;
                    }

                    {
                        #include "ElasticBounce.hcl"
                    }
                    j++;
                }
            }
        }
    }
}
