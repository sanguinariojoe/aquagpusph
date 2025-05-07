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

#ifndef DIMS
#ifdef HAVE_3D
#define DIMS 3
#else
#define DIMS 2
#endif
#endif

/** @brief Ensure that no particle can get closer to the boundary than its
 * volume allows to.
 * @param imove Moving flags.
 *   - imove > 0 for regular fluid particles.
 *   - imove = 0 for sensors.
 *   - imove < 0 for boundary elements/particles.
 * @param r Position \f$ \mathbf{r} \f$.
 * @param normal Normal \f$ \mathbf{n} \f$.
 * @param m Mass \f$ m \f$ or Area \f$ s \f$ (depending on @p imove).
 * @param rho Density \f$ \rho \f$.
 * @param N Number of particles.
 * @param icell Cell where each particle is located.
 * @param ihoc Head of chain for each cell (first particle found).
 * @param n_cells Number of cells in each direction
 */
__kernel void entry(const __global int* imove,
                    __global vec* r,
                    const __global vec* normal,
                    const __global float* m,
                    const __global float* rho,
                    usize N,
                    LINKLIST_LOCAL_PARAMS)
{
    const usize i = get_global_id(0);
    const usize it = get_local_id(0);
    if(i >= N)
        return;
    if(imove[i] != 1)
        return;

    const float Ri = 0.5f * pow(m[i] / rho[i], 1.f / DIMS);

    const usize c_i = icell[i];
    BEGIN_NEIGHS(c_i, N, n_cells, icell, ihoc){
        if (imove[j] != -3) {
            j++;
            continue;
        }
        const vec_xyz r_ij = r[j].XYZ - r[i].XYZ;
        const vec_xyz n_j = normal[j].XYZ;
        const float rn = dot(r_ij, n_j);
        if(fabs(rn) > Ri){
            // Either the particle is just OK, or it is too far away from the
            // boundary to consider it
            j++;
            continue;
        }

#ifdef HAVE_3D
        const float dr = sqrt(m[j]);
#else
        const float dr = m[j];
#endif
        const float Rj = __DR_FACTOR__ * dr;
        const vec_xyz rt = r_ij - rn * n_j;
        if(dot(rt, rt) >= Rj * Rj){
            // The particle is passing too far from the boundary element
            j++;
            continue;
        }

        r[i].XYZ += (rn - Ri) * n_j;
    }END_NEIGHS()
}
