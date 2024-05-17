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
 * @brief Boundary integral friction term.
 */

#if defined(LOCAL_MEM_SIZE) && defined(NO_LOCAL_MEM)
    #error NO_LOCAL_MEM has been set.
#endif

#include "resources/Scripts/types/types.h"
#include "resources/Scripts/KernelFunctions/Kernel.h"

#if __LAP_FORMULATION__ != __LAP_MORRIS__ && \
    __LAP_FORMULATION__ != __LAP_MONAGHAN__
    #error Unknown Laplacian formulation: __LAP_FORMULATION__
#endif

#if __LAP_FORMULATION__ == __LAP_MONAGHAN__
    #ifndef HAVE_3D
        #define __CLEARY__ 8.f
    #else
        #define __CLEARY__ 10.f
    #endif
#endif

/** @brief Performs the boundary friction effect on the fluid particles.
 * @param iset Set of particles index.
 * @param imove Moving flags.
 *   - imove > 0 for regular fluid particles.
 *   - imove = 0 for sensors.
 *   - imove < 0 for boundary elements/particles.
 * @param r Position \f$ \mathbf{r} \f$.
 * @param normal Normal \f$ \mathbf{n} \f$.
 * @param u Velocity \f$ \mathbf{u} \f$.
 * @param rho Density \f$ \rho \f$.
 * @param m Area of the boundary element \f$ s \f$.
 * @param lap_u Velocity laplacian \f$ \frac{\Delta \mathbf{u}}{rho} \f$.
 * @param icell Cell where each particle is located.
 * @param ihoc Head of chain for each cell (first particle found).
 * @param N Number of particles.
 * @param n_cells Number of cells in each direction
 * @param noslip_iset Particles set of the boundary terms which friction should
 * be taken into account.
 * @param dr Distance between particles \f$ \Delta r \f$.
 */
__kernel void entry(const __global uint* iset,
                    const __global int* imove,
                    const __global vec* r,
                    const __global vec* normal,
                    const __global vec* u,
                    const __global float* rho,
                    const __global float* m,
                    __global vec* lap_u,
                    usize N,
                    uint noslip_iset,
                    float dr,
                    LINKLIST_LOCAL_PARAMS)
{
    const usize i = get_global_id(0);
    const usize it = get_local_id(0);
    if(i >= N)
        return;
    if(imove[i] != 1)
        return;

    const vec_xyz r_i = r[i].XYZ;
    const vec_xyz u_i = u[i].XYZ;
    const float rho_i = rho[i];

    // Initialize the output
    #ifndef LOCAL_MEM_SIZE
        #define _LAPU_ lap_u[i].XYZ
    #else
        #define _LAPU_ lap_u_l[it]
        __local vec_xyz lap_u_l[LOCAL_MEM_SIZE];
        _LAPU_ = lap_u[i].XYZ;
    #endif

    const usize c_i = icell[i];
    BEGIN_NEIGHS(c_i, N, n_cells, icell, ihoc){
        if((imove[j] != -3) || (iset[j] != noslip_iset)){
            j++;
            continue;
        }
        const vec_xyz r_ij = r[j].XYZ - r_i;
        const float q = length(r_ij) / H;
        if(q >= SUPPORT)
        {
            j++;
            continue;
        }

        {
            const vec_xyz n_j = normal[j].XYZ;  // Assumed outwarding oriented
            const float area_j = m[j];

            const float w_ij = kernelW(q) * CONW * area_j;
            const vec_xyz du = u[j].XYZ - u_i;

            #if __LAP_FORMULATION__ == __LAP_MONAGHAN__
                const float r2 = (q * q + 0.01f) * H * H;
                _LAPU_ += __CLEARY__ * w_ij * dot(du, r_ij) / (r2 * rho_i) * n_j;
            #endif
            #if __LAP_FORMULATION__ == __LAP_MORRIS__ || \
                __LAP_FORMULATION__ == __LAP_MONAGHAN__
                const float dr_n = max(fabs(dot(r_ij, n_j)), dr);
                const vec_xyz du_t = du - dot(du, n_j) * n_j;
                _LAPU_ += 2.f * w_ij / (rho_i * dr_n) * du_t;
            #endif
        }
    }END_NEIGHS()

    #ifdef LOCAL_MEM_SIZE
        lap_u[i].XYZ = _LAPU_;
    #endif
}
