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
 * @brief Computation of the order one Morris correction term, due to the
 * interactions of fluid particles with boundary elements.
 */

#if defined(LOCAL_MEM_SIZE) && defined(NO_LOCAL_MEM)
    #error NO_LOCAL_MEM has been set.
#endif

#include "resources/Scripts/types/types.h"
#include "resources/Scripts/KernelFunctions/Kernel.h"

/** @brief Computation of the order one Morris correction term, due to the
 * interactions of fluid particles with boundary elements.
 *
 * The first order, divergent term, in the Morris expression is featured by the
 * following term:
 * \f$ \int_{\Omega} (\mathbf{y} - \mathbf{x})
 * \frac{\mathbf{y} - \mathbf{x}}{\vert \mathbf{y} - \mathbf{x} \vert^2}
 * \cdot \nabla W (\mathbf{y} - \mathbf{x}) d\mathbf{y}
 * \f$
 *
 * @param iset Set of particles index.
 * @param imove Moving flags.
 *   - imove > 0 for regular fluid particles.
 *   - imove = 0 for sensors.
 *   - imove < 0 for boundary elements/particles.
 * @param r Position \f$ \mathbf{r} \f$.
 * @param normal Normal \f$ \mathbf{n} \f$.
 * @param rho Density \f$ \rho \f$.
 * @param m Area of the boundary element \f$ s \f$.
 * @param morris_o1 Morris order 1 correction term.
 * @param icell Cell where each particle is located.
 * @param ihoc Head of chain for each cell (first particle found).
 * @param N Number of particles.
 * @param n_cells Number of cells in each direction
 * @param noslip_iset Particles set of the boundary terms which friction should
 * be taken into account.
 */
__kernel void order1(const __global uint* iset,
                     const __global int* imove,
                     const __global vec* r,
                     const __global vec* normal,
                     const __global float* rho,
                    const __global float* m,
                     __global vec* morris_o1,
                     // Link-list data
                     const __global uint *icell,
                     const __global uint *ihoc,
                     // Simulation data
                     uint N,
                     uivec4 n_cells,
                     uint noslip_iset)
{
    const uint i = get_global_id(0);
    const uint it = get_local_id(0);
    if(i >= N)
        return;
    if(imove[i] != 1){
        return;
    }

    const vec_xyz r_i = r[i].XYZ;

    // Initialize the output
    #ifndef LOCAL_MEM_SIZE
        #define _MORRISO1_ morris_o1[i]
    #else
        #define _MORRISO1_ morris_o1_l[it]
        __local vec_xyz morris_o1_l[LOCAL_MEM_SIZE];
        _MORRISO1_ = morris_o1[i].XYZ;
    #endif

    BEGIN_LOOP_OVER_NEIGHS(){
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

            const float r2 = (q * q + 0.01f) * H * H;
            _MORRISO1_ += 2.f * w_ij * dot(r_ij, n_j) / r2 * r_ij;
        }
    }END_LOOP_OVER_NEIGHS()

    #ifdef LOCAL_MEM_SIZE
        morris_o1[i].XYZ = _MORRISO1_;
    #endif
}
