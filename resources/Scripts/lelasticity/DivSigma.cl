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

/** @addtogroup lela
 * @{
 */

/** @file
 * @brief Solid particles interactions computation.
 */

#if defined(LOCAL_MEM_SIZE) && defined(NO_LOCAL_MEM)
    #error NO_LOCAL_MEM has been set.
#endif

#include "resources/Scripts/types/types.h"
#include "resources/Scripts/KernelFunctions/Kernel.h"

/** @brief Solid particles interactions computation.
 *
 * Compute the differential operators involved in the numerical scheme, taking
 * into account just the solid-solid interactions.
 *
 * @param imove Moving flags.
 *   - imove = 2 for regular solid particles.
 *   - imove = 0 for sensors (ignored by this module).
 *   - imove < 0 for boundary elements/particles.
 * @param r Position \f$ \mathbf{r} \f$.
 * @param rho Density \f$ \rho \f$.
 * @param m Mass \f$ m \f$.
 * @param sigma Stress tensor \f$ \sigma \f$.
 * @param div_sigma Divergence of the stress tensor
 * 	   \f$ \frac{\nabla \cdot \sigma}{rho} \f$.
 * @param icell Cell where each particle is located.
 * @param ihoc Head of chain for each cell (first particle found).
 * @param N Number of particles.
 * @param n_cells Number of cells in each direction.
 */
__kernel void entry(const __global int* imove,
                    const __global vec* r,
                    const __global float* rho,
                    const __global float* m,
                    const __global matrix* sigma,
                    __global vec* div_sigma,
                    // Link-list data
                    const __global uint *icell,
                    const __global uint *ihoc,
                    // Simulation data
                    uint N,
                    uivec4 n_cells)
{
    const uint i = get_global_id(0);
    const uint it = get_local_id(0);
    if(i >= N)
        return;
    if(imove[i] != 2){
        return;
    }

    const vec_xyz r_i = r[i].XYZ;
    const matrix s_i = sigma[i];
    const float rho_i = rho[i];

    // Initialize the output
    #ifndef LOCAL_MEM_SIZE
        #define _DIVS_ div_sigma[i].XYZ
    #else
        #define _DIVS_ div_sigma_l[it]
        __local vec_xyz div_sigma_l[LOCAL_MEM_SIZE];
        _DIVS_ = VEC_ZERO.XYZ;
    #endif

    BEGIN_LOOP_OVER_NEIGHS(){
        if((i == j) || (imove[j] != 2)){
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
            const float rho_j = rho[j];
            const matrix s_j = sigma[j];
            const float f_ij = kernelF(q) * CONF * m[j];

            _DIVS_ += MATRIX_DOT((s_i + s_j), f_ij / (rho_i * rho_j) * r_ij);
        }
    }END_LOOP_OVER_NEIGHS()

    #ifdef LOCAL_MEM_SIZE
        div_sigma[i].XYZ = _DIVS_;
    #endif
}

/*
 * @}
 */