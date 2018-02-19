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
 * @brief Generalized local coordinates renormalization factor, zeta.
 */

#if defined(LOCAL_MEM_SIZE) && defined(NO_LOCAL_MEM)
    #error NO_LOCAL_MEM has been set.
#endif

#include "resources/Scripts/types/types.h"
#include "resources/Scripts/KernelFunctions/Kernel.h"

/** @brief Zeta renormalization factor computation.
 *
 * \f[ \zeta(\mathbf{x}) = \int_{\partical \Omega}
 *     W(\mathbf{y} - \mathbf{x}) \mathrm{d}\mathbf{y} \f]
 *
 * Zeta renromalization factor is used just for the generalized local
 * coordinates integration close to the boundaries.
 *
 * @param imove Moving flags.
 *   - imove > 0 for regular fluid particles.
 *   - imove = 0 for sensors.
 *   - imove < 0 for boundary elements/particles.
 * @param r Position \f$ \mathbf{r} \f$.
 * @param m Area of the boundary element \f$ s \f$.
 * @param zeta Zeta term
 * \f$ \zeta(\mathbf{x}) = \int_{\partial \Omega}
 *     W(\mathbf{y} - \mathbf{x}) \mathrm{d}\mathbf{y} \f$.
 * @param icell Cell where each particle is located.
 * @param ihoc Head of chain for each cell (first particle found).
 * @param N Number of particles.
 * @param n_cells Number of cells in each direction
 */
__kernel void compute(const __global int* imove,
                      const __global vec* r,
                      const __global float* m,
                      __global float* zeta,
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
    if((imove[i] < -3) || (imove[i] > 1))
        return;

    const vec_xyz r_i = r[i].XYZ;

    // Initialize the output
    #ifndef LOCAL_MEM_SIZE
        #define _ZETA_ zeta[i]
    #else
        #define _ZETA_ zeta_l[it]
        __local float zeta_l[LOCAL_MEM_SIZE];
    #endif

    _ZETA_ = 0.f;

    BEGIN_LOOP_OVER_NEIGHS(){
        if(imove[j] != -3){
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

        const float area_j = m[j];
        _ZETA_ += kernelW(q) * CONW * m[j];
    }END_LOOP_OVER_NEIGHS()

    #ifdef LOCAL_MEM_SIZE
        zeta[i] = _ZETA_;
    #endif
}
