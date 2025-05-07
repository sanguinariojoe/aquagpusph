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
 * @brief Fluid particles interactions computation.
 */

#if defined(LOCAL_MEM_SIZE) && defined(NO_LOCAL_MEM)
    #error NO_LOCAL_MEM has been set.
#endif

#include "resources/Scripts/types/types.h"
#include "resources/Scripts/KernelFunctions/Kernel.h"

/** @brief Pressure interpolation at the boundary elements.
 *
 * The values are computed using just the fluid information. The resulting
 * interpolated values are not renormalized yet.
 *
 * Just the elements with the flag imove = -3 are considered boundary elements.
 *
 * @param iset Set of particles index.
 * @param imove Moving flags.
 *   - imove > 0 for regular fluid particles.
 *   - imove = 0 for sensors.
 *   - imove < 0 for boundary elements/particles.
 * @param r Position \f$ \mathbf{r} \f$.
 * @param m Mass \f$ m \f$.
 * @param rho Density \f$ \rho \f$.
 * @param grad_p Pressure gradient \f$ \frac{\nabla p}{rho} \f$.
 * @param p Pressure \f$ p \f$.
 * @param icell Cell where each particle is located.
 * @param ihoc Head of chain for each cell (first particle found).
 * @param refd Density of reference of the fluid \f$ \rho_0 \f$.
 * @param N Number of particles.
 * @param n_cells Number of cells in each direction
 * @see SensorsRenormalization.cl
 */
__kernel void entry(const __global uint* iset,
                    const __global int* imove,
                    const __global vec* r,
                    const __global float* m,
                    const __global float* rho,
                    const __global vec* grad_p,
                    __global float* p,
                    __constant float* refd,
                    usize N,
                    LINKLIST_LOCAL_PARAMS)
{
    const usize i = get_global_id(0);
    const usize it = get_local_id(0);
    if(i >= N)
        return;
    if(imove[i] != -3){
        return;
    }

    const vec_xyz r_i = r[i].XYZ;
    const vec_xyz gradp_i = rho[i] * grad_p[i].XYZ;

    // Initialize the output
    #ifndef LOCAL_MEM_SIZE
        #define _P_ p[i]
    #else
        #define _P_ p_l[it]
        __local float p_l[LOCAL_MEM_SIZE];
    #endif
    _P_ = 0.f;

    const usize c_i = icell[i];
    BEGIN_NEIGHS(c_i, N, n_cells, icell, ihoc){
        if(imove[j] != 1){
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
            const float w_ij = kernelW(q) * CONW * m[j] / rho[j];
            _P_ += (p[j] - dot(gradp_i, r_ij)) * w_ij;
        }
    }END_NEIGHS()

    #ifdef LOCAL_MEM_SIZE
        p[i] = _P_;
    #endif
}
