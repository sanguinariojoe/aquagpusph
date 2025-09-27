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
 * @brief Fixed ghost particles fields interpolation.
 */

#include "resources/Scripts/types/types.h"
#include "resources/Scripts/KernelFunctions/Kernel.h"

/** @brief Fixed ghost particles fields interpolation.
 *
 * In order to conveniently extend the flow fields we need to know the mirrored
 * values. Also the Shepard values should be recomputed taking into account the
 * mirrored position.
 *
 * @param imove Moving flags.
 *   - imove > 0 for regular fluid particles.
 *   - imove = 0 for sensors.
 *   - imove < 0 for boundary elements/particles.
 * @param r Position \f$ \mathbf{r} \f$.
 * @param normal Normal \f$ \mathbf{n} \f$.
 * @param m Mass \f$ m \f$.
 * @param rho Density \f$ \rho \f$.
 * @param p Pressure \f$ p \f$.
 * @param u Velocity \f$ \mathbf{u} \f$.
 * @param gp_rho Interpolated density in the mirrored position \f$ \rho \f$.
 * @param gp_p Interpolated pressure in the mirrored position \f$ p \f$.
 * @param gp_u Interpolated velocity in the mirrored position \f$ \mathbf{u} \f$.
 * @param shepard Shepard term
 * \f$ \gamma(\mathbf{x}) = \int_{\Omega}
 *     W(\mathbf{y} - \mathbf{x}) \mathrm{d}\mathbf{x} \f$.
 * @param icell Cell where each particle is located.
 * @param gp_icell Cell where each mirrored ghost particle is located.
 * @param ihoc Head of chain for each cell (first particle found).
 * @param N Number of particles.
 * @param n_cells Number of cells in each direction
 */
__kernel void entry(const __global int* restrict imove,
                    const __global vec* restrict r,
                    const __global vec* restrict normal,
                    const __global float* restrict m,
                    const __global float* restrict rho,
                    const __global float* restrict p,
                    const __global vec* restrict u,
                    __global float* restrict gp_rho,
                    __global float* restrict gp_p,
                    __global vec* restrict gp_u,
                    __global float* restrict shepard,
                    usize N,
                    const __global usize* restrict gp_icell,
                    LINKLIST_LOCAL_PARAMS)
{
    const usize i = get_global_id(0);
    const usize it = get_local_id(0);
    if(i >= N)
        return;
    if(imove[i] != -1)
        return;
    
    const vec_xyz r_i = r[i].XYZ;

    __private float __rho = 0.f;
    __private float __p = 0.f;
    __private vec_xyz __u = VEC_ZERO.XYZ;
    __private float __shepard = 0.f;

    #undef C_I()
    #define C_I() const usize c_i = gp_icell[i]
    BEGIN_LOOP_OVER_NEIGHS(){
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
            const float rho_j = rho[j];
            const float p_j = p[j];
            const vec_xyz u_j = u[j];
            const float m_j = m[j];

            {
                const float w_ij = kernelW(q) * CONW * m_j / rho_j;
                __rho += w_ij * rho_j;
                __p += w_ij * p_j; 
                __u += w_ij * u_j;
                __shepard += w_ij;
            }
        }
    }END_LOOP_OVER_NEIGHS()

    gp_rho[i] = __rho;
    gp_p[i] = __p;
    gp_u[i].XYZ = __u;
    shepard[i] = __shepard;
}
