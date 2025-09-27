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

/** @brief Fields interpolation at the sensors.
 *
 * The fields to be interpolated are:
 *    - Velocity \f$ \mathbf{u} \f$
 *    - Density \f$ \rho \f$
 *    - Pressure \f$ p \f$
 *
 * The values are computed using just the fluid information. The resulting
 * interpolated values are not renormalized yet.
 *
 * @param iset Set of particles index.
 * @param imove Moving flags.
 *   - imove > 0 for regular fluid particles.
 *   - imove = 0 for sensors.
 *   - imove < 0 for boundary elements/particles.
 * @param r Position \f$ \mathbf{r} \f$.
 * @param m Mass \f$ m \f$.
 * @param jhoc Head and tail of chains for each cell.
 * @param u Velocity \f$ \mathbf{u} \f$.
 * @param rho Density \f$ \rho \f$.
 * @param p Pressure \f$ p \f$.
 * @param N Number of particles.
 * @param g Gravity acceleration \f$ \mathbf{g} \f$.
 * @see SensorsRenormalization.cl
 */
__kernel void entry(const __global uint* restrict iset,
                    const __global int* restrict imove,
                    const __global vec* restrict r,
                    const __global float* restrict m,
                    const __global svec2* restrict jhoc,
                    __global vec* restrict u,
                    __global float* restrict rho,
                    __global float* restrict p,
                    usize N,
                    vec g)
{
    const usize i = get_global_id(0);
    const usize it = get_local_id(0);
    if(i >= N)
        return;
    if(imove[i] != 0){
        return;
    }

    const vec_xyz r_i = r[i].XYZ;

    __private vec_xyz __u = VEC_ZERO.XYZ;
    __private float __rho = 0.f;
    __private float __p = 0.f;

    FOR_NEIGHS(N, jhoc){
        if(i == j)
            continue;
        if(imove[j] != 1)
            continue;
        const vec_xyz r_ij = r[j].XYZ - r_i;
        const float q = length(r_ij) / H;
        if(q >= SUPPORT)
            continue;

        {
            const float rho_j = rho[j];
            const float m_j = m[j];
            const float p_j = p[j];
            const vec_xyz u_j = u[j].XYZ;
            const float w_ij = kernelW(q) * CONW * m_j / rho_j;

            __u += u_j * w_ij;
            __rho += rho_j * w_ij;
            __p += p_j * w_ij;
        }
    }END_FOR_NEIGHS()

    u[i].XYZ = __u;
    rho[i] = __rho;
    p[i] = __p;
}
