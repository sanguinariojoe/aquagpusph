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
 * @brief Boundary integral term for constant fields.
 */

#include "resources/Scripts/types/types.h"
#include "resources/Scripts/KernelFunctions/Kernel.h"

/** @brief Compute the boundary integrals.
 * @param imove Moving flags.
 *   - imove > 0 for regular fluid particles.
 *   - imove = 0 for sensors.
 *   - imove < 0 for boundary elements/particles.
 * @param r Position \f$ \mathbf{r} \f$.
 * @param normal Normal \f$ \mathbf{n} \f$.
 * @param u Velocity \f$ \mathbf{u} \f$.
 * @param m Area of the boundary element \f$ s \f$.
 * @param jhoc Head and tail of chains for each cell.
 * @param grad_w_bi Gradient of constant fields due to the boundary integral
 * \f$ \langle \nabla 1 \rangle^{\partial \Omega} \f$.
 * @param div_u Velocity divergence \f$ \nabla \cdot \mathbf{u} \f$. Actually
 * this is just the part that has to do with the boundary element velocity
 * @param N Number of particles.
 */
__kernel void entry(const __global int* restrict imove,
                    const __global vec* restrict r,
                    const __global vec* restrict normal,
                    const __global vec* restrict u,
                    const __global float* restrict m,
                    const __global svec2* restrict jhoc,
                    __global vec* restrict grad_w_bi,
                    __global float* restrict div_u_bi,
                    usize N)
{
    const usize i = get_global_id(0);
    const usize it = get_local_id(0);
    if(i >= N)
        return;
    if(imove[i] != 1)
        return;

    const vec_xyz r_i = r[i].XYZ;

    __private vec_xyz __grad_w = VEC_ZERO.XYZ;
    __private float __div_u = 0.f;

    FOR_NEIGHS(N, jhoc){
        if(imove[j] != -3)
            continue;
        const vec_xyz r_ij = r[j].XYZ - r_i;
        const float q = length(r_ij) / H;
        if(q >= SUPPORT)
            continue;

        {
            const vec_xyz n_j = normal[j].XYZ;  // Assumed outwarding oriented
            const vec_xyz u_j = u[j].XYZ;
            const float area_j = m[j];
            const vec_xyz grad_w = n_j.XYZ * (kernelW(q) * CONW * area_j);
            __grad_w += grad_w;
            __div_u -= dot(u_j, grad_w);
        }
    }END_FOR_NEIGHS()

    grad_w_bi[i].XYZ = __grad_w;
    div_u_bi[i] = __div_u;
}

/** @brief Compute the pressure on each boundary element.
 * @param imove Moving flags.
 *   - imove > 0 for regular fluid particles.
 *   - imove = 0 for sensors.
 *   - imove < 0 for boundary elements/particles.
 * @param r Position \f$ \mathbf{r} \f$.
 * @param m Particle mass \f$ m \f$.
 * @param rho Particle density \f$ \rho \f$.
 * @param jhoc Head and tail of chains for each cell.
 * @param p Particle pressure \f$ p \f$.
 * @param N Number of particles.
 */
__kernel void p_boundary(const __global int* restrict imove,
                         const __global vec* restrict r,
                         const __global float* restrict m,
                         const __global float* restrict rho,
                         const __global svec2* restrict jhoc,
                         __global float* restrict p,
                         usize N)
{
    const usize i = get_global_id(0);
    const usize it = get_local_id(0);
    if(i >= N)
        return;
    if(imove[i] != -3)
        return;

    const vec_xyz r_i = r[i].XYZ;

    __private float __p = 0.f;

    FOR_NEIGHS(N, jhoc){
        if(imove[j] != 1)
            continue;
        const vec_xyz r_ij = r[j].XYZ - r_i;
        const float q = length(r_ij) / H;
        if(q >= SUPPORT)
            continue;

        __p += 2.f * p[j] * kernelW(q) * CONW * m[j] / rho[j];
    }END_FOR_NEIGHS()

    p[i] = __p;
}
