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
 * @brief Fluid-Ghost particles interactions computation.
 */

#include "resources/Scripts/types/types.h"
#include "resources/Scripts/KernelFunctions/Kernel.h"

#if __LAP_FORMULATION__ == __LAP_MONAGHAN__
    #ifndef HAVE_3D
        #define __CLEARY__ 8.f
    #else
        #define __CLEARY__ 10.f
    #endif
#endif

/** @brief Fluid particles interactions with the ghost ones.
 *
 * Compute the differential operators except the Laplacian of the velocity
 * already computed, taking into account just the fluid-ghost interactions.
 *
 * @param iset Set of particles index.
 * @param imove Moving flags.
 *   - imove > 0 for regular fluid particles.
 *   - imove = 0 for sensors.
 *   - imove < 0 for boundary elements/particles.
 * @param r Position \f$ \mathbf{r} \f$.
 * @param u Velocity \f$ \mathbf{u} \f$.
 * @param rho Density \f$ \rho \f$.
 * @param m Mass \f$ m \f$.
 * @param p Pressure \f$ p \f$.
 * @param jhoc Head and tail of chains for each cell.
 * @param refd Density of reference \f$ \rho_0 \f$ (one per set of particles)
 * @param grad_p Pressure gradient \f$ \frac{\nabla p}{rho} \f$.
 * @param div_u Velocity divergence \f$ \rho \nabla \cdot \mathbf{u} \f$.
 * @param N Number of particles.
 * @param g Gravity acceleration \f$ \mathbf{g} \f$.
 */
__kernel void entry(const __global uint* restrict iset,
                    const __global int* restrict imove,
                    const __global vec* restrict r,
                    const __global vec* restrict u,
                    const __global float* restrict rho,
                    const __global float* restrict m,
                    const __global float* restrict p,
                    const __global svec2* restrict jhoc,
                    __constant float* restrict refd,
                    __global vec* restrict grad_p,
                    __global float* restrict div_u,
                    usize N,
                    vec g)
{
    const usize i = get_global_id(0);
    const usize it = get_local_id(0);
    if(i >= N)
        return;
    if(imove[i] != 1){
        return;
    }

    const vec_xyz r_i = r[i].XYZ;
    const vec_xyz u_i = u[i].XYZ;
    const float p_i = p[i];
    const float rho_i = rho[i];
    const float refd_i = refd[iset[i]];

    __private vec_xyz __grad_p = VEC_ZERO.XYZ;
    __private float __div_u = 0.f;

    FOR_NEIGHS(N, jhoc){
        if(imove[j] != -1)
            continue;
        const vec_xyz r_ij = r[j].XYZ - r_i;
        const float q = length(r_ij) / H;
        if(q >= SUPPORT)
            continue;
        {
            const float rho_j = rho[j];
            const float p_j = p[j];
            const float udr = dot(u[j].XYZ - u_i, r_ij);
            const float f_ij = kernelF(q) * CONF * m[j];

            __grad_p += (p_i + p_j) / (rho_i * rho_j) * f_ij * r_ij;
            __div_u += udr * f_ij * rho_i / rho_j;
        }
    }END_FOR_NEIGHS()

    grad_p[i].XYZ += __grad_p;
    div_u[i] += __div_u;
}
