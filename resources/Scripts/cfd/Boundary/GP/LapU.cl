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
 * @brief Ghost particles Laplacian of the velocity computation.
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

/** @brief Ghost particles Laplacian of the velocity computation.
 *
 * This method compute the Laplacian of the velocity also in the mirroring
 * boundary elements because such information will be required to extend the
 * pressure field.
 *
 * @param imove Moving flags.
 *   - imove > 0 for regular fluid particles.
 *   - imove = 0 for sensors.
 *   - imove < 0 for boundary elements/particles.
 * @param r Position \f$ \mathbf{r} \f$.
 * @param u Velocity \f$ \mathbf{u} \f$.
 * @param rho Density \f$ \rho \f$.
 * @param m Mass \f$ m \f$.
 * @param jhoc Head and tail of chains for each cell.
 * @param lap_u Velocity laplacian \f$ \frac{\Delta \mathbf{u}}{rho} \f$.
 * @param N Number of particles.
 */
__kernel void entry(const __global int* restrict imove,
                    const __global vec* restrict r,
                    const __global vec* restrict u,
                    const __global float* restrict rho,
                    const __global float* restrict m,
                    const __global svec2* restrict jhoc,
                    __global vec* restrict lap_u,
                    usize N)
{
    const usize i = get_global_id(0);
    const usize it = get_local_id(0);
    if(i >= N)
        return;
    const int imove_i = imove[i];
    if((imove_i != 1) && (imove_i != -2)){
        return;
    }

    const vec_xyz r_i = r[i].XYZ;
    const vec_xyz u_i = u[i].XYZ;
    const float rho_i = rho[i];

    __private vec_xyz __lap_u = 0.f;

    FOR_NEIGHS(N, jhoc){
        const int imove_j = imove[j];
        if((imove_j != 1) && (imove_j != -1))
            continue;
        // Don't take into account the fluid-fluid interactions
        if(imove_i == imove_j)
            continue;
        const vec_xyz r_ij = r[j].XYZ - r_i;
        const float q = length(r_ij) / H;
        if(q >= SUPPORT)
            continue;
        {
            const float rho_j = rho[j];
            const float udr = dot(u[j].XYZ - u_i, r_ij);
            const float f_ij = kernelF(q) * CONF * m[j];

            #if __LAP_FORMULATION__ == __LAP_MONAGHAN__
                const float r2 = (q * q + 0.01f) * H * H;
                __lap_u += f_ij * __CLEARY__ * udr / (r2 * rho_i * rho_j) * r_ij;
            #elif __LAP_FORMULATION__ == __LAP_MORRIS__
                __lap_u += f_ij * 2.f / (rho_i * rho_j) * (u[j].XYZ - u_i);
            #else
                #error Unknown Laplacian formulation: __LAP_FORMULATION__
            #endif
        }
    }END_FOR_NEIGHS()

    lap_u[i].XYZ += __lap_u;
}
