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
 * @brief Tool to compute the fluid pressure force and moment.
 */

#include "resources/Scripts/types/types.h"
#include "resources/Scripts/KernelFunctions/Kernel.h"

/** @brief Tool to compute the pressure force and moment for an especific body,
 * due to the Ghost particles.
 *
 * The pressure forces is computed using the same expression than in the
 * fluid-fluid interactions.
 *
 * @param pressureForces_f Force of each boundary element to be computed [N].
 * @param pressureForces_m Moment of each boundary element to be computed
 * [N \f$ \cdot \f$ m].
 * @param iset Set of particles index.
 * @param imove Moving flags.
 *   - imove > 0 for regular fluid particles.
 *   - imove = 0 for sensors.
 *   - imove < 0 for boundary elements/particles.
 * @param r Position \f$ \mathbf{r} \f$.
 * @param rho Density \f$ \rho \f$.
 * @param p Pressure \f$ p \f$.
 * @param m Mass \f$ m \f$.
 * @param jhoc Head and tail of chains for each cell.
 * @param N Number of particles.
 * @param pressureForces_iset Particles set to be computed.
 * @param pressureForces_r Point with respect the moments are computed
 * \f$ \mathbf{r}_0 \f$.
 */
__kernel void entry(__global vec* restrict pressureForces_f,
                    __global vec4* restrict pressureForces_m,
                    const __global uint* restrict iset,
                    const __global int* restrict imove,
                    const __global vec* restrict r,
                    const __global float* restrict rho,
                    const __global float* restrict p,
                    const __global float* restrict m,
                    const __global svec2* restrict jhoc,
                    usize N,
                    unsigned int pressureForces_iset,
                    vec pressureForces_r)
{
    // find position in global arrays
    const usize i = get_global_id(0);
    const usize it = get_local_id(0);
    if(i >= N)
        return;
    if((iset[i] != pressureForces_iset) || (imove[i] != -1)){
        pressureForces_f[i] = VEC_ZERO;
        pressureForces_m[i] = (vec4)(0.f, 0.f, 0.f, 0.f);
        return;
    }

    const vec_xyz r_i = r[i].XYZ;
    const float p_i = p[i];

    __private vec_xyz __f = VEC_ZERO.XYZ;

    FOR_NEIGHS(N, jhoc){
        if(imove[j] != 1)
            continue;
        const vec_xyz r_ij = r[j].XYZ - r_i;
        const float q = length(r_ij) / H;
        if(q >= SUPPORT)
            continue;

        {
            __f -= (p_i + p[j]) * kernelF(q) * CONF * m[j] / rho[j] * r_ij;
        }
    }END_FOR_NEIGHS()

    __f *= m[i] / rho[i];
    pressureForces_f[i].XYZ = __f;

    const vec arm = r[i] - pressureForces_r;
    pressureForces_m[i].z = arm.x * __f.y - arm.y * __f.x;
    pressureForces_m[i].w = 0.f;
    #ifdef HAVE_3D
        pressureForces_m[i].x = arm.y * __f.z - arm.z * __f.y;
        pressureForces_m[i].y = arm.z * __f.x - arm.x * __f.z;
    #else
        pressureForces_m[i].x = 0.f;
        pressureForces_m[i].y = 0.f;
    #endif
}
