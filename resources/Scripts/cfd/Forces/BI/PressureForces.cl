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

#if defined(LOCAL_MEM_SIZE) && defined(NO_LOCAL_MEM)
    #error NO_LOCAL_MEM has been set.
#endif

#include "resources/Scripts/types/types.h"
#include "resources/Scripts/KernelFunctions/Kernel.h"

/** @brief Tool to compute the pressure force and moment for an especific body.
 *
 * The force at each boundary element is just
 * \f$ \mathbf{f}_a = p_a \mathbf{n}_a s_a \f$
 * where \f$ s_a \f$ is the area of the element, stored in the masses array.
 * The moment is consequently computed as:
 * \f$ \mathbf{m}_a  = \mathbf{f}_a \times
 *     \left(\mathbf{r}_a - \mathbf{r}_0 \right) \f$
 * with \f$ \mathbf{r}_0 \f$ the reference point where the moment should be
 * computed.
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
 * @param normal Normal \f$ \mathbf{n} \f$.
 * @param p Pressure \f$ p \f$.
 * @param rho Density \f$ \rho \f$.
 * @param m Mass \f$ m \f$.
 * @param N Number of particles.
 * @param pressureForces_iset Particles set to be computed.
 * @param pressureForces_r Point with respect the moments are computed,
 * \f$ \mathbf{r}_0 \f$.
 */
__kernel void entry(__global vec* pressureForces_f,
                    __global vec4* pressureForces_m,
                    const __global uint* iset,
                    const __global int* imove,
                    const __global vec* r,
                    const __global vec* normal,
                    const __global float* p,
                    const __global float* rho,
                    const __global float* m,
                    uint N,
                    unsigned int pressureForces_iset,
                    vec pressureForces_r)
{
    // find position in global arrays
    const uint i = get_global_id(0);
    const uint it = get_local_id(0);
    if(i >= N)
        return;
    if((iset[i] != pressureForces_iset) || (imove[i] != -3)){
        pressureForces_f[i] = VEC_ZERO;
        pressureForces_m[i] = (vec4)(0.f, 0.f, 0.f, 0.f);
        return;
    }

    const vec_xyz force = p[i] * m[i] * normal[i].XYZ;

    pressureForces_f[i].XYZ = force;
    const vec_xyz arm = r[i].XYZ - pressureForces_r.XYZ;
    pressureForces_m[i].z = arm.x * force.y - arm.y * force.x;
    pressureForces_m[i].w = 0.f;
    #ifdef HAVE_3D
        pressureForces_f[i].w = 0.f;
        pressureForces_m[i].x = arm.y * force.z - arm.z * force.y;
        pressureForces_m[i].y = arm.z * force.x - arm.x * force.z;
    #else
        pressureForces_m[i].x = 0.f;
        pressureForces_m[i].y = 0.f;
    #endif
}
