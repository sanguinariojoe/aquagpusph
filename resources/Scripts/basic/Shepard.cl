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

/** @addtogroup basic
 * @{
 */

/** @file
 * @brief Shepard renormalization factor computation.
 */

#ifndef EXCLUDED_PARTICLE
    /** @brief Excluded particles from the Shepard renormalization factor
     * computation. 
     * 
     * By default all the particles are included. Therefore it is strongly
     * recommended to redefine this macro to specify whether the fluid
     * particles (imove != 1) or the solid particles (imove != 2) are used 
     * @note Redefining this macro this OpenCL script can be recicled
     * @remarks The Shepard renormalization factor is ever computed at the
     * boundary elements and sensors (imove <= 0)
     */
    #define EXCLUDED_PARTICLE(index) imove[index] >= 3
#endif

#if defined(LOCAL_MEM_SIZE) && defined(NO_LOCAL_MEM)
    #error NO_LOCAL_MEM has been set.
#endif

#include "resources/Scripts/types/types.h"
#include "resources/Scripts/KernelFunctions/Kernel.h"

/** @brief Shepard factor computation.
 *
 * \f[ \gamma(\mathbf{x}) = \int_{\Omega}
 *     W(\mathbf{y} - \mathbf{x}) \mathrm{d}\mathbf{y} \f]
 *
 * The shepard renormalization factor is applied for several purposes:
 *   - To interpolate values
 *   - To recover the consistency with the Boundary Integrals formulation
 *   - Debugging
 *
 * In the shepard factor computation the fluid extension particles are not taken
 * into account.
 *
 * @param imove Moving flags.
 *   - imove > 0 for regular fluid particles.
 *   - imove = 0 for sensors.
 *   - imove < 0 for boundary elements/particles.
 * @param r Position \f$ \mathbf{r} \f$.
 * @param rho Density \f$ \rho \f$.
 * @param m Mass \f$ m \f$.
 * @param jhoc Head and tail of chains for each cell.
 * @param shepard Shepard term
 * \f$ \gamma(\mathbf{x}) = \int_{\Omega}
 *     W(\mathbf{y} - \mathbf{x}) \mathrm{d}\mathbf{y} \f$.
 * @param N Number of particles.
 */
__kernel void entry(const __global int* imove,
                    const __global vec* r,
                    const __global float* rho,
                    const __global float* m,
                    const __global svec2* restrict jhoc,
                    __global float* shepard,
                    usize N)
{
    const usize i = get_global_id(0);
    const usize it = get_local_id(0);
    if(i >= N)
        return;
    if((imove[i] < -3) || ((imove[i] > 0) && (EXCLUDED_PARTICLE(i))))
        return;

    const vec_xyz r_i = r[i].XYZ;

    // Initialize the output
    #ifndef LOCAL_MEM_SIZE
        #define _SHEPARD_ shepard[i]
    #else
        #define _SHEPARD_ shepard_l[it]
        __local float shepard_l[LOCAL_MEM_SIZE];
        _SHEPARD_ = 0.f;
    #endif

    FOR_NEIGHS(N, jhoc){
        if(EXCLUDED_PARTICLE(j))
            continue;

        const vec_xyz r_ij = r[j].XYZ - r_i;
        const float q = length(r_ij) / H;
        if(q >= SUPPORT)
            continue;

        {
            _SHEPARD_ += kernelW(q) * CONW * m[j] / rho[j];
        }
    }END_FOR_NEIGHS()

    #ifdef LOCAL_MEM_SIZE
        shepard[i] = _SHEPARD_;
    #endif
}
