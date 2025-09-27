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

#include "resources/Scripts/types/types.h"
#include "resources/Scripts/KernelFunctions/Kernel.h"

/** @brief Shepard factor computation, due to the particles at the other portal
 * side.
 *
 * \f[ \gamma(\mathbf{x}) = \int_{\Omega}
 *     W(\mathbf{y} - \mathbf{x}) \mathrm{d}\mathbf{x} \f]
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
 * @param imirrored 0 if the particle has not been mirrored, 1 otherwise.
 * @param r Position \f$ \mathbf{r} \f$.
 * @param rho Density \f$ \rho \f$.
 * @param m Mass \f$ m \f$.
 * @param jhoc Head and tail of chains for each cell.
 * @param shepard Shepard term
 * \f$ \gamma(\mathbf{x}) = \int_{\Omega}
 *     W(\mathbf{y} - \mathbf{x}) \mathrm{d}\mathbf{x} \f$.
 * @param N Number of particles.
 */
__kernel void entry(const __global int* restrict imove,
                    const __global int* restrict imirrored,
                    const __global vec* restrict r,
                    const __global float* restrict rho,
                    const __global float* restrict m,
                    const __global svec2* restrict jhoc,
                    __global float* restrict shepard,
                    usize N)
{
    const usize i = get_global_id(0);
    const usize it = get_local_id(0);
    if(i >= N)
        return;
    if((imove[i] < -3) || (imove[i] > 1) || (!imirrored[i]))
        return;

    const vec_xyz r_i = r[i].XYZ;

    __private float __shepard = 0.f;

    FOR_NEIGHS(N, jhoc){
        if((imove[j] != 1) || (imirrored[j]))
            continue;

        const vec_xyz r_ij = r[j].XYZ - r_i;
        const float q = length(r_ij) / H;
        if(q >= SUPPORT)
            continue;

        {
            _SHEPARD_ += kernelW(q) * CONW * m[j] / rho[j];
        }
    }END_FOR_NEIGHS()

    shepard[i] += __shepard;
}
