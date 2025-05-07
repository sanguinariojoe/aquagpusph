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

#if defined(LOCAL_MEM_SIZE) && defined(NO_LOCAL_MEM)
    #error NO_LOCAL_MEM has been set.
#endif

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
 * @param h_var variable kernel lenght \f$ h \f$.
 * @param shepard Shepard term
 * \f$ \gamma(\mathbf{x}) = \int_{\Omega}
 *     W(\mathbf{y} - \mathbf{x}) \mathrm{d}\mathbf{x} \f$.
 * @param icell Cell where each particle is located.
 * @param ihoc Head of chain for each cell (first particle found).
 * @param N Number of particles.
 * @param n_cells Number of cells in each direction
 */
__kernel void entry(const __global int* imove,
                    const __global int* imirrored,
                    const __global vec* r,
                    const __global float* rho,
                    const __global float* m,
                    const __global float* h_var,
                    __global float* shepard,
                    usize N,
                    LINKLIST_LOCAL_PARAMS)
{
    const usize i = get_global_id(0);
    const usize it = get_local_id(0);
    if(i >= N)
        return;
    if((imove[i] < -3) || (imove[i] > 1) || (!imirrored[i]))
        return;

    const vec_xyz r_i = r[i].XYZ;
    const float h_i = h_var[i];
    #ifndef HAVE_3D
        const float conw = 1.f / (h_i * h_i);
    #else
        const float conw = 1.f / (h_i * h_i * h_i);
    #endif

    // Initialize the output
    #ifndef LOCAL_MEM_SIZE
        #define _SHEPARD_ shepard[i]
    #else
        #define _SHEPARD_ shepard_l[it]
        __local float shepard_l[LOCAL_MEM_SIZE];
        _SHEPARD_ = shepard[i];
    #endif

    const usize c_i = icell[i];
    BEGIN_NEIGHS(c_i, N, n_cells, icell, ihoc){
        if((imove[j] != 1) || (imirrored[j])){
            j++;
            continue;
        }

        const vec_xyz r_ij = r[j].XYZ - r_i;
        const float q = length(r_ij) / h_i;
        if(q >= SUPPORT)
        {
            j++;
            continue;
        }

        {
            _SHEPARD_ += conw * kernelW(q) * m[j] / rho[j];
        }
    }END_NEIGHS()

    #ifdef LOCAL_MEM_SIZE
        shepard[i] = _SHEPARD_;
    #endif
}

/*
 * @}
 */
