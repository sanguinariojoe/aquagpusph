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
 * @brief MLS kernel transformation matrix computation
 */

#if defined(LOCAL_MEM_SIZE) && defined(NO_LOCAL_MEM)
    #error NO_LOCAL_MEM has been set.
#endif

#include "resources/Scripts/types/types.h"
#include "resources/Scripts/KernelFunctions/Kernel.h"

/** @brief Compute the MLS transformation matrix inverse, \f$ L_i^{-1} \f$.
 * 
 * Such transformation matrix can be multiplied by the kernel gradient to
 * produce a new kernel gradient,
 * \f$ \nabla W^{L}_{ij} = L_i \cdot \nabla W_{ij} \f$, such that the lienar
 * fields differential operators are consistently computed.
 *
 * @param imove Moving flags.
 *   - imove > 0 for regular fluid/solid particles.
 *   - imove = 0 for sensors.
 *   - imove < 0 for boundary elements/particles.
 * @param r Position \f$ \mathbf{r}_{n+1} \f$.
 * @param rho Density \f$ \rho \f$.
 * @param m Mass \f$ m \f$.
 * @param jhoc Head of chain for each cell (first particle found).
 * @param mls Kernel MLS transformation matrix \f$ L \f$.
 * @param N Number of particles.
 * @param mls_imove Type of particles affected
 * @note The MLS kernel transformation will be computed just for the particles
 * with the moving flag mls_imove, and using just the information of the
 * particles with the moving flag mls_imove
 */
__kernel void entry(const __global int* imove,
                    const __global vec* r,
                    const __global float* rho,
                    const __global float* m,
                    const __global svec2* jhoc,
                    __global matrix* mls,
                    usize N,
                    uint mls_imove)
{
    const usize i = get_global_id(0);
    const usize it = get_local_id(0);
    if(i >= N)
        return;
    if(imove[i] != mls_imove){
        return;
    }

    const vec_xyz r_i = r[i].XYZ;

    // Initialize the output
    #ifndef LOCAL_MEM_SIZE
        #define _MLS_ mls[i]
    #else
        #define _MLS_ mls_l[it]
        __local matrix mls_l[LOCAL_MEM_SIZE];
    #endif
    _MLS_ = MAT_ZERO;

    FOR_NEIGHS(N, jhoc){
        if(i == j){
            continue;
        }
        if(imove[j] != mls_imove){
            continue;
        }
        const vec_xyz r_ij = r[j].XYZ - r_i;
        const float q = length(r_ij) / H;
        if(q >= SUPPORT)
        {
            continue;
        }
        {
            const float f_ij = kernelF(q) * CONF * m[j] / rho[j];
            _MLS_ += outer(r_ij, f_ij * r_ij);
        }
    }END_FOR_NEIGHS()

    #ifdef LOCAL_MEM_SIZE
        mls[i] = _MLS_;
    #endif
}

/** @brief Invert the matrix computed in entry() to get the final MLS
 * transformation matrix, \f$ L_i \f$.
 * 
 * Such transformation matrix can be multiplied by the kernel gradient to
 * produce a new kernel gradient,
 * \f$ \nabla W^{L}_{ij} = L_i \cdot \nabla W_{ij} \f$, such that the lienar
 * fields differential operators are consistently computed.
 *
 * @param imove Moving flags.
 *   - imove > 0 for regular fluid/solid particles.
 *   - imove = 0 for sensors.
 *   - imove < 0 for boundary elements/particles.
 * @param mls Kernel MLS transformation matrix \f$ L \f$.
 * @param N Number of particles.
 * @param mls_imove Type of particles affected
 */
__kernel void mls_inv(const __global int* imove,
                      __global matrix* mls,
                      usize N,
                      uint mls_imove)
{
    const usize i = get_global_id(0);
    if(i >= N)
        return;
    if(imove[i] != mls_imove){
        return;
    }

    mls[i] = MATRIX_INV(mls[i]);
}

/*
 * @}
 */
