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

#ifndef HAVE_3D
    #include "../../../types/2D.h"
    #include "../../../KernelFunctions/Wendland2D.hcl"
#else
    #include "../../../types/3D.h"
    #include "../../../KernelFunctions/Wendland3D.hcl"
#endif

/** @brief Compute the MLS transformation matrix inverse, \f$ L_i^{-1} \f$, due
 * to the particles at the other portal side.
 * 
 * Such transformation matrix can be multiplied by the kernel gradient to
 * produce a new kernel gradient,
 * \f$ \nabla W^{L}_{ij} = L_i \cdot \nabla W_{ij} \f$, such that the lienar
 * fields differential operators are consistently computed.
 *
 * @param imove Moving flags.
 *   - imove > 0 for regular fluid particles.
 *   - imove = 0 for sensors.
 *   - imove < 0 for boundary elements/particles.
 * @param imirrored 0 if the particle has not been mirrored, 1 otherwise.
 * @param r Position \f$ \mathbf{r} \f$.
 * @param rho Density \f$ \rho \f$.
 * @param m Mass \f$ m \f$.
 * @param mls Kernel MLS transformation matrix \f$ L \f$.
 * @param icell Cell where each particle is located.
 * @param ihoc Head of chain for each cell (first particle found).
 * @param N Number of particles.
 * @param n_cells Number of cells in each direction
 * @param mls_imove Type of particles affected
 * @note The MLS kernel transformation will be computed just for the particles
 * with the moving flag mls_imove, and using just the information of the
 * particles with the moving flag mls_imove
 */
__kernel void entry(const __global int* imove,
                    const __global int* imirrored,
                    const __global vec* r,
                    const __global float* rho,
                    const __global float* m,
                    __global matrix* mls,
                    // Link-list data
                    const __global uint *icell,
                    const __global uint *ihoc,
                    // Simulation data
                    uint N,
                    uivec4 n_cells,
                    uint mls_imove)
{
    const uint i = get_global_id(0);
    const uint it = get_local_id(0);
    if(i >= N)
        return;
    if((imove[i] != 1) || (!imirrored[i]))
        return;

    const vec_xyz r_i = r[i].XYZ;

    // Initialize the output
    #ifndef LOCAL_MEM_SIZE
        #define _MLS_ mls[i]
    #else
        #define _MLS_ mls_l[it]
        __local matrix mls_l[LOCAL_MEM_SIZE];
        _MLS_ = mls[i];
    #endif

    BEGIN_LOOP_OVER_NEIGHS(){
        if((imove[j] != mls_imove) || (imirrored[j])){
            j++;
            continue;
        }

        const vec_xyz r_ij = r[j].XYZ - r_i;
        const float q = length(r_ij) / H;
        if(q >= SUPPORT)
        {
            j++;
            continue;
        }

        {
            const float f_ij = kernelF(q) * CONF * m[j] / rho[j];
            _MLS_ += outer(r_ij, f_ij * r_ij);
        }
    }END_LOOP_OVER_NEIGHS()

    #ifdef LOCAL_MEM_SIZE
        mls[i] = _MLS_;
    #endif
}
