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

/** @addtogroup lela
 * @{
 */

/** @file
 * @brief delta-SPH methods, including the correction terms.
 */

#if defined(LOCAL_MEM_SIZE) && defined(NO_LOCAL_MEM)
    #error NO_LOCAL_MEM has been set.
#endif

#ifndef HAVE_3D
    #include "../types/2D.h"
    #include "../KernelFunctions/Wendland2D.hcl"
#else
    #include "../types/3D.h"
    #include "../KernelFunctions/Wendland3D.hcl"
#endif

/** @brief Velocity interpolation.
 *
 * @param imove Moving flags.
 *   - imove > 0 for regular fluid particles.
 *   - imove = 0 for sensors.
 *   - imove < 0 for boundary elements/particles.
 * @param r Position \f$ \mathbf{r} \f$.
 * @param u Velocity \f$ \mathbf{u} \f$.
 * @param rho Density \f$ \rho \f$.
 * @param m Mass \f$ m \f$.
 * @param xSPH_u Interpolated velocity (Shepard renormalization is not applied
 * yet).
 * @param icell Cell where each particle is located.
 * @param ihoc Head of chain for each cell (first particle found).
 * @param N Number of particles.
 * @param n_cells Number of cells in each direction
 */
__kernel void interpolation(const __global int* imove,
                            const __global vec* r,
                            const __global vec* u,
                            const __global float* rho,
                            const __global float* m,
                            __global vec* xSPH_u,
                            const __global uint *icell,
                            const __global uint *ihoc,
                            uint N,
                            uivec4 n_cells)
{
    const uint i = get_global_id(0);
    const uint it = get_local_id(0);
    if(i >= N)
        return;
    if(imove[i] != 2){
        return;
    }

    const vec_xyz r_i = r[i].XYZ;
    const vec_xyz u_i = u[i].XYZ;

    // Initialize the output
    #ifndef LOCAL_MEM_SIZE
        #define _U_ xSPH_u[i].XYZ
    #else
        #define _U_ xSPH_u_l[it]
        __local vec_xyz xSPH_u_l[LOCAL_MEM_SIZE];
    #endif
    _U_ = VEC_ZERO.XYZ;

    BEGIN_LOOP_OVER_NEIGHS(){
        if( (i == j) || (imove[j] != 2)){
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
            const float w_ij = kernelW(q) * CONW * m[j] / rho[j];
            _U_ += (u[j].XYZ - u_i) * w_ij;
        }
    }END_LOOP_OVER_NEIGHS()

    #ifdef LOCAL_MEM_SIZE
        xSPH_u[i] = _U_;
    #endif
}

/** @brief xSPH velocity replacement.
 *
 * @param imove Moving flags.
 *   - imove > 0 for regular fluid particles.
 *   - imove = 0 for sensors.
 *   - imove < 0 for boundary elements/particles.
 * @param xSPH_u Interpolated velocity (Shepard renormalization is not applied
 * yet).
 * @param shepard Shepard term
 * \f$ \gamma(\mathbf{x}) = \int_{\Omega}
 *     W(\mathbf{y} - \mathbf{x}) \mathrm{d}\mathbf{x} \f$.
 * @param u Velocity \f$ \mathbf{u} \f$.
 * @param N Number of particles.
 * @param xSPH_factor Ratio between the instantaneous velocity, and the
 * interpolated one.
 */
__kernel void xSPH(const __global int* imove,
                   const __global vec* xSPH_u,
                   const __global float* shepard,
                   __global vec* u,
                   unsigned int N,
                   float xSPH_factor)
{
    unsigned int i = get_global_id(0);
    if(i >= N)
        return;
    if(imove[i] != 2)
        return;

    u[i] += xSPH_factor / shepard[i] * xSPH_u[i];
}

/*
 * @}
 */