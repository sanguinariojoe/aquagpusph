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
 * @brief Fixed boundary elements methods.
 */

#if defined(LOCAL_MEM_SIZE) && defined(NO_LOCAL_MEM)
    #error NO_LOCAL_MEM has been set.
#endif

#include "resources/Scripts/types/types.h"
#include "resources/Scripts/KernelFunctions/Kernel.h"

/** @brief Pressure and stress deviation interpolation at the boundary elements.
 *
 * The values are computed using just the fluid information. The resulting
 * interpolated values are not renormalized yet.
 *
 * Just the elements with the flag imove = -3 are considered boundary elements.
 *
 * @param iset Set of particles index.
 * @param imove Moving flags.
 *   - imove > 0 for regular fluid particles.
 *   - imove = 0 for sensors.
 *   - imove < 0 for boundary elements/particles.
 * @param r Position \f$ \mathbf{r} \f$.
 * @param m Mass \f$ m \f$.
 * @param rho Density \f$ \rho \f$.
 * @param p Pressure \f$ p \f$.
 * @param S Deviatory stress \f$ S \f$.
 * @param icell Cell where each particle is located.
 * @param ihoc Head of chain for each cell (first particle found).
 * @param refd Density of reference of the fluid \f$ \rho_0 \f$.
 * @param N Number of particles.
 * @param n_cells Number of cells in each direction
 * @param g Gravity acceleration \f$ \mathbf{g} \f$.
 * @param BIstress_iset Set of particles affected
 */
__kernel void interpolation(const __global uint* iset,
                            const __global int* imove,
                            const __global vec* r,
                            const __global float* m,
                            const __global float* rho,
                            __global float* p,
                            __global matrix* S,
                            // Link-list data
                            const __global uint *icell,
                            const __global uint *ihoc,
                            // Simulation data
                            __constant float* refd,
                            uint N,
                            uivec4 n_cells,
                            vec g,
                            uint BIstress_iset)
{
    const uint i = get_global_id(0);
    const uint it = get_local_id(0);
    if(i >= N)
        return;
    if((imove[i] != -3) || (iset[i] != BIstress_iset)){
        return;
    }

    const vec_xyz r_i = r[i].XYZ;
    const float rdenf = refd[iset[i]];

    // Initialize the output
    #ifndef LOCAL_MEM_SIZE
        #define _P_ p[i]
        #define _S_ S[i]
    #else
        #define _P_ p_l[it]
        __local float p_l[LOCAL_MEM_SIZE];
        #define _S_ S_l[it]
        __local matrix S_l[LOCAL_MEM_SIZE];
    #endif
    _P_ = 0.f;
    _S_ = MAT_ZERO;

    BEGIN_LOOP_OVER_NEIGHS(){
        if(imove[j] != 2){
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
            _P_ += (p[j] - rdenf * dot(g.XYZ, r_ij)) * w_ij;
            _S_ += S[j] * w_ij;
        }
    }END_LOOP_OVER_NEIGHS()

    #ifdef LOCAL_MEM_SIZE
        p[i] = _P_;
        S[i] = _S_;
    #endif
}

/** @brief Pressure and stress deviation renormalization.
 *
 * @param iset Set of particles index.
 * @param imove Moving flags.
 *   - imove = 2 for regular solid particles.
 *   - imove = 0 for sensors.
 *   - imove < 0 for boundary elements/particles.
 * @param shepard Shepard term
 * \f$ \gamma(\mathbf{x}) = \int_{\Omega}
 *     W(\mathbf{y} - \mathbf{x}) \mathrm{d}\mathbf{x} \f$.
 * @param p Pressure \f$ p \f$.
 * @param S Deviatory stress \f$ S \f$.
 * @param N Total number of particles and boundary elements.
 * @param BIstress_iset Set of particles affected
 */
__kernel void shepard(const __global uint* iset,
                      const __global int* imove,
                      const __global float* shepard,
                      __global float* p,
                      __global matrix* S,
                      uint N,
                      uint BIstress_iset)
{
    uint i = get_global_id(0);
    if(i >= N)
        return;
    if((imove[i] != -3) || (iset[i] != BIstress_iset)){
        return;
    }
    float shepard_i = shepard[i];
    if(shepard_i < 1.0E-6f){
        // It will be considered that there are not enough
        // particles to interpolate
        shepard_i = 1.f;
    }

    p[i] /= shepard_i;
    S[i] /= shepard_i;
}

/*
 * @}
 */