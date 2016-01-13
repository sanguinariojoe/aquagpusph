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

#ifndef HAVE_3D
    #include "../../../types/2D.h"
    #include "../../../KernelFunctions/Wendland2D.hcl"
#else
    #include "../../../types/3D.h"
    #include "../../../KernelFunctions/Wendland3D.hcl"
#endif

/** @brief Velocity interpolation at the boundary elements.
 *
 * The values are computed using just the solid information. The resulting
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
 * @param u Velocity \f$ \mathbf{u} \f$.
 * @param icell Cell where each particle is located.
 * @param ihoc Head of chain for each cell (first particle found).
 * @param N Number of particles.
 * @param n_cells Number of cells in each direction
 * @param BImotion_iset Set of particles affected
 */
__kernel void interpolation(const __global uint* iset,
                            const __global int* imove,
                            const __global vec* r,
                            const __global float* m,
                            const __global float* rho,
                            __global vec* u,
                            // Link-list data
                            const __global uint *icell,
                            const __global uint *ihoc,
                            // Simulation data
                            uint N,
                            uivec4 n_cells,
                            uint BImotion_iset)
{
    const uint i = get_global_id(0);
    const uint it = get_local_id(0);
    if(i >= N)
        return;
    if((imove[i] != -3) || (iset[i] != BImotion_iset)){
        return;
    }

    const vec_xyz r_i = r[i].XYZ;

    // Initialize the output
    #ifndef LOCAL_MEM_SIZE
        #define _U_ u[i].XYZ
    #else
        #define _U_ u_l[it]
        __local vec_xyz u_l[LOCAL_MEM_SIZE];
    #endif
    _U_ = VEC_ZERO.XYZ;

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
            _U_ += u[j].XYZ * w_ij;
        }
    }END_LOOP_OVER_NEIGHS()

    #ifdef LOCAL_MEM_SIZE
        u[i].XYZ = _U_;
    #endif
}

/** @brief Velocity renormalization.
 *
 * @param iset Set of particles index.
 * @param imove Moving flags.
 *   - imove = 2 for regular solid particles.
 *   - imove = 0 for sensors.
 *   - imove < 0 for boundary elements/particles.
 * @param shepard Shepard term
 * \f$ \gamma(\mathbf{x}) = \int_{\Omega}
 *     W(\mathbf{y} - \mathbf{x}) \mathrm{d}\mathbf{x} \f$.
 * @param u Velocity \f$ \mathbf{u} \f$.
 * @param N Total number of particles and boundary elements.
 * @param BImotion_iset Set of particles affected
 */
__kernel void shepard(const __global uint* iset,
                      const __global int* imove,
                      const __global float* shepard,
                      __global vec* u,
                      uint N,
                      uint BImotion_iset)
{
    uint i = get_global_id(0);
    if(i >= N)
        return;
    if((imove[i] != -3) || (iset[i] != BImotion_iset)){
        return;
    }
    const float shepard_i = shepard[i];
    if(shepard_i < 1.0E-6f){
        // It will be considered that there are not enough
        // particles to interpolate
        u[i] = VEC_ZERO;
        return;
    }

    u[i] /= shepard_i;
}

/** @brief Simple Euler time integration of the velocity
 *
 * Since the velocity is resulting from the interpolation of the solid particle,
 * it does not make any sense to integrate it with an improved Euler scheme.
 * 
 *   \f$ \mathbf{u}_{n+1} = \mathbf{u}_{n} + \Delta t
 *      \left. \frac{\mathrm{d} \mathbf{u}}{\mathrm{d}t} \right\vert_{n+1/2}
 *   \right) \f$
 *
 * @param iset Set of particles index.
 * @param imove Moving flags.
 *   - imove = 2 for regular solid particles.
 *   - imove = 0 for sensors.
 *   - imove < 0 for boundary elements/particles.
 * @param r Position \f$ \mathbf{r} \f$.
 * @param u Velocity \f$ \mathbf{u} \f$.
 * @param N Number of particles.
 * @param dt Time step \f$ \Delta t \f$.
 * @param BImotion_iset Set of particles affected
 */
__kernel void euler(const __global uint* iset,
                    const __global int* imove,
                    __global vec* r,
                    const __global vec* u,
                    unsigned int N,
                    float dt,
                    uint BImotion_iset)
{
    unsigned int i = get_global_id(0);
    if(i >= N)
        return;
    if((imove[i] != -3) || (iset[i] != BImotion_iset)){
        return;
    }

    r[i] += dt * u[i];
}

/*
 * @}
 */