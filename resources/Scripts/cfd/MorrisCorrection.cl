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
 * @brief Computation of \f$ \nabla \mathbf{u} \f$ due to fluid particles
 * interactions.
 */

#if defined(LOCAL_MEM_SIZE) && defined(NO_LOCAL_MEM)
    #error NO_LOCAL_MEM has been set.
#endif

#include "resources/Scripts/types/types.h"
#include "resources/Scripts/KernelFunctions/Kernel.h"

/** @brief Computation of \f$ \nabla \mathbf{u} \f$ due to fluid particles
 * interactions.
 *
 * Compute the gradienht of the velocity matrix, taking into account just the
 * fluid-fluid interactions.
 *
 * @param imove Moving flags.
 *   - imove > 0 for regular fluid particles.
 *   - imove = 0 for sensors.
 *   - imove < 0 for boundary elements/particles.
 * @param r Position \f$ \mathbf{r} \f$.
 * @param rho Density \f$ \rho \f$.
 * @param m Mass \f$ m \f$.
 * @param morris_o1 Morris order 1 correction term.
 * @param icell Cell where each particle is located.
 * @param ihoc Head of chain for each cell (first particle found).
 * @param N Number of particles.
 * @param n_cells Number of cells in each direction
 */
__kernel void order1(const __global int* imove,
                     const __global vec* r,
                     const __global float* rho,
                    const __global float* m,
                     __global vec* morris_o1,
                     // Link-list data
                     const __global uint *icell,
                     const __global uint *ihoc,
                     // Simulation data
                     uint N,
                     uivec4 n_cells)
{
    const uint i = get_global_id(0);
    const uint it = get_local_id(0);
    if(i >= N)
        return;
    if(imove[i] != 1){
        return;
    }

    const vec_xyz r_i = r[i].XYZ;

    // Initialize the output
    #ifndef LOCAL_MEM_SIZE
        #define _MORRISO1_ morris_o1[i]
    #else
        #define _MORRISO1_ morris_o1_l[it]
        __local vec_xyz morris_o1_l[LOCAL_MEM_SIZE];
    #endif
    _MORRISO1_ = VEC_ZERO.XYZ;

    BEGIN_LOOP_OVER_NEIGHS(){
        if(i == j){
            j++;
            continue;
        }
        if(imove[j] != 1){
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

            _MORRISO1_ += 2.f * f_ij * r_ij;
        }
    }END_LOOP_OVER_NEIGHS()

    #ifdef LOCAL_MEM_SIZE
        morris_o1[i].XYZ = _MORRISO1_;
    #endif
}

/** @brief Morris velocity Laplacian term correction.
 *
 * @param imove Moving flags.
 *   - imove > 0 for regular fluid particles.
 *   - imove = 0 for sensors.
 *   - imove < 0 for boundary elements/particles.
 * @param rho Density \f$ \rho_{n+1} \f$.
 * @param grad_u Velocity gradient \f$ \nabla \mathbf{u} \f$.
 * @param morris_o1 Morris order 1 correction term.
 * @param lap_u Velocity laplacian \f$ \frac{\Delta \mathbf{u}}{rho} \f$.
 * @param N Number of particles.
 */
__kernel void entry(const __global int* imove,
                    const __global float* rho,
                    const __global matrix* grad_u,
                    const __global vec* morris_o1,
                    __global vec* lap_u,
                    unsigned int N)
{
    unsigned int i = get_global_id(0);
    if(i >= N)
        return;
    if(imove[i] != 1)
        return;

    lap_u[i].XYZ -= MATRIX_DOT(grad_u[i], morris_o1[i]).XYZ / rho[i];
}
