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
 * @brief Boundary Integrals for Particles packing algorithm
 */

#if defined(LOCAL_MEM_SIZE) && defined(NO_LOCAL_MEM)
    #error NO_LOCAL_MEM has been set.
#endif

#include "resources/Scripts/types/types.h"
#include "resources/Scripts/KernelFunctions/Kernel.h"

/** @brief Performs the boundary effect on the fluid particles.
 * @param imove Moving flags.
 *   - imove > 0 for regular fluid particles.
 *   - imove = 0 for sensors.
 *   - imove < 0 for boundary elements/particles.
 * @param r Position \f$ \mathbf{r} \f$.
 * @param normal Normal \f$ \mathbf{n} \f$.
 * @param m Area of the boundary element \f$ s \f$.
 * @param grad_shepard Shepard gradient \f$ \nabla \gamma \f$.
 * @param icell Cell where each particle is located.
 * @param ihoc Head of chain for each cell (first particle found).
 * @param N Number of particles.
 * @param n_cells Number of cells in each direction
 */
__kernel void interactions(const __global int* imove,
                           const __global vec* r,
                           const __global vec* normal,
                           const __global float* m,
                           __global vec* grad_shepard,
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
    if(imove[i] != 1)
        return;

    const vec_xyz r_i = r[i].XYZ;

    // Initialize the output
    #ifndef LOCAL_MEM_SIZE
        #define _GRADS_ grad_shepard[i].XYZ
    #else
        #define _GRADS_ grad_shepard_l[it]
        __local vec_xyz grad_shepard_l[LOCAL_MEM_SIZE];
        _GRADS_ = grad_shepard[i].XYZ;
    #endif

    BEGIN_LOOP_OVER_NEIGHS(){
        if(imove[j] != -3){
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
            _GRADS_ += kernelW(q) * CONW * m[j] * normal[j].XYZ;
        }
    }END_LOOP_OVER_NEIGHS()

    #ifdef LOCAL_MEM_SIZE
        grad_shepard[i].XYZ = _GRADS_;
    #endif
}

/** @brief Renormalize the differential operators.
 *
 * The main drawback of the boundary integrals formulation is the requirement
 * of the renormalization of the computed differentiqal operators, which is
 * destroying several conservation properties.
 *
 * @param imove Moving flags.
 *   - imove > 0 for regular fluid particles.
 *   - imove = 0 for sensors.
 *   - imove < 0 for boundary elements/particles.
 * @param shepard Shepard term
 * \f$ \gamma(\mathbf{x}) = \int_{\Omega}
 *     W(\mathbf{y} - \mathbf{x}) \mathrm{d}\mathbf{x} \f$.
 * @param grad_shepard Shepard gradient \f$ \nabla \gamma \f$.
 * @param N Total number of particles and boundary elements.
 */
__kernel void apply(const __global int* imove,
                    const __global float* shepard,
                    __global vec* grad_shepard,
                    uint N)
{
    uint i = get_global_id(0);
    if(i >= N)
        return;
    if(imove[i] != 1){
        return;
    }

    grad_shepard[i] /= shepard[i];
}