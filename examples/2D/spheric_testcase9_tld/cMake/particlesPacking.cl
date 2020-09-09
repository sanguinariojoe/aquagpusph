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
 * @brief Particles packing algorithm
 */

#if defined(LOCAL_MEM_SIZE) && defined(NO_LOCAL_MEM)
    #error NO_LOCAL_MEM has been set.
#endif

#include "resources/Scripts/types/types.h"
#include "resources/Scripts/KernelFunctions/Kernel.h"

/** @brief Fluid particles interactions computation.
 *
 * Compute the differential operators involved in the numerical scheme, taking
 * into account just the fluid-fluid and fluid-ghost interactions.
 *
 * @param imove Moving flags.
 *   - imove > 0 for regular fluid particles.
 *   - imove = 0 for sensors.
 *   - imove < 0 for boundary elements/particles.
 * @param r Position \f$ \mathbf{r} \f$.
 * @param rho Density \f$ \rho \f$.
 * @param m Mass \f$ m \f$.
 * @param grad_shepard Shepard gradient \f$ \nabla \gamma \f$.
 * @param icell Cell where each particle is located.
 * @param ihoc Head of chain for each cell (first particle found).
 * @param N Number of particles.
 * @param n_cells Number of cells in each direction
 */
__kernel void interactions(const __global int* imove,
                           const __global vec* r,
                           const __global float* rho,
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
        _GRADS_ = VEC_ZERO.XYZ;
    #endif

    BEGIN_LOOP_OVER_NEIGHS(){
        if(i == j){
            j++;
            continue;
        }
        if((imove[j] != 1) && (imove[j] != -1)){
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
            const float rho_j = rho[j];
            const float f_ij = kernelF(q) * CONF * m[j];

            _GRADS_ += kernelF(q) * CONF * m[j] / rho[j] * r_ij;
        }
    }END_LOOP_OVER_NEIGHS()

    #ifdef LOCAL_MEM_SIZE
        grad_shepard[i].XYZ = _GRADS_;
    #endif
}

/** @brief Carry out the particles shifting.
 *
 * The particles are slightly moved to minimize the \f$ \nabla \gamma \f$ term.
 *
 * @param imove Moving flags.
 *   - imove > 0 for regular fluid particles.
 *   - imove = 0 for sensors.
 *   - imove < 0 for boundary elements/particles.
 * @param rho Density \f$ \rho \f$.
 * @param m Mass \f$ m \f$.
 * @param grad_shepard Shepard gradient \f$ \nabla \gamma \f$.
 * @param r Position \f$ \mathbf{r} \f$.
 * @param relax_factor Relaxation factor to grant the convergence \f$ f \f$.
 * @param N Number of particles.
 */
__kernel void shifting(const __global int* imove,
                       const __global float* rho,
                       const __global float* m,
                       const __global vec* grad_shepard,
                       __global vec* r,
                       unsigned int N,
                       float relax_factor)
{
    const uint i = get_global_id(0);
    if(i >= N)
        return;
    if(imove[i] != 1)
        return;

    const float dr = pow(m[i] / rho[i], 1.f / (float)(DIMS));
    r[i].XYZ -= relax_factor * dr * H * grad_shepard[i].XYZ;
}

/** @brief Setup an hydrostatic pressure field
 *
 * To this end, a Stiffness Equation Of State (EOS) that relates the pressure
 * and the density fields is considered,
 * \f$ p = c_s^2 \left(\rho - \rho_0 \right) \f$
 *
 * @param iset Set of particles index.
 * @param imove Moving flags.
 *   - imove > 0 for regular fluid particles.
 *   - imove = 0 for sensors.
 *   - imove < 0 for boundary elements/particles.
 * @param r Position \f$ \mathbf{r} \f$.
 * @param rho Density \f$ \rho_{n+1/2} \f$.
 * @param m Mass \f$ m \f$.
 * @param refd Density of reference of the fluid \f$ \rho_0 \f$.
 * @param N Number of particles.
 * @param cs Speed of sound \f$ c_s \f$.
 */
__kernel void hydrostatic(const __global unsigned int* iset,
                          const __global int* imove,
                          const __global vec* r,
                          __global float* rho,
                          __global float* m,
                          __constant float* refd,
                          unsigned int N,
                          vec r0,
                          vec g,
                          float cs)
{
    unsigned int i = get_global_id(0);
    if(i >= N)
        return;
    if(imove[i] < 1)
        return;

    const float v = m[i] / rho[i];
    const float p = refd[iset[i]] * dot(g.XYZ, r[i].XYZ - r0.XYZ);
    rho[i] = refd[iset[i]] + p / (cs * cs);
    m[i] = v * rho[i];
}
