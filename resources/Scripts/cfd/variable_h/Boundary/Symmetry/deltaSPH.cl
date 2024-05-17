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
 * @brief Particles interactions computation.
 */

#if defined(LOCAL_MEM_SIZE) && defined(NO_LOCAL_MEM)
    #error NO_LOCAL_MEM has been set.
#endif

#include "resources/Scripts/types/types.h"
#include "resources/Scripts/KernelFunctions/Kernel.h"

#if __LAP_FORMULATION__ == __LAP_MONAGHAN__
    #ifndef HAVE_3D
        #define __CLEARY__ 8.f
    #else
        #define __CLEARY__ 10.f
    #endif
#endif

/** @brief Laplacian of the pressure computation.
 *
 * @param imove Moving flags.
 *   - imove > 0 for regular fluid particles.
 *   - imove = 0 for sensors.
 *   - imove < 0 for boundary elements/particles.
 * @param imirrored 0 if the particle has not been mirrored, 1 otherwise.
 * @param r Position \f$ \mathbf{r} \f$.
 * @param rmirrored Mirrored position of the particle, \a r if \a imirrored is
 * false (0).
 * @param rho Density \f$ \rho \f$.
 * @param m Mass \f$ m \f$.
 * @param p Pressure \f$ p \f$.
 * @param h_var variable kernel lenght \f$ h \f$.
 * @param lap_p Pressure laplacian \f$ \Delta p \f$.
 * @param icell Cell where each particle is located.
 * @param ihoc Head of chain for each cell (first particle found).
 * @param N Number of particles.
 * @param n_cells Number of cells in each direction
 */
__kernel void lapp(const __global int* imove,
                   const __global int* imirrored,
                   const __global vec* r,
                   const __global vec* rmirrored,
                   const __global float* rho,
                   const __global float* m,
                   const __global float* p,
                   const __global float* h_var,
                   __global float* lap_p,
                   // Link-list data
                   __global uint *icell,
                   __global uint *ihoc,
                   // Simulation data
                   usize N,
                   LINKLIST_LOCAL_PARAMS)
{
    const usize i = get_global_id(0);
    const usize it = get_local_id(0);
    if(i >= N)
        return;
    if((!imirrored[i]) || (imove[i] != 1))
        return;

    const vec_xyz r_i = r[i].XYZ;
    const float p_i = p[i];
    const float h_i = h_var[i];
    #ifndef HAVE_3D
        const float conf_i = 1.f / (h_i * h_i * h_i * h_i);
    #else
        const float conf_i = 1.f / (h_i * h_i * h_i * h_i * h_i);
    #endif

    // Initialize the output
    #ifndef LOCAL_MEM_SIZE
        #define _LAPP_ lap_p[i]
    #else
        #define _LAPP_ lap_p_l[it]
        __local float lap_p_l[LOCAL_MEM_SIZE];
        _LAPP_ = lap_p[i];
    #endif

    const usize c_i = icell[i];
    BEGIN_NEIGHS(c_i, N, n_cells, icell, ihoc){
        if(!imirrored[j] || (imove[j] != 1)){
            j++;
            continue;
        }
        const vec_xyz r_ij = rmirrored[j].XYZ - r_i;
        const float h_j = h_var[j];
        const float l_ij = length(r_ij);
        const float q_i = min(l_ij / h_i, SUPPORT);
        const float q_j = min(l_ij / h_j, SUPPORT);
        if((q_i == SUPPORT) && (q_j == SUPPORT))
        {
            j++;
            continue;
        }

        {
            #ifndef HAVE_3D
                const float conf_j = 1.f / (h_j * h_j * h_j * h_j);
            #else
                const float conf_j = 1.f / (h_j * h_j * h_j * h_j * h_j);
            #endif
            const float fi_ij = conf_i * kernelF(q_i);
            const float fj_ij = conf_j * kernelF(q_j);
            _LAPP_ += (p[j] - p_i) * 0.5f * (fi_ij + fj_ij) * m[j] / rho[j];
        }
    }END_NEIGHS()

    #ifdef LOCAL_MEM_SIZE
        lap_p[i] = _LAPP_;
    #endif
}
