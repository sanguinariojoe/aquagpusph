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
 * @brief delta-SPH methods, including the correction terms
 */

#ifndef EXCLUDED_PARTICLE
    /** @brief Condition to exclude a particle from the delta-SPH model
     * 
     * By default all the boundary elements are excluded. Even though it is
     * enough for simulation where fluid and solid mechanics are not combined,
     * it is strongly recommended to conveniently overload this macro. 
     * @note Redefining this macro this OpenCL script can be recicled
     */
    #define EXCLUDED_PARTICLE(index) imove[index] <= 0
#endif

#if defined(LOCAL_MEM_SIZE) && defined(NO_LOCAL_MEM)
    #error NO_LOCAL_MEM has been set.
#endif

#ifndef HAVE_3D
    #include "../../types/2D.h"
    #include "../../KernelFunctions/Wendland2D.hcl"
#else
    #include "../../types/3D.h"
    #include "../../KernelFunctions/Wendland3D.hcl"
#endif

/** @brief Laplacian of the pressure computation.
 *
 * @param imove Moving flags.
 *   - imove > 0 for regular fluid particles.
 *   - imove = 0 for sensors.
 *   - imove < 0 for boundary elements/particles.
 * @param r Position \f$ \mathbf{r} \f$.
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
                   const __global vec* r,
                   const __global float* rho,
                   const __global float* m,
                   const __global float* p,
                   const __global float* h_var,
                   __global float* lap_p,
                   const __global uint *icell,
                   const __global uint *ihoc,
                   uint N,
                   uivec4 n_cells)
{
    const uint i = get_global_id(0);
    const uint it = get_local_id(0);
    if(i >= N)
        return;
    if(EXCLUDED_PARTICLE(i)){
        return;
    }

    const vec_xyz r_i = r[i].XYZ;
    const float h_i = h_var[i];
    const float p_i = p[i];
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
        _LAPP_ = 0.f;
    #endif

    BEGIN_LOOP_OVER_NEIGHS(){
        if( (i == j) || (EXCLUDED_PARTICLE(j))){
            j++;
            continue;
        }
        const vec_xyz r_ij = r[j].XYZ - r_i;
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
    }END_LOOP_OVER_NEIGHS()

    #ifdef LOCAL_MEM_SIZE
        lap_p[i] = _LAPP_;
    #endif
}

/*
 * @}
 */