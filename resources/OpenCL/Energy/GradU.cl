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
 * @brief Compute all the spatial derivatives of the velocity.
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

/** @brief Compute all the spatial derivatives of the velocity.
 *
 * @param imove Moving flags.
 *   - imove > 0 for regular fluid particles.
 *   - imove = 0 for sensors.
 *   - imove < 0 for boundary elements/particles.
 * @param r Position \f$ \mathbf{r} \f$.
 * @param u Velocity \f$ \mathbf{u} \f$.
 * @param rho Density \f$ \rho \f$.
 * @param m Mass \f$ m \f$.
 * @param grad_ux Gradient of the first component of the velocity:
 *   \f$ \nabla \left(\mathbf{u} \cdot \mathbf{e_1}\right) \f$
 * @param grad_uy Gradient of the second component of the velocity:
 *   \f$ \nabla \left(\mathbf{u} \cdot \mathbf{e_2}\right) \f$
 * @param grad_uz Gradient of the third component of the velocity:
 *   \f$ \nabla \left(\mathbf{u} \cdot \mathbf{e_3}\right) \f$
 * @param icell Cell where each particle is located.
 * @param ihoc Head of chain for each cell (first particle found).
 * @param N Number of particles.
 * @param n_cells Number of cells in each direction
 * @param g Gravity acceleration \f$ \mathbf{g} \f$.
 */
__kernel void main(const __global int* imove,
                   const __global vec* r,
                   const __global vec* u,
                   const __global float* rho,
                   const __global float* m,
                   __global vec4* grad_ux,
                   __global vec4* grad_uy,
                   __global vec4* grad_uz,
                   // Link-list data
                   const __global uint *icell,
                   const __global uint *ihoc,
                   // Simulation data
                   uint N,
                   uivec4 n_cells,
                   vec g)
{
    const uint i = get_global_id(0);
    const uint it = get_local_id(0);
    if(i >= N)
        return;

    grad_ux[i] = (vec4)(0.f, 0.f, 0.f, 0.f);
    grad_uy[i] = (vec4)(0.f, 0.f, 0.f, 0.f);
    grad_uz[i] = (vec4)(0.f, 0.f, 0.f, 0.f);

    const int move_i = imove[i];
    if(move_i <= 0){
        // Not a fluid particle
        return;
    }

    const uint c_i = icell[i];
    const vec_xyz r_i = r[i].XYZ;
    const vec_xyz u_i = u[i].XYZ;
    const float rho_i = rho[i];

    // Initialize the output
    #ifndef LOCAL_MEM_SIZE
        #define _GRADUX_ grad_ux[i].XYZ
        #define _GRADUY_ grad_uy[i].XYZ
        #define _GRADUZ_ grad_uz[i].XYZ
    #else
        #define _GRADUX_ grad_ux_l[i]
        #define _GRADUY_ grad_uy_l[i]
        #define _GRADUZ_ grad_uz_l[i]
        __local vec_xyz grad_ux_l[LOCAL_MEM_SIZE];
        __local vec_xyz grad_uy_l[LOCAL_MEM_SIZE];
        __local vec_xyz grad_uz_l[LOCAL_MEM_SIZE];
        _GRADUX_ = VEC_ZERO;
        _GRADUY_ = VEC_ZERO;
        _GRADUZ_ = VEC_ZERO;
    #endif

    // Loop over neighs
    // ================
    for(int ci = -1; ci <= 1; ci++) {
        for(int cj = -1; cj <= 1; cj++) {
            #ifdef HAVE_3D
            for(int ck = -1; ck <= 1; ck++) {
            #else
            const int ck = 0; {
            #endif
                const uint c_j = c_i +
                                ci +
                                cj * n_cells.x +
                                ck * n_cells.x * n_cells.y;
                uint j = ihoc[c_j];
                while((j < N) && (icell[j] == c_j)) {
                    if(i == j){
                        j++;
                        continue;
                    }
                    const int move_j = imove[j];
                    if((move_j != 1) && (move_j != -1)){
                        j++;
                        continue;
                    }
                    const vec_xyz r_ij = r[j].XYZ - r_i;
                    const float q = fast_length(r_ij) / H;
                    if(q >= SUPPORT)
                    {
                        j++;
                        continue;
                    }
                    {
                        #include "GradU.hcl"
                    }
                    j++;
                }
            }
        }
    }

    #ifdef LOCAL_MEM_SIZE
        grad_ux[i].XYZ = _GRADUX_;
        grad_uy[i].XYZ = _GRADUY_;
        grad_uz[i].XYZ = _GRADUZ_;
    #endif
}
