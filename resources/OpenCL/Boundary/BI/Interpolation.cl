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
 * @brief Fluid particles interactions computation.
 */

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

/** @brief Pressure interpolation at the boundary elements.
 *
 * The values are computed using just the fluid information. The resulting
 * interpolated values are not renormalized yet.
 *
 * Just the elements with the flag imove = -3 are considered boundary elements.
 *
 * @param imove Moving flags.
 *   - imove > 0 for regular fluid particles.
 *   - imove = 0 for sensors.
 *   - imove < 0 for boundary elements/particles.
 * @param r Position \f$ \mathbf{r} \f$.
 * @param m Mass \f$ m \f$.
 * @param rho Density \f$ \rho \f$.
 * @param p Pressure \f$ p \f$.
 * @param icell Cell where each particle is located.
 * @param ihoc Head of chain for each cell (first particle found).
 * @param N Number of particles.
 * @param n_cells Number of cells in each direction
 * @param g Gravity acceleration \f$ \mathbf{g} \f$.
 * @see SensorsRenormalization.cl
 */
__kernel void main(const __global int* imove,
                   const __global vec* r,
                   const __global float* m,
                   const __global float* rho,
                   __global float* p,
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
    if(imove[i] != -3){
        return;
    }

    const uint c_i = icell[i];
    const vec_xyz r_i = r[i].XYZ;

    // Initialize the output
    #ifndef LOCAL_MEM_SIZE
        #define _P_ p[i]
    #else
        #define _P_ p_l[it]
        __local float p_l[LOCAL_MEM_SIZE];
    #endif
    _P_ = 0.f;

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
                    if(imove[j] != 1){
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
                        const float w_ij = kernelW(q) * CONW * m[j] / rho[j];
                        _P_ += (p[j] - dot(g.XYZ, r_ij)) * w_ij;
                    }
                    j++;
                }
            }
        }
    }

    #ifdef LOCAL_MEM_SIZE
        p[i] = _P_;
    #endif
}
