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
 * @brief Boundary integral friction term.
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

/** @brief Compute the viscous term of the energy conservation equation
 * (inserting it in the external work components):
 *
 * \f$ \left. \frac{dw_a}{dt} \right\vert_\mu =
 * + \frac{1}{\rho_a} \mathbf{u}_a \cdot \Delta \mathbf{u}
 * + \frac{1}{\rho_a} \nabla \mathbf{u}_a \cdot T \f$
 *
 * This tool is computing just the part associated to a no-slip boundary
 * condition.
 *
 * @see EnergyVisc.cl
 *
 * @param iset Set of particles index.
 * @param imove Moving flags.
 *   - imove > 0 for regular fluid particles.
 *   - imove = 0 for sensors.
 *   - imove < 0 for boundary elements/particles.
 * @param r Position \f$ \mathbf{r} \f$.
 * @param normal Normal \f$ \mathbf{n} \f$.
 * @param u Velocity \f$ \mathbf{u} \f$.
 * @param rho Density \f$ \rho \f$.
 * @param m Mass \f$ m \f$.
 * @param energy_dwdt External work (due to the viscous term).
 * @param icell Cell where each particle is located.
 * @param ihoc Head of chain for each cell (first particle found).
 * @param N Number of particles.
 * @param n_cells Number of cells in each direction
 * @param noslip_iset Particles set of the boundary terms which friction should
 * be taken into account.
 * @param dr Distance between particles \f$ \Delta r \f$.
 */
__kernel void main(const __global uint* iset,
                   const __global int* imove,
                   const __global vec* r,
                   const __global vec* normal,
                   const __global vec* u,
                   const __global float* rho,
                   const __global float* m,
                   __global float* energy_dwdt,
                   // Link-list data
                   __global uint *icell,
                   __global uint *ihoc,
                   // Simulation data
                   uint N,
                   uivec4 n_cells,
                   uint noslip_iset,
                   float dr)
{
    const uint i = get_global_id(0);
    const uint it = get_local_id(0);
    if(i >= N)
        return;
    if(imove[i] <= 0)
        return;

    const uint c_i = icell[i];
    const vec_xyz r_i = r[i].XYZ;
    const vec_xyz u_i = u[i].XYZ;
    const float rho_i = rho[i];

    // Initialize the output
    #ifndef LOCAL_MEM_SIZE
        #define _DWDT_ energy_dwdt[i]
    #else
        #define _DWDT_ energy_dwdt_l[it]
        __local float energy_dwdt_l[LOCAL_MEM_SIZE];
        _DWDT_ = energy_dwdt[i];
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
                    if((imove[j] != -3) || (iset[j] != noslip_iset)){
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
                        #include "EnergyViscNoSlipBI.hcl"
                    }
                    j++;
                }
            }
        }
    }

    #ifdef LOCAL_MEM_SIZE
        energy_dwdt[i] = _DWDT_;
    #endif
}
