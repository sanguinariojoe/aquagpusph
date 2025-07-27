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
 * @brief Carry out operations regarding the neightbour chains of each
 * particle
 */

#if defined(LOCAL_MEM_SIZE) && defined(NO_LOCAL_MEM)
    #error NO_LOCAL_MEM has been set.
#endif

#include "resources/Scripts/types/types.h"
#include "resources/Scripts/KernelFunctions/Kernel.h"

/** @brief Move to a cell based neighbour chain to a particle-by-particle one.
 *
 * While traversing the cell based link-list requires reading 2 arrays (icell
 * and ihoc) as well as a scalar (n_cells), the particle-by-particle one only
 * requires reading jhoc. On top of that, the cell-by-cell link-list is making
 * more complex to get an efficient memory banking, while on the
 * particle-by-particle the information is read on a quite nice, with boosted
 * coalescence.
 *
 * Thus, spending a bit of computation here will easy the subsequent
 * interactions computations.
 *
 * @param icell Cell where each particle is located.
 * @param ihoc Head and tail of chain for each cell.
 * @param jhoc First and last particles on each neighbour chain.
 * @param N Number of particles.
 * @param n_cells Number of cells in each direction
 */
__kernel void neigh_chains(const __global usize *icell,
                           const __global svec2 * ihoc,
                           __global svec2 * jhoc,
                           usize N,
                           svec4 n_cells)
{
    const usize i = get_global_id(0);
    if(i >= N)
        return;

    const usize c_i = icell[i];
    for(int cy = -1; cy <= 1; cy++) {
#ifdef HAVE_3D
        for(int cz = -1; cz <= 1; cz++) {
            const usize iout = i + N * ((cz + 1) + (cy + 1) * 3);
#else
        const int cz = 0; {
            const usize iout = i + N * (cy + 1);
#endif
            jhoc[iout] = (svec2)(N);
            // Look for candidates on the same row of cells
            for(int cx = -1; cx <= 1; cx++) {
                const uint c_j = c_i +
                                 cx +
                                 cy * n_cells.x +
                                 cz * n_cells.x * n_cells.y;
                if (ihoc[c_j].x >= N)
                    continue;
                // We are always taking the last possible tail of chain
                jhoc[iout].y = ihoc[c_j].y;
                // But we only want the first hit for head of chain
                if (jhoc[iout].x < N) {
                    jhoc[iout].x = ihoc[c_j].x;
                }
            }
        }
    }
}

/** @brief Number of neighbours of each particle.
 *
 * One of the main targets of this kernel is checking that the number of
 * neighbours is not excessively large. Along this line, if #neighs_limit
 * neighbours are reached, the kernel will stop the execution.
 *
 * @param imove Moving flags.
 *   - imove > 0 for regular fluid particles.
 *   - imove = 0 for sensors.
 *   - imove < 0 for boundary elements/particles.
 * @param jhoc First and last particles on each neighbour chain.
 * @param n_neighs Number of neighbours per particle.
 * @param N Number of particles.
 */
__kernel void entry(const __global int* imove,
                    const __global svec2* jhoc,
                    __global uint* n_neighs,
                    uint neighs_limit,
                    usize N,
                    LINKLIST_LOCAL_PARAMS)
{
    const usize i = get_global_id(0);
    const usize it = get_local_id(0);
    if(i >= N)
        return;
    if(imove[i] <= -255) {
        n_neighs[i] = 0;
        return;
    }

    // Initialize the output
    #ifndef LOCAL_MEM_SIZE
        #define _NEIGHS_ n_neighs[i]
    #else
        #define _NEIGHS_ n_neighs_l[it]
        __local uint n_neighs_l[LOCAL_MEM_SIZE];
    #endif
    _NEIGHS_ = 0;

    for(unsigned int row = 0; row < NNC; row++) {
        const unsigned int jhoc_id = i + row * N;
        _NEIGHS_ += jhoc[jhoc_id].y - jhoc[jhoc_id].x;
    }

    #ifdef LOCAL_MEM_SIZE
        n_neighs[i] = _NEIGHS_;
    #endif
}

/*
 * @}
 */
