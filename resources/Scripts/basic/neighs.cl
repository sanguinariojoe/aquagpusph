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
 * @brief Compute the number of neighbours of each particle
 */

#if defined(LOCAL_MEM_SIZE) && defined(NO_LOCAL_MEM)
    #error NO_LOCAL_MEM has been set.
#endif

#include "resources/Scripts/types/types.h"
#include "resources/Scripts/KernelFunctions/Kernel.h"

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
 * @param r Position \f$ \mathbf{r} \f$.
 * @param n_neighs Number of neighbours per particle.
 * @param icell Cell where each particle is located.
 * @param ihoc Head of chain for each cell (first particle found).
 * @param neighs_limit The largest number of neighbours accepted.
 * @param N Number of particles.
 * @param n_cells Number of cells in each direction
 */
__kernel void entry(const __global int* imove,
                    __global uint* n_neighs,
                    uint neighs_limit,
                    uint N,
                    LINKLIST_LOCAL_PARAMS)
{
    const uint i = get_global_id(0);
    const uint it = get_local_id(0);
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

    const unsigned int c_i = icell[i];
    BEGIN_NEIGHS(c_i, N, n_cells, icell, ihoc){
        _NEIGHS_ += 1;
        if(_NEIGHS_ >= neighs_limit){
            // Ops! Too much neighbours! Stop right now!
            #ifdef LOCAL_MEM_SIZE
                n_neighs[i] = _NEIGHS_;
            #endif
            return;
        }
    }END_NEIGHS()

    #ifdef LOCAL_MEM_SIZE
        n_neighs[i] = _NEIGHS_;
    #endif
}

/*
 * @}
 */
