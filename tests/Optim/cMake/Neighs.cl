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

#include "resources/Scripts/types/types.h"
#include "resources/Scripts/KernelFunctions/Kernel.h"

__kernel void legacy(const __global int* imove,
                     __global uint* n_neighs_legacy,
                     usize N,
                     LINKLIST_LOCAL_PARAMS)
{
    const usize i = get_global_id(0);
    const usize it = get_local_id(0);
    if(i >= N)
        return;
    if(imove[i] <= -255) {
        n_neighs_legacy[i] = 0;
        return;
    }

    __private uint __n_neighs = 0;

    const usize c_i = icell[i];
    BEGIN_NEIGHS(c_i, N, n_cells, icell, ihoc){
        __n_neighs += 1;
    }END_NEIGHS()

    n_neighs_legacy[i] = __n_neighs;
}

__kernel void optim(const __global int* imove,
                    const __global svec2* jhoc,
                    __global uint* n_neighs_optim,
                    usize N)
{
    const usize i = get_global_id(0);
    const usize it = get_local_id(0);
    if(i >= N)
        return;
    if(imove[i] <= -255) {
        n_neighs_optim[i] = 0;
        return;
    }

    __private uint __n_neighs = 0;

    for(unsigned int row = 0; row < NNC; row++) {
        const unsigned int jhoc_id = i + row * N;
        __n_neighs += jhoc[jhoc_id].y - jhoc[jhoc_id].x;
    }

    n_neighs_optim[i] = __n_neighs;
}

__kernel void err(const __global uint* n_neighs_legacy,
                  const __global uint* n_neighs_optim,
                  __global uint* n_neighs_err,
                  usize N)
{
    const usize i = get_global_id(0);
    if(i >= N)
        return;

    n_neighs_err[i] = abs(n_neighs_optim[i] - n_neighs_legacy[i]);
}
