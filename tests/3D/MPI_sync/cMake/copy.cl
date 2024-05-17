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

__kernel void copy_in(const __global vec* r,
                      const __global unsigned int* mask,
                      __global vec* mpi_r,
                      __global unsigned int* mpi_mask,
                      usize N)
{
    const usize i = get_global_id(0);
    if(i >= N)
        return;

    mpi_r[i] = r[i];
    mpi_mask[i] = mask[i];
}

__kernel void copy_out(__global vec* r,
                       __global unsigned int* mask,
                       const __global vec* mpi_r,
                       const __global unsigned int* mpi_mask,
                       usize N)
{
    const usize i = get_global_id(0);
    if(i >= N)
        return;

    r[i] = mpi_r[i];
    mask[i] = mpi_mask[i];
}
