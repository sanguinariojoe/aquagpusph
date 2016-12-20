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
 *  @brief Sort all the particle variables by the cell indexes
 */

#include "resources/Scripts/types/types.h"

/** @brief Sort all the particle variables by the cell indexes
 *
 * Due to the large number of registers consumed (21 global memory arrays are
 * loaded), it is safer carrying out the sorting process in to stages. This
 * is the first stage.
 *
 * @param h_var_in Unsorted variable kernel lenght \f$ h \f$.
 * @param h_var Sorted variable kernel lenght \f$ h \f$.
 * @param id_sorted Permutations list from the unsorted space to the sorted
 * one.
 * @param N Number of particles.
 */
__kernel void entry(const __global float *h_var_in, __global float *h_var,
                    const __global unit *id_sorted,
                    unsigned int N)
{
    uint i = get_global_id(0);
    if(i >= N)
        return;

    const uint i_out = id_sorted[i];

    h_var[i_out] = h_var_in[i];
}

/*
 * @}
 */
