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

/** @addtogroup lela
 * @{
 */

/** @file
 *  @brief Sort all the particle variables by the cell indexes
 */

#include "resources/Scripts/types/types.h"

/** @brief Sort all the particle variables by the cell indexes.
 *
 * @param p_in Unsorted pressure \f$ p \f$.
 * @param p Sorted pressure \f$ p \f$.
 * @param S_in Unsorted Deviatory stress \f$ S \f$.
 * @param S Sorted Deviatory stress \f$ S \f$.
 * @param dSdt_in Unsorted Deviatory stress rate of change
 * \f$ \frac{d S}{d t} \f$.
 * @param dSdt Sorted Deviatory stress rate of change
 * \f$ \frac{d S}{d t} \f$.
 * \f$ \left. \frac{d S}{d t} \right\vert_{n+1} \f$.
 * @param id_sorted Permutations list from the unsorted space to the sorted
 * one.
 * @param N Number of particles.
 */
__kernel void entry(const __global float *p_in, __global float *p,
                    const __global matrix *S_in, __global matrix *S,
                    const __global matrix *dSdt_in, __global matrix *dSdt,
		            const __global unit *id_sorted,
                    unsigned int N)
{
    uint i = get_global_id(0);
    if(i >= N)
        return;

    const uint i_out = id_sorted[i];

    p[i_out] = p_in[i];
    S[i_out] = S_in[i];
    dSdt[i_out] = dSdt_in[i];
}

/*
 * @}
 */