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

#ifndef HAVE_3D
    #include "../types/2D.h"
#else
    #include "../types/3D.h"
#endif

/** @brief Sort all the particle variables by the cell indexes.
 *
 * @param r0_in Unsorted equilibrium position \f$ \mathbf{r}_0 \f$.
 * @param r0 Sorted equilibrium position \f$ \mathbf{r}_0 \f$.
 * @param r_r0_in Unsorted deformation \f$ \mathbf{r}^{*} - \mathbf{r}_0 \f$.
 * @param r_r0 Sorted deformation \f$ \mathbf{r}^{*} - \mathbf{r}_0 \f$.
 * @param id_sorted Permutations list from the unsorted space to the sorted
 * one.
 * @param N Number of particles.
 */
__kernel void entry(const __global vec *r0_in, __global vec *r0,
		            const __global vec *r_r0_in, __global vec *r_r0,
		            const __global unit *id_sorted,
                    unsigned int N)
{
    uint i = get_global_id(0);
    if(i >= N)
        return;

    const uint i_out = id_sorted[i];

    r0[i_out] = r0_in[i];
    r_r0[i_out] = r_r0_in[i];
}

/*
 * @}
 */