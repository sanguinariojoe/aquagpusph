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
 * @param id_in Unsorted particle indexes
 * @param id Sorted particle indexes
 * @param iset_in Unsorted set of particles indexes.
 * @param iset Sorted set of particles indexes.
 * @param imove_in Unsorted moving flags.
 *   - imove > 0 for regular fluid particles.
 *   - imove = 0 for sensors.
 *   - imove < 0 for boundary elements/particles.
 * @param imove Sorted moving flags.
 *   - imove > 0 for regular fluid particles.
 *   - imove = 0 for sensors.
 *   - imove < 0 for boundary elements/particles.
 * @param r_in Unsorted position \f$ \mathbf{r} \f$.
 * @param r Sorted position \f$ \mathbf{r} \f$.
 * @param normal_in Unsorted normal \f$ \mathbf{n} \f$.
 * @param normal Sorted normal \f$ \mathbf{n} \f$.
 * @param tangent_in Unsorted tangent \f$ \mathbf{t} \f$.
 * @param tangent Sorted tangent \f$ \mathbf{t} \f$.
 * @param id_sorted Permutations list from the unsorted space to the sorted
 * one.
 * @param N Number of particles.
 */
__kernel void entry(const __global float *eint_in, __global float *eint,
                    const __global float *deintdt, __global float *deintdt_in,
                    const __global usize *id_sorted,
                    usize N)
{
    usize i = get_global_id(0);
    if(i >= N)
        return;

    const usize i_out = id_sorted[i];

    eint[i_out] = eint_in[i];
    deintdt_in[i_out] = deintdt[i];
}

/*
 * @}
 */
