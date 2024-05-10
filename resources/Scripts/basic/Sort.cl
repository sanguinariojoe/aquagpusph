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
__kernel void stage1(const __global uint *id_in, __global uint *id,
                     const __global uint *iset_in, __global uint *iset,
                     const __global int *imove_in, __global int *imove,
                     const __global vec *r_in, __global vec *r,
                     const __global vec *normal_in, __global vec *normal,
                     const __global vec *tangent_in, __global vec *tangent,
                     const __global unit *id_sorted,
                     unsigned int N)
{
    uint i = get_global_id(0);
    if(i >= N)
        return;

    const uint i_out = id_sorted[i];

    id[i_out] = id_in[i];
    iset[i_out] = iset_in[i];
    imove[i_out] = imove_in[i];
    r[i_out] = r_in[i];
    normal[i_out] = normal_in[i];
    tangent[i_out] = tangent_in[i];
}

/** @brief Sort all the particle variables by the cell indexes.
 *
 * Due to the large number of registers consumed (21 global memory arrays are
 * loaded), it is safer carrying out the sorting process in to stages. This
 * is the first stage.
 *
 * @param rho_in Unsorted density \f$ \rho \f$.
 * @param rho Sorted density \f$ \rho \f$.
 * @param m_in Unsorted mass \f$ m \f$.
 * @param m Sorted mass \f$ m \f$.
 * @param u_in Unsorted velocity \f$ \mathbf{u} \f$.
 * @param u Sorted velocity \f$ \mathbf{u} \f$.
 * @param dudt Unsorted velocity rate of change
 * \f$ \frac{d \mathbf{u}}{d t} \f$.
 * @param dudt_in Sorted velocity rate of change
 * \f$ \frac{d \mathbf{u}}{d t} \f$.
 * @param drhodt Unsorted density rate of change \f$ \frac{d \rho}{d t} \f$.
 * @param drhodt_in Sorted density rate of change \f$ \frac{d \rho}{d t} \f$.
 * @param id_sorted Permutations list from the unsorted space to the sorted
 * one.
 * @param N Number of particles.
 */
__kernel void stage2(const __global float *rho_in, __global float *rho,
                     const __global float *m_in, __global float *m,
                     const __global vec *u_in, __global vec *u,
                     const __global vec *dudt, __global vec *dudt_in,
                     const __global float *drhodt, __global float *drhodt_in,
                     const __global unit *id_sorted,
                     unsigned int N)
{
    uint i = get_global_id(0);
    if(i >= N)
        return;

    const uint i_out = id_sorted[i];

    rho[i_out] = rho_in[i];
    m[i_out] = m_in[i];
    u[i_out] = u_in[i];
    // Take care with the variation rates, since they are epheremeral. Thus,
    // the output is actually the _in variable, whilst the other one will be
    // overwritten
    dudt_in[i_out] = dudt[i];
    drhodt_in[i_out] = drhodt[i];
}

/*
 * @}
 */
