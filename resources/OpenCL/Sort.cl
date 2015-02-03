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

#ifndef HAVE_3D
    #include "types/2D.h"
#else
    #include "types/3D.h"
#endif

/** 

@brief Sort all the particle variables by the cell indexes.
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
 * @param u_in Unsorted velocity \f$ \mathbf{u} \f$.
 * @param u Sorted velocity \f$ \mathbf{u} \f$.
 * @param dudt_in Unsorted velocity rate of change
 * \f$ \frac{d \mathbf{u}}{d t} \f$.
 * @param dudt Sorted velocity rate of change
 * \f$ \frac{d \mathbf{u}}{d t} \f$.
 * @param rho_in Unsorted density \f$ \rho \f$.
 * @param rho Sorted density \f$ \rho \f$.
 * @param drhodt_in Unsorted density rate of change \f$ \frac{d \rho}{d t} \f$.
 * @param drhodt Sorted density rate of change \f$ \frac{d \rho}{d t} \f$.
 * @param p_in Unsorted pressure \f$ p \f$.
 * @param p Sorted pressure \f$ p \f$.
 * @param m_in Unsorted mass \f$ m \f$.
 * @param m Sorted mass \f$ m \f$.
 * @param id_sorted Permutations list from the unsorted space to the sorted
 * one.
 * @param N Number of particles.
 */
__kernel void main(__global uint *id_in, __global uint *id,
                   __global uint *iset_in, __global uint *iset,
                   __global int *imove_in, __global int *imove,
                   __global vec *r_in, __global vec *r,
                   __global vec *normal_in, __global vec *normal,
                   __global vec *u_in, __global vec *u,
                   __global vec *dudt_in, __global vec *dudt,
                   __global float *rho_in, __global float *rho,
                   __global float *drhodt_in, __global float *drhodt,
                   __global float *p_in, __global float *p,
                   __global float *m_in, __global float *m,
                   __global unit *id_sorted,
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
    u[i_out] = u_in[i];
    dudt[i_out] = dudt_in[i];
    rho[i_out] = rho_in[i];
    drhodt[i_out] = drhodt_in[i];
    p[i_out] = p_in[i];
    m[i_out] = m_in[i];
}
