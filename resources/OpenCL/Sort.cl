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

/** Sort all the variables in order to be sorted by the cell index of each
 * particle.
 */
__kernel void main(__global uint *id_in, __global uint *id,
                   __global uint *iset_in, __global uint *iset,
                   __global int *imove_in, __global int *imove,
                   __global vec *pos_in, __global vec *pos,
                   __global vec *normal_in, __global vec *normal,
                   __global vec *v_in, __global vec *v,
                   __global vec *dvdt_in, __global vec *dvdt,
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
    pos[i_out] = pos_in[i];
    normal[i_out] = normal_in[i];
    v[i_out] = v_in[i];
    dvdt[i_out] = dvdt_in[i];
    rho[i_out] = rho_in[i];
    drhodt[i_out] = drhodt_in[i];
    p[i_out] = p_in[i];
    m[i_out] = m_in[i];
}
