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

/** @file
 * @brief Mirroring process for the symmetry boundary condition.
 */

#include "resources/Scripts/types/types.h"

/** @brief Compute the cell where the particle is allocated.
 *
 * @param r Position \f$ \mathbf{r} \f$
 * @param r_min Minimum of r (considering all the particles).
 * @param n_cells Number of cells at each direction, and the total number of
 * allocated cells.
 */
unsigned int cell(vec r, vec r_min, uivec4 n_cells)
{
    uivec cell;

    const float idist = 1.f / (SUPPORT * H);
    cell.x = (unsigned int)((r.x - r_min.x) * idist) + 3u;
    cell.y = (unsigned int)((r.y - r_min.y) * idist) + 3u;
    #ifdef HAVE_3D
        cell.z = (unsigned int)((r.z - r_min.z) * idist) + 3u;
        return cell.x - 1u +
               (cell.y - 1u) * n_cells.x +
               (cell.z - 1u) * n_cells.x * n_cells.y;
    #else
        return cell.x - 1u +
               (cell.y - 1u) * n_cells.x;
    #endif
}

/** @brief Compute the mirrored position of the fluid particles.
 *
 * The mirrored particles (the ones close enough to the out portal plane) will
 * be marked with \a imirrored = 1.
 * 
 * @param r Position \f$ \mathbf{r} \f$
 * @param imirrored 0 if the particle has not been mirrored, 1 otherwise.
 * @param icell Cell where each particle is allocated.
 * @param N Number of particles
 * @param portal_in_r In portal infinite plane position
 * @param portal_out_r Out portal infinite plane position
 * @param portal_n Portal infinite planes normal
 * @param r_min Minimum of r.
 * @param n_cells Number of cells at each direction, and the total number of
 * allocated cells.
 */
__kernel void mirror(__global vec* r,
                     __global int* imirrored,
                     __global unsigned int *icell,
                     unsigned int N,
                     vec portal_in_r,
                     vec portal_out_r,
                     vec portal_n,
                     vec r_min,
                     uivec4 n_cells)
{
    unsigned int i = get_global_id(0);
    if(i >= N)
        return;

    // Discard the particles far away of the out portal
    const vec r_ij = r[i] - portal_out_r;
    const float dist = dot(r_ij, portal_n);
    if(fabs(dist) > SUPPORT * H){
        imirrored[i] = 0;
        return;
    }

    imirrored[i] = 1;
    r[i] = portal_in_r + r_ij;
    icell[i] = cell(r[i], r_min, n_cells);
}

/** @brief Unmirror the affected fluid particles.
 *
 * The mirrored particles (marked with \a imirrored = 1) should be unmirrored,
 * including the particles that will be teleported later.
 * Otherwise, the following particles interactions will missconsider the
 * non-unmirrored particles, like subsequent portal instances.
 * 
 * @param r Position \f$ \mathbf{r} \f$
 * @param imirrored 0 if the particle has not been mirrored, 1 otherwise.
 * @param N Number of particles
 * @param portal_in_r In portal infinite plane position
 * @param portal_out_r Out portal infinite plane position
 * @param portal_n Portal infinite planes normal
 */
__kernel void unmirror(__global vec* r,
                       const __global int* imirrored,
                       unsigned int N,
                       vec portal_in_r,
                       vec portal_out_r,
                       vec portal_n)
{
    unsigned int i = get_global_id(0);
    if(i >= N)
        return;

    if(!imirrored[i])
        return;

    const vec r_ij = r[i] - portal_in_r;
    const float dist = dot(r_ij, portal_n);
    r[i] = portal_out_r + r_ij;
}

/** @brief Teleport the particles tresspassing the out portal.
 *
 * The interactions with the other side of the portal should be computed before
 * carring out this operation.
 * 
 * @param r Position \f$ \mathbf{r} \f$
 * @param N Number of particles
 * @param portal_in_r In portal infinite plane position
 * @param portal_out_r Out portal infinite plane position
 * @param portal_n Portal infinite planes normal
 */
__kernel void teleport(__global vec* r,
                       unsigned int N,
                       vec portal_in_r,
                       vec portal_out_r,
                       vec portal_n)
{
    unsigned int i = get_global_id(0);
    if(i >= N)
        return;

    // Discard the particles far away of the out portal
    const vec r_ij = r[i] - portal_out_r;
    const float dist = dot(r_ij, portal_n);
    if(dist < 0.f)
        return;

    r[i] = portal_in_r + r_ij;
}