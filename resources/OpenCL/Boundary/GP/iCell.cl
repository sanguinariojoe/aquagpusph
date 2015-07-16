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
 * @brief Compute the cell for the mirrored dummy particle.
 */

#ifndef HAVE_3D
    #include "../../types/2D.h"
#else
    #include "../../types/3D.h"
#endif

/** @brief Compute the cell for the mirrored dummy particle.
 *
 * Since the dummy particle will be mirrored with respect to another associated
 * particle in order to ionterpolate the fields, it is eventually moving to
 * another cell, and it should be taken into account.
 *
 * @param mirrored_cell Cell where each mirrored particle is allocated.
 * @param icell Cell where each particle is located.
 * @param imove Moving flags.
 *   - imove > 0 for regular fluid particles.
 *   - imove = 0 for sensors.
 *   - imove < 0 for boundary elements/particles.
 * @param associations Mirroring particles.
 * @param r Position \f$ \mathbf{r} \f$.
 * @param normal Normal \f$ \mathbf{n} \f$.
 * @param N Number of particles.
 * @param r_min Minimum of r.
 * @param n_cells Number of cells at each direction, and the total number of
 * allocated cells.
 */
__kernel void main(__global unsigned int *mirrored_cell,
                   const __global unsigned int *icell,
                   const __global int* imove,
                   const __global uint *associations,
                   const __global vec *r,
                   const __global vec* normal,
                   unsigned int N,
                   vec r_min,
                   uivec4 n_cells)
{
    // find position in global arrays
    unsigned int i = get_global_id(0);
    if(i >= N)
        return;
    mirrored_cell[i] = icell[i];
    if(imove[i] != -1)
        return;
    const uint iref = associations[i];
    if(iref >= N)
        return;

    // Get the mirrored position
    const vec_xyz deltar = dot(r[iref].XYZ - r[i].XYZ, normal[iref].XYZ)
                           * normal[iref].XYZ;
    const vec_xyz rr = r[i].XYZ + 2.f * deltar;

    // Compute the cell
    uivec cell;
    unsigned int cell_id;
    const float idist = 1.f / (SUPPORT * H);
    cell.x = (unsigned int)((rr.x - r_min.x) * idist) + 3u;
    cell.y = (unsigned int)((rr.y - r_min.y) * idist) + 3u;
    #ifdef HAVE_3D
        cell.z = (unsigned int)((rr.z - r_min.z) * idist) + 3u;
        cell_id = cell.x - 1u +
                  (cell.y - 1u) * n_cells.x +
                  (cell.z - 1u) * n_cells.x * n_cells.y;
        mirrored_cell[i] = cell_id;
    #else
        cell_id = cell.x - 1u +
                  (cell.y - 1u) * n_cells.x;
        mirrored_cell[i] = cell_id;
    #endif
}
