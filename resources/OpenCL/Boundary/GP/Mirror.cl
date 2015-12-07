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
 * @brief Fixed ghost particles mirroring process.
 */

#ifndef HAVE_3D
    #include "../../types/2D.h"
    #include "../../KernelFunctions/Wendland2D.hcl"
#else
    #include "../../types/3D.h"
    #include "../../KernelFunctions/Wendland3D.hcl"
#endif

/** @brief Fixed ghost particles mirroring.
 *
 * The fixed ghost particles require to mirror the particles in order to can
 * interpolate the field values, mirroring back later to can compute the
 * interactions with the fluid particles.
 *
 * In the mirroring process, both the particle and its cell information should
 * be changed.
 *
 * @param imove Moving flags.
 *   - imove > 0 for regular fluid particles.
 *   - imove = 0 for sensors.
 *   - imove < 0 for boundary elements/particles.
 * @param associations Mirroring particles.
 * @param normal Normal \f$ \mathbf{n} \f$.
 * @param r Position \f$ \mathbf{r} \f$.
 * @param gp_icell Cell where each particle is located.
 * @param N Number of particles.
 * @param r_min Minimum position of a particle
 * @param n_cells Number of cells in each direction
 */
__kernel void main(const __global int* imove,
                   const __global uint* associations,
                   const __global vec* normal,
                   __global vec* r,
                   __global uint *gp_icell,
                   uint N,
                   vec r_min,
                   uivec4 n_cells)
{
    const uint i = get_global_id(0);
    if(i >= N)
        return;
    if(imove[i] != -1)
        return;

    // Let's get the associated boundary element, and check it is valid
    const uint iref = associations[i];
    if(iref >= N)
        return;

    // Mirror the particle
    const vec_xyz r_iref = r[iref].XYZ;
    const vec_xyz n_iref = normal[iref].XYZ;
    const vec_xyz dr_i = dot(r_iref - r[i].XYZ, n_iref) * n_iref;
    r[i].XYZ += 2.f * dr_i;

    // Compute the new cell
    uivec cell;
    unsigned int cell_id;
    const float idist = 1.f / (SUPPORT * H);
    cell.x = (unsigned int)((r[i].x - r_min.x) * idist) + 3u;
    cell.y = (unsigned int)((r[i].y - r_min.y) * idist) + 3u;
    #ifdef HAVE_3D
        cell.z = (unsigned int)((r[i].z - r_min.z) * idist) + 3u;
        cell_id = cell.x - 1u +
                  (cell.y - 1u) * n_cells.x +
                  (cell.z - 1u) * n_cells.x * n_cells.y;
        gp_icell[i] = cell_id;
    #else
        cell_id = cell.x - 1u +
                  (cell.y - 1u) * n_cells.x;
        gp_icell[i] = cell_id;
    #endif
}
