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
 * @brief Bounding boxes for processes
 */

#include "resources/Scripts/types/types.h"

/** @brief Set the masks to let AQUAgpusph to which process it shall send the
 * particle
 *
 * Just fluid particles are exchanged
 *
 * @param imove Moving flags
 *   - imove > 0 for regular fluid particles
 *   - imove = 0 for sensors
 *   - imove < 0 for boundary elements/particles
 * @param mpi_local_mask Particles that shall be handled by another process
 * from now on
 * @param mpi_neigh_mask Particles that shall be sent to another process so it
 * can consider them as neighbours
 * @param mpi_plane_r Bounding plane position
 * @param mpi_plane_n Outward bounding plane normal
 * @param mpi_plane_proc Neighbour process
 * @param N Number of particles
 */
__kernel void entry(const __global int* imove,
                    const __global vec* r,
                    __global unsigned int* mpi_local_mask,
                    __global unsigned int* mpi_neigh_mask,
                    vec mpi_plane_r,
                    vec mpi_plane_n,
                    unsigned int mpi_plane_proc,
                    unsigned int N)
{
    unsigned int i = get_global_id(0);
    if(i >= N)
        return;
    if(imove[i] <= 0)
        return;

    const float d = dot(r[i] - mpi_plane_r, mpi_plane_n);
    if (d > 0.f) {
        // The particle has crossed the boundary, so it is not our problem
        // anymore
        mpi_local_mask[i] = mpi_plane_proc;
    } else if (d >= -SUPPORT * H) {
        // The particle is needed by the neighboring process to compute
        mpi_neigh_mask[i] = mpi_plane_proc;
    }
}
