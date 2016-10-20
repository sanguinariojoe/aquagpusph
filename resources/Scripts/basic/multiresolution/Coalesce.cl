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
 *  @brief Splitting particles methods
 */

#ifndef HAVE_3D
    #include "../../types/2D.h"
#else
    #include "../../types/3D.h"
#endif

#ifndef M_ITERS
    #define M_ITERS 10
#endif


/** @brief Look for all the seed candidates, i.e. all the particles which have a
 * refinement level target lower than their current value.
 *
 * This tool is also placing the particles into cells with the volume of the
 * partner particle to become generated 
 *
 * @param iset Set of particles index.
 * @param imove Moving flags.
 *   - imove > 0 for regular fluid/solid particles.
 *   - imove = 0 for sensors.
 *   - imove < 0 for boundary elements/particles.
 * @param ilevel Current refinement level of the particle.
 * @param level Target refinement level of the particle.
 * @param m0 Mass \f$ m_0 \f$.
 * @param miter Mass transfer iteration (Positive for growing particles,
 * negative for shrinking particles).
 * @param r Position \f$ \mathbf{r} \f$.
 * @param isplit 0 if the particle should not become coalesced, 1 for the
 * coalescing particles, 2 for the seeds.
 * @param split_cell Seed cell where the particle is located.
 * @param split_dist Distance to the center of the cell. Just the closest
 * particle to the cell center will be kept as seed
 * @param dr_level0 Theoretical distance between particles at the lowest
 * refinement level.
 * @param N Number of particles.
 */
__kernel void seed_candidates(__global const unsigned int* iset,
                              __global const int* imove,
                              __global const unsigned int* ilevel,
                              __global const unsigned int* level,
                              __global const float* m0,
                              __global const vec* r,
                              __global int* miter,
                              __global unsigned int* isplit,
                              __global ivec* split_cell,
                              __global float* split_dist,
                              __constant float* dr_level0,
                              unsigned int N)
{
    unsigned int i = get_global_id(0);
    if(i >= N)
        return;

    isplit[i] = 0;
    if(imove[i] <= 0){
        // Neglect boundary elements/particles
        return;
    }
    else if(miter[i] <= M_ITERS){
        // The particle is growing or shrinking... Let's wait
        return;
    }
    else if(level[i] >= ilevel[i]){
        // The particle is not a candidate to become coalesced
        return;
    }

    isplit[i] = 2;
    miter[i] = -1;
    const float dr = dr_level0[iset[i]] / (ilevel[i] - 1);
    split_cell[i] = CONVERT(ivec, r[i] / dr);
    const vec r_cell = r[i] - (CONVERT(vec, split_cell[i]) + VEC_ONE * 0.5f * dr);
    split_dist[i] = length(r_cell);
}

/** @brief Get only one seed per cell.
 *
 * To do that, all the seed candidates will inspect their neighbours, and in
 * case they find another seed candidate closer to the cell center, they are
 * desisting to become a seed.
 *
 * @param iset Set of particles index.
 * @param isplit 0 if the particle should not become coalesced, 1 for the
 * coalescing particles, 2 for the seeds.
 * @param split_cell Seed cell where the particle is located.
 * @param split_dist Distance to the center of the cell. Just the closest
 * particle to the cell center will be kept as seed
 * @param N Number of particles.
 */
__kernel void seeds(__global const unsigned int* iset,
                    __global unsigned int* isplit,
                    __global const ivec* split_cell,
                    __global const float* split_dist,
                    __global const uint *icell,
                    __global const uint *ihoc,
                    unsigned int N,
                    uivec4 n_cells)
{
    unsigned int i = get_global_id(0);
    if(i >= N)
        return;
    if(isplit[i] != 2)
        return;

    const ivec_xyz isplit_cell = split_cell[i].XYZ;
    const float isplit_dist = split_dist[i];

    BEGIN_LOOP_OVER_NEIGHS(){
        if(i == j){
            j++;
            continue;
        }
        if((isplit[j] != 2) ||
           (iset[i] != iset[j]) ||
           (any(isplit_cell != split_cell[j].XYZ))
        ){
            j++;
            continue;
        }

        {
            const float jsplit_dist = split_dist[i];
            if(isplit_dist > jsplit_dist){
                // I'm not a seed, because there is another better candidate
                isplit[i] = 1;
                return;
            }
            else if(isplit_dist == jsplit_dist){
                // A very unlikely situation... Let use just the last particle
                if(i < j){
                    isplit[i] = 1;
                    return;
                }
            }
        }
    }END_LOOP_OVER_NEIGHS()
}

/*
 * @}
 */