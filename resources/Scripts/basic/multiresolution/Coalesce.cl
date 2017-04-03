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

#include "resources/Scripts/types/types.h"

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
                              __global int* miter,
                              __global const vec* r,
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
    split_cell[i] = CONVERT(ivec, VEC_ZERO);
    split_dist[i] = MAXFLOAT;

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
    const float dr = dr_level0[iset[i]] / ilevel[i];
    // We add a "virtual" cell to clearly split negative and positive
    // (otherwise it maybe rounded and mixed)
    split_cell[i] = CONVERT(ivec, r[i] / dr + sign(r[i]));
    const vec r_cell = (CONVERT(vec, split_cell[i]) + 0.5f * VEC_ONE) * dr;
    split_dist[i] = length(r[i] - r_cell);
}

/** @brief Get only one seed per cell.
 *
 * To do that, all the seed candidates will inspect their neighbours, and in
 * case they find another seed candidate closer to the cell center, they are
 * desisting to become a seed.
 *
 * @param iset Set of particles index.
 * @param isplit_in 0 if the particle should not become coalesced, 1 for the
 * coalescing particles, 2 for the seeds. This array is used to read.
 * @param isplit 0 if the particle should not become coalesced, 1 for the
 * coalescing particles, 2 for the seeds. This array is used to write.
 * @param miter Mass transfer iteration (Positive for growing particles,
 * negative for shrinking particles).
 * @param split_cell Seed cell where the particle is located.
 * @param split_dist Distance to the center of the cell. Just the closest
 * particle to the cell center will be kept as seed
 * @param icell Cell where each particle is located.
 * @param ihoc Head of chain for each cell (first particle found).
 * @param N Number of particles.
 * @param n_cells Number of cells in each direction
 */
__kernel void seeds(__global const unsigned int* iset,
                    __global const unsigned int* isplit_in,
                    __global unsigned int* isplit,
                    __global int* miter,
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
    if(isplit_in[i] != 2)
        return;

    const ivec_xyz isplit_cell = split_cell[i].XYZ;
    const float isplit_dist = split_dist[i];

    BEGIN_LOOP_OVER_NEIGHS(){
        if(i == j){
            j++;
            continue;
        }

        if((isplit_in[j] != 2) ||
           (iset[i] != iset[j]) ||
           (any(isplit_cell != split_cell[j].XYZ))
        ){
            j++;
            continue;
        }

        {
            const float jsplit_dist = split_dist[j];
            if(isplit_dist > jsplit_dist){
                // I'm not a seed, because there is another better candidate
                isplit[i] = 0;
                miter[i] = M_ITERS + 1;
            }
            else if(isplit_dist == jsplit_dist){
                // A very unlikely situation... Let use just the last particle
                if(i < j){
                    isplit[i] = 0;
                    miter[i] = M_ITERS + 1;
                }
            }
        }
    }END_LOOP_OVER_NEIGHS()
}

/** @brief Create a copy of isplit, where everything is 0 except the seeds,
 * which take the value 1. Such array can be used to count the number of new
 * particles to become generated.
 *
 * Since isplit is used to sort the particles, it should have "n_radix" items,
 * which is bigger than "N". But in order to conveniently sort isplit, you must
 * ensure that the values "out of bounds" are bigger than the other ones (and
 * therefore kept at the end of the list). in this case a value of 3 is
 * selected.
 *
 * @param isplit 0 if the particle should not become split, 1 otherwise
 * @param isplit_in 0 if the particle should not become split, 1 otherwise
 * @param N Number of particles.
 */
__kernel void set_isplit_in(__global const unsigned int* isplit,
                            __global unsigned int* isplit_in,
                            unsigned int N)
{
    unsigned int i = get_global_id(0);
    if(i >= N)
        return;

    if(isplit[i] == 2)
        isplit_in[i] = 1;
}


/** @brief Look for all children, close enough to the seed.
 *
 * When a particle is not already splitting or coalescing, and is close enough
 * to one or more seeds (which should has the same refinement level and belongs
 * to the same particles set), then it is marked as coalescing child.
 *
 * Hence, since a particle may be close enough to several partner particles,
 * later we should compute contribution weights.
 *
 * @param imove Moving flags.
 *   - imove > 0 for regular fluid/solid particles.
 *   - imove = 0 for sensors.
 *   - imove < 0 for boundary elements/particles.
 * @param iset Set of particles index.
 * @param r Position \f$ \mathbf{r} \f$.
 * @param ilevel Current refinement level of the particle.
 * @param isplit 0 if the particle should not become coalesced, 1 for the
 * coalescing particles, 2 for the seeds.
 * @param miter Mass transfer iteration (Positive for growing particles,
 * negative for shrinking particles).
 * @param icell Cell where each particle is located.
 * @param ihoc Head of chain for each cell (first particle found).
 * @param dr_level0 Theoretical distance between particles at the lowest
 * refinement level.
 * @param N Number of particles.
 * @param n_cells Number of cells in each direction
 */
__kernel void children(__global const int* imove,
                       __global const unsigned int* iset,
                       __global const vec* r,
                       __global const unsigned int* ilevel,
                       __global unsigned int* isplit,
                       __global int* miter,
                       __global const uint *icell,
                       __global const uint *ihoc,
                       __constant float* dr_level0,
                       unsigned int N,
                       uivec4 n_cells)
{
    unsigned int i = get_global_id(0);
    if(i >= N)
        return;
    if(isplit[i] != 2)
        return;

    const vec_xyz r_i = r[i].XYZ;
    const float dr = dr_level0[iset[i]] / ilevel[i];


    BEGIN_LOOP_OVER_NEIGHS(){
        if((isplit[j] != 0) ||        // It is a seed or already a child
           (imove[j] <= 0) ||         // Neglect boundaries/sensors
           (miter[j] <= M_ITERS) ||   // It is already splitting/coalescing
           (iset[i] != iset[j]) ||    // Different set of particles
           (ilevel[i] != ilevel[j])   // Different level of refinement
        ){
            j++;
            continue;
        }

        if(all(islessequal(fabs(r[j].XYZ - r_i),
                           0.6f * dr * VEC_ONE.XYZ))
        ){
            // It does not matter if several threads write this data at the same
            // time, because are constant values
            isplit[j] = 1;
            miter[j] = -1;
        }
    }END_LOOP_OVER_NEIGHS()
}


/** @brief Compute the contribution weight of each children particle.
 *
 * The children particles (including the seed itself) may be close enough to
 * several partner particles, such that, at the time of computing the partner
 * particles properties, it is eventually contributing to several partners.
 * Hence, the contributions should be conveniently weighted
 *
 * @see children()
 * @param iset Set of particles index.
 * @param ilevel Current refinement level of the particle.
 * @param r Position \f$ \mathbf{r} \f$.
 * @param isplit 0 if the particle should not become coalesced, 1 for the
 * coalescing particles, 2 for the seeds.
 * @param split_weight Coalescing contribution weight.
 * @param icell Cell where each particle is located.
 * @param ihoc Head of chain for each cell (first particle found).
 * @param dr_level0 Theoretical distance between particles at the lowest
 * refinement level.
 * @param N Number of particles.
 * @param n_cells Number of cells in each direction
 */
__kernel void weights(__global const unsigned int* iset,
                      __global const unsigned int* ilevel,
                      __global const vec* r,
                      __global const unsigned int* isplit,
                      __global float* split_weight,
                      __global const uint *icell,
                      __global const uint *ihoc,
                      __constant float* dr_level0,
                      unsigned int N,
                      uivec4 n_cells)
{
    unsigned int i = get_global_id(0);
    if(i >= N)
        return;
    if((isplit[i] != 1) && (isplit[i] != 2))
        return;

    const float dr = dr_level0[iset[i]] / ilevel[i];
    const vec_xyz r_i = r[i].XYZ;
    // split_weight[i] = 1.f;

    BEGIN_LOOP_OVER_NEIGHS(){
        if((isplit[j] != 2) ||             // Not a seed
           (iset[i] != iset[j]) ||         // Different set of particles
           (ilevel[i] != ilevel[j])        // Different level of refinement
        ){
            j++;
            continue;
        }

        if(all(islessequal(fabs(r[j].XYZ - r_i),
                           0.6f * dr * VEC_ONE.XYZ))
        ){
            split_weight[i] += 1.f;
        }
    }END_LOOP_OVER_NEIGHS()

    split_weight[i] = 1.f / split_weight[i];
}


/** @brief Associate a buffer particle to a seed.
 *
 * Later we are looking for the rest of children, used to compute the variable
 * values of this buffer/partner particle
 *
 * @param imove Moving flags.
 *   - imove > 0 for regular fluid/solid particles.
 *   - imove = 0 for sensors.
 *   - imove < 0 for boundary elements/particles.
 * @param iset Index of the set of particles.
 * @param isplit 0 if the particle should not become coalesced, 1 for the
 * coalescing particles, 2 for the seeds.
 * @param split_invperm Permutation to find the index of the particle in the
 * list of particles to become split.
 * @param mybuffer Index of the partner buffer particle.
 * @param ilevel Current refinement level of the particle.
 * @param level Target refinement level of the particle.
 * @param miter Mass transfer iteration (Positive for growing particles,
 * negative for shrinking particles).
 * @param mybuffer Partner particle associated to the seed
 * @param N Number of particles.
 * @param nbuffer Number of available buffer particles.
 */
__kernel void generate(__global int* imove,
                       __global int* iset,
                       __global const unsigned int* isplit,
                       __global unsigned int* split_invperm,
                       __global unsigned int* ilevel,
                       __global unsigned int* level,
                       __global int* miter,
                       __global uint *mybuffer,
                       unsigned int N,
                       unsigned int nbuffer)
{
    unsigned int i = get_global_id(0);
    if(i >= N)
        return;
    mybuffer[i] = N;

    // Neglect boundary elements/particles
    if(imove[i] <= 0)
        return;

    // Check whether the particle is a seed or not
    const unsigned int j = split_invperm[i];
    if(isplit[j] != 2){
        return;
    }

    // Compute the index of the buffer particle to steal. Take care, the radix
    // sort is storing the particles to become split at the end of the list
    const unsigned int i0 = N - nbuffer;
    unsigned int ii = i0 + (N - j - 1);
    // Check that there are buffer particles enough
    if(ii >= N){
        // PROBLEMS! This seed cannot generate a partner
        return;
    }
    mybuffer[i] = ii;

    // Set the already known variables
    imove[ii] = imove[i];
    iset[ii] = iset[i];
    ilevel[ii] = ilevel[i] - 1;
    level[ii] = ilevel[i] - 1;
    miter[ii] = 1;
}


/** @brief Collect the children, and the seed itself, in order to compute the
 * field values of the buffer partner particle.
 *
 * Since each children may contribute to several partners, we are weighting
 * their field values.
 *
 * @param isplit 0 if the particle should not become coalesced, 1 for the
 * coalescing particles, 2 for the seeds.
 * @param split_weight Coalescing contribution weight.
 * @param m0 Target mass, \f$ m_0 \f$.
 * @param m Current mass, \f$ m \f$.
 * @param r Position \f$ \mathbf{r} \f$.
 * @param u Velocity \f$ \mathbf{u} \f$.
 * @param dudt Velocity rate of change \f$ \frac{d \mathbf{u}}{d t} \f$.
 * @param rho Density \f$ \rho \f$.
 * @param drhodt Density rate of change \f$ \frac{d \rho}{d t} \f$.
 * @param icell Cell where each particle is located.
 * @param ihoc Head of chain for each cell (first particle found).
 * @param N Number of particles.
 * @param n_cells Number of cells in each direction
 */
__kernel void fields(__global const unsigned int* iset,
                     __global const uint* isplit,
                     __global const uint* mybuffer,
                     __global const unsigned int* ilevel,
                     __global const float* split_weight,
                     __global float* m0,
                     __global float* m,
                     __global vec* r,
                     __global vec* u,
                     __global vec* dudt,
                     __global float* rho,
                     __global float* drhodt,
                     __global const uint *icell,
                     __global const uint *ihoc,
                     __constant float* dr_level0,
                     unsigned int N,
                     uivec4 n_cells)
{
    unsigned int i = get_global_id(0);
    if(i >= N)
        return;

    const uint ii = mybuffer[i];
    if(ii >= N){
        return;
    }
    const vec_xyz r_i = r[i].XYZ;
    const float dr = dr_level0[iset[i]] / ilevel[i];


    float w = 0.f;
    m0[ii] = 0.f;
    r[ii] = VEC_ZERO;
    u[ii] = VEC_ZERO;
    dudt[ii] = VEC_ZERO;
    rho[ii] = 0.f;
    drhodt[ii] = 0.f;

    BEGIN_LOOP_OVER_NEIGHS(){
        if(((isplit[j] != 1) && (isplit[j] != 2)) ||  // Not a child
           (iset[i] != iset[j]) ||                    // Different set of particles
           (ilevel[i] != ilevel[j])                   // Different level of refinement
        ){
            j++;
            continue;
        }

        if(all(islessequal(fabs(r[j].XYZ - r_i),
                           0.6f * dr * VEC_ONE.XYZ))
        ){
            const float w_j = split_weight[j];
            w += w_j;
            m0[ii] += w_j * m0[j];
            r[ii] += w_j * r[j];
            u[ii] += w_j * u[j];
            dudt[ii] += w_j * dudt[j];
            rho[ii] += w_j * rho[j];
            drhodt[ii] += w_j * drhodt[j];
        }
    }END_LOOP_OVER_NEIGHS()

    m[ii] = 0.f;
    // m0[ii] /= 1.f;                 // The mass is integrated, not averaged
    r[ii] /= w;
    u[ii] /= w;
    dudt[ii] /= w;
    rho[ii] /= w;
    drhodt[ii] /= w;    
}

/*
 * @}
 */