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

#ifndef HAVE_3D
    #define N_DAUGHTER 4u
    #define DAUGHTER_ID_TYPE uint4
    #define DAUGHTER_DIST_TYPE float4
#else
    #define N_DAUGHTER 8u
    #define DAUGHTER_ID_TYPE uint8
    #define DAUGHTER_DIST_TYPE float8
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
    split_cell[i] = CONVERT(ivec, r[i] / dr);
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
 * @param isplit 0 if the particle should not become coalesced, 1 for the
 * coalescing particles, 2 for the seeds.
 * @param split_cell Seed cell where the particle is located.
 * @param split_dist Distance to the center of the cell. Just the closest
 * particle to the cell center will be kept as seed
 * @param icell Cell where each particle is located.
 * @param ihoc Head of chain for each cell (first particle found).
 * @param N Number of particles.
 * @param n_cells Number of cells in each direction
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
            const float jsplit_dist = split_dist[j];
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
 * @param isplit 0 if the particle should not become split, 1 otherwise.
 * @param split_invperm Permutation to find the index of the particle in the
 * list of particles to become split.
 * @param mybuffer Index of the partner buffer particle.
 * @param ilevel Current refinement level of the particle.
 * @param level Target refinement level of the particle.
 * @param miter Mass transfer iteration (Positive for growing particles,
 * negative for shrinking particles).
 * @param N Number of particles.
 * @param nbuffer Number of available buffer particles.
 */
__kernel void generate(__global int* imove,
                       __global int* iset,
                       __global unsigned int* isplit,
                       __global unsigned int* split_invperm,
                       __global unsigned int* mybuffer,
                       __global unsigned int* ilevel,
                       __global unsigned int* level,
                       __global int* miter,
                       unsigned int N,
                       unsigned int nbuffer)
{
    unsigned int i = get_global_id(0);
    if(i >= N)
        return;

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
        mybuffer[i] = N;
        return;
    }

    // Associate the buffer particle to the seed
    mybuffer[i] = ii;

    // Set the already known variables
    imove[ii] = imove[i];
    iset[ii] = iset[i];                
    ilevel[ii] = ilevel[i] - 1;
    level[ii] = ilevel[i] - 1;
    miter[ii] = 1;
}

/** @brief Look for all the candidates to become children, close enough to the
 * seed.
 *
 * When a particle is not already splitting or coalescing, is close enough to
 * the seed, has the same refinement level, and belongs to the same particles
 * set, it is a children candidates.
 *
 * A particle marked as children candidate is invariably coalescing. However,
 * we don't know yet which particle will be its partner. That's why it is just a
 * candidate. Later we are looking all the neighbours of the children candidates
 * to determine the best seed.
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
__kernel void children_candidates(__global const int* imove,
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
           (iset[i] != iset[j]) ||    // Another set of particles
           (ilevel[i] != ilevel[j])   // Another level of refinement
        ){
            j++;
            continue;
        }

        if(all(islessequal(fabs(r[j].XYZ - r_i),
                           0.75f * dr * VEC_ONE.XYZ))
        ){
            // It does not matter if several threads write this data at the same
            // time, because are constant values
            isplit[j] = 1;
            miter[j] = -1;
        }
    }END_LOOP_OVER_NEIGHS()
}

/** @brief Associate each children to the closest seed.
 *
 * We are actually not associating it to the seed, but with the buffer partner
 * particle.
 *
 * @param iset Set of particles index.
 * @param ilevel Current refinement level of the particle.
 * @param r Position \f$ \mathbf{r} \f$.
 * @param isplit 0 if the particle should not become coalesced, 1 for the
 * coalescing particles, 2 for the seeds.
 * @param miter Mass transfer iteration (Positive for growing particles,
 * negative for shrinking particles).
 * @param mybuffer Index of the partner buffer particle associated.
 * @param icell Cell where each particle is located.
 * @param ihoc Head of chain for each cell (first particle found).
 * @param dr_level0 Theoretical distance between particles at the lowest
 * refinement level.
 * @param N Number of particles.
 * @param n_cells Number of cells in each direction
 */
__kernel void children(__global const unsigned int* iset,
                       __global const unsigned int* ilevel,
                       __global const vec* r,
                       __global unsigned int* isplit,
                       __global int* miter,
                       __global unsigned int* mybuffer,
                       __global const uint *icell,
                       __global const uint *ihoc,
                       __constant float* dr_level0,
                       unsigned int N,
                       uivec4 n_cells)
{
    unsigned int i = get_global_id(0);
    if(i >= N)
        return;
    if((isplit[i] != 1) || (mybuffer[i] != N))
        return;

    // Let's set the particle as not eventually coalescing anymore (e.g. it
    // cannot found a non-complete seed, see count_children())
    isplit[i] = 0;
    miter[i] = M_ITERS + 1;

    const float dr = dr_level0[iset[i]] / ilevel[i];
    float dist = MAXFLOAT;
    const vec_xyz r_i = r[i].XYZ;

    BEGIN_LOOP_OVER_NEIGHS(){
        if((isplit[j] != 2) ||             // Not a seed/the seed is complete
           (iset[i] != iset[j])            // Another set of particles
        ){
            j++;
            continue;
        }

        const vec_xyz r_ij = r[j].XYZ - r_i;
        if(any(isgreater(fabs(r_ij), 0.75f * dr * VEC_ONE.XYZ))){
            // Too far from him
            j++;
            continue;
        }
        {
            const float l = length(r_ij);
            if(l < dist){
                // This is a better seed candidate
                dist = l;
                mybuffer[i] = mybuffer[j];
                isplit[i] = 1;
                miter[i] = -1;
            }
        }
    }END_LOOP_OVER_NEIGHS()
}

#ifdef LOCAL_MEM_SIZE
    #define MEM_ADDRESS __local
#else
    #define MEM_ADDRESS __private
#endif

/** @brief Check if a child can be inserted into the list.
 *
 * If the child is closer than any other particle in the list, it will be
 * inserted into the list, pushing out the farthest one.
 *
 * @param dist Distance between the particle and the seed.
 * @param id ID of the child candidate
 * @param dists List of distances between selected children and seed
 * @param ids List of selected children ids
 * @return true if the particle has been inserted, false otherwise.
 */
bool check_child(float dist,
                 uint id,
                 MEM_ADDRESS float *dists,
                 MEM_ADDRESS uint *ids)
{
    unsigned int i, j;
    for(i = 0; i < N_DAUGHTER; i++){
        if(dist >= dists[i]){
            continue;
        }
        // Push out the last particle (the farthest one)
        for(j = N_DAUGHTER - 1; j > i; j--){
            dists[j] = dists[j - 1];
            ids[j] = ids[j - 1];
        }
        // And insert the new particle
        dists[i] = dist;
        ids[i] = id;
        return true;
    }
    return false;
}


/** @brief Collect the children to check if the seed is complete.
 *
 * In order to avoid that a seed may generate a partner particle with too much
 * children, which may create noisy fields, we are counting the number of
 * children associated to the same partner, keeping just the N_DAUGHTER closest
 * ones.
 *
 * The drop children will be set with a mybuffer
 *
 * @param r Position \f$ \mathbf{r} \f$.
 * @param isplit 0 if the particle should not become coalesced, 1 for the
 * coalescing particles, 2 for the seeds.
 * @param mybuffer Index of the partner buffer particle associated.
 * @param icell Cell where each particle is located.
 * @param ihoc Head of chain for each cell (first particle found).
 * @param N Number of particles.
 * @param n_cells Number of cells in each direction
 */
__kernel void count_children(__global const vec* r,
                             __global unsigned int* isplit,
                             __global unsigned int* mybuffer,
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

    const vec_xyz r_i = r[i].XYZ;

    const unsigned int ii = mybuffer[i];
    if(ii == N){
        // A problematic particle which has not found available buffer...
        return;
    }
    unsigned int n_children = 1;
    // Create children id and distance lists
    MEM_ADDRESS uint id_children[N_DAUGHTER];
    MEM_ADDRESS float dist_children[N_DAUGHTER];
    for(uint j = 0; j < N_DAUGHTER; j++){
        id_children[j] = N;
        dist_children[j] = MAXFLOAT;
    }
    // Add the seed as a child
    id_children[0] = i;
    dist_children[0] = 0.f;

    BEGIN_LOOP_OVER_NEIGHS(){
        if((isplit[j] != 1) ||    // Not a child
           (mybuffer[j] == N) ||  // Not partner associated?
           (mybuffer[j] != ii)    // Not the same partner associated
        ){
            j++;
            continue;
        }

        {
            // Set the particle as non-selected. Later we are checking if it is
            // actually present in the list of the selected ones
            mybuffer[j] = N;

            const float dist = length(r[j].XYZ - r_i);
            if(check_child(dist, j, dist_children, id_children))
            {
                n_children++;
            }
        }
    }END_LOOP_OVER_NEIGHS()

    n_children = min(n_children, N_DAUGHTER);
    if(n_children == N_DAUGHTER){
        // Set the seed as complete
        isplit[i] = 3;
    }

    // Reselect the children (all the children candidates out of id_children
    // list are therefore kept as unselected, such that they should look for
    // another partner, or eventually give up of trying to coalesce)
    for(uint j = 0; j < n_children; j++){
        mybuffer[id_children[j]] = ii;
    }
}

/** @brief Collect the children, and the seed itself, in order to compute the
 * field values of the buffer partner particle.
 *
 * Now we have a ssociated the children to the same partner than the seed, we
 * can make a neighbours loop on the seed, integrating the field values to
 * get the average value at the end of the kernel.
 *
 * @param isplit 0 if the particle should not become coalesced, 1 for the
 * coalescing particles, 2 for the seeds.
 * @param mybuffer Index of the partner buffer particle.
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
__kernel void fields(__global const unsigned int* isplit,
                     __global const unsigned int* mybuffer,
                     __global float* m0,
                     __global float* m,
                     __global vec* r,
                     __global vec* u,
                     __global vec* dudt,
                     __global float* rho,
                     __global float* drhodt,
                     __global const uint *icell,
                     __global const uint *ihoc,
                     unsigned int N,
                     uivec4 n_cells)
{
    unsigned int i = get_global_id(0);
    if(i >= N)
        return;
    if((isplit[i] != 2) && (isplit[i] != 3))
        return;

    const unsigned int ii = mybuffer[i];
    if(ii == N){
        // A problematic particle which has not found available buffer...
        return;
    }
    unsigned int n_children = 1;
    m0[ii] = m0[i];
    r[ii] = r[i];
    u[ii] = u[i];
    dudt[ii] = dudt[i];
    rho[ii] = rho[i];
    drhodt[ii] = drhodt[i];

    BEGIN_LOOP_OVER_NEIGHS(){
        if((isplit[j] != 1) ||   // Not a child
           (mybuffer[j] != ii)   // Not the same partner
        ){
            j++;
            continue;
        }

        {
            n_children++;
            m0[ii] += m0[j];
            r[ii] += r[j];
            u[ii] += u[j];
            dudt[ii] += dudt[j];
            rho[ii] += rho[j];
            drhodt[ii] += drhodt[j];
        }
    }END_LOOP_OVER_NEIGHS()

    m[ii] = m0[ii];            // The mass is integrated, not averaged
    r[ii] /= n_children;
    u[ii] /= n_children;
    dudt[ii] /= n_children;
    rho[ii] /= n_children;
    drhodt[ii] /= n_children;    
}

/*
 * @}
 */