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

/** @addtogroup cfd
 * @{
 */

/** @file
 * @brief Inflow boundary condition functions.
 */

#include "resources/Scripts/types/types.h"

/** @brief Mark the particles at the inflow side of the boundary. Just fluid
 * particles (imove > 0) are considered
 *
 * @param imove Moving flags.
 *   - imove > 0 for regular fluid particles.
 *   - imove = 0 for sensors.
 *   - imove < 0 for boundary elements/particles.
 * @param r Position \f$ \mathbf{r} \f$.
 * @param iinflow_in 1 if the particle is a fluid particle at the inflow side of
 * the boundary, 0 otherwise.
 * @param N Number of particles.
 * @param inlet_r Lower corner of the inlet square.
 * @param inlet_n Inflow boundary plane normal (outwards fluid domain).
 */
__kernel void marker(const __global int* imove,
                     const __global vec* r,
                     __global uint* iinflow_in,
                     unsigned int N,
                     vec inlet_r,
                     vec inlet_n)
{
    const unsigned int i = get_global_id(0);
    if(i >= N){
        return;
    }
    iinflow_in[i] = 0;
    if(imove[i] <= 0) {
        return;
    }

    if(dot(r[i].XYZ - inlet_r.XYZ, inlet_n.XYZ) >= 0) {
        iinflow_in[i] = 1;
    }
}

/** @brief Sort the inflow particle marks
 *
 * @param iinflow_in Unsorted original mass \f$ m_0 \f$.
 * @param iinflow Sorted original mass \f$ m_0 \f$.
 * @param id_sorted Permutations list from the unsorted space to the sorted
 * one.
 * @param N Number of particles.
 */
__kernel void sort(const __global float *iinflow_in,
                   __global float *iinflow,
                   const __global unit *id_sorted,
                   unsigned int N)
{
    uint i = get_global_id(0);
    if(i >= N)
        return;

    const uint i_out = id_sorted[i];

    iinflow[i_out] = iinflow_in[i];
}

/** @brief Mark the inflow seed particles, i.e. those particles that used to be
 * inflow ones, but now shall become regular fluid particles
 *
 * @param imove Moving flags.
 *   - imove > 0 for regular fluid particles.
 *   - imove = 0 for sensors.
 *   - imove < 0 for boundary elements/particles.
 * @param r Position \f$ \mathbf{r} \f$.
 * @param iinflow 1 if the particle is a fluid particle at the inflow side of
 * the boundary, 0 otherwise.
 * @param iinflowseed 1 if a particle marked as #iinflow moves inside the fluid
 * domain (becoming a regular fluid particle), 0 otherwise.
 * @param N Number of particles.
 * @param inlet_r Lower corner of the inlet square.
 * @param inlet_n Inflow boundary plane normal (outwards fluid domain).
 */
__kernel void seeds(const __global int* imove,
                    const __global vec* r,
                    const __global uint* iinflow,
                    __global uint* iinflowseed,
                    unsigned int N,
                    vec inlet_r,
                    vec inlet_n)
{
    const unsigned int i = get_global_id(0);
    if(i >= N){
        return;
    }
    iinflowseed[i] = 0;
    if(iinflow[i] == 0) {
        return;
    }

    if(dot(r[i].XYZ - inlet_r.XYZ, inlet_n.XYZ) < 0) {
        iinflowseed[i] = 1;
    }
}

/** @brief Particles generation at the inflow boundary condition.
 *
 * Particles are generated just when another particle which used to be at the
 * inflow side of the boundary, suddently moves in the fluid domain, becoming
 * a regular fluid particle. Then the particle is duplicated, moving the new
 * particle in the inflow region a distance of 2*h
 *
 * @param iinflowseed 1 if a particle marked as #iinflow moves inside the fluid
 * domain (becoming a regular fluid particle), 0 otherwise.
 * @param inflowseed_invperm Permutation to find the index of the particle in
 * the list of particles to become duplicated.
 * @param imove Moving flags.
 *   - imove > 0 for regular fluid particles.
 *   - imove = 0 for sensors.
 *   - imove < 0 for boundary elements/particles.
 * @param iset Set of particles index.
 * @param r Position \f$ \mathbf{r} \f$.
 * @param u Velocity \f$ \mathbf{u} \f$.
 * @param dudt Velocity rate of change \f$ \frac{d \mathbf{u}}{d t} \f$.
 * @param rho Density \f$ \rho \f$.
 * @param drhodt Density rate of change \f$ \frac{d \rho}{d t} \f$.
 * @param m Mass \f$ m \f$.
 * @param p Pressure \f$ p \f$.
 * @param N Number of particles.
 * @param nbuffer Number of buffer particles.
 * @param inlet_n = Inflow boundary plane normal (outwards fluid domain).
 */
__kernel void feed(const __global uint* iinflowseed,
                   const __global unsigned int* inflowseed_invperm,
                   __global int* imove,
                   __global vec* r,
                   __global vec* u,
                   __global vec* dudt,
                   __global float* rho,
                   __global float* drhodt,
                   __global float* m,
                   __global float* p,
                   unsigned int N,
                   unsigned int nbuffer,
                   vec inlet_n)
{
    const unsigned int i = get_global_id(0);
    if(i >= N){
        return;
    }

    // Check whether the particle should become duplicated or not
    const unsigned int j = inflowseed_invperm[i];
    if(iinflowseed[j] == 0){
        return;
    }

    // Compute the index of the first buffer particle to steal. Take care, the
    // radix sort is storing the particles to become split at the end of the
    // list
    const unsigned int i0 = N - nbuffer;
    unsigned int ii = i0 + (N - j - 1);
    // Check that there are buffer particles enough
    if(ii >= N){
        // PROBLEMS! This particle cannot be duplicated because we have not
        // buffer particles enough to create the children. That's should be
        // theoretically catched up by "cfd check buffer particles" tool
        return;
    }

    imove[ii] = 1;
    r[ii] = r[i] + 2 * inlet_n;
    u[ii] = u[i];
    dudt[ii] = VEC_ZERO;
    rho[ii] = rho[i];
    drhodt[ii] = 0.f;
    m[ii] = m[i];
    p[ii] = p[i];
}

/** @brief Vanish the velocity and density rates of variation for the dummy
 * particles at the inflow.
 *
 * @param iinflow 1 if the particle is a fluid particle at the inflow side of
 * the boundary, 0 otherwise.
 * @param dudt Velocity rate of change \f$ \frac{d \mathbf{u}}{d t} \f$.
 * @param drhodt Density rate of change \f$ \frac{d \rho}{d t} \f$.
 * @param N Number of particles.
 */
__kernel void rates(const __global uint* iinflow,
                    __global vec* dudt,
                    __global float* drhodt,
                    unsigned int N)
{
    // find position in global arrays
    unsigned int i = get_global_id(0);
    if(i >= N)
        return;
    if(iinflow[i] == 0) {
        return;
    }

    dudt[i] = VEC_ZERO;
    drhodt[i] = 0.f;
}

/*
 * @}
 */
