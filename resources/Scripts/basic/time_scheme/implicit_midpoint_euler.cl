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

/** @defgroup basic Basic preset
 *
 * @brief Basic preset of tools to build more complex sets of tools later
 * 
 * @{
 */

/** @file
 *  @brief Implicit Midpoint Euler time integration scheme
 */

#include "resources/Scripts/types/types.h"

/** @brief Implicit Midpoint Euler time integration scheme predictor stage
 * @param imove Moving flags.
 *   - imove > 0 for regular fluid/solid particles.
 *   - imove = 0 for sensors.
 *   - imove < 0 for boundary elements/particles.
 * @param r Position \f$ \mathbf{r}_{n} \f$.
 * @param u_0 Velocity \f$ \mathbf{u}_{n} \f$.
 * @param dudt Velocity rate of change
 * \f$ \left. \frac{d \mathbf{u}}{d t} \right\vert_{n+1/2} \f$.
 * @param rho_0 Density \f$ \rho_{n} \f$.
 * @param drhodt Density rate of change
 * \f$ \left. \frac{d \rho}{d t} \right\vert_{n+1/2} \f$.
 * @param r_in Position backup \f$ \mathbf{r}_{n} \f$.
 * @param u_in Midpoint velocity \f$ \mathbf{u}_{n+1/2} \f$.
 * @param dudt_in Velocity rate of change backup
 * \f$ \left. \frac{d \mathbf{u}}{d t} \right\vert_{n+1/2} \f$.
 * @param rho_in Midpoint density \f$ \rho_{n+1/2} \f$.
 * @param drhodt_in Density rate of change backup
 * \f$ \left. \frac{d \rho}{d t} \right\vert_{n+1/2} \f$.
 * @param N Number of particles.
 * @param dt Time step \f$ \Delta t \f$.
 */
__kernel void predictor(const __global int* imove,
                        const __global vec* r,
                        const __global vec* u_0,
                        const __global vec* dudt,
                        const __global float* rho_0,
                        const __global float* drhodt,
                        __global vec* r_in,
                        __global vec* u_in,
                        __global vec* dudt_in,
                        __global float* rho_in,
                        __global float* drhodt_in,
                        unsigned int N,
                        float dt)
{
    unsigned int i = get_global_id(0);
    if(i >= N)
        return;

    // We want to move to the mid point of the time step
    float hdt = 0.f;
    if(imove[i] > 0)
        hdt = 0.5f * dt;

    r_in[i] = r[i];
    u_in[i] = u_0[i] + hdt * dudt[i];
    rho_in[i] = rho_0[i] + hdt * drhodt[i];
    
    dudt_in[i] = dudt[i];
    drhodt_in[i] = drhodt[i];
}

/** @brief Implicit Midpoint Euler time integration scheme variables sorting
 * @param u_0 Sorted velocity \f$ \mathbf{u}_{n} \f$.
 * @param rho_0 Sorted density \f$ \rho_{n} \f$.
 * @param u_0_in Unsorted velocity \f$ \mathbf{u}_{n} \f$.
 * @param rho_0_in Unsorted density \f$ \rho_{n+1/2} \f$.
 * @param id_sorted Permutations list from the unsorted space to the sorted
 * one.
 * @param N Number of particles.
 */
__kernel void sort(const __global float *rho_0_in, __global float *rho_0,
                   const __global vec *u_0_in, __global vec *u_0,
                   const __global unit *id_sorted,
                   unsigned int N)
{
    uint i = get_global_id(0);
    if(i >= N)
        return;

    const uint i_out = id_sorted[i];

    rho_0[i_out] = rho_0_in[i];
    u_0[i_out] = u_0_in[i];
}

/** @brief Relax the variation rates modification at the fixed point iterative
 * solution.
 * @param imove Moving flags.
 *   - imove > 0 for regular fluid/solid particles.
 *   - imove = 0 for sensors.
 *   - imove < 0 for boundary elements/particles.
 * @param dudt_in Velocity rate of change from the previous subiteration
 * \f$ \left. \frac{d \mathbf{u}}{d t} \right\vert_{n+1/2} \f$.
 * @param drhodt_in Density rate of change from the previous subiteration
 * \f$ \left. \frac{d \rho}{d t} \right\vert_{n+1/2} \f$.
 * @param dudt New velocity rate of change, to be relaxed
 * \f$ \left. \frac{d \mathbf{u}}{d t} \right\vert_{n+1/2} \f$.
 * @param drhodt New density rate of change, to be relaxed
 * \f$ \left. \frac{d \rho}{d t} \right\vert_{n+1/2} \f$.
 * @param N Number of particles.
 */
__kernel void corrector(const __global int* imove,
                        const __global vec* dudt_in,
                        const __global float* drhodt_in,
                        __global vec* dudt,
                        __global float* drhodt,
                        unsigned int N,
                        float subiter_relax)
{
    unsigned int i = get_global_id(0);
    if(i >= N)
        return;
    if(imove[i] <= 0)
        return;

    float f = subiter_relax;
    const vec ddudt = dudt[i] - dudt_in[i];
    const float ddrhodt = drhodt[i] - drhodt_in[i];

    dudt[i] = dudt_in[i] + f * ddudt;
    drhodt[i] = drhodt_in[i] + f * ddrhodt;
}

/** @brief Integrate the system to the next time step, when the implicit
 * subiterator already finished his job.
 * @param imove Moving flags.
 *   - imove > 0 for regular fluid/solid particles.
 *   - imove = 0 for sensors.
 *   - imove < 0 for boundary elements/particles.
 * @param r_in Position \f$ \mathbf{r}_{n} \f$.
 * @param u_0 Velocity \f$ \mathbf{u}_{n} \f$.
 * @param dudt Velocity rate of change
 * \f$ \left. \frac{d \mathbf{u}}{d t} \right\vert_{n+1/2} \f$.
 * @param rho_0 Density \f$ \rho_{n} \f$.
 * @param drhodt Density rate of change
 * \f$ \left. \frac{d \rho}{d t} \right\vert_{n+1/2} \f$.
 * @param r Position \f$ \mathbf{r}_{n+1} \f$.
 * @param u Velocity \f$ \mathbf{u}_{n+1} \f$.
 * @param rho Density \f$ \rho_{n+1} \f$.
 * @param N Number of particles.
 * @param dt Time step \f$ \Delta t \f$.
 */
__kernel void integrate(const __global int* imove,
                        const __global vec* r_in,
                        const __global vec* u_0,
                        const __global vec* dudt,
                        const __global float* rho_0,
                        const __global float* drhodt,
                        __global vec* r,
                        __global vec* u,
                        __global float* rho,
                        unsigned int N,
                        float dt)
{
    unsigned int i = get_global_id(0);
    if(i >= N)
        return;
    if(imove[i] <= 0)
        return;

    r[i] += dt * u_0[i] + 0.5f * dt * dt * dudt[i];
    u[i] = u_0[i] + dt * dudt[i];
    rho[i] = rho_0[i] + dt * drhodt[i];
}

/*
 * @}
 */
 
