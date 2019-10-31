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
 *  @brief Improved Euler time integration scheme
 *
 * Time integration is based in the following quasi-second order
 * Predictor-Corrector integration scheme:
 *   - \f$ \mathbf{u}_{n+1} = \mathbf{u}_{n}
     + \frac{\Delta t}{2} \left(
        \left. \frac{\mathrm{d}\mathbf{u}}{\mathrm{d}t} \right\vert_{n+1/2} +
        \left. \frac{\mathrm{d}\mathbf{u}}{\mathrm{d}t} \right\vert_{n-1/2}
     \right)
     \f$
 *   - \f$ \mathbf{r}_{n+1} = \mathbf{r}_{n} + \Delta t \, \mathbf{u}_{n}
     + \frac{\Delta t^2}{4} \left(
        \left. \frac{\mathrm{d}\mathbf{u}}{\mathrm{d}t} \right\vert_{n+1/2} +
        \left. \frac{\mathrm{d}\mathbf{u}}{\mathrm{d}t} \right\vert_{n-1/2}
     \right)
     \f$
 *   - \f$ \rho_{n+1} = \rho_{n}
     + \frac{\Delta t}{2} \left(
        \left. \frac{\mathrm{d}\rho}{\mathrm{d}t} \right\vert_{n+1/2} +
        \left. \frac{\mathrm{d}\rho}{\mathrm{d}t} \right\vert_{n-1/2}
     \right)
     \f$
 */

#include "resources/Scripts/types/types.h"

#ifndef TSCHEME_ADAMS_BASHFORTH_STEPS
    #define TSCHEME_ADAMS_BASHFORTH_STEPS 5u
#endif

/** @brief Adams-Bashforth time integration scheme predictor stage
 * @param imove Moving flags.
 *   - imove > 0 for regular fluid/solid particles.
 *   - imove = 0 for sensors.
 *   - imove < 0 for boundary elements/particles.
 * @param r Position \f$ \mathbf{r}_{n+1} \f$.
 * @param u Velocity \f$ \mathbf{u}_{n+1} \f$.
 * @param dudt Velocity rate of change
 * \f$ \left. \frac{d \mathbf{u}}{d t} \right\vert_{n+1} \f$.
 * @param rho Density \f$ \rho_{n+1} \f$.
 * @param drhodt Density rate of change
 * \f$ \left. \frac{d \rho}{d t} \right\vert_{n+1} \f$.
 * @param r_in Position \f$ \mathbf{r}_{n+1/2} \f$.
 * @param u_in Velocity \f$ \mathbf{u}_{n+1/2} \f$.
 * @param dudt_in Velocity rate of change
 * \f$ \left. \frac{d \mathbf{u}}{d t} \right\vert_{n+1/2} \f$.
 * @param rho_in Density \f$ \rho_{n+1/2} \f$.
 * @param drhodt_in Density rate of change
 * \f$ \left. \frac{d \rho}{d t} \right\vert_{n+1/2} \f$.
 * @param N Number of particles.
 */
__kernel void predictor(__global vec* r,
                        __global vec* u,
                        __global vec* dudt,
                        __global float* rho,
                        __global float* drhodt,
                        __global vec* r_in,
                        __global vec* u_in,
                        __global vec* dudt_in,
                        __global float* rho_in,
                        __global float* drhodt_in,
                        unsigned int N)
{
    unsigned int i = get_global_id(0);
    if(i >= N)
        return;

    dudt_in[i] = dudt[i];
    u_in[i] = u[i];
    r_in[i] = r[i];
    
    drhodt_in[i] = drhodt[i];
    rho_in[i] = rho[i];
}

/** @brief Sort the stored Adams-Bashforth stored variation rates
 *
 * @param id_in Unsorted particle indexes
 * @param id Sorted particle indexes
 * @param iset_in Unsorted set of particles indexes.
 * @param iset Sorted set of particles indexes.
 * @param imove_in Unsorted moving flags.
 *   - imove > 0 for regular fluid particles.
 *   - imove = 0 for sensors.
 *   - imove < 0 for boundary elements/particles.
 * @param imove Sorted moving flags.
 *   - imove > 0 for regular fluid particles.
 *   - imove = 0 for sensors.
 *   - imove < 0 for boundary elements/particles.
 * @param r_in Unsorted position \f$ \mathbf{r} \f$.
 * @param r Sorted position \f$ \mathbf{r} \f$.
 * @param normal_in Unsorted normal \f$ \mathbf{n} \f$.
 * @param normal Sorted normal \f$ \mathbf{n} \f$.
 * @param u_in Unsorted velocity \f$ \mathbf{u} \f$.
 * @param u Sorted velocity \f$ \mathbf{u} \f$.
 * @param id_sorted Permutations list from the unsorted space to the sorted
 * one.
 * @param N Number of particles.
 */
__kernel void sort(const __global vec *dudt_as1_in, __global vec *dudt_as1,
                   const __global vec *dudt_as2_in, __global vec *dudt_as2,
                   const __global vec *dudt_as3_in, __global vec *dudt_as3,
                   const __global vec *dudt_as4_in, __global vec *dudt_as4,
                   const __global float *drhodt_as1_in, __global float *drhodt_as1,
                   const __global float *drhodt_as2_in, __global float *drhodt_as2,
                   const __global float *drhodt_as3_in, __global float *drhodt_as3,
                   const __global float *drhodt_as4_in, __global float *drhodt_as4,
                   const __global unit *id_sorted,
                   unsigned int N)
{
    uint i = get_global_id(0);
    if(i >= N)
        return;

    const uint i_out = id_sorted[i];

    dudt_as1[i_out] = dudt_as1_in[i];
    dudt_as2[i_out] = dudt_as2_in[i];
    dudt_as3[i_out] = dudt_as3_in[i];
    dudt_as4[i_out] = dudt_as4_in[i];
    drhodt_as1[i_out] = drhodt_as1_in[i];
    drhodt_as2[i_out] = drhodt_as2_in[i];
    drhodt_as3[i_out] = drhodt_as3_in[i];
    drhodt_as4[i_out] = drhodt_as4_in[i];
}


#define DYDT_1(dydt_0) dydt_0
#define DYDT_2(dydt_0, dydt_1) 1.5f * dydt_0 - 0.5f * dydt_1
#define DYDT_3(dydt_0, dydt_1, dydt_2) \
    23.f / 12.f * dydt_0 - 4.f / 3.f * dydt_1 + 5.f / 12.f * dydt_2
#define DYDT_4(dydt_0, dydt_1, dydt_2, dydt_3) \
    55.f / 24.f * dydt_0 - 59.f / 24.f * dydt_1 + 37.f / 24.f * dydt_2 - \
    3.f / 8.f * dydt_3
#define DYDT_5(dydt_0, dydt_1, dydt_2, dydt_3, dydt_4) \
    1901.f / 720.f * dydt_0 - 1387.f / 360.f * dydt_1 + 109.f / 30.f * dydt_2 - \
    637.f / 360.f * dydt_3 + 251.f / 720.f * dydt_4


/** @brief Improved Euler time integration scheme corrector stage
 * @param imove Moving flags.
 *   - imove > 0 for regular fluid particles.
 *   - imove = 0 for sensors.
 *   - imove < 0 for boundary elements/particles.
 * @param iset Set of particles index.
 * @param r Position \f$ \mathbf{r}_{n+1/2} \f$.
 * @param u Velocity \f$ \mathbf{u}_{n+1/2} \f$.
 * @param dudt Velocity rate of change
 * \f$ \left. \frac{d \mathbf{u}}{d t} \right\vert_{n+1/2} \f$.
 * @param rho Density \f$ \rho_{n+1/2} \f$.
 * @param drhodt Density rate of change
 * \f$ \left. \frac{d \rho}{d t} \right\vert_{n+1/2} \f$.
 * @param dudt_as1 Velocity rate of change of previous step
 * \f$ \left. \frac{d \mathbf{u}}{d t} \right\vert_{n-1/2} \f$.
 * @param drhodt_as1 Density rate of change of previous step
 * \f$ \left. \frac{d \rho}{d t} \right\vert_{n-1/2} \f$.
 * @param N Number of particles.
 * @param dt Time step \f$ \Delta t \f$.
 * @param iter Current iteration.
 */
__kernel void corrector(__global int* imove,
                        __global unsigned int* iset,
                        __global vec* r,
                        __global vec* u,
                        __global vec* dudt,
                        __global float* rho,
                        __global float* drhodt,
                        __global vec* dudt_as1,
                        __global float* drhodt_as1,
                        __global vec* dudt_as2,
                        __global float* drhodt_as2,
                        __global vec* dudt_as3,
                        __global float* drhodt_as3,
                        __global vec* dudt_as4,
                        __global float* drhodt_as4,
                        unsigned int N,
                        float dt,
                        unsigned int iter)
{
    unsigned int i = get_global_id(0);
    if(i >= N)
        return;

    if(imove[i] > 0) {
        const uint local_iter = min(iter, TSCHEME_ADAMS_BASHFORTH_STEPS);
        vec dudt_as;
        float drhodt_as;
        if(local_iter < 1) {
            dudt_as = DYDT_1(dudt[i]);
            drhodt_as = DYDT_1(drhodt[i]);
        } else if(local_iter < 2) {
            dudt_as = DYDT_2(dudt[i], dudt_as1[i]);
            drhodt_as = DYDT_2(drhodt[i], drhodt_as1[i]);
        } else if(local_iter < 3) {
            dudt_as = DYDT_3(dudt[i], dudt_as1[i], dudt_as2[i]);
            drhodt_as = DYDT_3(drhodt[i], drhodt_as1[i], drhodt_as2[i]);
        } else if(local_iter < 4) {
            dudt_as = DYDT_4(dudt[i], dudt_as1[i], dudt_as2[i], dudt_as3[i]);
            drhodt_as = DYDT_4(drhodt[i], drhodt_as1[i], drhodt_as2[i],
                               drhodt_as3[i]);
        } else {
            dudt_as = DYDT_5(dudt[i], dudt_as1[i], dudt_as2[i], dudt_as3[i],
                             dudt_as4[i]);
            drhodt_as = DYDT_5(drhodt[i], drhodt_as1[i], drhodt_as2[i],
                               drhodt_as3[i], drhodt_as4[i]);
        }
        r[i] += dt * u[i] + 0.5f * dt * dt * dudt_as;
        u[i] += dt * dudt_as;
        rho[i] += dt * drhodt_as;
    }
}

/** @brief Backup de data for the sorting algorithm
 * @param dudt_as1 Velocity rate of change of previous step
 * \f$ \left. \frac{d \mathbf{u}}{d t} \right\vert_{n-1/2} \f$.
 * @param drhodt_as1 Density rate of change of previous step
 * \f$ \left. \frac{d \rho}{d t} \right\vert_{n-1/2} \f$.
 * @param dudt_as2 Velocity rate of change of 2 steps ago
 * \f$ \left. \frac{d \mathbf{u}}{d t} \right\vert_{n-1/2} \f$.
 * @param drhodt_as2 Density rate of change of 2 steps ago
 * \f$ \left. \frac{d \rho}{d t} \right\vert_{n-1/2} \f$.
 * @param dudt_as3 Velocity rate of change of 3 steps ago
 * \f$ \left. \frac{d \mathbf{u}}{d t} \right\vert_{n-1/2} \f$.
 * @param drhodt_as3 Density rate of change of 3 steps ago
 * \f$ \left. \frac{d \rho}{d t} \right\vert_{n-1/2} \f$.
 * @param dudt Velocity rate of change
 * \f$ \left. \frac{d \mathbf{u}}{d t} \right\vert_{n+1/2} \f$.
 * @param drhodt Density rate of change
 * \f$ \left. \frac{d \rho}{d t} \right\vert_{n+1/2} \f$.
 * @param dudt_as1_in Velocity rate of change of previous step
 * \f$ \left. \frac{d \mathbf{u}}{d t} \right\vert_{n-1/2} \f$.
 * @param drhodt_as1_in Density rate of change of previous step
 * \f$ \left. \frac{d \rho}{d t} \right\vert_{n-1/2} \f$.
 * @param dudt_as2_in Velocity rate of change of 2 steps ago
 * \f$ \left. \frac{d \mathbf{u}}{d t} \right\vert_{n-1/2} \f$.
 * @param drhodt_as2_in Density rate of change of 2 steps ago
 * \f$ \left. \frac{d \rho}{d t} \right\vert_{n-1/2} \f$.
 * @param dudt_as3_in Velocity rate of change of 3 steps ago
 * \f$ \left. \frac{d \mathbf{u}}{d t} \right\vert_{n-1/2} \f$.
 * @param drhodt_as3_in Density rate of change of 3 steps ago
 * \f$ \left. \frac{d \rho}{d t} \right\vert_{n-1/2} \f$.
 * @param dudt_as4_in Velocity rate of change of 4 steps ago
 * \f$ \left. \frac{d \mathbf{u}}{d t} \right\vert_{n-1/2} \f$.
 * @param drhodt_as4_in Density rate of change of 4 steps ago
 * \f$ \left. \frac{d \rho}{d t} \right\vert_{n-1/2} \f$.
 * @param N Number of particles.
 */
__kernel void postcorrector(const __global vec* dudt_as1,
                            const __global float* drhodt_as1,
                            const __global vec* dudt_as2,
                            const __global float* drhodt_as2,
                            const __global vec* dudt_as3,
                            const __global float* drhodt_as3,
                            const __global vec* dudt,
                            const __global float* drhodt,
                            __global vec* dudt_as1_in,
                            __global float* drhodt_as1_in,
                            __global vec* dudt_as2_in,
                            __global float* drhodt_as2_in,
                            __global vec* dudt_as3_in,
                            __global float* drhodt_as3_in,
                            __global vec* dudt_as4_in,
                            __global float* drhodt_as4_in,
                            unsigned int N)
{
    unsigned int i = get_global_id(0);
    if(i >= N)
        return;

    dudt_as1_in[i] = dudt[i];
    dudt_as2_in[i] = dudt_as1[i];
    dudt_as3_in[i] = dudt_as2[i];
    dudt_as4_in[i] = dudt_as3[i];
    drhodt_as1_in[i] = drhodt[i];
    drhodt_as2_in[i] = drhodt_as1[i];
    drhodt_as3_in[i] = drhodt_as2[i];
    drhodt_as4_in[i] = drhodt_as3[i];
}

/*
 * @}
 */
 
