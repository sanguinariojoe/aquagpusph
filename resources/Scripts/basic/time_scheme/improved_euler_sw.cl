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

/** @brief Improved Euler time integration scheme predictor stage
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
 * @param dt Time step \f$ \Delta t \f$.
 */
__kernel void predictor(const __global int* imove,
                        const __global vec* r,
                        const __global vec* u,
                        const __global vec* dudt,
                        const __global float* rho,
                        const __global float* drhodt, 
                        const __global float* eee,
                        const __global float* dedt,
                        __global vec* r_in,
                        __global vec* u_in,
                        __global vec* dudt_in,
                        __global float* rho_in,
                        __global float* drhodt_in,   
                        __global float* e_in,
                        __global float* dedt_in,
                        const usize N,
                        const float dt)
{
    usize i = get_global_id(0);
    if(i >= N)
        return;

    float DT = dt;
    if(imove[i] <= 0)
        DT = 0.f;

    dudt_in[i] = dudt[i];
    u_in[i] = u[i] + DT * dudt[i];
    r_in[i] = r[i] + DT * u[i] + 0.5f * DT * DT * dudt[i];
    
    drhodt_in[i] = drhodt[i];
    rho_in[i] = rho[i] + DT * drhodt[i];
    
    dedt_in[i] = dedt[i];
    e_in[i] = eee[i] + DT * dedt[i];
}

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
 * @param dudt_in Velocity rate of change
 * \f$ \left. \frac{d \mathbf{u}}{d t} \right\vert_{n-1/2} \f$.
 * @param drhodt_in Density rate of change
 * \f$ \left. \frac{d \rho}{d t} \right\vert_{n-1/2} \f$.
 * @param N Number of particles.
 * @param dt Time step \f$ \Delta t \f$.
 */
__kernel void corrector(const __global int* imove,
                        const __global unsigned int* iset,
                        __global vec* r,
                        __global vec* u,
                        const __global vec* dudt,
                        __global float* rho,
                        const __global float* drhodt,
                        __global float* eee,
                        const __global float* dedt,
                        const __global vec* dudt_in,
                        const __global float* drhodt_in,
                        const __global float* dedt_in,
                        const usize N,
                        const float dt)
{
    usize i = get_global_id(0);
    if(i >= N)
        return;

    if(imove[i] > 0) {
        const float DT = 0.5f * dt;
        u[i] += DT * (dudt[i] - dudt_in[i]);
        r[i] += DT * DT * (dudt[i] - dudt_in[i]);
        rho[i] += DT * (drhodt[i] - drhodt_in[i]);
        eee[i] += DT * (dedt[i] - dedt_in[i]);
    }
}

/*
 * @}
 */
 
