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
 *  @brief 1st order Euler time integration scheme
 *
 * Time integration is based in the following 1st order integration scheme:
 *   - \f$ \mathbf{u}_{n+1} = \mathbf{u}_{n} + \Delta t
        \left. \frac{\mathrm{d}\mathbf{u}}{\mathrm{d}t} \right\vert_{n}
     \f$
 *   - \f$ \mathbf{r}_{n+1} = \mathbf{r}_{n} + \Delta t \, \mathbf{u}_{n}
     + \frac{\Delta t^2}{2}
        \left. \frac{\mathrm{d}\mathbf{u}}{\mathrm{d}t} \right\vert_{n}
     \f$
 *   - \f$ \rho_{n+1} = \rho_{n} + \Delta t
        \left. \frac{\mathrm{d}\rho}{\mathrm{d}t} \right\vert_{n}
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

/** @brief Euler time integration scheme corrector stage
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
 * @param r_in Position \f$ \mathbf{r}_{n} \f$.
 * @param u_in Velocity \f$ \mathbf{u}_{n} \f$.
 * @param dudt_in Velocity rate of change
 * \f$ \left. \frac{d \mathbf{u}}{d t} \right\vert_{n-1/2} \f$.
 * @param rho_in Density \f$ \rho_{n} \f$.
 * @param drhodt_in Density rate of change
 * \f$ \left. \frac{d \rho}{d t} \right\vert_{n-1/2} \f$.
 * @param N Number of particles.
 * @param dt Time step \f$ \Delta t \f$.
 */
__kernel void corrector(__global int* imove,
                        __global unsigned int* iset,
                        __global vec* r,
                        __global vec* u,
                        __global vec* dudt,
                        __global float* rho,
                        __global float* drhodt,
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

    if(imove[i] > 0) {
        r[i] += dt * u[i] + 0.5f * dt * dt * dudt[i];
        u[i] += dt * dudt[i];
        rho[i] += dt * drhodt[i];
    }

    r_in[i] = r[i];
    u_in[i] = u[i];
    rho_in[i] = rho[i];
    dudt_in[i] = dudt[i];
    drhodt_in[i] = drhodt[i];
}

/*
 * @}
 */
 
