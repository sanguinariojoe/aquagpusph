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
 *  @brief Semi-implicit Midpoint Euler time integration scheme
 */

#include "resources/Scripts/types/types.h"

/** @brief Semi-implicit Midpoint Euler time integration scheme predictor stage
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

/** @brief Advance to the tims step midpoint
 * @param imove Moving flags.
 *   - imove > 0 for regular fluid particles.
 *   - imove = 0 for sensors.
 *   - imove < 0 for boundary elements/particles.
 * @param u_in Velocity \f$ \mathbf{u}_{n} \f$.
 * @param u Velocity \f$ \mathbf{u}_{n+1/2} \f$.
 * @param dudt Velocity rate of change
 * \f$ \left. \frac{d \mathbf{u}}{d t} \right\vert_{n+1/2} \f$.
 * @param rho_in Density \f$ \rho_{n} \f$.
 * @param rho Density \f$ \rho_{n+1/2} \f$.
 * @param drhodt Density rate of change
 * \f$ \left. \frac{d \rho}{d t} \right\vert_{n+1/2} \f$.
 * @param N Number of particles.
 * @param dt Time step \f$ \Delta t \f$.
 */
__kernel void midpoint(__global int* imove,
                       __global vec* u_in,
                       __global vec* u,
                       __global vec* dudt,
                       __global float* rho_in,
                       __global float* rho,
                       __global float* drhodt,
                       unsigned int N,
                       float dt)
{
    unsigned int i = get_global_id(0);
    if(i >= N)
        return;
    if(imove[i] <= 0)
        return;

    u[i] = u_in[i] + 0.5f * dt * dudt[i];
    rho[i] = rho_in[i] + 0.5f * dt * drhodt[i];
}

__kernel void relax(const __global int* imove,
                    __global vec* dudt_in,
                    __global vec* dudt,
                    __global float* drhodt_in,
                    __global float* drhodt,
                    unsigned int N,
                    float relax_midpoint)
{
    unsigned int i = get_global_id(0);
    if(i >= N)
        return;
    if(imove[i] <= 0)
        return;

    dudt[i] = (1.f - relax_midpoint) * dudt_in[i] + relax_midpoint * dudt[i];
    drhodt[i] = (1.f - relax_midpoint) * drhodt_in[i] + relax_midpoint * drhodt[i];
}

__kernel void residuals(const __global int* imove,
                        const __global float* m,
                        const __global vec* u,
                        const __global vec* dudt_in,
                        const __global vec* dudt,
                        const __global float* rho,
                        const __global float* p,
                        const __global float* drhodt_in,
                        const __global float* drhodt,
                        __global float* residual_midpoint,
                        unsigned int N)
{
    unsigned int i = get_global_id(0);
    if(i >= N)
        return;
    if(imove[i] <= 0) {
        residual_midpoint[i] = 0.f;
        return;
    }

    const float rho2 = rho[i] * rho[i];
    residual_midpoint[i] = m[i] * (
        fabs(dot(u[i], dudt[i] - dudt_in[i])) +
        fabs(p[i] / rho2 * (drhodt[i] - drhodt_in[i]))
    );
}


/** @brief 1st order Semi-implicit Euler time integration scheme corrector
 * stage
 * @param imove Moving flags.
 *   - imove > 0 for regular fluid particles.
 *   - imove = 0 for sensors.
 *   - imove < 0 for boundary elements/particles.
 * @param r Position \f$ \mathbf{r}_{n+1/2} \f$.
 * @param u_in Velocity \f$ \mathbf{u}_{n} \f$.
 * @param u Velocity \f$ \mathbf{u}_{n+1/2} \f$.
 * @param dudt Velocity rate of change
 * \f$ \left. \frac{d \mathbf{u}}{d t} \right\vert_{n+1/2} \f$.
 * @param rho_in Density \f$ \rho_{n} \f$.
 * @param rho Density \f$ \rho_{n+1/2} \f$.
 * @param drhodt Density rate of change
 * \f$ \left. \frac{d \rho}{d t} \right\vert_{n+1/2} \f$.
 * @param N Number of particles.
 * @param dt Time step \f$ \Delta t \f$.
 */
__kernel void corrector(__global int* imove,
                        __global vec* r,
                        __global vec* u_in,
                        __global vec* u,
                        __global vec* dudt,
                        __global float* rho_in,
                        __global float* rho,
                        __global float* drhodt,
                        unsigned int N,
                        float dt)
{
    unsigned int i = get_global_id(0);
    if(i >= N)
        return;
    if(imove[i] <= 0)
        return;

    r[i] += dt * u_in[i] + 0.5f * dt * dt * dudt[i];
    u[i] = u_in[i] + dt * dudt[i];
    rho[i] = rho_in[i] + dt * drhodt[i];
}

/*
 * @}
 */
 
