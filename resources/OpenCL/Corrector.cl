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
 * @brief Leap-frog time integration scheme corrector stage.
 */

#ifndef HAVE_3D
    #include "types/2D.h"
#else
    #include "types/3D.h"
#endif

/** @brief Leap-frog time integration scheme corrector stage.
 *
 * Time integration is based in the following quasi-second order
 * Predictor-Corrector integration scheme:
 *   - \f$ \mathbf{u}_{n+1} = \mathbf{u}_{n} + \Delta t \left(
        \mathbf{g} +
        \left. \frac{\mathrm{d}\mathbf{u}}{\mathrm{d}t} \right\vert_{n+1/2}
     \right)
     + \frac{\Delta t}{2} \left(
        \left. \frac{\mathrm{d}\mathbf{u}}{\mathrm{d}t} \right\vert_{n + 1/2} -
        \left. \frac{\mathrm{d}\mathbf{u}}{\mathrm{d}t} \right\vert_{n - 1/2}
     \right)
     \f$
 *   - \f$ \mathbf{r}_{n+1} = \mathbf{r}_{n} + \Delta t \, \mathbf{u}_{n}
     + \frac{\Delta t^2}{2} \left(
        \mathbf{g} +
        \left. \frac{\mathrm{d}\mathbf{u}}{\mathrm{d}t} \right\vert_{n+1/2}
     \right)
     \f$
 *   - \f$ \rho_{n+1} = \rho_{n} + \Delta t
        \left. \frac{\mathrm{d}\rho}{\mathrm{d}t} \right\vert_{n+1/2}
     + \frac{\Delta t}{2} \left(
        \left. \frac{\mathrm{d}\rho}{\mathrm{d}t} \right\vert_{n + 1/2} -
        \left. \frac{\mathrm{d}\rho}{\mathrm{d}t} \right\vert_{n - 1/2}
     \right)
     \f$
 *
 * @param imove Moving flags.
 *   - imove > 0 for regular fluid particles.
 *   - imove = 0 for sensors.
 *   - imove < 0 for boundary elements/particles.
 * @param r Position \f$ \mathbf{r}_{n+1/2} \f$.
 * @param v Velocity \f$ \mathbf{u}_{n+1/2} \f$.
 * @param dvdt Velocity rate of change
 * \f$ \left. \frac{d \mathbf{u}}{d t} \right\vert_{n+1/2} \f$.
 * @param rho Density \f$ \rho_{n+1/2} \f$.
 * @param drhodt Density rate of change
 * \f$ \left. \frac{d \rho}{d t} \right\vert_{n+1/2} \f$.
 * @param r_in Position \f$ \mathbf{r}_{n} \f$.
 * @param v_in Velocity \f$ \mathbf{u}_{n} \f$.
 * @param dvdt_in Velocity rate of change
 * \f$ \left. \frac{d \mathbf{u}}{d t} \right\vert_{n-1/2} \f$.
 * @param rho_in Density \f$ \rho_{n} \f$.
 * @param drhodt_in Density rate of change
 * \f$ \left. \frac{d \rho}{d t} \right\vert_{n-1/2} \f$.
 * @param N Number of particles.
 * @param dt Time step \f$ \Delta t \f$.
 * @see Predictor.cl
 */
__kernel void main(__global int* imove,
                   __global unsigned int* iset,
                   __global vec* r,
                   __global vec* v,
                   __global vec* dvdt,
                   __global float* rho,
                   __global float* drhodt,
                   __global vec* r_in,
                   __global vec* v_in,
                   __global vec* dvdt_in,
                   __global float* rho_in,
                   __global float* drhodt_in,
                   unsigned int N,
                   float dt)
{
    // find position in global arrays
    unsigned int i = get_global_id(0);
    if(i >= N)
        return;

    float DT, HDT;
    DT = dt;
    if(imove[i] <= 0)
        DT = 0.f;
    HDT = 0.5f * DT;

    // Corrector step for the fluid
    v[i] += HDT * (dvdt[i] - dvdt_in[i]);
    if(imove[i] == -1){
        // Continuity equation must be solved for fixed particles too
        HDT = 0.5f * dt;
    }
    rho[i] += HDT * (drhodt[i] - drhodt_in[i]);

    r_in[i] = r[i];
    v_in[i] = v[i];
    rho_in[i] = rho[i];
    dvdt_in[i] = dvdt[i];
    drhodt_in[i] = drhodt[i];
}
