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
 * @brief Leap-frog time integration scheme predictor stage.
 * (See Aqua::CalcServer::Predictor for details)
 */

#ifndef HAVE_3D
    #include "types/2D.h"
#else
    #include "types/3D.h"
#endif

/** @brief Leap-frog time integration scheme predictor stage.
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
 * @param iset Set of particles index.
 * @param pos Position \f$ \mathbf{r}_{n+1} \f$.
 * @param v Velocity \f$ \mathbf{u}_{n+1} \f$.
 * @param dvdt Velocity rate of change
 * \f$ \left. \frac{d \mathbf{u}}{d t} \right\vert_{n+1} \f$.
 * @param rho Density \f$ \rho_{n+1} \f$.
 * @param drhodt Density rate of change
 * \f$ \left. \frac{d \rho}{d t} \right\vert_{n+1} \f$.
 * @param pos_in Position \f$ \mathbf{r}_{n+1/2} \f$.
 * @param v_in Velocity \f$ \mathbf{u}_{n+1/2} \f$.
 * @param dvdt_in Velocity rate of change
 * \f$ \left. \frac{d \mathbf{u}}{d t} \right\vert_{n+1/2} \f$.
 * @param rho_in Density \f$ \rho_{n+1/2} \f$.
 * @param drhodt_in Density rate of change
 * \f$ \left. \frac{d \rho}{d t} \right\vert_{n+1/2} \f$.
 * @param mass Mass \f$ m_{n+1} \f$.
 * @param press Pressure \f$ p \f$.
 * @param gamma Eq. of state exponent \f$ \gamma \f$.
 * @param refd Density of reference of the fluid \f$ \rho_0 \f$.
 * @param N Number of particles.
 * @param dt Time step \f$ \Delta t \f$.
 * @param cs Speed of sound \f$ c_s \f$.
 * @param g Gravity acceleration \f$ \mathbf{g} \f$.
 * @param rho_min Minimum tolerated density value \f$ \rho_{min} \f$.
 * @param rho_max Maximum tolerated density value \f$ \rho_{max} \f$.
 * @see Corrector.cl
 */
__kernel void main(__global int* imove,
                   __global unsigned int* iset,
                   __global vec* pos,
                   __global vec* v,
                   __global vec* dvdt,
                   __global float* rho,
                   __global float* drhodt,
                   __global vec* pos_in,
                   __global vec* v_in,
                   __global vec* dvdt_in,
                   __global float* rho_in,
                   __global float* drhodt_in,
                   __global float* p_in,
                   __constant float* gamma,
                   __constant float* refd,
                   unsigned int N,
                   float dt,
                   float cs,
                   vec g,
                   float rho_min,
                   float rho_max)
{
    // find position in global arrays
    unsigned int i = get_global_id(0);
    if(i >= N)
        return;

    // Momentum equation is solved just for the fluid particles
    float DT = dt;
    if(imove[i] <= 0)
        DT = 0.f;
    dvdt_in[i] = dvdt[i];
    v_in[i] = v[i] + DT * dvdt[i];
    pos_in[i] = pos[i] + DT * v[i] + 0.5f * DT * DT * dvdt[i];
    
    // Continuity equation must be solved for the fixed particles too
    if(imove[i] == -1){
        DT = dt;
    }
    drhodt_in[i] = drhodt[i];
    rho_in[i] = rho[i] + DT * drhodt[i];
    if(rho_in[i] < rho_min) rho_in[i] = rho_min;
    if(rho_in[i] > rho_max) rho_in[i] = rho_max;

    // Batchelor 1967
    {
        const float ddenf = rho_in[i] / refd[iset[i]];
        const float prb = cs * cs * refd[iset[i]] / gamma[iset[i]];
        p_in[i] = prb * (pow(ddenf, gamma[iset[i]]) - 1.f);
    }
}
