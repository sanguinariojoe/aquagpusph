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

#include EOS_MODEL
#include SOUND_MODEL
#include "resources/Scripts/types/types.h"


/** @brief Compute the characteristics of the outwards waves.
 * @param imove Moving flags.
 *   - imove > 0 for regular fluid particles.
 *   - imove = 0 for sensors.
 *   - imove < 0 for boundary elements/particles.
 * @param iset Set of particles index.
 * @param r Position \f$ \mathbf{r} \f$.
 * @param u Velocity \f$ \mathbf{u} \f$.
 * @param rho Density \f$ \rho \f$.
 * @param p Pressure \f$ p \f$.
 * @param j1 First characteristic \f$ J_1 \f$.
 * @param j2 Second characteristic \f$ J_2 \f$.
 * @param j3 Third characteristic \f$ J_3 \f$.
 * @param refd Density of reference of the fluid \f$ \rho_0 \f$.
 * @param N Number of particles.
 * @param dt Time step \f$ \Delta t \f$.
 * @param cs Speed of sound \f$ c_s \f$.
 * @param p0 Background pressure \f$ p_0 \f$.
 * @param g Gravity acceleration \f$ \mathbf{g} \f$.
 * @param io_r Lower corner of the inlet/outlet square.
 * @param io_n = Velocity direction.
 * @param io_U = Constant velocity magnitude.
 * @param io_rFS The point where the pressure is the reference one
 * (\f$ p_0 \f$).
 */
__kernel void characteristics(const __global int* restrict imove,
                              const __global unsigned int* restrict iset,
                              const __global vec* restrict r,
                              const __global vec* restrict u,
                              const __global float* restrict rho,  
                              const __global float* restrict eint,
                              const __global float* restrict gamma,
                              __global float* restrict j1,
                              __global float* restrict j2,
                              __global float* restrict j3,
                              usize N,
                              usize nbuffer,
                              float dt,
                              vec g,
                              vec io_r,
                              vec io_n,
                              float io_U,
                              float io_rho,
                              float io_eint,
                              float io_gamma)
{
    const usize i = get_global_id(0);
    if(i >= N)
        return;
    if(imove[i] != 1)
        return;

    // Discard the particles at the inlet/outlet
    // left of the door
    if(dot(r[i] - io_r, INWARD_NORMAL_SIGN * io_n) < 0.f)
        return;
    const float p_i = p_from_rho_eint(gamma[iset[i]], rho[i], eint[i]);
    const float cs_i = sound_speed_perfect_gas(gamma[iset[i]], p_i, rho[i]);

    const float cs_i2 = cs_i * cs_i;
    const float un = dot(u[i], io_n);

    // Get the reference values
    const float uref = io_U;
    const float pref = p_from_rho_eint(io_gamma, io_rho, io_eint);
    const float rhoref = io_rho;

    j1[i] = -cs_i2 * (rho[i] - rhoref) + p_i - pref;
    j2[i] = rho[i] * cs_i * (un - uref) + p_i - pref;
    j3[i] = -rho[i] * cs_i * (un - uref) + p_i - pref;
}



__kernel void values(const __global int* restrict imove,
                     const __global unsigned int* restrict iset,
                     const __global vec* restrict r,
                     __global vec* restrict u,
                     __global float* restrict rho,  
                     __global float* restrict eint,
                     const __global float* restrict gamma,
                     __global float* restrict p,
                     const __global float* restrict j1,
                     const __global float* restrict j2,
                     const __global float* restrict j3,
                     usize N,
                     vec g, //maybe keep it here for further treatment
                     vec io_r,
                     vec io_n,
                     float io_U,
                     float io_rho,
                     float io_eint, 
                     float io_gamma)
{
    //printf("reaches inlet-outlet\n");
    const usize i = get_global_id(0);
    if(i >= N)
        return;
    if(imove[i] != 1)
        return;

    // Discard the particles that already passed through the inlet/outlet
    if(dot(r[i] - io_r, INWARD_NORMAL_SIGN * io_n) > 0.f)
        return;

    const float p_i = p_from_rho_eint(gamma[iset[i]], rho[i], eint[i]);
    const float cs_i = sound_speed_perfect_gas(gamma[iset[i]], p_i, rho[i]);

    const float cs_i2 = cs_i * cs_i;

    // Get the reference values
    const float uref = io_U;
    const float pref = p_from_rho_eint(io_gamma, io_rho, io_eint);
    const float rhoref = io_rho;

    // const float uref = io_U;
    // const float pref = refd[iset[i]] * dot(g, r[i] - io_rFS) + p0;
    // const float rhoref = refd[iset[i]] + (p[i] - p0) / cs2;

    rho[i] = rhoref + 1.f / cs_i2 * (-j1[i] + 0.5f * j2[i] + 0.5f * j3[i]);
    u[i] = (uref + 1.f / (2.f * rho[i] * cs_i) * (j2[i] - j3[i])) * io_n;
    p[i] = pref + 0.5f * (j2[i] + j3[i]);
    eint[i] = eint_from_rho_p(gamma[iset[i]], rho[i], p[i]);
    
    printf("\n");    
    printf("Routine values:\n");
    printf("Inlet-outlet rho: %g\n", rho[i]);
    printf("Inlet-outlet p: %g\n", p[i]);
    printf("Inlet-outlet eint: %g\n", eint[i]);
    printf("Inlet-outlet velocity: %g, %g, %g\n", u[i].x, u[i].y, u[i].z);
    printf("\n");
}