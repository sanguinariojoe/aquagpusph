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
 * @brief Tool to compute the power energy components.
 */

#ifndef HAVE_3D
    #include "../../types/2D.h"
#else
    #include "../../types/3D.h"
#endif

/** @brief Tool to compute the power components due to the interactions
 * between fluid particles, i.e. excluding the effect of the boundaries.
 *
 * @param dekdt Variation of the kinetic energy:
 * \f$ \frac{dE^{kin}_a}{dt} =
 *   m_a \mathbf{u}_a \cdot \frac{d \mathbf{u}_a}{dt}\f$
 * @param depdt Variation of the potential energy:
 * \f$ \frac{dE^{pot}_a}{dt} =
 *   - m_a \mathbf{g} \cdot \mathbf{u}_a\f$
 * @param devdt Viscous dissipation function:
 * \f$ \frac{dE^{\mu}_a}{dt} =
 *   - \mu m_a \mathbf{u}_a \cdot \frac{\Delta \mathbf{u}}{rho} \f$
 * @param decdt Variation of the compressibility energy:
 * \f$ \frac{dE^{com}_a}{dt} =
 *   \frac{m_a}{\rho_a} \frac{p_a}{\rho_a} \frac{d \rho_a}{dt} \f$
 * @param deddt Energy due tot he compressibility dissipated by the
 * \f$\delta\f$-SPH term:
 * \f$ \frac{dE^{\delta}_a}{dt} =
 *   - \delta \Delta t \frac{m_a}{\rho_0} \frac{p_a}{\rho_a^2} \Delta p \f$
 * @param iset Set of particles index.
 * @param imove Moving flags.
 *   - imove > 0 for regular fluid particles.
 *   - imove = 0 for sensors.
 *   - imove < 0 for boundary elements/particles.
 * @param u Velocity \f$ \mathbf{u} \f$.
 * @param rho Density \f$ \rho \f$.
 * @param m Mass \f$ m \f$.
 * @param p Pressure \f$ p \f$.
 * @param dudt Velocity rate of change
 * \f$ \frac{d \mathbf{u}}{d t} \f$.
 * @param drhodt Density rate of change
 * \f$ \frac{d \rho}{d t} \f$.
 * @param lap_u Velocity laplacian \f$ \frac{\Delta \mathbf{u}}{rho} \f$.
 * @param lap_p Pressure laplacian \f$ \Delta p \f$.
 * @param refd Density of reference of the fluid \f$ \rho_0 \f$.
 * @param visc_dyn Dynamic viscosity \f$ \mu \f$.
 * @param delta Diffusive term \f$ \delta \f$ multiplier.
 * @param N Number of particles.
 * @param g Gravity acceleration \f$ \mathbf{g} \f$.
 * @param dt Time step \f$ \Delta t \f$.
 */
__kernel void fluid(__global float* dekdt,
                    __global float* depdt,
                    __global float* devdt,
                    __global float* decdt,
                    __global float* deddt,
                    const __global uint* iset,
                    const __global int* imove,
                    const __global vec* u,
                    const __global float* rho,
                    const __global float* m,
                    const __global float* p,
                    const __global vec* dudt,
                    const __global float* drhodt,
                    const __global vec* lap_u,
                    const __global float* lap_p,
                    __constant float* refd,
                    __constant float* visc_dyn,
                    __constant float* delta,
                    unsigned int N,
                    vec g,
                    float dt)
{
    // find position in global arrays
    unsigned int i = get_global_id(0);
    if(i >= N)
        return;
    if(imove[i] != 1){
        dekdt[i] = 0.f;
        depdt[i] = 0.f;
        devdt[i] = 0.f;
        decdt[i] = 0.f;
        deddt[i] = 0.f;
        return;
    }

    depdt[i] = -m[i] * dot(u[i], g);
    dekdt[i] = m[i] * dot(u[i], dudt[i]);
    decdt[i] = m[i] * p[i] / (rho[i] * rho[i]) * drhodt[i];

    const uint set_i = iset[i];
    const float delta_f = delta[set_i] * dt * rho[i] / refd[set_i];
    devdt[i] = - m[i] * dot(u[i], visc_dyn[iset[i]] * lap_u[i]);
    deddt[i] = - m[i] * p[i] / (rho[i] * rho[i]) * delta_f * lap_p[i];
}

/** @brief Tool to compute the power due to the interactions of fluid particles
 * with the boundaries.
 *
 * @param desdt power due to the interactions of fluid particles with the
 * boundaries.
 * @param iset Set of particles index.
 * @param imove Moving flags.
 *   - imove > 0 for regular fluid particles.
 *   - imove = 0 for sensors.
 *   - imove < 0 for boundary elements/particles.
 * @param u Velocity \f$ \mathbf{u} \f$.
 * @param rho Density \f$ \rho \f$.
 * @param m Mass \f$ m \f$.
 * @param p Pressure \f$ p \f$.
 * @param dudt Velocity rate of change
 * \f$ \frac{d \mathbf{u}}{d t} \f$.
 * @param drhodt Density rate of change
 * \f$ \frac{d \rho}{d t} \f$.
 * @param dekdt Variation of the kinetic energy due to the interactions between
 * fluid particles:
 * \f$ \frac{dE^{kin}_a}{dt} =
 *   m_a \mathbf{u}_a \cdot \frac{d \mathbf{u}_a}{dt}\f$
 * @param decdt Variation of the compressibility energy due to the interactions
 * between fluid particles:
 * \f$ \frac{dE^{com}_a}{dt} =
 *   \frac{m_a}{\rho_a} \frac{p_a}{\rho_a} \frac{d \rho_a}{dt} \f$
 * @param N Number of particles.
 */
__kernel void bound(__global float* desdt,
                    const __global uint* iset,
                    const __global int* imove,
                    const __global vec* u,
                    const __global float* rho,
                    const __global float* m,
                    const __global float* p,
                    const __global vec* dudt,
                    const __global float* drhodt,
                    const __global float* dekdt,
                    const __global float* decdt,
                    unsigned int N)
{
    // find position in global arrays
    unsigned int i = get_global_id(0);
    if(i >= N)
        return;
    if(imove[i] != 1){
        desdt[i] = 0.f;
        return;
    }

    const float pk = m[i] * dot(u[i], dudt[i]);
    const float pc = m[i] * p[i] / (rho[i] * rho[i]) * drhodt[i];
    
    desdt[i] = (dekdt[i] - pk) + (decdt[i] - pc);
}