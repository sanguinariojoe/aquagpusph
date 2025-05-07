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
 * @brief Tool to compute the fluid global energy components.
 */

#include "resources/Scripts/types/types.h"

/** @brief Tool to compute the fluid energy variation rate components.
 *
 * Actually, in this kernel the energy componets variation are computed per
 * particle.
 *
 * @param energy_dekdt Variation of the kinetic energy:
 * \f$ \frac{dE^{kin}_a}{dt} =
 *   m_a \mathbf{u}_a \cdot \frac{d \mathbf{u}_a}{dt}\f$
 * @param energy_depdt Variation of the potential energy:
 * \f$ \frac{dE^{pot}_a}{dt} =
 *   - m_a \mathbf{g} \cdot \mathbf{u}_a\f$
 * @param energy_decdt Variation of the compressibility energy:
 * \f$ \frac{dE^{com}_a}{dt} =
 *   \frac{m_a}{\rho_a} \frac{p_a}{\rho_a} \frac{d \rho_a}{dt} \f$
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
 * @param N Number of particles.
 * @param g Gravity acceleration \f$ \mathbf{g} \f$.
 */
__kernel void power(__global float* energy_dekdt,
                    __global float* energy_depdt,
                    __global float* energy_decdt,
                    const __global int* imove,
                    const __global vec* u,
                    const __global float* rho,
                    const __global float* m,
                    const __global float* p,
                    const __global vec* dudt,
                    const __global float* drhodt,
                    usize N,
                    vec g)
{
    // find position in global arrays
    const usize i = get_global_id(0);
    if(i >= N)
        return;
    if(imove[i] != 1){
        energy_dekdt[i] = 0.f;
        energy_depdt[i] = 0.f;
        energy_decdt[i] = 0.f;
        return;
    }

    energy_depdt[i] = -m[i] * dot(g, u[i]);
    energy_dekdt[i] = m[i] * dot(u[i], dudt[i]);
    energy_decdt[i] = m[i] * p[i] / (rho[i] * rho[i]) * drhodt[i];
}

/** @brief Tool to compute the fluid energy components.
 *
 * Actually, in this kernel the energy componets are computed per particle.
 *
 * @param energy_ek Kinetic energy:
 * \f$ E^{kin}_a =
 *   \frac{1}{2} m_a \vert \mathbf{u}_a \vert^2\f$
 * @param energy_ep Potential energy:
 * \f$ E^{pot}_a =
 *   - m_a \mathbf{g} \cdot \mathbf{r}_a\f$
 * @param energy_ec Compressibility energy:
 * \f$ E^{com}_a =
 *   m_a c_0^2 \left( \frac{\rho_0}{\rho_a} + \log(\rho_a) \right) \f$
 * @param iset Set of particles index.
 * @param imove Moving flags.
 *   - imove > 0 for regular fluid particles.
 *   - imove = 0 for sensors.
 *   - imove < 0 for boundary elements/particles.
 * @param r Position \f$ \mathbf{r} \f$.
 * @param u Velocity \f$ \mathbf{u} \f$.
 * @param rho Density \f$ \rho \f$.
 * @param m Mass \f$ m \f$.
 * @param N Number of particles.
 * @param g Gravity acceleration \f$ \mathbf{g} \f$.
 * @param cs Speed of sound \f$ c_s \f$.
 */
__kernel void energy(__global float* energy_ek,
                     __global float* energy_ep,
                     __global float* energy_ec,
                     const __global uint* iset,
                     const __global int* imove,
                     const __global vec* r,
                     const __global vec* u,
                     const __global float* rho,
                     const __global float* m,
                     __constant float* refd,
                     usize N,
                     vec g,
                     float cs)
{
    // find position in global arrays
    const usize i = get_global_id(0);
    if(i >= N)
        return;
    if(imove[i] != 1){
        energy_ek[i] = 0.f;
        energy_ep[i] = 0.f;
        energy_ec[i] = 0.f;
        return;
    }

    energy_ek[i] = 0.5f * m[i] * dot(u[i], u[i]);
    energy_ep[i] = -m[i] * dot(g, r[i]);
    const float rho0 = refd[iset[i]];
    energy_ec[i] = m[i] * cs * cs * (rho0 / rho[i] + log(rho[i] / rho0) - 1.f);
}
