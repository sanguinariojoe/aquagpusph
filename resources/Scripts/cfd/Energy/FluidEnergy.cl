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

/** @brief Tool to compute the fluid energy components.
 *
 * Actually, in this kernel the energy componets variation are computed per each
 * particle.
 *
 * @param fluid_energy_dekindt Variation of the kinetic energy:
 * \f$ \frac{dE^{kin}_a}{dt} =
 *   m_a \mathbf{u}_a \cdot \frac{d \mathbf{u}_a}{dt}\f$
 * @param fluid_energy_depotdt Variation of the potential energy:
 * \f$ \frac{dE^{pot}_a}{dt} =
 *   - m_a \mathbf{g} \cdot \mathbf{u}_a\f$
 * @param fluid_energy_decomdt Variation of the compressibility energy:
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
__kernel void entry(__global float* fluid_energy_dekindt,
                    __global float* fluid_energy_depotdt,
                    __global float* fluid_energy_decomdt,
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
        fluid_energy_dekindt[i] = 0.f;
        fluid_energy_depotdt[i] = 0.f;
        fluid_energy_decomdt[i] = 0.f;
        return;
    }

    fluid_energy_depotdt[i] = -m[i] * dot(g, u[i]);
    fluid_energy_dekindt[i] = m[i] * dot(u[i], dudt[i]);
    fluid_energy_decomdt[i] = m[i] * p[i] / (rho[i] * rho[i]) * drhodt[i];
}
