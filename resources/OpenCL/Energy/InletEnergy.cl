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
 * @brief Tool to compute the energy injected into the system by the inlet.
 */

#ifndef HAVE_3D
    #include "../types/2D.h"
#else
    #include "../types/3D.h"
#endif

/** @brief Tool to compute the energy injected into the system by the inlet.
 *
 * Actually, in this kernel the energy componets variation are computed per each
 * particle.
 *
 * @param inlet_energy_dekindt Variation of the kinetic energy due to the
 * inflow.
 * @param inlet_energy_depotdt Variation of the potential energy due to the
 * inflow.
 * @param inlet_energy_decomdt Variation of the compressibility energy due to
 * the inflow.
 * @param fluid_energy_dekindt Variation of the kinetic energy (restricted to
 * the fluid particles interactions): \f$ \frac{dE^{kin}_a}{dt} =
 *   m_a \mathbf{u}_a \cdot \frac{d \mathbf{u}_a}{dt}\f$
 * @param fluid_energy_depotdt Variation of the potential energy (restricted to
 * the fluid particles interactions): \f$ \frac{dE^{pot}_a}{dt} =
 *   - m_a \mathbf{g} \cdot \mathbf{u}_a\f$
 * @param fluid_energy_decomdt Variation of the compressibility energy
 * (restricted to the fluid particles interactions): \f$ \frac{dE^{com}_a}{dt} =
 *   \frac{m_a}{\rho_a} \frac{p_a}{\rho_a} \frac{d \rho_a}{dt} \f$
 * @param imove Moving flags.
 *   - imove > 0 for regular fluid particles.
 *   - imove = 0 for sensors.
 *   - imove < 0 for boundary elements/particles.
 * @param r Position \f$ \mathbf{r} \f$.
 * @param inlet_r Lower corner of the inlet square.
 * @param inlet_n Velocity direction of the generated particles.
 * @param N Number of particles.
 */
__kernel void main(__global float* inlet_energy_dekindt,
                   __global float* inlet_energy_depotdt,
                   __global float* inlet_energy_decomdt,
                   __global float* fluid_energy_dekindt,
                   __global float* fluid_energy_depotdt,
                   __global float* fluid_energy_decomdt,
                   const __global int* imove,
                   const __global vec* r,
                   vec inlet_r,
                   vec inlet_n,
                   unsigned int N)
{
    // find position in global arrays
    unsigned int i = get_global_id(0);

    // Discard the particles already passed through the inlet
    if(dot(r[i] - inlet_r, inlet_n) > 0.f)
        return;


    if(i >= N)
        return;
    if((imove[i] != 1) || (dot(r[i] - inlet_r, inlet_n) > 0.f)){
        inlet_energy_dekindt[i] = 0.f;
        inlet_energy_depotdt[i] = 0.f;
        inlet_energy_decomdt[i] = 0.f;
        return;
    }

    inlet_energy_depotdt[i] = fluid_energy_depotdt[i];
    inlet_energy_dekindt[i] = fluid_energy_dekindt[i];
    inlet_energy_decomdt[i] = fluid_energy_decomdt[i];
    fluid_energy_depotdt[i] = 0.f;
    fluid_energy_dekindt[i] = 0.f;
    fluid_energy_decomdt[i] = 0.f;
}
