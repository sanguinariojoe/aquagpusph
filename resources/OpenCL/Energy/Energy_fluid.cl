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

#ifndef HAVE_3D
    #include "../types/2D.h"
#else
    #include "../types/3D.h"
#endif

/** @brief Tool to compute the fluid energy components.
 *
 * Actually, in this kernel the energy componets variation are computed per each
 * particle.
 *
 * @see EnergyVisc.cl
 * 
 * @param energy_fluid_deintdt Variation of the internal energy:
 * \f$ \frac{dU_a}{dt} =
 *   \frac{dW_a}{dt} - \frac{dE^{pot}_a}{dt} - \frac{dE^{kin}_a}{dt} \f$
 * @param energy_fluid_dsdt Variation of the entropy:
 * \f$ T \frac{d S_a}{d t} =
 *   \frac{dU_a}{dt}
 *   - m_a \frac{p_a}{\rho_a^2} \frac{d \rho_a}{d t} \f$
 * @param energy_fluid_dekindt Variation of the kinetic energy:
 * \f$ \frac{dE^{kin}_a}{dt} =
 *   m_a \mathbf{u}_a \cdot \frac{d \mathbf{u}_a}{dt}\f$
 * @param energy_fluid_depotdt Variation of the potential energy:
 * \f$ \frac{dE^{pot}_a}{dt} =
 *   - m_a \mathbf{g} \cdot \mathbf{u}_a\f$
 * @param energy_fluid_dwdt Variation of the potential energy:
 * \f$ \frac{dW_a}{dt} =
 *   - \frac{m_a}{\rho_a} \nabla p_a \cdot \mathbf{u}_a
 *   - \frac{m_a}{\rho_a} p_a \nabla \cdot \mathbf{u}_a 
 *   + \mu \frac{m_a}{\rho_a} \left. \frac{dw_a}{dt} \right\vert_\mu \f$
 * The last term was previously computed in EnergyVisc.cl.
 * @param imove Moving flags.
 *   - imove > 0 for regular fluid particles.
 *   - imove = 0 for sensors.
 *   - imove < 0 for boundary elements/particles.
 * @param iset Set of particles index.
 * @param u Velocity \f$ \mathbf{u} \f$.
 * @param rho Density \f$ \rho \f$.
 * @param m Mass \f$ m \f$.
 * @param p Pressure \f$ p \f$.
 * @param grad_p Pressure gradient \f$ \frac{\nabla p}{rho} \f$.
 * @param lap_u Velocity laplacian \f$ \frac{\Delta \mathbf{u}}{rho} \f$.
 * @param div_u Velocity divergence \f$ \rho \nabla \cdot \mathbf{u} \f$.
 * @param lap_p Pressure laplacian \f$ \Delta p \f$.
 * @param visc_dyn Dynamic viscosity \f$ \mu \f$.
 * @param delta Diffusive term \f$ \delta \f$ multiplier.
 * @param refd Density of reference of the fluid \f$ \rho_0 \f$.
 * @param N Number of particles.
 * @param dt Time step \f$ \Delta t \f$.
 * @param g Gravity acceleration \f$ \mathbf{g} \f$.
 */
__kernel void main(__global float* energy_fluid_deintdt,
                   __global float* energy_fluid_dsdt,
                   __global float* energy_fluid_dekindt,
                   __global float* energy_fluid_depotdt,
                   __global float* energy_fluid_dwdt,
                   const __global int* imove,
                   const __global int* iset,
                   const __global vec* u,
                   const __global float* rho,
                   const __global float* m,
                   const __global float* p,
                   const __global vec* grad_p,
                   const __global vec* lap_u,
                   const __global float* div_u,
                   const __global float* lap_p,
                   __constant float* visc_dyn,
                   __constant float* delta,
                   __constant float* refd,
                   unsigned int N,
                   float dt,
                   vec g)
{
    // find position in global arrays
    unsigned int i = get_global_id(0);
    if(i >= N)
        return;
    if(imove[i] <= 0){
        energy_fluid_deintdt[i] = 0.f;
        energy_fluid_dsdt[i] = 0.f;
        energy_fluid_dekindt[i] = 0.f;
        energy_fluid_depotdt[i] = 0.f;
        energy_fluid_dwdt[i] = 0.f;
        return;
    }

    const uint set_i = iset[i];
    const float mass = m[i];
    const float dens = rho[i];
    const float prfac = p[i] / (dens * dens);
    const float mu = visc_dyn[iset[i]];
    const float delta_f = delta[set_i] * dt * dens / refd[set_i];

    // Since we can be computing the energy before the rates computation, it is
    // convenient to recompute the acceleration and density variation rate
    const vec dudt = (-grad_p[i] + visc_dyn[set_i] * lap_u[i]);
    const float drhodt = -div_u[i] + delta_f * lap_p[i];

    // External work
    energy_fluid_depotdt[i] = 0.f;
    // Since its not possible to make a local analysis of the internal energy
    // exchanged due to viscous terms, it is assumed symmetric, such that the
    // global energy is conserved
    energy_fluid_dwdt[i] = mass * (- dot(grad_p[i], u[i])
                                   - prfac * div_u[i]);

    // Fluid particle energy
    energy_fluid_dekindt[i] = mass * dot(u[i], dudt);
    energy_fluid_deintdt[i] =   energy_fluid_dwdt[i]
                              - energy_fluid_dekindt[i]
                              - energy_fluid_depotdt[i];

    // Entropy part of the energy
    energy_fluid_dsdt[i] = energy_fluid_deintdt[i] - mass * prfac * drhodt;
}
