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
 * @param energy_deintdt Variation of the internal energy:
 * \f$ \frac{dU_a}{dt} =
 *   \frac{dW_a}{dt} - \frac{dE^{pot}_a}{dt} - \frac{dE^{kin}_a}{dt} \f$
 * @param energy_dsdt Variation of the entropy:
 * \f$ T \frac{d S_a}{d t} =
 *   \frac{dU_a}{dt}
 *   - m_a \frac{p_a}{\rho_a^2} \frac{d \rho_a}{d t} \f$
 * @param energy_dekindt Variation of the kinetic energy:
 * \f$ \frac{dE^{kin}_a}{dt} =
 *   m_a \mathbf{u}_a \cdot \frac{d \mathbf{u}_a}{dt}\f$
 * @param energy_depotdt Variation of the potential energy:
 * \f$ \frac{dE^{pot}_a}{dt} =
 *   - m_a \mathbf{g} \cdot \mathbf{u}_a\f$
 * @param energy_dwdt Variation of the potential energy:
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
 * @param shepard Shepard term
 * \f$ \gamma(\mathbf{x}) = \int_{\Omega}
 *     W(\mathbf{y} - \mathbf{x}) \mathrm{d}\mathbf{x} \f$.
 * @param dudt Velocity rate of change
 * \f$ \frac{d \mathbf{u}}{d t} \f$.
 * @param drhodt Density rate of change
 * \f$ \frac{d \rho}{d t} \f$.
 * @param grad_ux Gradient of the first component of the velocity:
 *   \f$ \nabla \left(\mathbf{u} \cdot \mathbf{e_1}\right) \f$
 * @param grad_uy Gradient of the second component of the velocity:
 *   \f$ \nabla \left(\mathbf{u} \cdot \mathbf{e_2}\right) \f$
 * @param grad_uz Gradient of the third component of the velocity:
 *   \f$ \nabla \left(\mathbf{u} \cdot \mathbf{e_3}\right) \f$
 * @param visc_dyn Dynamic viscosity \f$ \mu \f$.
 * @param N Number of particles.
 * @param dt Time step \f$ \Delta t \f$.
 * @param g Gravity acceleration \f$ \mathbf{g} \f$.
 */
__kernel void main(__global float* energy_deintdt,
                   __global float* energy_dsdt,
                   __global float* energy_dekindt,
                   __global float* energy_depotdt,
                   __global float* energy_dwdt,
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
                   const __global float* shepard,
                   const __global vec* dudt,
                   const __global float* drhodt,
                   const __global vec4* grad_ux,
                   const __global vec4* grad_uy,
                   const __global vec4* grad_uz,
                   __constant float* visc_dyn,
                   unsigned int N,
                   float dt,
                   vec g)
{
    // find position in global arrays
    unsigned int i = get_global_id(0);
    if(i >= N)
        return;
    if(imove[i] <= 0){
        energy_deintdt[i] = 0.f;
        energy_dsdt[i] = 0.f;
        energy_dekindt[i] = 0.f;
        energy_depotdt[i] = 0.f;
        energy_dwdt[i] = 0.f;
        return;
    }

    const float mass = m[i];
    const float dens = rho[i];
    const float prfac = p[i] / (dens * dens);
    const float mu = visc_dyn[iset[i]];
    const float gradu_visc_term = 2.f * (grad_ux[i].x * grad_ux[i].x +
                                         grad_uy[i].y * grad_uy[i].y +
                                         grad_uz[i].z * grad_uz[i].z)
        + (grad_uy[i].x + grad_ux[i].y) * (grad_uy[i].x + grad_ux[i].y)
        + (grad_uz[i].x + grad_ux[i].z) * (grad_uz[i].x + grad_uz[i].y)
        + (grad_uz[i].y + grad_uy[i].z) * (grad_uz[i].y + grad_uy[i].z);


    // External work
    energy_depotdt[i] = -mass * dot(g, u[i]);
    energy_dwdt[i] = mass / shepard[i] * (mu * dot(lap_u[i], u[i])
                                          + mu * gradu_visc_term
                                          - dot(grad_p[i], u[i])
                                          - prfac * div_u[i]);

    // Fluid particle energy
    energy_dekindt[i] = mass * dot(u[i], dudt[i]);
    energy_deintdt[i] = energy_dwdt[i] - energy_dekindt[i] - energy_depotdt[i];

    // Entropy part of the energy
    energy_dsdt[i] = energy_deintdt[i] - mass * prfac * drhodt[i];
}
