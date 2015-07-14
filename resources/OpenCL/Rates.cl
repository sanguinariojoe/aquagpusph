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
 * @brief Velocity and density variation rates computation.
 */

#ifndef HAVE_3D
    #include "types/2D.h"
#else
    #include "types/3D.h"
#endif

/** @brief Velocity and density variation rates computation.
 *
 * The mass conservation and momentum equations are applied from the already
 * computed differential operators:
 *
 *   - \f$ \frac{\mathrm{d} \mathbf{u}}{\mathrm{d} t} =
 *     - \frac{\nabla p}{rho}
 *     + \frac{\mu}{rho} \Delta \mathbf{u}
 *     + \mathbf{g}\f$
 *   - \f$ \frac{\mathrm{d} \rho}{\mathrm{d} t} =
 *     - \rho \nabla \cdot \mathbf{u}
 *     + \delta \Delta t \frac{\rho_a}{\rho_0} \Delta p\f$
 *
 * @param iset Set of particles index.
 * @param imove Moving flags.
 *   - imove > 0 for regular fluid particles.
 *   - imove = 0 for sensors.
 *   - imove < 0 for boundary elements/particles.
 * @param rho Density \f$ \rho_{n+1} \f$.
 * @param grad_p Pressure gradient \f$ \frac{\nabla p}{rho} \f$.
 * @param lap_u Velocity laplacian \f$ \frac{\Delta \mathbf{u}}{rho} \f$.
 * @param div_u Velocity divergence \f$ \rho \nabla \cdot \mathbf{u} \f$.
 * @param lap_p Pressure laplacian \f$ \Delta p \f$.
 * @param shepard Shepard term
 * \f$ \gamma(\mathbf{x}) = \int_{\Omega}
 *     W(\mathbf{y} - \mathbf{x}) \mathrm{d}\mathbf{x} \f$.
 * @param dudt Velocity rate of change
 * \f$ \left. \frac{d \mathbf{u}}{d t} \right\vert_{n+1} \f$.
 * @param drhodt Density rate of change
 * \f$ \left. \frac{d \rho}{d t} \right\vert_{n+1} \f$.
 * @param visc_dyn Dynamic viscosity \f$ \mu \f$.
 * @param delta Diffusive term \f$ \delta \f$ multiplier.
 * @param refd Density of reference of the fluid \f$ \rho_0 \f$.
 * @param N Number of particles.
 * @param dt Time step \f$ \Delta t \f$.
 * @param g Gravity acceleration \f$ \mathbf{g} \f$.
 */
__kernel void main(const __global uint* iset,
                   const __global int* imove,
                   const __global float* rho,
                   const __global vec* grad_p,
                   const __global vec* lap_u,
                   const __global float* div_u,
                   const __global float* lap_p,
                   const __global float* shepard,
                   __global vec* dudt,
                   __global float* drhodt,
                   __constant float* visc_dyn,
                   __constant float* delta,
                   __constant float* refd,
                   unsigned int N,
                   float dt,
                   vec g)
{
    unsigned int i = get_global_id(0);
    if(i >= N)
        return;
    if(imove[i] <= 0)
        return;

    const uint set_i = iset[i];
    const float delta_f = delta[set_i] * dt * rho[i] / refd[set_i];

    // Momentum equation
    dudt[i] = (-grad_p[i] + visc_dyn[set_i] * lap_u[i]) / shepard[i] + g;
    // Conservation of mass equation
    drhodt[i] = -div_u[i] + delta_f * lap_p[i];
}
