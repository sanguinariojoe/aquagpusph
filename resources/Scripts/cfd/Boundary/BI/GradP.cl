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
 * @brief Compute the gradient of the pressure at the boundary to impose the
 * Neumann BC.
 */

#include "resources/Scripts/types/types.h"

/** @brief Compute the gradient of the pressure at the boundary to impose the
 * free-slip BC, i.e.
 *
 * \f$ \left. \frac{d\mathbf{u}}{dt} \right\vert_{\partial \Omega}
 *       \cdot \mathbf{n} =
 *     \left. \frac{d\mathbf{u}}{dt} \right\vert_{\partial \Omega}
 *       \cdot \mathbf{n}\f$.
 *
 * Such that the pressure can be computed from the Momentum Equation.
 *
 * @param iset Set of particles index.
 * @param imove Moving flags.
 *   - imove > 0 for regular fluid particles.
 *   - imove = 0 for sensors.
 *   - imove < 0 for boundary elements/particles.
 * @param shepard Shepard term
 * \f$ \gamma(\mathbf{x}) = \int_{\Omega}
 *     W(\mathbf{y} - \mathbf{x}) \mathrm{d}\mathbf{x} \f$.
 * @param rho Density \f$ \rho \f$.
 * @param p Pressure \f$ p \f$.
 * @param rho Density \f$ \rho_{n+1} \f$.
 * @param lap_u Velocity laplacian \f$ \Delta \mathbf{u} \f$.
 * @param dudt Velocity rate of change \f$ \frac{d \mathbf{u}}{d t} \f$.
 * @param grad_p Pressure gradient \f$ \frac{\nabla p}{rho} \f$.
 * @param visc_dyn Dynamic viscosity \f$ \mu \f$.
 * @param refd Density of reference \f$ \rho_0 \f$.
 * @param N Number of particles.
 * @param g Gravity acceleration \f$ \mathbf{g} \f$.
 */
__kernel void freeslip(const __global uint* iset,
                       const __global int* imove,
                       const __global float* shepard,
                       const __global float* rho,
                       const __global vec* lap_u,
                       const __global vec* dudt,
                       __global vec* grad_p,
                       __constant float* visc_dyn,
                       __constant float* refd,
                       unsigned int N,
                       vec g)
{
    uint i = get_global_id(0);
    if(i >= N)
        return;

    if(imove[i] != -3){
        return;
    }

    float shepard_i = shepard[i];
    if(shepard_i < 1.0E-6f){
        // It will be considered that there are not enough
        // particles to interpolate
        shepard_i = 1.f;
    }

    grad_p[i] = visc_dyn[iset[i]] / refd[iset[i]] * lap_u[i] / shepard_i + g - dudt[i];
}
