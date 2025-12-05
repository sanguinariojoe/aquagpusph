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
 * @brief Pressure gradient and velocity divergence updating with the boundary
 * integrals.
 */

#include "resources/Scripts/types/types.h"

/** @brief Pressure gradient and velocity divergence updating with the boundary
 * integrals.
 *
 * @param imove Moving flags.
 *   - imove > 0 for regular fluid particles.
 *   - imove = 0 for sensors.
 *   - imove < 0 for boundary elements/particles.
 * @param rho Density \f$ \rho_{n+1} \f$.
 * @param p Pressure \f$ p \f$.
 * @param u Velocity \f$ \mathbf{u} \f$.
 * @param grad_w_bi Gradient of constant fields due to the boundary integral
 * \f$ \langle \nabla 1 \rangle^{\partial \Omega} \f$.
 * @param div_u Velocity divergence \f$ \nabla \cdot \mathbf{u} \f$. Actually
 * this is just the part that has to do with the boundary element velocity
 * @param work_density Work density.
 * @param N Number of particles.
 */
__kernel void entry(const __global int* restrict imove,
                    const __global float* restrict rho,
                    const __global float* restrict p,
                    const __global vec* restrict u,
                    const __global vec* restrict grad_w_bi,
                    const __global float* restrict div_u_bi,
                    __global float* restrict work_density,
                    usize N)
{
    const usize i = get_global_id(0);
    if(i >= N)
        return;
    if(imove[i] != 1)
        return;

    work_density[i] -= p[i] / rho[i] *(2.f * dot(u[i], grad_w_bi[i]) + div_u_bi[i]);
   
}

