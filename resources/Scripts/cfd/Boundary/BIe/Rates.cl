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
 * @param grad_p Pressure gradient \f$ \frac{\nabla p}{rho} \f$.
 * @param div_u Velocity divergence \f$ \rho \nabla \cdot \mathbf{u} \f$.
 * @param N Number of particles.
 */
__kernel void entry(const __global int* imove,
                    const __global float* rho,
                    const __global float* p,
                    const __global vec* u,
                    const __global vec* grad_w_bi,
                    const __global float* div_u_bi,
                    __global vec* grad_p,
                    __global float* div_u,
                    usize N)
{
    const usize i = get_global_id(0);
    if(i >= N)
        return;
    if(imove[i] != 1)
        return;

    grad_p[i] += 2.f * p[i] / rho[i] * grad_w_bi[i];
    div_u[i] -= 2.f * rho[i] * (dot(u[i], grad_w_bi[i]) + div_u_bi[i]);
}

/** @brief Filter out the fluid particles from the forces.
 *
 * @param imove Moving flags.
 *   - imove > 0 for regular fluid particles.
 *   - imove = 0 for sensors.
 *   - imove < 0 for boundary elements/particles.
 * @param r Position \f$ \mathbf{r} \f$.
 * @param normal Normal \f$ \mathbf{n} \f$.
 * @param m Area of the boundary element \f$ s \f$.
 * @param p Pressure of the boundary element \f$ p \f$.
 * @param force_p Pressure force on the boundary element
 * @param moment_p Pressure moment on the boundary element
 * @param forces_r Point with respect the moments are computed
 * \f$ \mathbf{r}_0 \f$.
 * @param N Number of particles.
 */
__kernel void force_press(const __global int* imove,
                          const __global vec* r,
                          const __global vec* normal,
                          const __global float* m,
                          const __global float* p,
                          __global vec* force_p,
                          __global vec4* moment_p,
                          vec forces_r,
                          usize N)
{
    const usize i = get_global_id(0);
    if(i >= N)
        return;
    if(imove[i] != -3) {
        force_p[i] = VEC_ZERO;
        return;
    }

    // Expand the force and the arm as 3D variables
    vec4 F = (vec4)(0.f);
    vec4 R = (vec4)(0.f);
    F.XYZ = p[i] * m[i] * normal[i].XYZ;
    R.XYZ = r[i].XYZ - forces_r.XYZ;
    force_p[i] = F.XYZ;
    moment_p[i] = cross(R, F);    
}
