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
 * @brief Fixed ghost particles mirroring process.
 */

#include "resources/Scripts/types/types.h"
#include "resources/Scripts/KernelFunctions/Kernel.h"

/** @brief Pressure and velocity extensions for the computation of the
 * differential operators, except the Laplacian of the velocity already
 * computed.
 *
 * For the Divergence of the velocity, the normmal velocity is mirrored (ASM),
 * while the mirrored interpolated value is used for the tangential velocity
 * (SSM).
 *
 * For the pressure the following expression is used:
 *
 * \f$ p(\mathbf{r}) = p(\hat{\mathbf{r}}) +
 * \frac{\partial{p}}{\partial{\mathbf{n_b}}} \cdot \left( \hat{\mathbf{r}} -
 * \mathbf{r}\right) \f$
 *
 * where the derivative of the pressure is computed as:
 *
 * \f$ \frac{\partial{p}}{\partial{\mathbf{n_b}}} = \rho \left[
 * \mathbf{g} - \frac{\mathrm{d}\mathbf{u_b}}{\mathrm{d}t}
 * + \frac{\mu}{\rho} \Delta \mathbf{u} \right] \cdot \mathbf{n_b} \f$
 *
 * @param iset Set of particles index.
 * @param imove Moving flags.
 *   - imove > 0 for regular fluid particles.
 *   - imove = 0 for sensors.
 *   - imove < 0 for boundary elements/particles.
 * @param associations Mirroring particles.
 * @param r Position \f$ \mathbf{r} \f$.
 * @param normal Normal \f$ \mathbf{n} \f$.
 * @param dudt Velocity rate of change
 * \f$ \left. \frac{d \mathbf{u}}{d t} \right\vert_{n+1} \f$.
 * @param lap_u Velocity laplacian \f$ \frac{\Delta \mathbf{u}}{rho} \f$.
 * @param gp_p Interpolated pressure in the mirrored position \f$ p \f$.
 * @param gp_u Interpolated velocity in the mirrored position \f$ \mathbf{u} \f$.
 * @param rho Density \f$ \rho \f$.
 * @param p Pressure \f$ p \f$.
 * @param u Velocity \f$ \mathbf{u} \f$.
 * @param visc_dyn Dynamic viscosity \f$ \mu \f$.
 * @param N Number of particles.
 * @param g Gravity acceleration \f$ \mathbf{g} \f$.
 */
__kernel void entry(const __global uint* iset,
                    const __global int* imove,
                    const __global usize* associations,
                    const __global vec* r,
                    const __global vec* normal,
                    const __global vec* dudt,
                    const __global vec* lap_u,
                    const __global float* gp_p,
                    const __global vec* gp_u,
                    __global float* rho,
                    __global float* p,
                    __global vec* u,
                    __constant float* visc_dyn,
                    usize N,
                    vec g)
{
    const usize i = get_global_id(0);
    if(i >= N)
        return;
    if(imove[i] != -1)
        return;

    // Let's get the associated boundary element, and check it is valid
    const usize iref = associations[i];
    if(iref >= N)
        return;
    const vec_xyz n_i = normal[iref].XYZ;

    // Mirror the normal velocity interpolated (ASM), and keep the interpolated
    // tangential velocity
    const vec_xyz u_n = dot(gp_u[i].XYZ, n_i) * n_i;
    const vec_xyz u_t = gp_u[i].XYZ - u_n;
    u[i].XYZ = 2.f * dot(u[iref].XYZ, n_i) * n_i - u_n + u_t;

    // Extend the pressure through the gradient
    const float dpdn = rho[iref] * dot((g.XYZ
        - dudt[iref].XYZ + visc_dyn[iset[iref]] * lap_u[iref].XYZ), n_i);
    p[i] = gp_p[i] + 2.f * dpdn * dot(r[i].XYZ - r[iref].XYZ, n_i);
}
