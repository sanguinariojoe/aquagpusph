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
 * @brief Velocity extension for the computation of the Laplacian of the
 * velocity.
 */

#ifndef HAVE_3D
    #include "../../types/2D.h"
    #include "../../KernelFunctions/Wendland2D.hcl"
#else
    #include "../../types/3D.h"
    #include "../../KernelFunctions/Wendland3D.hcl"
#endif

/** @brief Velocity extension for the computation of the Laplacian of the
 * velocity.
 *
 * For the Laplacian of the velocity, the normmal velocity is taken from the
 * mirroring interpolated value (SSM), while the tangential velocity is taken
 * directly from the body motion (U0M).
 *
 * In the thesis of Benjamin Bouscasse the tangential velocity is asymetrically
 * extended (ASM). However it is producing unphysical results.
 *
 * Herein the density is also got from the mirroring interpolated value (SSM)
 *
 * @param imove Moving flags.
 *   - imove > 0 for regular fluid particles.
 *   - imove = 0 for sensors.
 *   - imove < 0 for boundary elements/particles.
 * @param associations Mirroring particles.
 * @param normal Normal \f$ \mathbf{n} \f$.
 * @param gp_rho Interpolated density in the mirrored position \f$ \rho \f$.
 * @param gp_u Interpolated velocity in the mirrored position \f$ \mathbf{u} \f$.
 * @param rho Density \f$ \rho \f$.
 * @param u Velocity \f$ \mathbf{u} \f$.
 * @param N Number of particles.
 */
__kernel void main(const __global int* imove,
                   const __global uint* associations,
                   const __global vec* normal,
                   const __global float* gp_rho,
                   const __global vec* gp_u,
                   __global float* rho,
                   __global vec* u,
                   uint N)
{
    const uint i = get_global_id(0);
    if(i >= N)
        return;
    if(imove[i] != -1)
        return;

    // Let's get the associated boundary element, and check it is valid
    const uint iref = associations[i];
    if(iref >= N)
        return;

    // Use the interpolated density (SSM)
    rho[i] = gp_rho[i];

    #if __GP_LAPU_MODEL__ == __GP_LAPU_TAKEDA__
        // Benjamin Bouscasse approach: Normal velocity extended by an SSM
        // model, and tangential velocity extended by an ASM model
        const vec_xyz n_iref = normal[iref].XYZ;
        const vec_xyz u_iref = u[iref].XYZ - dot(u[iref].XYZ, n_iref) * n_iref;

        const vec_xyz u_n = dot(gp_u[i].XYZ, n_iref) * n_iref;
        const vec_xyz u_t = 2.f * u_iref - (gp_u[i].XYZ - u_n);

        u[i].XYZ = u_n + u_t;
    #elif __GP_LAPU_MODEL__ == __GP_LAPU_U0M__
        // Forzen velocity approach: The body velocity is kept
        return;
    #else
        #error Unknown extension model: __GP_LAPU_MODEL__
    #endif
}
