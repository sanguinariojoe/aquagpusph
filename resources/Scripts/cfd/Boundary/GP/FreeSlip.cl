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
 * velocity, when free-slip should be enforced.
 */

#ifndef HAVE_3D
    #include "../../../types/2D.h"
    #include "../../../KernelFunctions/Wendland2D.hcl"
#else
    #include "../../../types/3D.h"
    #include "../../../KernelFunctions/Wendland3D.hcl"
#endif

/** @brief Velocity extension for the computation of the Laplacian of the
 * velocity, when free-slip should be enforced.
 *
 * In this case, both the normal and the tangential velocities are taken from
 * the mirroring interpolated value (SSM).
 *
 * Herein the density is also got from the mirroring interpolated value (SSM)
 *
 * @param iset Set of particles index.
 * @param imove Moving flags.
 *   - imove > 0 for regular fluid particles.
 *   - imove = 0 for sensors.
 *   - imove < 0 for boundary elements/particles.
 * @param associations Mirroring particles.
 * @param gp_u Interpolated velocity in the mirrored position \f$ \mathbf{u} \f$.
 * @param u Velocity \f$ \mathbf{u} \f$.
 * @param N Number of particles.
 * @param freeslip_iset Particles set which should be considered free-slip.
 */
__kernel void entry(const __global uint* iset,
                    const __global int* imove,
                    const __global uint* associations,
                    const __global vec* gp_u,
                    __global vec* u,
                    uint N,
                    uint freeslip_iset)
{
    const uint i = get_global_id(0);
    if(i >= N)
        return;
    if((imove[i] != -1) || (iset[i] != freeslip_iset))
        return;

    // Let's get the associated boundary element, and check it is valid
    const uint iref = associations[i];
    if(iref >= N)
        return;

    u[i].XYZ = gp_u[i].XYZ;
}
