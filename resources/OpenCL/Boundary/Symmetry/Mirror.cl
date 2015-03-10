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
 * @brief Mirroring process for the symmetry boundary condition.
 */

#ifndef HAVE_3D
    #include "../../types/2D.h"
#else
    #include "../../types/3D.h"
#endif

/** @brief Compute the mirrored position of the fluid particles.
 *
 * The mirrored particles (the ones close enough to the symmetry plane) will be
 * marked with \a imirrored = 1.
 * 
 * @param r Position \f$ \mathbf{r} \f$.
 * @param u Velocity \f$ \mathbf{u} \f$.
 * @param normal Normal \f$ \mathbf{n} \f$.
 * @param imirrored 0 if the particle has not been mirrored, 1 otherwise.
 * @param rmirrored Mirrored position of the particle, \a r if \a imirrored is
 * false (0).
 * @param umirrored Mirrored velocity of the particle, \a u if \a imirrored is
 * false (0).
 * @param nmirrored Mirrored normal of the particle, \a normal if \a imirrored
 * is false (0).
 * @param N Number of particles.
 * @param symmetry_r Position of the symmetry plane.
 * @param symmetry_n Normal of the symmetry plane. It is assumed as normalized.
 */
__kernel void main(const __global vec* r,
                   const __global vec* u,
                   const __global vec* normal,
                   __global int* imirrored,
                   __global vec* rmirrored,
                   __global vec* umirrored,
                   __global vec* nmirrored,
                   unsigned int N,
                   vec symmetry_r,
                   vec symmetry_n)
{
    // find position in global arrays
    unsigned int i = get_global_id(0);
    if(i >= N)
        return;

    // Discard non fluid particles
    if(imove[i] <= 0){
        imirrored[i] = 0;
        rmirrored[i] = r[i];
        umirrored[i] = u[i];
        nmirrored[i] = normal[i];
        return;
    }

    // Get the minimum distance to the plane
    const float dr_n = dot(symmetry_r - r[i], symmetry_n);
    // Discard the particles outside the plane, or far away
    if((dr_n < 0.f) || (dr_n > SUPPORT * H)){
        imirrored[i] = 0;        
        rmirrored[i] = r[i];
        nmirrored[i] = normal[i];
        return;
    }

    // Compute the mirrored particle
    const float du_n = dot(-u[i], symmetry_n);
    const float dn_n = dot(-normal[i], symmetry_n);
    imirrored[i] = 0;        
    rmirrored[i] = r[i] + 2.f * dr_n * symmetry_n;
    umirrored[i] = u[i] + 2.f * du_n * symmetry_n;
    nmirrored[i] = normal[i] + 2.f * dn_n * symmetry_n;
}
