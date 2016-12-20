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

#include "resources/Scripts/types/types.h"

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
__kernel void entry(const __global vec* r,
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
    unsigned int i = get_global_id(0);
    if(i >= N)
        return;

    // Get the minimum distance to the plane
    const float dr_n = dot(symmetry_r - r[i], symmetry_n);
    // Discard the particles far away from the plane
    if(fabs(dr_n) > SUPPORT * H){
        imirrored[i] = 0;        
        rmirrored[i] = r[i];
        nmirrored[i] = normal[i];
        return;
    }

    // Compute the mirrored particle
    const float du_n = dot(-u[i], symmetry_n);
    const float dn_n = dot(-normal[i], symmetry_n);
    imirrored[i] = 1;        
    rmirrored[i] = r[i] + 2.f * dr_n * symmetry_n;
    umirrored[i] = u[i] + 2.f * du_n * symmetry_n;
    nmirrored[i] = normal[i] + 2.f * dn_n * symmetry_n;
}


/** Deflection vector.
 *
 * The deflection vector is defined as follows:
 *
 * \f$ \mathbf{v}(\mathbf{u}, \mathbf{n}) = -2 \left(
 *     \mathbf{u} \cdot \mathbf{n}
 * \right) \mathbf{n} \f$
 *
 * where \f$\mathbf{u}\f$ is the vector to become deflected, \f$\mathbf{n}\f$
 * is the reflection plane normal, and \f$\mathbf{v}\f$ is the deflection
 * vector, which added to the original vector returns its reflected version.
 *
 * @note The input vector should be relative to the symmetry plane. That's
 * important in case of position vectors, from which an arbitrary point of the
 * plane should be substracted
 */
vec_xyz reflection(vec_xyz u, vec_xyz n)
{
    return -2.f * dot(u, n) * n;
}


/** @brief Reflect the particles tresspassing the symmetry plane.
 *
 * The particles tresspassing the symmetry plane should be inserted again inside
 * the fluid domain.
 * 
 * @param r Position \f$ \mathbf{r}_{n+1/2} \f$.
 * @param u Velocity \f$ \mathbf{u}_{n+1/2} \f$.
 * @param dudt Velocity rate of change
 * \f$ \left. \frac{d \mathbf{u}}{d t} \right\vert_{n+1/2} \f$.
 * @param r_in Position \f$ \mathbf{r}_{n} \f$.
 * @param u_in Velocity \f$ \mathbf{u}_{n} \f$.
 * @param dudt_in Velocity rate of change
 * \f$ \left. \frac{d \mathbf{u}}{d t} \right\vert_{n-1/2} \f$.
 * @param N Number of particles.
 * @param symmetry_r Position of the symmetry plane.
 * @param symmetry_n Normal of the symmetry plane. It is assumed as normalized.
 */
__kernel void teleport(__global vec* r,
                       __global vec* normal,
                       __global vec* u,
                       __global vec* dudt,
                       __global vec* r_in,
                       __global vec* u_in,
                       __global vec* dudt_in,
                       unsigned int N,
                       vec symmetry_r,
                       vec symmetry_n)
{
    unsigned int i = get_global_id(0);
    if(i >= N)
        return;

    // Get the minimum distance to the plane
    const float dr_n = dot(symmetry_r - r[i], symmetry_n);
    // Discard the particles at the good side of the plane, or far away from the
    // symmetry plane
    if((dr_n >= 0.f) || (fabs(dr_n) > SUPPORT * H)){
        return;
    }

    // Compute the mirrored particle
    r[i].XYZ += reflection(r[i].XYZ - symmetry_r.XYZ, symmetry_n.XYZ);
    normal[i].XYZ += reflection(normal[i].XYZ, symmetry_n.XYZ);
    u[i].XYZ += reflection(u[i].XYZ, symmetry_n.XYZ);
    dudt[i].XYZ += reflection(dudt[i].XYZ, symmetry_n.XYZ);
    r_in[i].XYZ += reflection(r_in[i].XYZ - symmetry_r.XYZ, symmetry_n.XYZ);
    u_in[i].XYZ += reflection(u_in[i].XYZ, symmetry_n.XYZ);
    dudt_in[i].XYZ += reflection(dudt_in[i].XYZ, symmetry_n.XYZ);
}