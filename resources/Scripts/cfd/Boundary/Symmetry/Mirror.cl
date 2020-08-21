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

/** @brief Remove the particles at the other side of the mirror.
 * 
 * @param imove Moving flags.
 *   - imove > 0 for regular fluid particles.
 *   - imove = 0 for sensors.
 *   - imove < 0 for boundary elements/particles.
 * @param r Position \f$ \mathbf{r} \f$.
 * @param N Number of particles.
 * @param symmetry_r Position of the symmetry plane.
 * @param symmetry_n Normal of the symmetry plane. It is assumed as normalized.
 * @param domain_max Top-right-back corner of the computational domain.
 */
__kernel void drop(__global int* imove,
                   __global vec* r,
                   unsigned int N,
                   vec symmetry_r,
                   vec symmetry_n,
                   vec domain_max)
{
    unsigned int i = get_global_id(0);
    if(i >= N)
        return;
    if(imove[i] <= -255)
        return;

    const float dr_n = dot(r[i] - symmetry_r, symmetry_n);
    if(dr_n >= 0.f){
        r[i] = domain_max + VEC_ONE;
        imove[i] = -256;
    }
}

/** @brief Detect the particles to be mirrored.
 *
 * The mirroring particles (the ones close enough to the symmetry plane) will be
 * marked with \a imirror = 1.
 * 
 * @param imove Moving flags.
 *   - imove > 0 for regular fluid particles.
 *   - imove = 0 for sensors.
 *   - imove < 0 for boundary elements/particles.
 * @param r_in Position \f$ \mathbf{r} \f$.
 * @param imirror 0 if the particle has not been mirrored, 1 otherwise.
 * @param N Number of particles.
 * @param symmetry_r Position of the symmetry plane.
 * @param symmetry_n Normal of the symmetry plane. It is assumed as normalized.
 */
__kernel void detect(const __global int* imove,
                     const __global vec* r_in,
                     __global unsigned int* imirror,
                     unsigned int N,
                     vec symmetry_r,
                     vec symmetry_n)
{
    unsigned int i = get_global_id(0);
    if(i >= N)
        return;
    if(imove[i] <= -255) {
        imirror[i] = 0;
        return;
    }

    // Get the minimum distance to the plane
    const float dr_n = dot(symmetry_r - r_in[i], symmetry_n);
    // Discard the particles far away from the plane
    if(fabs(dr_n) <= SUPPORT * H) {
        imirror[i] = 1;
        return;
    }
    imirror[i] = 0;
}

/** Reflection vector.
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

/** @brief Mirror the particles marked with a flag \a imirror = 1.
 *
 * @param imove Moving flags.
 *   - imove > 0 for regular fluid/solid particles.
 *   - imove = 0 for sensors.
 *   - imove < 0 for boundary elements/particles.
 * @param iset Index of the set of particles.
 * @param imirror 0 if the particle should not be mirrored, 1 otherwise.
 * @param imirror_invperm Permutation to find the index of the particle in the
 * list of particles to become split.
 * @param m Mass, \f$ m \f$.
 * @param normal Normal, \f$ \mathbf{n} \f$.
 * @param tangent Tangent, \f$ \mathbf{t} \f$.
 * @param r_in Position \f$ \mathbf{r} \f$.
 * @param u_in Velocity \f$ \mathbf{u} \f$.
 * @param dudt_in Velocity rate of change \f$ \frac{d \mathbf{u}}{d t} \f$.
 * @param rho_in Density \f$ \rho \f$.
 * @param drhodt_in Density rate of change \f$ \frac{d \rho}{d t} \f$.
 * @param N Number of particles.
 * @param nbuffer Number of available buffer particles.
 * @param symmetry_r Position of the symmetry plane.
 * @param symmetry_n Normal of the symmetry plane. It is assumed as normalized.
 */
__kernel void feed(__global int* imove,
                   __global int* iset,
                   __global unsigned int* imirror,
                   __global unsigned int* imirror_invperm,
                   __global float* m,
                   __global vec* normal,
                   __global vec* tangent,
                   __global vec* r_in,
                   __global vec* u_in,
                   __global vec* dudt_in,
                   __global float* rho_in,
                   __global float* drhodt_in,
                   unsigned int N,
                   unsigned int nbuffer,
                   vec symmetry_r,
                   vec symmetry_n)
{
    unsigned int i = get_global_id(0);
    if(i >= N)
        return;
    if(imove[i] <= -255)
        return;

    // Check whether the particle should become mirrored or not
    const unsigned int j = imirror_invperm[i];
    if(imirror[j] != 1)
        return;

    // Compute the index of the first buffer particle to steal. Take care, the
    // radix sort is storing the particles to become split at the end of the
    // list
    const unsigned int i0 = N - nbuffer;
    unsigned int ii = i0 + (N - j - 1);
    // Check that there are buffer particles enough
    if(ii >= N){
        // PROBLEMS! This particle cannot be mirrored because we have not buffer
        // particles enough to create the children.
        // How the hell this happened??
        return;
    }

    // Create the mirrored particle
    imove[ii] = imove[i];
    iset[ii] = iset[i];
    m[ii] = m[i];
    rho_in[ii] = rho_in[i];
    drhodt_in[ii] = drhodt_in[i];
    normal[ii].XYZ = normal[i].XYZ + reflection(normal[i].XYZ,
                                                symmetry_n.XYZ);
    tangent[ii].XYZ = tangent[i].XYZ + reflection(tangent[i].XYZ,
                                                  symmetry_n.XYZ);
    r_in[ii].XYZ = r_in[i].XYZ + reflection(r_in[i].XYZ - symmetry_r.XYZ,
                                            symmetry_n.XYZ);
    u_in[ii].XYZ = u_in[i].XYZ + reflection(u_in[i].XYZ,
                                            symmetry_n.XYZ);
    dudt_in[ii].XYZ = dudt_in[i].XYZ + reflection(dudt_in[i].XYZ,
                                                  symmetry_n.XYZ);
}
