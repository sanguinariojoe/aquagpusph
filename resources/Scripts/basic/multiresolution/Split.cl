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

/** @addtogroup basic
 * @{
 */

/** @file
 *  @brief Splitting particles methods
 */

#ifndef HAVE_3D
    #include "../../types/2D.h"
#else
    #include "../../types/3D.h"
#endif

#ifndef HAVE_3D
    #define N_DAUGHTER 4
#else
    #define N_DAUGHTER 8
#endif

/** @brief Check and store wether a particle should become split or not.
 *
 * A particle should be split if the target refinement level is bigger than
 * its current one.
 *
 * @param imove Moving flags.
 *   - imove > 0 for regular fluid/solid particles.
 *   - imove = 0 for sensors.
 *   - imove < 0 for boundary elements/particles.
 * @param ilevel0 Level of refinement of the particle, by construction.
 * @param ilevel Current refinement level of the particle.
 * @param level Target refinement level of the particle.
 * @param isplit 0 if the particle should not become split, 1 otherwise
 * @param N Number of particles.
 * @see P.N. Sun, A. Colagrossi, S. Marrone, A.M. Zhang. Multi-resolution
 * delta-SPH model. 2016
 */
__kernel void check_split(__global const int* imove,
                          __global const unsigned int* ilevel0,
                          __global const unsigned int* ilevel,
                          __global const unsigned int* level,
                          __global unsigned int* isplit,
                          unsigned int N)
{
    unsigned int i = get_global_id(0);
    if(i >= N)
        return;

    isplit[i] = 0;
    if((imove[i] <= 0) ||           // Neglect boundary elements/particles
       (ilevel[i] != ilevel0[i]) || // It is already a mother particle
       (ilevel[i] >= level[i])) {   // Not asked to split
        return;
    }

    isplit[i] = 1;
}


/** @brief Split a particle in a set of daughter particles.
 *
 * The daughter particles will be a duplicate of the mother split particle, but
 * conveniently modifying the position and the mass (m0 and gamma_m).
 *
 * @param imove Moving flags.
 *   - imove > 0 for regular fluid/solid particles.
 *   - imove = 0 for sensors.
 *   - imove < 0 for boundary elements/particles.
 * @param iset Index of the set of particles.
 * @param isplit 0 if the particle should not become split, 1 otherwise.
 * @param split_invperm Permutation to find the index of the particle in the
 * list of particles to become split.
 * @param ilevel0 Level of refinement of the particle, by construction.
 * @param ilevel Current refinement level of the particle.
 * @param level Target level.
 * @param m0 Original mass, before applying the gamma factor, \f$ m_0 \f$
 * @param gamma_m Mass multiplier \f$ \gamma_m \f$.
 * @param r Position \f$ \mathbf{r} \f$.
 * @param u Velocity \f$ \mathbf{u} \f$.
 * @param dudt Velocity rate of change \f$ \frac{d \mathbf{u}}{d t} \f$.
 * @param rho Density \f$ \rho \f$.
 * @param drhodt Density rate of change \f$ \frac{d \rho}{d t} \f$.
 * @param N Number of particles.
 * @see P.N. Sun, A. Colagrossi, S. Marrone, A.M. Zhang. Multi-resolution
 * delta-SPH model. 2016
 */
__kernel void generate(__global int* imove,
                       __global int* iset,
                       __global unsigned int* isplit,
                       __global unsigned int* split_invperm,
                       __global unsigned int* ilevel0,
                       __global unsigned int* ilevel,
                       __global unsigned int* level,
                       __global float* m0,
                       __global float* gamma_m,
                       __global vec* r,
                       __global vec* u,
                       __global vec* dudt,
                       __global float* rho,
                       __global float* drhodt,
                       unsigned int N)
{
    unsigned int i = get_global_id(0);
    if(i >= N)
        return;

    // Neglect boundary elements/particles
    if(imove[i] <= 0)
        return;

    // Check whether the particle should become split or not
    unsigned int j = split_invperm[i];
    if(isplit[j] != 1)
        return;

    // Compute the index of the first buffer particle to steal. Take care, the
    // radix sort is storing the particles to become split at the end of the
    // list 
    j = N - (N - j) * N_DAUGHTER;
    #ifdef __MULTIRESOLUTION_N_SKIP__
        j -= __MULTIRESOLUTION_N_SKIP__;
    #endif

    // Insert the daughters
    ilevel[i]++;
    const float dr = 0.25f * pow(m0[i] / rho[i], 1 / DIMS);
    for(int ci = -1; ci <= 1; ci += 2) {
        for(int cj = -1; cj <= 1; cj += 2) {
#ifdef HAVE_3D
            for(int ck = -1; ck <= 1; ck += 2) {
                vec_xyz r_j = r[i].XYZ + dr * (float3)(ci, cj, ck);
#else
                vec_xyz r_j = r[i].XYZ + dr * (float2)(ci, cj);
#endif
                // Set the new particle properties
                ilevel0[j] = ilevel[i];
                ilevel[j] = ilevel[i];
                level[j] = ilevel[i];  // Needed to avoid coalescing
                imove[j] = imove[i];
                iset[j] = iset[i];                
                m0[j] = m0[i] / N_DAUGHTER;
                gamma_m[j] = 1.f - gamma_m[i];
                r[j] = r_j;
                u[j] = u[i];
                dudt[j] = dudt[i];
                rho[j] = rho[i];
                drhodt[j] = drhodt[i];
                // Move to the next buffer particle
                j--;
#ifdef HAVE_3D
            }
#endif
        }
    }
}

/*
 * @}
 */