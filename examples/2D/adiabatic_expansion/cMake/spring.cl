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
 *  @brief Spring related tools
 */

#include "resources/Scripts/types/types.h"

/** @brief Sort all the spring related arrays
 * @param r0_in Unsorted initial position \f$ \mathbf{r}(t=0) \f$.
 * @param r0 Sorted initial position \f$ \mathbf{r}(t=0) \f$.
 * @param id_sorted Permutations list from the unsorted space to the sorted
 * one.
 * @param N Number of particles.
 */
__kernel void sort(const __global vec *r0_in, __global vec *r0,
                   const __global unit *id_sorted,
                   unsigned int N)
{
    uint i = get_global_id(0);
    if(i >= N)
        return;

    const uint i_out = id_sorted[i];
    r0[i_out] = r0_in[i];
}

/** @brief Drop the forces on all the boundaries but the piston
 * @param imove Moving flags.
 *   - imove > 0 for regular fluid particles.
 *   - imove = 0 for sensors.
 *   - imove < 0 for boundary elements/particles.
 * @param normal Sorted normal \f$ \mathbf{n} \f$.
 * @param force_p Force on the boundary elements.
 * @param N Number of particles.
 */
__kernel void force(const __global int* imove,
                    const __global vec* normal,
                    __global vec *force_p,
                    unsigned int N)
{
    uint i = get_global_id(0);
    if(i >= N)
        return;
    if(imove[i] != -3)
        return;

    if (dot(normal[i], (vec2)(1.f, 0.f)) < 0.9f)
        force_p[i] = VEC_ZERO;
}

/** @brief Set the velocity of the piston particles
 *
 * This method is used within the semi-implicit midpoint iterator
 * @param imove Moving flags.
 *   - imove > 0 for regular fluid particles.
 *   - imove = 0 for sensors.
 *   - imove < 0 for boundary elements/particles.
 * @param normal Sorted normal \f$ \mathbf{n} \f$.
 * @param u Velocity \f$ \mathbf{u}_{n+1/2} \f$.
 * @param dxdt Piston velocity \f$ \frac{d x}{d t} \f$.
 * @param N Number of particles.
 */
__kernel void piston_u(const __global int* imove,
                       const __global vec* normal,
                       __global vec* u,
                       float dxdt,
                       unsigned int N)
{
    uint i = get_global_id(0);
    if(i >= N)
        return;
    if(imove[i] != -3)
        return;

    if (dot(normal[i], (vec2)(1.f, 0.f)) < 0.9f)
        return;

    u[i].x = dxdt;
}


/** @brief Set the position and velocity of the piston particles
 * @param imove Moving flags.
 *   - imove > 0 for regular fluid particles.
 *   - imove = 0 for sensors.
 *   - imove < 0 for boundary elements/particles.
 * @param normal Sorted normal \f$ \mathbf{n} \f$.
 * @param r Position \f$ \mathbf{r}_{n+1/2} \f$.
 * @param u Velocity \f$ \mathbf{u}_{n+1/2} \f$.
 * @param x Piston position \f$ x \f$.
 * @param dxdt Piston velocity \f$ \frac{d x}{d t} \f$.
 * @param N Number of particles.
 */
__kernel void piston(const __global int* imove,
                     const __global vec* normal,
                     __global vec* r,
                     __global vec* u,
                     float x,
                     float dxdt,
                     unsigned int N)
{
    uint i = get_global_id(0);
    if(i >= N)
        return;
    if(imove[i] != -3)
        return;

    if (dot(normal[i], (vec2)(1.f, 0.f)) < 0.9f)
        return;

    r[i].x = x;
    u[i].x = dxdt;
}

/** @brief Set the position and mass of the top and bottom walls
 * @param imove Moving flags.
 *   - imove > 0 for regular fluid particles.
 *   - imove = 0 for sensors.
 *   - imove < 0 for boundary elements/particles.
 * @param normal Sorted normal \f$ \mathbf{n} \f$.
 * @param r0 Initial Position \f$ \mathbf{r}(t=0) \f$.
 * @param r Position \f$ \mathbf{r}_{n+1/2} \f$.
 * @param m Area \f$ s \f$.
 * @param x Piston position \f$ x \f$.
 * @param x0 Initial position of the piston
 * @param L Length of the box at equilibrium
 * @param nx Number of particles on the top or bottom boundaries.
 * @param N Number of particles.
 */
__kernel void bounds(const __global int* imove,
                     const __global vec* normal,
                     const __global vec* r0,
                     __global vec* r,
                     __global float* m,
                     float x,
                     float x0,
                     float L,
                     unsigned int nx,
                     unsigned int N)
{
    uint i = get_global_id(0);
    if(i >= N)
        return;
    if(imove[i] != -3)
        return;

    if (fabs(dot(normal[i], (vec2)(0.f, 1.f))) < 0.9f)
        return;

    const float l0 = L + x0;
    const float l = L + x;
    r[i].x = -L + (L + r0[i].x) / l0 * l;
    m[i] = l / nx;
}
