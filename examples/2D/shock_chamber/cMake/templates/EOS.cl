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
 * @brief Set the velocity according to the initialization stage.
 */

#ifndef HAVE_3D
    #include "@RESOURCES_DIR@/Scripts/types/2D.h"
#else
    #include "@RESOURCES_DIR@/Scripts/types/3D.h"
#endif

/** @brief Improved Euler time integration scheme predictor stage
 *
 * Here the internal energy is integrated
 *
 * @param imove Moving flags.
 *   - imove > 0 for regular fluid/solid particles.
 *   - imove = 0 for sensors.
 *   - imove < 0 for boundary elements/particles.
 * @param e Internal energy \f$ e_{n+1} \f$.
 * @param dedt Internal energy rate of change
 * \f$ \left. \frac{d e}{d t} \right\vert_{n+1} \f$.
 * @param e_in Internal energy \f$ e_{n+1/2} \f$.
 * @param dedt_in Internal energy rate of change
 * \f$ \left. \frac{d e}{d t} \right\vert_{n+1/2} \f$.
 * @param N Number of particles.
 * @param dt Time step \f$ \Delta t \f$.
 */
__kernel void predictor(__global const int* imove,
                        __global const float* e,
                        __global const float* dedt,
                        __global float* e_in,
                        __global float* dedt_in,
                        unsigned int N,
                        float dt)
{
    unsigned int i = get_global_id(0);
    if(i >= N)
        return;

    float DT = dt;
    if(imove[i] <= 0)
        DT = 0.f;

    dedt_in[i] = dedt[i];
    e_in[i] = e[i] + DT * dedt[i];
}

/** @brief Sort the internal energy, and its variation rate.
 *
 * @param e_in Unsorted internal energy \f$ e \f$.
 * @param e Sorted internal energy \f$ e \f$.
 * @param dedt_in Unsorted internal energy rate of change \f$ \frac{d e}{d t} \f$.
 * @param dedt Sorted internal energy rate of change \f$ \frac{d e}{d t} \f$.
 * @param id_sorted Permutations list from the unsorted space to the sorted
 * one.
 * @param N Number of particles.
 */
__kernel void sort(const __global float *e_in, __global float *e,
                   const __global float *dedt_in, __global float *dedt,
                   const __global unit *id_sorted,
                   unsigned int N)
{
    uint i = get_global_id(0);
    if(i >= N)
        return;

    const uint i_out = id_sorted[i];

    e[i_out] = e_in[i];
    dedt[i_out] = dedt_in[i];
}

/** @brief Equation Of State (EOS) computation
 *
 * The equation of state relates the pressure, the density, and the internal
 * energy,
 * \f$ p = \left( \gamma - 1 \right) \rho e \f$
 *
 * @param rho Density \f$ \rho_{n+1/2} \f$.
 * @param e Internal energy \f$ e_{n+1/2} \f$.
 * @param p Pressure \f$ p_{n+1/2} \f$.
 * @param N Number of particles.
 * @param gamma EOS \f$ \gamma \f$ factor.
 */
__kernel void eos(__global const float* rho,
                  __global const float* e,
                  __global float* p,
                  unsigned int N,
                  float gamma)
{
    unsigned int i = get_global_id(0);
    if(i >= N)
        return;

    p[i] = (gamma - 1.f) * rho[i] * e[i];
}

/** @brief Internal energy variation rate computation.
 *
 * The variation rate of the internal energy can be deduced from the EOS, and
 * viceversa
 *
 * @param imove Moving flags.
 *   - imove > 0 for regular fluid particles.
 *   - imove = 0 for sensors.
 *   - imove < 0 for boundary elements/particles.
 * @param p Pressure \f$ p_{n+1/2} \f$.
 * @param rho Density \f$ \rho_{n+1/2} \f$.
 * @param drhodt Density rate of change
 * \f$ \left. \frac{d \rho}{d t} \right\vert_{n+1/2} \f$.
 * @param dedt Internal energy rate of change
 * \f$ \left. \frac{d e}{d t} \right\vert_{n+1/2} \f$.
 * @param N Number of particles.
 */
__kernel void rates(__global const int* imove,
                    __global const float* p,
                    __global const float* rho,
                    __global const float* drhodt,
                    __global float* dedt,
                    unsigned int N)
{
    unsigned int i = get_global_id(0);
    if(i >= N)
        return;
    if(imove[i] != 1)
        return;

    dedt[i] = p[i] / (rho[i] * rho[i]) * drhodt[i];
}

/** @brief Improved Euler time integration scheme corrector stage
 *
 * Here the internal energy is integrated
 *
 * @param imove Moving flags.
 *   - imove > 0 for regular fluid/solid particles.
 *   - imove = 0 for sensors.
 *   - imove < 0 for boundary elements/particles.
 * @param e Internal energy \f$ e_{n+1/2} \f$.
 * @param dedt Internal energy rate of change
 * \f$ \left. \frac{d e}{d t} \right\vert_{n+1/2} \f$.
 * @param e_in Internal energy \f$ e_{n-1/2} \f$.
 * @param dedt_in Internal energy rate of change
 * \f$ \left. \frac{d e}{d t} \right\vert_{n-1/2} \f$.
 * @param N Number of particles.
 * @param dt Time step \f$ \Delta t \f$.
 */
__kernel void corrector(__global int* imove,
                        __global float* e,
                        __global const float* dedt,
                        __global float* e_in,
                        __global float* dedt_in,
                        unsigned int N,
                        float dt)
{
    unsigned int i = get_global_id(0);
    if(i >= N)
        return;

    float DT = 0.5f * dt;
    if(imove[i] <= 0)
        DT = 0.f;

    // Integrate it
    e[i] += DT * (dedt[i] - dedt_in[i]);
    e_in[i] = e[i];
    dedt_in[i] = dedt[i];
}