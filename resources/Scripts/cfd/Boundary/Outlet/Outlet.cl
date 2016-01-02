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
 * @brief Vanish the density varaition rate of the particles crossing the outlet
 * infinite plane, destroying (moving out the domain) the particles far away.
 */

#ifndef HAVE_3D
    #include "../../../types/2D.h"
#else
    #include "../../../types/3D.h"
#endif

/** @brief Outlet boundary condition, which is consisting into vanishing the
 * density variation rate of the particles tresspassing the outlet plane, as
 * well as moving out of the computational domain the particles far away, such
 * that they can be used as a buffer.
 *
 * @param imove Moving flags.
 *   - imove > 0 for regular fluid particles.
 *   - imove = 0 for sensors.
 *   - imove < 0 for boundary elements/particles.
 * @param iset Set of particles index.
 * @param r Position \f$ \mathbf{r} \f$.
 * @param u Velocity \f$ \mathbf{u} \f$.
 * @param rho Density \f$ \rho \f$.
 * @param p Pressure \f$ p \f$.
 * @param dudt Velocity rate of change \f$ \frac{d \mathbf{u}}{d t} \f$.
 * @param dudt_in Prev velocity rate of change \f$ \frac{d \mathbf{u}}{d t} \f$.
 * @param drhodt Density rate of change \f$ \frac{d \rho}{d t} \f$.
 * @param drhodt_in Prev density rate of change \f$ \frac{d \rho}{d t} \f$.
 * @param gamma Eq. of state exponent \f$ \gamma \f$.
 * @param refd Density of reference of the fluid \f$ \rho_0 \f$.
 * @param N Number of particles.
 * @param cs Speed of sound \f$ c_s \f$.
 * @param g Gravity acceleration \f$ \mathbf{g} \f$.
 * @param p0 Background pressure \f$ p_0 \f$.
 * @param domain_max Maximum point of the computational domain.
 * @param outlet_r Lower corner of the outlet square.
 * @param outlet_n = Velocity direction of the generated particles.
 * @param outlet_U = Constant outlet velocity magnitude
 * @param outlet_rFS The point where the pressure is the reference one (0 Pa).
 */
__kernel void entry(__global int* imove,
                    __global unsigned int* iset,
                    __global vec* r,
                    __global vec* u,
                    __global float* rho,
                    __global float* p,
                    __global vec* dudt,
                    __global vec* dudt_in,
                    __global float* drhodt,
                    __global float* drhodt_in,
                    __constant float* gamma,
                    __constant float* refd,
                    unsigned int N,
                    float cs,
                    vec g,
                    float p0,
                    vec domain_max,
                    vec outlet_r,
                    vec outlet_n,
                    float outlet_U,
                    vec outlet_rFS)
{
    // find position in global arrays
    unsigned int i = get_global_id(0);
    if(i >= N)
        return;
    if(imove[i] != 1)
        return;

    // Compute the distance to the outlet plane
    const float dist = dot(r[i] - outlet_r, outlet_n);
    if(dist < 0.f)
        return;

    drhodt[i] = 0.f;
    drhodt_in[i] = 0.f;
    dudt[i] = VEC_ZERO;
    dudt_in[i] = VEC_ZERO;
    u[i] = outlet_U * outlet_n;
    p[i] = refd[iset[i]] * dot(g, r[i] - outlet_rFS);
    // Batchelor 1967
    const float prb = cs * cs * refd[iset[i]] / gamma[iset[i]];
    rho[i] = refd[iset[i]] * pow(p[i] / prb + 1.f, 1.f / gamma[iset[i]]);
    p[i] += p0;

    // Destroy the particles far away from the outlet plane
    if(dist > SUPPORT * H)
        r[i] = domain_max + VEC_ONE;
}
