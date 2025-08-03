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

/** @addtogroup cfd
 * @{
 */

/** @file
 * @brief Vanish the velocity and desnity rates of variation of the velocity
 * and density for the dummy particles of the inflow.
 */

#include "resources/Scripts/types/types.h"

/** @brief Particles generation at the inflow boundary condition.
 *
 * Particles are generated just when the inflow is starving, i.e. the
 * previously generated layer of particles have moved more than dr. To do that
 * /outlet is extracting the particles from the "buffer", which are the last
 * particles in the sorted list.
 *
 * @param imove Moving flags.
 *   - imove > 0 for regular fluid particles.
 *   - imove = 0 for sensors.
 *   - imove < 0 for boundary elements/particles.
 * @param iset Set of particles index.
 * @param r Position \f$ \mathbf{r} \f$.
 * @param u Velocity \f$ \mathbf{u} \f$.
 * @param dudt Velocity rate of change \f$ \frac{d \mathbf{u}}{d t} \f$.
 * @param rho Density \f$ \rho \f$.
 * @param drhodt Density rate of change \f$ \frac{d \rho}{d t} \f$.
 * @param m Mass \f$ m \f$.
 * @param p Pressure \f$ p \f$.
 * @param refd Density of reference of the fluid \f$ \rho_0 \f$.
 * @param N Number of particles.
 * @param nbuffer Number of buffer particles.
 * @param dt Time step \f$ \Delta t \f$.
 * @param cs Speed of sound \f$ c_s \f$.
 * @param p0 Background pressure \f$ p_0 \f$.
 * @param g Gravity acceleration \f$ \mathbf{g} \f$.
 * @param dr Distance between particles \f$ \Delta r \f$.
 * @param inflow_r Lower corner of the inflow square.
 * @param inflow_ru Square U vector.
 * @param inflow_rv Square V vector.
 * @param inflow_N Number of particles to be generated in each direction.
 * @param inflow_n = Velocity direction of the generated particles.
 * @param inflow_U = Constant inflow velocity magnitude
 * @param inflow_rFS The point where the pressure is the reference one (0 Pa).
 * @param inflow_R Accumulated displacement (to be added to the generation point)
 * @param inflow_starving Is the inflow starving, so we need to feed it?
 */
__kernel void feed(__global int* restrict imove,
                   __global unsigned int* restrict iset,
                   __global vec* restrict r,
                   __global vec* restrict u,
                   __global vec* restrict dudt,
                   __global float* restrict rho,
                   __global float* restrict drhodt,
                   __global float* restrict m,
                   __global float* restrict p,
                   __constant float* restrict refd,
                   usize N,
                   usize nbuffer,
                   float dt,
                   float cs,
                   float p0,
                   vec g,
                   float dr,
                   vec inflow_r,
                   vec inflow_ru,
                   vec inflow_rv,
                   svec2 inflow_N,
                   vec inflow_n,
                   float inflow_U,
                   vec inflow_rFS,
                   float inflow_R,
                   int inflow_starving)
{
    // find position in global arrays
    const usize i = get_global_id(0);
    if(inflow_starving == 0)
        return;
    if((i >= nbuffer) || (i >= (inflow_N.x * inflow_N.y))){
        // Either the thread has not a buffer particle to consume or such buffer
        // particle is not required
        return;
    }
    const usize ii = i + N - nbuffer;

    // Compute the generation point
    #ifndef HAVE_3D
        const float u_fac = ((float)i + 0.5f) / inflow_N.x;
        const float v_fac = 0.f;
    #else
        const usize u_id = i % inflow_N.x;
        const usize v_id = i / inflow_N.x;
        const float u_fac = ((float)u_id + 0.5f) / inflow_N.x;
        const float v_fac = ((float)v_id + 0.5f) / inflow_N.y;
    #endif
    r[ii] = inflow_r + u_fac * inflow_ru + v_fac * inflow_rv
            + (inflow_R - SUPPORT * H - 0.5f * dr) * inflow_n;

    // Set the particle data
    imove[ii] = 1;
    dudt[ii] = VEC_ZERO;
    drhodt[ii] = 0.f;
    u[ii] = inflow_U * inflow_n;
    p[ii] = refd[iset[ii]] * dot(g, r[ii] - inflow_rFS);
    #ifdef HAVE_3D
        m[ii] = refd[iset[ii]] * dr * dr * dr;
    #else
        m[ii] = refd[iset[ii]] * dr * dr;
    #endif
    // reversed EOS
    rho[ii] = refd[iset[ii]] + p[ii] / (cs * cs);
    p[ii] += p0;
}

/** @brief Vanish the velocity and desnity rates of variation of the velocity
 * and density for the dummy particles of the inflow.
 *
 * @param imove Moving flags.
 *   - imove > 0 for regular fluid particles.
 *   - imove = 0 for sensors.
 *   - imove < 0 for boundary elements/particles.
 * @param r Position \f$ \mathbf{r} \f$.
 * @param u Velocity \f$ \mathbf{u} \f$.
 * @param dudt Velocity rate of change \f$ \frac{d \mathbf{u}}{d t} \f$.
 * @param drhodt Density rate of change \f$ \frac{d \rho}{d t} \f$.
 * @param N Number of particles.
 * @param inflow_r Lower corner of the inflow square.
 * @param inflow_U Velocity magnitude of the generated particles.
 * @param inflow_n Velocity direction of the generated particles.
 */
__kernel void rates(__global int* restrict imove,
                    __global vec* restrict r,
                    __global vec* restrict u,
                    __global vec* restrict dudt,
                    __global float* restrict drhodt,
                    usize N,
                    vec inflow_r,
                    float inflow_U,
                    vec inflow_n)
{
    // find position in global arrays
    const usize i = get_global_id(0);
    if(i >= N)
        return;
    if(imove[i] != 1)
        return;

    // Discard the particles already passed through the inflow
    if(dot(r[i] - inflow_r, inflow_n) > 0.f)
        return;

    u[i] = inflow_U * inflow_n;
    dudt[i] = VEC_ZERO;
    drhodt[i] = 0.f;
}

/*
 * @}
 */
