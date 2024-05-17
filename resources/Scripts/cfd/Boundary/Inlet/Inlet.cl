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
 * and density for the dummy particles of the inlet.
 */

#include "resources/Scripts/types/types.h"

/** @brief Particles generation at the inlet (i.e. inflow) boundary condition.
 *
 * Particles are generated just when the inlet is starving, i.e. the previously
 * generated layer of particles have moved more than dr. To do that inlet is
 * extracting the particles from the "buffer", which are the last particles in
 * the sorted list.
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
 * @param inlet_r Lower corner of the inlet square.
 * @param inlet_ru Square U vector.
 * @param inlet_rv Square V vector.
 * @param inlet_N Number of particles to be generated in each direction.
 * @param inlet_n = Velocity direction of the generated particles.
 * @param inlet_U = Constant inlet velocity magnitude
 * @param inlet_rFS The point where the pressure is the reference one (0 Pa).
 * @param inlet_R Accumulated displacement (to be added to the generation point)
 * @param inlet_starving Is the inlet starving, so we need to feed it?
 */
__kernel void feed(__global int* imove,
                   __global unsigned int* iset,
                   __global vec* r,
                   __global vec* u,
                   __global vec* dudt,
                   __global float* rho,
                   __global float* drhodt,
                   __global float* m,
                   __global float* p,
                   __constant float* refd,
                   usize N,
                   usize nbuffer,
                   float dt,
                   float cs,
                   float p0,
                   vec g,
                   float dr,
                   vec inlet_r,
                   vec inlet_ru,
                   vec inlet_rv,
                   svec2 inlet_N,
                   vec inlet_n,
                   float inlet_U,
                   vec inlet_rFS,
                   float inlet_R,
                   int inlet_starving)
{
    // find position in global arrays
    const usize i = get_global_id(0);
    if(inlet_starving == 0)
        return;
    if((i >= nbuffer) || (i >= (inlet_N.x * inlet_N.y))){
        // Either the thread has not a buffer particle to consume or such buffer
        // particle is not required
        return;
    }
    const usize ii = N - nbuffer + i;

    // Compute the generation point
    #ifndef HAVE_3D
        const float u_fac = ((float)i + 0.5f) / inlet_N.x;
        const float v_fac = 0.f;
    #else
        const usize u_id = i % inlet_N.x;
        const usize v_id = i / inlet_N.x;
        const float u_fac = ((float)u_id + 0.5f) / inlet_N.x;
        const float v_fac = ((float)v_id + 0.5f) / inlet_N.y;
    #endif
    r[ii] = inlet_r + u_fac * inlet_ru + v_fac * inlet_rv
            + (inlet_R - SUPPORT * H - 0.5f * dr) * inlet_n;

    // Set the particle data
    imove[ii] = 1;
    dudt[ii] = VEC_ZERO;
    drhodt[ii] = 0.f;
    u[ii] = inlet_U * inlet_n;
    p[ii] = refd[iset[ii]] * dot(g, r[ii] - inlet_rFS);
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
 * and density for the dummy particles of the inlet.
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
 * @param inlet_r Lower corner of the inlet square.
 * @param inlet_U Velocity magnitude of the generated particles.
 * @param inlet_n Velocity direction of the generated particles.
 */
__kernel void rates(__global int* imove,
                    __global vec* r,
                    __global vec* u,
                    __global vec* dudt,
                    __global float* drhodt,
                    usize N,
                    vec inlet_r,
                    float inlet_U,
                    vec inlet_n)
{
    // find position in global arrays
    const usize i = get_global_id(0);
    if(i >= N)
        return;
    if(imove[i] != 1)
        return;

    // Discard the particles already passed through the inlet
    if(dot(r[i] - inlet_r, inlet_n) > 0.f)
        return;

    u[i] = inlet_U * inlet_n;
    dudt[i] = VEC_ZERO;
    drhodt[i] = 0.f;
}

/*
 * @}
 */
