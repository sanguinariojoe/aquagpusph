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
 * @brief Particles generation at the Inlet (i.e. inflow) boundary condition.
 */

#ifndef HAVE_3D
    #include "../../types/2D.h"
#else
    #include "../../types/3D.h"
#endif

/** @brief Particles generation at the inlet.
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
 * @param gamma Eq. of state exponent \f$ \gamma \f$.
 * @param refd Density of reference of the fluid \f$ \rho_0 \f$.
 * @param N Number of particles.
 * @param dt Time step \f$ \Delta t \f$.
 * @param cs Speed of sound \f$ c_s \f$.
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
__kernel void main(__global int* imove,
                   __global unsigned int* iset,
                   __global vec* r,
                   __global vec* u,
                   __global vec* dudt,
                   __global float* rho,
                   __global float* drhodt,
                   __global float* m,
                   __global float* p,
                   __constant float* gamma,
                   __constant float* refd,
                   unsigned int N,
                   float dt,
                   float cs,
                   vec g,
                   float dr,
                   vec inlet_r,
                   vec inlet_ru,
                   vec inlet_rv,
                   uivec2 inlet_N,
                   vec inlet_n,
                   float inlet_U,
                   vec inlet_rFS,
                   float inlet_R,
                   int inlet_starving)
{
    // find position in global arrays
    const unsigned int i = get_global_id(0);
    if(inlet_starving == 0)
        return;
    if(i >= N)
        return;
    const unsigned int i0 = N - (inlet_N.x * inlet_N.y);
    if(i < i0)
        return;

    // Compute the generation point
    const unsigned int j = i - i0;
    #ifndef HAVE_3D
        const float u_fac = ((float)j + 0.5f) / inlet_N.x;
        const float v_fac = 0.f;
    #else
        const unsigned int u_id = j % inlet_N.x;
        const unsigned int v_id = j / inlet_N.x;
        const float u_fac = ((float)u_id + 0.5f) / inlet_N.x;
        const float v_fac = ((float)v_id + 0.5f) / inlet_N.y;
    #endif
    r[i] = inlet_r + u_fac * inlet_ru + v_fac * inlet_rv
           + (inlet_R - SUPPORT * H + 0.5f * dr) * inlet_n;
    
    // Set the particle data
    imove[i] = 1;
    dudt[i] = VEC_ZERO;
    drhodt[i] = 0.f;
    u[i] = inlet_U * inlet_n;
    p[i] = refd[iset[i]] * dot(g, inlet_rFS - r[i]);
    // Batchelor 1967
    const float prb = cs * cs * refd[iset[i]] / gamma[iset[i]];
    rho[i] = refd[iset[i]] * pow(p[i] / prb + 1.f, 1.f / gamma[iset[i]]);
}
