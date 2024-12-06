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
 * @brief Equation Of State (EOS) computation
 */

#ifndef EXCLUDED_PARTICLE
    /** @brief Condition to exclude a particle from the EOS computation
     * 
     * By default all the boundary elements are excluded, excepting the fixed
     * particles (imove = -1)
     * @note Redefining this macro this OpenCL script can be recicled for
     * different particle types.
     */
    #define EXCLUDED_PARTICLE(index) (imove[index] <= 0) && (imove[index] != -1)
#endif

#include "resources/Scripts/types/types.h"

/** @brief Stiffness Equation Of State (EOS) computation
 *
 * The equation of state relates the pressure and the density fields,
 * \f$ p = p_0 + c_s^2 \left(\rho - \rho_0 \right) \f$
 *
 * @param iset Set of particles index.
 * @param imove Moving flags.
 *   - imove > 0 for regular fluid particles.
 *   - imove = 0 for sensors.
 *   - imove < 0 for boundary elements/particles.
 * @param rho Density \f$ \rho_{n+1/2} \f$.
 * @param p Pressure \f$ p_{n+1/2} \f$.
 * @param refd Density of reference of the fluid \f$ \rho_0 \f$.
 * @param N Number of particles.
 * @param cs Speed of sound \f$ c_s \f$.
 * @param p0 Background pressure \f$ p_0 \f$.
 */

#define gamma 1.44f

__kernel void entry(const __global unsigned int* iset,
                    const __global int* imove,
                    const __global vec* r,
                    const __global float* rho,
                    __global float* p,
                    const __global float* eee,
                    const __global vec* u,
                    unsigned int N)
{
    unsigned int i = get_global_id(0);
    if(i >= N)
        return;
    if(EXCLUDED_PARTICLE(i))
        return;

//     const float epsilon = e[i] - 0.5f * rho[i] * MODULE(u[i]);
//
//     p[i] = (gamma-1.0f)*rho[i]*epsilon;

    //const float epsilon = e[i] - 0.5f * rho[i] * MODULE(u[i]);
    //if ((r[i].x>0) && (r[i].y>0) &&  (r[i].x<0.001f) && (r[i].y<0.001f))
    //    printf("Density(%f, %f)= %f Energy(%f, %f)= %f\n",r[i].x, r[i].y, rho[i], r[i].x, r[i].y, eee[i]);
    //exit(0);
    p[i] = (gamma-1.0f)*rho[i]*eee[i];
    //p[i] = (gamma-1.0f)*rho[i];
}

/*
 * @}
 */
