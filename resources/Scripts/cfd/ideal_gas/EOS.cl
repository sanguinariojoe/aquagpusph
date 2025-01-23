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

/** @defgroup ideal_gas Preset to model ideal gases
 *
 * @brief A preset to model ideal gases within @ref cfd preset
 * 
 * @{
 */

/** @file
 * @brief Equation Of State (EOS) for ideal gases.
 *
 * This is an replacement for resources/Scripts/basic/EOS.cl
 */

#ifndef EXCLUDED_PARTICLE
    #define EXCLUDED_PARTICLE(index) (imove[index] <= 0) && (imove[index] != -1)
#endif

#include "resources/Scripts/types/types.h"

/** @brief Ideal gas Equation Of State (EOS) computation
 *
 * The equation of state relates the pressure, density and internal energy
 * fields,
 * \f$ p = \rho \left( \gamma - 1 \right) e \f$
 *
 * @param iset Set of particles index.
 * @param imove Moving flags.
 *   - imove > 0 for regular fluid particles.
 *   - imove = 0 for sensors.
 *   - imove < 0 for boundary elements/particles.
 * @param rho Density \f$ \rho_{n+1/2} \f$.
 * @param eint Internal energy \f$ e_{n+1/2} \f$.
 * @param p Pressure \f$ p_{n+1/2} \f$.
 * @param gamma Heat capacity ratio \f$ \gamma \f$.
 * @param N Number of particles.
 */
__kernel void entry(const __global unsigned int* iset,
                    const __global int* imove,
                    const __global float* rho,
                    const __global float* eint,
                    __global float* p,
                    __global float* T,
                    __constant float* gamma,
                    __constant float* cv,
                    usize N)
{
    usize i = get_global_id(0);
    if(i >= N)
        return;
    if(EXCLUDED_PARTICLE(i))
        return;

    p[i] = (gamma[iset[i]] - 1.0f) * rho[i] * eint[i];
    T[i] = eint[i]/cv[iset[i]];
}

/*
 * @}
 */
