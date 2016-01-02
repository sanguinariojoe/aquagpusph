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

/** @defgroup cfd CFD simulations preset
 *
 * @brief Preset of tools to perform CFD simulations
 *
 * @see @ref basic
 *
 * @{
 */

/** @file
 * @brief Equation Of State (EOS) computation
 */

#ifndef HAVE_3D
    #include "../types/2D.h"
#else
    #include "../types/3D.h"
#endif

/** @brief Equation Of State (EOS) computation
 *
 * Pressure computation using the density data,
 * \f$ p = p_0 + \frac{c_s^2 \rho_0}{\gamma} \left(
 *     \left( \frac{\rho}{\rho_0} \right)^\gamma - 1
 * \right) \f$
 *
 * @param iset Set of particles index.
 * @param rho_in Density \f$ \rho_{n+1/2} \f$.
 * @param p_in Pressure \f$ \left. p \right\vert_{n+1/2} \f$.
 * @param gamma Eq. of state exponent \f$ \gamma \f$.
 * @param refd Density of reference of the fluid \f$ \rho_0 \f$.
 * @param N Number of particles.
 * @param cs Speed of sound \f$ c_s \f$.
 * @param p0 Background pressure \f$ p_0 \f$.
 */
__kernel void entry(__global unsigned int* iset,
                    __global float* rho_in,
                    __global float* p_in,
                    __constant float* gamma,
                    __constant float* refd,
                    unsigned int N,
                    float cs,
                    float p0)
{
    unsigned int i = get_global_id(0);
    if(i >= N)
        return;

    // Batchelor 1967
	const float ddenf = rho_in[i] / refd[iset[i]];
	const float prb = cs * cs * refd[iset[i]] / gamma[iset[i]];
	p_in[i] = p0 + prb * (pow(ddenf, gamma[iset[i]]) - 1.f);
}

/*
 * @}
 */