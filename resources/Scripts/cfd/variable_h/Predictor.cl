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
 *  @brief Kernel length computation
 */

#include "resources/Scripts/types/types.h"

/** @brief Recompute the kernel length for each particle
 *
 * \f$ h_i = \nu \left( \frac{m_i}{\rho_i} \right)^{1 / d} \f$
 *
 * where \f$ \nu \f$ is the hfac ratio, and \f$ d \f$ is the dimension, 2 and 3
 * for 2D and 3D respectively.
 *
 * @param imove Moving flags.
 *   - imove > 0 for regular fluid/solid particles.
 *   - imove = 0 for sensors.
 *   - imove < 0 for boundary elements/particles.
 * @param rho_in Density \f$ \rho \f$.
 * @param m Mass \f$ m \f$.
 * @param h_var_in variable kernel lenght \f$ h \f$.
 * @param N Number of particles.
 * @param hfac Ratio between the kernel length and the distance between
 * particles.
 * @see Iason Zisis, Bas van der Linden, Christina Giannopapa, Barry Koren. On
 * the derivation of SPH schemes for shocks through inhomogeneous media. Int.
 * Jnl. of Multiphysics. Vol. 9, Number 2. 2015
 */
__kernel void entry(__global const int* imove,
                    __global const float* rho_in,
                    __global const float* m,
                    __global float* h_var_in,
                    unsigned int N,
                    float hfac)
{
    unsigned int i = get_global_id(0);
    if(i >= N)
        return;

    if((imove[i] == -2) || (imove[i] == -3)){
        // It is a boundary element, the mass is actually the area
        h_var_in[i] = hfac * pow(m[i], 1.f / (DIMS - 1));
        return;
    }

    h_var_in[i] = hfac * pow(m[i] / rho_in[i], 1.f / DIMS);
}

/*
 * @}
 */
