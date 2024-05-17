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
 * @brief Compute the tangent vector in 2D sims, and the binormal in 3D sims
 */

#include "resources/Scripts/types/types.h"

/** @brief Compute the tanget vector in 2D simulations, and the binormal vector
 * in 3D simulations.
 *
 * @param normal Normal vector \f$ \mathbf{n} \f$.
 * @param tangent Tangent vector \f$ \mathbf{t} \f$.
 * @param binormal Binormal vector \f$ \mathbf{t} \f$.
 * @param N Number of particles.
 */
__kernel void entry(const __global vec* normal,
                    __global vec* tangent,
                    __global vec* binormal,
                    usize N)
{
    const usize i = get_global_id(0);
    if(i >= N)
        return;

    #ifndef HAVE_3D
        binormal[i] = VEC_ZERO;
        tangent[i].x = normal[i].y;
        tangent[i].y = -normal[i].x;
    #else
        binormal[i] = cross(normal[i], tangent[i]);
    #endif
}

/*
 * @}
 */
