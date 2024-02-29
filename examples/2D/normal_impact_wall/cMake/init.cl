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
 * @brief Init the boundary elements normal.
 */

#include "resources/Scripts/types/types.h"

/** @brief Initialize the boundary elements normal direction
 * @param imove Moving flags.
 *   - imove > 0 for regular fluid particles.
 *   - imove = 0 for sensors.
 *   - imove < 0 for boundary elements/particles.
 * @param normal Normal \f$ \mathbf{n} \f$.
 * @param normal Tangent \f$ \mathbf{t} \f$.
 */
__kernel void entry(const __global int* imove,
                    __global vec* normal,
                    __global vec* tangent,
                    uint N)
{
    const uint i = get_global_id(0);
    if(i >= N)
        return;

    normal[i] = VEC_ZERO;
    tangent[i] = VEC_ZERO;

    if(imove[i] != -3)
        return;

    normal[i].y = -1.f;
    tangent[i].x = -1.f;
}
