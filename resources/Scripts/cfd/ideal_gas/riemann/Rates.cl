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
 * @brief Velocity and density variation rates computation.
 */

#include "resources/Scripts/types/types.h"

/** @brief Energy variation rates computation.
 *
 * The energy conservation are applied from the already
 * computed differential operators:
 *
 * @param iset Set of particles index.
 * @param imove Moving flags.
 *   - imove > 0 for regular fluid particles.
 *   - imove = 0 for sensors.
 *   - imove < 0 for boundary elements/particles.
 * @param work_density compression work
 * @param deintdt Energy rate of change
 * @param N Number of particles.
 */
__kernel void entry(const __global uint* iset,
                    const __global int* imove,
                    const __global float* work_density,
                    __global float* deintdt,
                    const usize N)
{
    const usize i = get_global_id(0);
    if(i >= N)
        return;
    if(imove[i] != 1)
        return;

    // Conservation of energy equation
    deintdt[i] = -work_density[i];
}
