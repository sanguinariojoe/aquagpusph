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

#ifndef EXCLUDED_PARTICLE
    #define EXCLUDED_PARTICLE(index) (imove[index] <= 0) && (imove[index] != -1)
#endif

#include "resources/Scripts/types/types.h"

#define gamma 1.44f

__kernel void entry(const __global unsigned int* iset,
                    const __global int* imove,
                    const __global float* rho,
                    __global float* p,
                    const __global float* eint,
                    usize N)
{
    usize i = get_global_id(0);
    if(i >= N)
        return;
    if(EXCLUDED_PARTICLE(i))
        return;

    p[i] = (gamma - 1.0f) * rho[i] * eint[i];
}
