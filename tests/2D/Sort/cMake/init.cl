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

#include "resources/Scripts/types/types.h"

__kernel void init(__global int* imove,
                   __global float* rho,
                   __global float* m,
                   __global vec* u,
                   __global vec* dudt,
                   uint N)
{
    const uint i = get_global_id(0);
    if(i >= N)
        return;
    
    imove[i] = 1;
    rho[i] = 1.0;
    m[i] = 1.0 / N;
    u[i] = VEC_ZERO;
    dudt[i] = VEC_ZERO;
}
