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

__kernel void predictor(const __global float* eint,
                        const __global float* deintdt,
                        __global float* eint_in,
                        __global float* deintdt_in,
                        const usize N)
{
    const usize i = get_global_id(0);
    if(i >= N)
        return;

    eint_in[i] = eint[i];
    deintdt_in[i] = deintdt[i];
}

__kernel void corrector(const __global int* imove,
                        __global float* eint,
                        const __global float* deintdt,
                        const unsigned int N,
                        const float dt)
{
    usize i = get_global_id(0);
    if(i >= N)
        return;

    if(imove[i] > 0) {
        eint[i] += dt * deintdt[i];
    }
}
