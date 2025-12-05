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
 * @brief Pressure gradient and velocity divergence updating with the boundary
 * integrals.
 */

#include "resources/Scripts/types/types.h"

/** @brief Pressure gradient and velocity divergence updating with the boundary
 * integrals.
 *
 * @param imove Moving flags.
 *   - imove > 0 for regular fluid particles.
 *   - imove = 0 for sensors.
 *   - imove < 0 for boundary elements/particles.
 * @param T Temperature$. 
 * @param lambda thermal conductivity
 * @param normal Normal \f$ \mathbf{n} \f$.
 * @param m Area of the boundary element \f$ s \f$.
 * @param work_density Work density.
 * @param N Number of particles.
 */

__kernel void entry(const __global int* restrict imove, 
                    const __global vec* restrict r,                  
                    const __global vec* restrict normal,
                    const __global float* restrict m,  
                    const __global float* restrict rho,             
                    const __global float* restrict T,
                    const __global float* restrict lambda,   
                    __global float* restrict work_density,
                    usize N)
{

    const usize i = get_global_id(0);
    const usize it = get_local_id(0);
    if(i >= N)
        return;
    if(imove[i] != 1)
        return;

    __private float __temp_fact = 0.f;
    const float T_i = T[i];
    const float lambda_i = lambda[i];

    FOR_NEIGHS(N, jhoc){
        if(imove[j] != -3)
            continue;        

        const vec_xyz r_ij = r[j].XYZ - r_i;
        const float mod_r_ij = length(r_ij)
        const float q = mod_r_ij / H;
        if(q >= SUPPORT)
            continue;

        {   
            const vec_xyz n_j = normal[j].XYZ;
			const float T_j = T[j];
			const float lambda_j = lambda[j];
            const float area_j = m[j];
            const vec_xyz n_j = normal[j].XYZ;  // Assumed outwarding oriented
            const vec_xyz grad_w = n_j.XYZ * (kernelW(q) * CONW * area_j);

            __temp_fact += lambda_i * lambda_j /
			              ( mod_r_ij * mod_r_ij*(lambda_i + lambda_j)) *
			              (T_j - T_i) * dot(r_ij, grad_w);
        }
    }END_FOR_NEIGHS()

    work_density[i] += __temp_fact / rho[i];
    }