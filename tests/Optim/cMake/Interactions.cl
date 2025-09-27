/*
 *  This file is part of AQUAgpusph, a free CFD program based on SPH.
 *  Copyright (C) 2025  Jose Luis Cercos Pita <jl.cercos@upm.es>
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

#if defined(LOCAL_MEM_SIZE) && defined(NO_LOCAL_MEM)
    #error NO_LOCAL_MEM has been set.
#endif

#include "resources/Scripts/types/types.h"
#include "resources/Scripts/KernelFunctions/Kernel.h"

__kernel void legacy(const __global int* imove,
                     const __global vec* r,
                     const __global vec* u,
                     const __global float* rho,
                     const __global float* m,
                     const __global float* p,
                     __global vec* grad_p,
                     __global float* div_u,
                     usize N,
                     LINKLIST_LOCAL_PARAMS)
{
    const usize i = get_global_id(0);
    const usize it = get_local_id(0);
    if(i >= N)
        return;
    if(imove[i] != 1){
        return;
    }

    const vec_xyz r_i = r[i].XYZ;
    const vec_xyz u_i = u[i].XYZ;
    const float p_i = p[i];
    const float rho_i = rho[i];

    __private vec_xyz __grad_p = VEC_ZERO.XYZ;
    __private float __div_u = 0.f;

    const usize c_i = icell[i];
    BEGIN_NEIGHS(c_i, N, n_cells, icell, ihoc){
        if(i == j){
            j++;
            continue;
        }
        if(imove[j] != 1){
            j++;
            continue;
        }
        const vec_xyz r_ij = r[j].XYZ - r_i;
        const float q = length(r_ij) / H;
        if(q >= SUPPORT)
        {
            j++;
            continue;
        }
        {
            const float rho_j = rho[j];
            const float p_j = p[j];
            const float udr = dot(u[j].XYZ - u_i, r_ij);
            const float f_ij = kernelF(q) * CONF * m[j];

            __grad_p += (p_i + p_j) / (rho_i * rho_j) * f_ij * r_ij;
            __div_u += udr * f_ij * rho_i / rho_j;
        }
    }END_NEIGHS()

    grad_p[i].XYZ = __grad_p;
    div_u[i] = __div_u;
}


__kernel void optim(const __global int* imove,
                    const __global vec* r,
                    const __global vec* u,
                    const __global float* rho,
                    const __global float* m,
                    const __global float* p,
                    const __global svec2* jhoc,
                    __global vec* grad_p,
                    __global float* div_u,
                    usize N)
{
    const usize i = get_global_id(0);
    const usize it = get_local_id(0);
    if(i >= N)
        return;
    if(imove[i] != 1){
        return;
    }

    const vec_xyz r_i = r[i].XYZ;
    const vec_xyz u_i = u[i].XYZ;
    const float p_i = p[i];
    const float rho_i = rho[i];

    __private vec_xyz __grad_p = VEC_ZERO.XYZ;
    __private float __div_u = 0.f;

    FOR_NEIGHS(N, jhoc){
        if(i == j){
            continue;
        }
        if(imove[j] != 1){
            continue;
        }
        const vec_xyz r_ij = r[j].XYZ - r_i;
        const float q = length(r_ij) / H;
        if(q >= SUPPORT)
        {
            continue;
        }
        {
            const float rho_j = rho[j];
            const float p_j = p[j];
            const float udr = dot(u[j].XYZ - u_i, r_ij);
            const float f_ij = kernelF(q) * CONF * m[j];

            __grad_p += (p_i + p_j) / (rho_i * rho_j) * f_ij * r_ij;
            __div_u += udr * f_ij * rho_i / rho_j;
        }
    }END_FOR_NEIGHS()

    grad_p[i].XYZ = __grad_p;
    div_u[i] = __div_u;
}
