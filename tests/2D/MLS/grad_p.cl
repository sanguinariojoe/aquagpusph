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
#include "resources/Scripts/KernelFunctions/Kernel.h"

__kernel void entry(const __global vec* r,
                    const __global float* m,
                    const __global float* p,
                    const __global matrix* mls,
                    __global vec* grad_p,
                    const __global uint *icell,
                    const __global uint *ihoc,
                    uint N,
                    uivec4 n_cells)
{
    const uint i = get_global_id(0);
    const uint it = get_local_id(0);
    if(i >= N)
        return;

    const vec_xyz r_i = r[i].XYZ;
    const float p_i = p[i];

    // Initialize the output
    #ifndef LOCAL_MEM_SIZE
        #define _GRADP_ grad_p[i].XYZ
    #else
        #define _GRADP_ grad_p_l[it]
        __local vec_xyz grad_p_l[LOCAL_MEM_SIZE];
        _GRADP_ = VEC_ZERO.XYZ;
    #endif

    BEGIN_LOOP_OVER_NEIGHS(){
        if(i == j) {
            j++;
            continue;
        }
        const vec_xyz r_ij = r[j].XYZ - r_i;
        const float q = length(r_ij) / H;
        if(q >= SUPPORT) {
            j++;
            continue;
        }
        {
            const float f_ij = kernelF(q) * CONF * m[j];
            _GRADP_ += (p[j] - p_i) * f_ij * r_ij;
        }
    }END_LOOP_OVER_NEIGHS()

    grad_p[i].XYZ = MATRIX_DOT(mls[i], _GRADP_);
}

__kernel void squared_error(const __global vec* grad_p,
                            __global float* se,
                            uint N)
{
    const uint i = get_global_id(0);
    if(i >= N)
        return;

    vec_xyz ref = VEC_ZERO.XYZ;
    ref.y = 1.0;
    const vec_xyz error = grad_p[i].XYZ - ref;
    
    se[i] = dot(error, error);
}
