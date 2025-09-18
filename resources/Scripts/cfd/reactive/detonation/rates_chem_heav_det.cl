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
#include "resources/Scripts/cfd/reactive/reaction_generic.hcl"
#include "resources/Scripts/cfd/reactive/detonation/heav_detonation.hcl"



/** @brief Compute the reaction sources sinks for each particle.
 *
 * @param rho density.
 * @param eint internal energy.
 * @param p pressure.
 * @param T temperature.
 * @param z tracer.
 * @param ys chemical species.
 * @param trigger virtual plug.
 * @param deintdt change of internal energy.
 * @param dysdt change of species.
 * @param zeta_dot rate of change of the reaction. 
 */
__kernel void entry(const __global unsigned int* iset,
                    const __global int* imove,
                    const __global float* rho,
                    const __global float* eint,
                    const __global float* p,
                    const __global float* T,
                    const __global float* z,
                    const __global species_t* ys,
                    const __global int* trigger,
                    __global float* deintdt,
                    __global species_t* dysdt,
                    __global float* zeta_dot,                    
                    usize N,
                    LINKLIST_LOCAL_PARAMS)
{
    const usize i = get_global_id(0);
    
    if(i >= N)
        return;
    if(imove[i] != 1){
        return;
    }

    w_rhos_heav_det(z[i], T[i], ys[i], trigger[i], dysdt+i, deintdt+i, zeta_dot+i);
}