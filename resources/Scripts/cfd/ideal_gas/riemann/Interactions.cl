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
 * @brief Fluid particles interactions computation.
 */

#if defined(LOCAL_MEM_SIZE) && defined(NO_LOCAL_MEM)
    #error NO_LOCAL_MEM has been set.
#endif

#include "resources/Scripts/types/types.h"
#include "resources/Scripts/KernelFunctions/Kernel.h"
#include "resources/Scripts/cfd/ideal_gas/sound_speed.hcl"


/** @brief Fluid particles interactions computation.
 *
 * Compute the differential operators involved in the numerical scheme, taking
 * into account just the fluid-fluid interactions.
 *
 * @param imove Moving flags.
 *   - imove > 0 for regular fluid particles.
 *   - imove = 0 for sensors.
 *   - imove < 0 for boundary elements/particles.
 * @param r Position \f$ \mathbf{r} \f$.
 * @param u Velocity \f$ \mathbf{u} \f$.
 * @param rho Density \f$ \rho \f$.
 * @param m Mass \f$ m \f$.
 * @param p Pressure \f$ p \f$.
 * @param grad_p Pressure gradient \f$ \frac{\nabla p}{rho} \f$.
 * @param lap_u Velocity laplacian \f$ \frac{\Delta \mathbf{u}}{rho} \f$.
 * @param div_u Velocity divergence \f$ \rho \nabla \cdot \mathbf{u} \f$.
 * @param N Number of particles.
 * @param icell Cell where each particle is located.
 * @param ihoc Head of chain for each cell (first particle found).
 * @param n_cells Number of cells in each direction
 */



__kernel void entry(const __global unsigned int* iset,
                    const __global int* imove,
                    const __global vec* r,
                    const __global vec* u,
                    const __global float* rho,                    
                    const __global float* m,
                    const __global float* p,
                    __global vec* grad_p,
                    __global float* div_u,
                    __global float* work_density,
                    __constant float* gamma,
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
    const float gamma_i = gamma[iset[i]];    
    const float s_i = sound_speed_perfect_gas(gamma_i, p_i, rho_i);    
    const float m_i = m[i];


    //const float D_i = sqrt(m_i/rho_i);

    // Initialize the output
    #ifndef LOCAL_MEM_SIZE
        #define _GRADP_ grad_p[i].XYZ
        #define _W_DEN_ work_density[i]
        #define _DIVU_ div_u[i]
    #else
        #define _GRADP_ grad_p_l[it]
        #define _W_DEN_ work_density_l[it]
        #define _DIVU_ div_u_l[it]
        __local vec_xyz grad_p_l[LOCAL_MEM_SIZE];
        __local float work_density_l[LOCAL_MEM_SIZE];
        __local float div_u_l[LOCAL_MEM_SIZE];
        _GRADP_ = VEC_ZERO.XYZ;
        _W_DEN_ = 0.f;
        _DIVU_ = 0.f;
    #endif


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
            const float m_j = m[j];
            const float gamma_j = gamma[iset[j]];            
            const float s_j = sound_speed_perfect_gas(gamma_j, p_j, rho_j); 
            
            const vec_xyz l_ij = r_ij / length(r_ij);
            const float u_R_i = dot(u[i].XYZ, l_ij);
            const float u_R_j = dot(u[j].XYZ, l_ij);

            const float u_star = (u_R_j * rho_j * s_j + u_R_i * rho_i * s_i - p_j + p_i)/(rho_j * s_j + rho_i * s_i);
            const float p_star = (p_j * rho_i * s_i + p_i * rho_j * s_j - rho_j * s_j * rho_i * s_i * (u_R_j- u_R_i))/(rho_j * s_j + rho_i * s_i);

            const float Wij_prima = -q * kernelF(q) * CONW;

            const float aux = 2.0f * m_j / (rho_j * H) * (u_R_i - u_star) * Wij_prima;

            //_DIVU_ +=  2.0f * m_j * rho_i / (rho_j * H) * (u_R_i - u_star) * kernelF(q) * CONW;
            _DIVU_ += rho_i * aux;
            
            _GRADP_ -= 2.0f * m_j * p_star /(rho_j * rho_i * H) * Wij_prima * l_ij;
            
            //_W_DEN_ += 2.0f * m_j * p_star /(rho_j * rho_i * H) * (u_R_i - u_star) * Wij_prima
            _W_DEN_ += p_star / rho_i * aux;
           
        }
    }END_NEIGHS()

    #ifdef LOCAL_MEM_SIZE
        grad_p[i].XYZ = _GRADP_;
        work_density[i].XYZ = _W_DEN_;
        div_u[i] = _DIVU_;
    #endif
}
