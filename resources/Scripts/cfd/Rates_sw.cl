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
#include "resources/Scripts/KernelFunctions/Kernel.h"

/** @brief Velocity and density variation rates computation.
 *
 * The mass conservation and momentum equations are applied from the already
 * computed differential operators:
 *
 *   - \f$ \frac{\mathrm{d} \mathbf{u}}{\mathrm{d} t} =
 *     - \frac{\nabla p}{rho}
 *     + \frac{\mu}{rho} \Delta \mathbf{u}
 *     + \mathbf{g}\f$
 *   - \f$ \frac{\mathrm{d} \rho}{\mathrm{d} t} =
 *     - \rho \nabla \cdot \mathbf{u}
 *     + \delta \Delta t \frac{\rho_a}{\rho_0} \Delta p\f$
 *
 * @param iset Set of particles index.
 * @param imove Moving flags.
 *   - imove > 0 for regular fluid particles.
 *   - imove = 0 for sensors.
 *   - imove < 0 for boundary elements/particles.
 * @param rho Density \f$ \rho_{n+1} \f$.
 * @param grad_p Pressure gradient \f$ \frac{\nabla p}{rho} \f$.
 * @param lap_u Velocity laplacian \f$ \frac{\Delta \mathbf{u}}{rho} \f$.
 * @param div_u Velocity divergence \f$ \rho \nabla \cdot \mathbf{u} \f$.
 * @param dudt Velocity rate of change
 * \f$ \left. \frac{d \mathbf{u}}{d t} \right\vert_{n+1} \f$.
 * @param drhodt Density rate of change
 * \f$ \left. \frac{d \rho}{d t} \right\vert_{n+1} \f$.
 * @param visc_dyn Dynamic viscosity \f$ \mu \f$.
 * @param N Number of particles.
 * @param g Gravity acceleration \f$ \mathbf{g} \f$.
 */

#define gamma 1.4f

__kernel void entry(//const __global uint* iset,
                    const __global int* imove,
                    const __global vec* r,
                    const __global vec* u,
                    const __global float* rho,
                    const __global float* m,
                    const __global float* p,
                    const __global float* eee,
                    __global vec* dudt,
                    __global float* drhodt,
                    __global float* dedt,
                    __global float* div_u,
                    __global vec* grad_p,
/*                     const __global uint *icell,
                    const __global uint *ihoc,
                    const uivec4 n_cells, */
                    // Simulation data
                    const uint N,
                    LINKLIST_LOCAL_PARAMS
                    )
{
    unsigned int i = get_global_id(0);
    if(i >= N)
        return;
    if(imove[i] != 1)
        return;

    const vec_xyz r_i = r[i].XYZ;
    const vec_xyz u_i = u[i].XYZ;
    const float p_i = p[i];
    const float rho_i = rho[i];
    const float e_i = eee[i];
    const float s_i = sqrt(gamma * p_i / rho_i);
    const float m_i = m[i];
    const float D_i = sqrt(4.0f*m_i*M_1_PI_F/rho_i);
    //const float D_i = sqrt(m_i/rho_i);
    
    drhodt[i] = 0.0f;
    dudt[i] = VEC_ZERO.XYZ;
    dedt[i] = 0.0f;

    const usize c_i = icell[i];
    BEGIN_NEIGHS(c_i, N, n_cells, icell, ihoc){
    //BEGIN_LOOP_OVER_NEIGHS(){
        if(i == j){
            j++;
            continue;
        }
        if(imove[j] != 1){
            j++;
            continue;
        }

        vec_xyz r_ij = r[j].XYZ - r_i;

        const float tyni = 1.0e-12f;
        float auxval = 1.0f / (length(r_ij)+tyni);
        vec_xyz l_ij = r_ij * auxval;

        float rho_j = rho[j];
        float p_j = p[j];
        float e_j = eee[j];
        float m_j = m[j];
        float D_j = sqrt(4.0f*m_j*M_1_PI_F/rho_j);
        //float D_j = sqrt(m_j/rho_j);
        float s_j = sqrt(gamma * p_j / rho_j);

        //float hh = 0.5 * (D_i + D_j);
        float hh = 2.0f*(D_i + D_j);
        float beta = 0.7f * M_PI_F * hh * hh;
        float im_beta = 1.0f / (beta + tyni);
        float q = length(r_ij) / (hh + tyni);
        
//        printf("SUPPORT=%f\n",SUPPORT);
        if(q >= SUPPORT)
        {
            j++;
            continue;
        }

        //float Wij_prima = 0.0f;
        //if(q<=2.0f)
        //{
        float Wij_prima = -0.75f * ( 2.0f - q) * ( 2.0f - q) * im_beta;
        //}
        if(q<=1.0f)
        {
            Wij_prima = (2.25f* q * q - 3.0f * q) *  im_beta;
        }


     //   const float u_R_i = 1.0f/H*dot(u[i].XYZ, l_ij);
     //   const float u_R_j = 1.0f/H*dot(u[j].XYZ, l_ij);
     
        float u_R_i = dot(u[i].XYZ, l_ij);
        float u_R_j = dot(u[j].XYZ, l_ij);
        
        float u_star = (u_R_j * rho_j * s_j + u_R_i * rho_i * s_i - p_j + p_i)/(rho_j * s_j + rho_i * s_i);
        float p_star = (p_j * rho_i * s_i + p_i * rho_j * s_j - rho_j * s_j * rho_i * s_i * (u_R_j- u_R_i))/(rho_j * s_j + rho_i * s_i);

        //const vec_xyz auxu = 2.0f * m_j * p_star /(rho_j * rho_i * H) * kernelF(q) * l_ij;
        //vec_xyz auxu = 2.0f * m_j * p_star /(rho_j * rho_i * hh) * Wij_prima * l_ij;

        //drhodt[i] -= 2.0f * m_j * rho_i / (rho_j * H) * (u_R_i - u_star) * kernelF(q); 
    
        drhodt[i] -= 2.0f * m_j * rho_i / (rho_j * hh) * (u_R_i - u_star) * Wij_prima;
        //dudt[i]   += 2.0f * m_j * p_star /(rho_j * rho_i * H) * kernelF(q) * l_ij;
        //dudt[i]   += auxu;
        dudt[i]   += 2.0f * m_j * p_star /(rho_j * rho_i * hh) * Wij_prima * l_ij;
        //dedt[i]   -= 2.0f * m_j * p_star /(rho_j * rho_i * H) * (u_R_i - u_star) * kernelF(q);
        dedt[i]   -= 2.0f * m_j * p_star /(rho_j * rho_i * hh) * (u_R_i - u_star) * Wij_prima;

    //}END_LOOP_OVER_NEIGHS()
    }END_NEIGHS()

    div_u[i] = drhodt[i] / rho_i;
    grad_p[i] = dudt[i] / rho_i;

}
