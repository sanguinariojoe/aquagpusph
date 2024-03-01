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
 * @brief Boundary integral term for constant fields.
 */

#if defined(LOCAL_MEM_SIZE) && defined(NO_LOCAL_MEM)
    #error NO_LOCAL_MEM has been set.
#endif

#include "resources/Scripts/types/types.h"
#include "resources/Scripts/KernelFunctions/Kernel.h"

/** @brief Compute the boundary integrals.
 * @param imove Moving flags.
 *   - imove > 0 for regular fluid particles.
 *   - imove = 0 for sensors.
 *   - imove < 0 for boundary elements/particles.
 * @param r Position \f$ \mathbf{r} \f$.
 * @param normal Normal \f$ \mathbf{n} \f$.
 * @param u Velocity \f$ \mathbf{u} \f$.
 * @param m Area of the boundary element \f$ s \f$.
 * @param grad_w_bi Gradient of constant fields due to the boundary integral
 * \f$ \langle \nabla 1 \rangle^{\partial \Omega} \f$.
 * @param div_u Velocity divergence \f$ \nabla \cdot \mathbf{u} \f$. Actually
 * this is just the part that has to do with the boundary element velocity
 * @param N Number of particles.
 * @param icell Cell where each particle is located.
 * @param ihoc Head of chain for each cell (first particle found).
 * @param n_cells Number of cells in each direction
 */
__kernel void entry(const __global int* imove,
                    const __global vec* r,
                    const __global vec* normal,
                    const __global vec* u,
                    const __global float* m,
                    __global vec* grad_w_bi,
                    __global float* div_u_bi,
                    uint N,
                    LINKLIST_LOCAL_PARAMS)
{
    const uint i = get_global_id(0);
    const uint it = get_local_id(0);
    if(i >= N)
        return;
    if(imove[i] != 1)
        return;

    const vec_xyz r_i = r[i].XYZ;

    // Initialize the output
    #ifndef LOCAL_MEM_SIZE
        #define _GRADW_ grad_w_bi[i].XYZ
        #define _DIVU_ div_u_bi[i]
    #else
        #define _GRADW_ grad_w_bi_l[it]
        #define _DIVU_ div_u_l[it]
        __local vec_xyz grad_w_bi_l[LOCAL_MEM_SIZE];
        __local float div_u_l[LOCAL_MEM_SIZE];
        _GRADW_ = VEC_ZERO.XYZ;
        _DIVU_ = 0.f;
    #endif

    const unsigned int c_i = icell[i];
    BEGIN_NEIGHS(c_i, N, n_cells, icell, ihoc){
        if(imove[j] != -3){
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
            const vec_xyz n_j = normal[j].XYZ;  // Assumed outwarding oriented
            const vec_xyz u_j = u[j].XYZ;
            const float area_j = m[j];
            const vec_xyz grad_w = n_j.XYZ * (kernelW(q) * CONW * area_j);
            _GRADW_ += grad_w;
            _DIVU_ -= dot(u_j, grad_w);
        }
    }END_NEIGHS()

    #ifdef LOCAL_MEM_SIZE
        grad_w_bi[i].XYZ = _GRADW_;
        div_u_bi[i] = _DIVU_;
    #endif
}

/** @brief Compute the pressure on each boundary element.
 * @param imove Moving flags.
 *   - imove > 0 for regular fluid particles.
 *   - imove = 0 for sensors.
 *   - imove < 0 for boundary elements/particles.
 * @param r Position \f$ \mathbf{r} \f$.
 * @param m Particle mass \f$ m \f$.
 * @param rho Particle density \f$ \rho \f$.
 * @param p Particle pressure \f$ p \f$.
 * @param N Number of particles.
 * @param icell Cell where each particle is located.
 * @param ihoc Head of chain for each cell (first particle found).
 * @param n_cells Number of cells in each direction
 */
__kernel void p_boundary(const __global int* imove,
                         const __global vec* r,
                         const __global float* m,
                         const __global float* rho,
                         __global float* p,
                         uint N,
                         LINKLIST_LOCAL_PARAMS)
{
    const uint i = get_global_id(0);
    const uint it = get_local_id(0);
    if(i >= N)
        return;
    if(imove[i] != -3)
        return;

    const vec_xyz r_i = r[i].XYZ;

    // Initialize the output
    #ifndef LOCAL_MEM_SIZE
        #define _P_ p[i]
    #else
        #define _P_ p_l[it]
        __local float p_l[LOCAL_MEM_SIZE];
    #endif
    _P_ = 0.0;

    const unsigned int c_i = icell[i];
    BEGIN_NEIGHS(c_i, N, n_cells, icell, ihoc){
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

        _P_ += 2.f * p[j] * kernelW(q) * CONW * m[j] / rho[j];
    }END_NEIGHS()

    #ifdef LOCAL_MEM_SIZE
        p[i] = _P_;
    #endif
}
