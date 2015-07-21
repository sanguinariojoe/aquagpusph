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
 * @brief Fixed ghost particles fields interpolation.
 */

#if defined(LOCAL_MEM_SIZE) && defined(NO_LOCAL_MEM)
    #error NO_LOCAL_MEM has been set.
#endif

#ifndef HAVE_3D
    #include "../../types/2D.h"
    #include "../../KernelFunctions/Wendland2D.hcl"
#else
    #include "../../types/3D.h"
    #include "../../KernelFunctions/Wendland3D.hcl"
#endif

/** @brief Fixed ghost particles fields interpolation.
 *
 * In order to conveniently extend the flow fields we need to know the mirrored
 * values. Also the Shepard values should be recomputed taking into account the
 * mirrored position.
 *
 * @param imove Moving flags.
 *   - imove > 0 for regular fluid particles.
 *   - imove = 0 for sensors.
 *   - imove < 0 for boundary elements/particles.
 * @param r Position \f$ \mathbf{r} \f$.
 * @param normal Normal \f$ \mathbf{n} \f$.
 * @param m Mass \f$ m \f$.
 * @param rho Density \f$ \rho \f$.
 * @param p Pressure \f$ p \f$.
 * @param u Velocity \f$ \mathbf{u} \f$.
 * @param gp_rho Interpolated density in the mirrored position \f$ \rho \f$.
 * @param gp_p Interpolated pressure in the mirrored position \f$ p \f$.
 * @param gp_u Interpolated velocity in the mirrored position \f$ \mathbf{u} \f$.
 * @param shepard Shepard term
 * \f$ \gamma(\mathbf{x}) = \int_{\Omega}
 *     W(\mathbf{y} - \mathbf{x}) \mathrm{d}\mathbf{x} \f$.
 * @param icell Cell where each particle is located.
 * @param gp_icell Cell where each mirrored ghost particle is located.
 * @param ihoc Head of chain for each cell (first particle found).
 * @param N Number of particles.
 * @param n_cells Number of cells in each direction
 */
__kernel void main(const __global int* imove,
                   const __global vec* r,
                   const __global vec* normal,
                   const __global float* m,
                   const __global float* rho,
                   const __global float* p,
                   const __global vec* u,
                   __global float* gp_rho,
                   __global float* gp_p,
                   __global vec* gp_u,
                   __global float* shepard,
                   // Link-list data
                   const __global uint *icell,
                   const __global uint *gp_icell,
                   const __global uint *ihoc,
                   // Simulation data
                   uint N,
                   uivec4 n_cells)
{
    const uint i = get_global_id(0);
    const uint it = get_local_id(0);
    if(i >= N)
        return;
    if(imove[i] != -1)
        return;
    
    const vec_xyz r_i = r[i].XYZ;

    // Initialize the output
    #ifndef LOCAL_MEM_SIZE
        #define _RHO_ gp_rho[i]
        #define _P_ gp_p[i]
        #define _U_ gp_u[i].XYZ
        #define _SHEPARD_ shepard[i]
    #else
        #define _RHO_ rho_l[it]
        #define _P_ p_l[it]
        #define _U_ u_l[it]
        #define _SHEPARD_ shepard_l[it]
        __local float rho_l[LOCAL_MEM_SIZE];
        __local float p_l[LOCAL_MEM_SIZE];
        __local vec_xyz u_l[LOCAL_MEM_SIZE];
        __local float shepard_l[LOCAL_MEM_SIZE];
        _RHO_ = 0.f;
        _P_ = 0.f;
        _U_ = VEC_ZERO.XYZ;
    #endif
    _SHEPARD_ = 0.f;

    #undef C_I()
    #define C_I() const uint c_i = gp_icell[i]
    BEGIN_LOOP_OVER_NEIGHS(){
        if(imove[j] != 1){
            j++;
            continue;
        }
        const vec_xyz r_ij = r[j].XYZ - r_i;
        const float q = fast_length(r_ij) / H;
        if(q >= SUPPORT)
        {
            j++;
            continue;
        }

        {
            const float rho_j = rho[j];
            const float p_j = p[j];
            const vec_xyz u_j = u[j];
            const float m_j = m[j];

            {
                const float w_ij = kernelW(q) * CONW * m_j / rho_j;
                _RHO_ += w_ij * rho_j;
                _P_ += w_ij * p_j; 
                _U_ += w_ij * u_j;
                _SHEPARD_ += w_ij;
            }
        }
    }END_LOOP_OVER_NEIGHS()

    #ifdef LOCAL_MEM_SIZE
        gp_rho[i] = _RHO_;
        gp_p[i] = _P_;
        gp_u[i].XYZ = _U_;
        shepard[i] = _SHEPARD_;
    #endif
}
