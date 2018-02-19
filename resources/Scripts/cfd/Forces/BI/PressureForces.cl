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
 * @brief Tool to compute the fluid pressure force and moment.
 */

#if defined(LOCAL_MEM_SIZE) && defined(NO_LOCAL_MEM)
    #error NO_LOCAL_MEM has been set.
#endif

#include "resources/Scripts/types/types.h"
#include "resources/Scripts/KernelFunctions/Kernel.h"

/** @brief Tool to compute the pressure force and moment for an especific body.
 *
 * .
 *
 * @param pressureForces_f Force of each boundary element to be computed [N].
 * @param pressureForces_m Moment of each boundary element to be computed
 * [N \f$ \cdot \f$ m].
 * @param iset Set of particles index.
 * @param imove Moving flags.
 *   - imove > 0 for regular fluid particles.
 *   - imove = 0 for sensors.
 *   - imove < 0 for boundary elements/particles.
 * @param r Position \f$ \mathbf{r} \f$.
 * @param grad_p Pressure gradient \f$ \frac{\nabla p}{rho} \f$.
 * @param grad_p_bulk Pressure gradient, due to the fluid-fluid interactions,
 * \f$ \left. \frac{\nabla p}{rho} \right\vert_{\Omega \rightarrow \Omega} \f$.
 * @param zeta Zeta term
 * \f$ \zeta(\mathbf{x}) = \int_{\partial \Omega}
 *     W(\mathbf{y} - \mathbf{x}) \mathrm{d}\mathbf{y} \f$.
 * @param m Mass \f$ m \f$.
 * @param icell Cell where each particle is located.
 * @param ihoc Head of chain for each cell (first particle found).
 * @param N Number of particles.
 * @param n_cells Number of cells in each direction
 * @param pressureForces_iset Particles set to be computed.
 * @param pressureForces_r Point with respect the moments are computed
 * \f$ \mathbf{r}_0 \f$.
 */
__kernel void entry(__global vec* pressureForces_f,
                    __global vec4* pressureForces_m,
                    const __global uint* iset,
                    const __global int* imove,
                    const __global vec* r,
                    const __global vec* grad_p,
                    const __global vec* grad_p_bulk,
                    const __global float* zeta,
                    const __global float* m,
                    const __global uint *icell,
                    const __global uint *ihoc,
                    uint N,
                    uivec4 n_cells,
                    unsigned int pressureForces_iset,
                    vec pressureForces_r)
{
    // find position in global arrays
    const uint i = get_global_id(0);
    const uint it = get_local_id(0);
    if(i >= N)
        return;
    if((iset[i] != pressureForces_iset) || (imove[i] != -3)){
        pressureForces_f[i] = VEC_ZERO;
        pressureForces_m[i] = (vec4)(0.f, 0.f, 0.f, 0.f);
        return;
    }

    const vec_xyz r_i = r[i].XYZ;
    const float area_i = m[i];

    // Initialize the output
    #ifndef LOCAL_MEM_SIZE
        #define _F_ pressureForces_f[i].XYZ
    #else
        #define _F_ f_l[it]
        __local vec_xyz f_l[LOCAL_MEM_SIZE];
    #endif
    _F_ = VEC_ZERO.XYZ;

    BEGIN_LOOP_OVER_NEIGHS(){
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
            const float m_j = m[j];
            const float zeta_j = zeta[j];
            // grad_p is already divided by the density
            const vec_xyz gradp_j = grad_p[j].XYZ - grad_p_bulk[j].XYZ;
            const float w_ij = kernelW(q) * CONW * m_j;

            _F_ += area_i * w_ij / zeta_j * gradp_j;
        }
    }END_LOOP_OVER_NEIGHS()

    #ifdef LOCAL_MEM_SIZE
        pressureForces_f[i].XYZ = _F_;
    #endif

    const vec_xyz arm = r_i - pressureForces_r.XYZ;
    pressureForces_m[i].z = arm.x * _F_.y - arm.y * _F_.x;
    pressureForces_m[i].w = 0.f;
    #ifdef HAVE_3D
        pressureForces_m[i].x = arm.y * _F_.z - arm.z * _F_.y;
        pressureForces_m[i].y = arm.z * _F_.x - arm.x * _F_.z;
    #else
        pressureForces_m[i].x = 0.f;
        pressureForces_m[i].y = 0.f;
    #endif
}
