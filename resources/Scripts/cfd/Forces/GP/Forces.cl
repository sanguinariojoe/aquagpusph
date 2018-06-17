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

/** @brief Tool to compute the pressure force and moment for an especific body,
 * due to the Ghost particles.
 *
 * The pressure forces is computed using the same expression than in the
 * fluid-fluid interactions.
 *
 * @param forces_f Force of each boundary element to be computed [N].
 * @param forces_m Moment of each boundary element to be computed
 * [N \f$ \cdot \f$ m].
 * @param iset Set of particles index.
 * @param imove Moving flags.
 *   - imove > 0 for regular fluid particles.
 *   - imove = 0 for sensors.
 *   - imove < 0 for boundary elements/particles.
 * @param r Position \f$ \mathbf{r} \f$.
 * @param u Velocity \f$ \mathbf{u} \f$.
 * @param rho Density \f$ \rho \f$.
 * @param p Pressure \f$ p \f$.
 * @param m Mass \f$ m \f$.
 * @param visc_dyn Dynamic viscosity \f$ \mu \f$.
 * @param icell Cell where each particle is located.
 * @param ihoc Head of chain for each cell (first particle found).
 * @param N Number of particles.
 * @param n_cells Number of cells in each direction
 * @param forces_iset Particles set to be computed.
 * @param forces_r Point with respect the moments are computed
 * \f$ \mathbf{r}_0 \f$.
 */
__kernel void entry(__global vec* forces_f,
                    __global vec4* forces_m,
                    const __global uint* iset,
                    const __global int* imove,
                    const __global vec* r,
                    const __global vec* u,
                    const __global float* rho,
                    const __global float* p,
                    const __global float* m,
                    __constant float* visc_dyn,
                    // Link-list data
                    const __global uint *icell,
                    const __global uint *ihoc,
                    // Simulation data
                    uint N,
                    uivec4 n_cells,
                    unsigned int forces_iset,
                    vec forces_r)
{
    // find position in global arrays
    const uint i = get_global_id(0);
    const uint it = get_local_id(0);
    if(i >= N)
        return;
    if((iset[i] != forces_iset) || (imove[i] != -1)){
        forces_f[i] = VEC_ZERO;
        forces_m[i] = (vec4)(0.f, 0.f, 0.f, 0.f);
        return;
    }

    const vec_xyz r_i = r[i].XYZ;
    const float p_i = p[i];
    const float mu = visc_dyn[iset[i]];

    // Initialize the output
    #ifndef LOCAL_MEM_SIZE
        #define _F_ forces_f[i].XYZ
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
            const float f_ij = kernelF(q) * CONF * m[j] / rho[j];
            const float p_ij = p_i + p[j];
            #if __LAP_FORMULATION__ == __LAP_MONAGHAN__
                const float r2 = (q * q + 0.01f) * H * H;
                const float udr = dot(u[j].XYZ - u_i, r_ij);
                _F_ += f_ij * (mu * __CLEARY__ * udr / r2 - p_ij) * r_ij;
            #elif __LAP_FORMULATION__ == __LAP_MORRIS__
                _F_ += f_ij * (mu * 2.f * (u[j].XYZ - u_i) - p_ij * r_ij);
            #else
                #error Unknown Laplacian formulation: __LAP_FORMULATION__
            #endif

        }
    }END_LOOP_OVER_NEIGHS()

    _F_ *= m[i] / rho[i];

    #ifdef LOCAL_MEM_SIZE
        forces_f[i].XYZ = _F_;
    #endif

    const vec arm = r[i] - forces_r;
    forces_m[i].z = arm.x * _F_.y - arm.y * _F_.x;
    forces_m[i].w = 0.f;
    #ifdef HAVE_3D
        forces_m[i].x = arm.y * _F_.z - arm.z * _F_.y;
        forces_m[i].y = arm.z * _F_.x - arm.x * _F_.z;
    #else
        forces_m[i].x = 0.f;
        forces_m[i].y = 0.f;
    #endif
}
