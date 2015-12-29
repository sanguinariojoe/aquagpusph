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
 * @brief Tool to compute the fluid viscous force and moment.
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

#if __LAP_FORMULATION__ == __LAP_MONAGHAN__
    #ifndef HAVE_3D
        #define __CLEARY__ 8.f
    #else
        #define __CLEARY__ 10.f
    #endif
#endif

/** @brief Tool to compute the viscous force and moment for an especific body.
 *
 * In this approach the following operation is performed for the boundary
 * elements:
 * \f$ \mathbf{f}_a = \mu \, \sum_b -\frac{
 *     \left(\mathbf{u}_b - \mathbf{u}_a\right) - 
 *     \left(\left(\mathbf{u}_b - \mathbf{u}_a\right) \cdot \mathbf{n}_a \right)
 *     \mathbf{n}_a}{
 *     \left(\mathbf{r}_b - \mathbf{r}_a\right) \cdot \mathbf{n}_a}s_a
 *     \, W\left(\mathbf{u}_b - \mathbf{u}_a\right) \frac{m_b}{\rho_b}\f$
 * where \f$ s_a \f$ is the area of the element, stored in the masses array.
 * The moment is computed therefore as:
 * \f$ \mathbf{m}_a  = \mathbf{f}_a \times
 * \left(\mathbf{r}_a - \mathbf{r}_0 \right) \f$
 * becoming \f$ \mathbf{r}_0 \f$ the reference point where the moment should be
 * computed.
 *
 * @param viscousForces_f Force of each boundary element to be computed [N].
 * @param viscousForces_m Moment of each boundary element to be computed
 * @param iset Set of particles index.
 * @param imove Moving flags.
 *   - imove > 0 for regular fluid particles.
 *   - imove = 0 for sensors.
 *   - imove < 0 for boundary elements/particles.
 * @param r Position \f$ \mathbf{r} \f$.
 * @param u Velocity \f$ \mathbf{u} \f$.
 * @param rho Density \f$ \rho \f$.
 * @param m Mass \f$ m \f$.
 * @param visc_dyn Dynamic viscosity \f$ \mu \f$.
 * @param icell Cell where each particle is located.
 * @param ihoc Head of chain for each cell (first particle found).
 * @param N Number of particles.
 * @param n_cells Number of cells in each direction
 * @param viscousForces_iset Particles set to be computed.
 * @param viscousForces_r Point with respect the moments are computed
 * \f$ \mathbf{r}_0 \f$.
 */
__kernel void entry(__global vec* viscousForces_f,
                    __global vec4* viscousForces_m,
                    const __global uint* iset,
                    const __global int* imove,
                    const __global vec* r,
                    const __global vec* u,
                    const __global float* rho,
                    const __global float* m,
                    __constant float* visc_dyn,
                    // Link-list data
                    const __global uint *icell,
                    const __global uint *ihoc,
                    // Simulation data
                    uint N,
                    uivec4 n_cells,
                    unsigned int viscousForces_iset,
                    vec viscousForces_r)
{
    const uint i = get_global_id(0);
    const uint it = get_local_id(0);
    if(i >= N)
        return;
    if((iset[i] != viscousForces_iset) || (imove[i] != -1)){
        viscousForces_f[i] = VEC_ZERO;
        viscousForces_m[i] = (vec4)(0.f, 0.f, 0.f, 0.f);
        return;
    }
    
    const vec_xyz r_i = r[i].XYZ;
    const vec_xyz u_i = u[i].XYZ;
    const float visc_dyn_i = visc_dyn[iset[i]];

    // Initialize the output
    #ifndef LOCAL_MEM_SIZE
        #define _F_ viscousForces_f[i].XYZ
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
            #if __LAP_FORMULATION__ == __LAP_MONAGHAN__
                const float r2 = (q * q + 0.01f) * H * H;
                const float udr = dot(u[j].XYZ - u_i, r_ij);
                _F_ += f_ij * __CLEARY__ * udr / r2 * r_ij;
            #elif __LAP_FORMULATION__ == __LAP_MORRIS__
                _F_ += f_ij * 2.f * (u[j].XYZ - u_i);
            #else
                #error Unknown Laplacian formulation: __LAP_FORMULATION__
            #endif
        }
    }END_LOOP_OVER_NEIGHS()

    _F_ *= visc_dyn_i * m[i] / rho[i];

    #ifdef LOCAL_MEM_SIZE
        viscousForces_f[i].XYZ = _F_;
    #endif

    const vec_xyz arm = r_i - viscousForces_r;
    viscousForces_m[i].z = arm.x * _F_.y - arm.y * _F_.x;
    viscousForces_m[i].w = 0.f;
    #ifdef HAVE_3D
        viscousForces_f[i].w = 0.f;
        viscousForces_m[i].x = arm.y * _F_.z - arm.z * _F_.y;
        viscousForces_m[i].y = arm.z * _F_.x - arm.x * _F_.z;
    #else
        viscousForces_m[i].x = 0.f;
        viscousForces_m[i].y = 0.f;
    #endif
}