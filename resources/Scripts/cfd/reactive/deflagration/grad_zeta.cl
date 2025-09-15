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
#include "resources/Scripts/cfd/reaction/reaction_generic.hcl"

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
__kernel void entry(const __global int* imove,
                    const __global vec* r,
                    const __global float* rho,
                    const __global float* m,
                    __global vec* grad_zeta,
                    const __global float* z,
                    const __global float* ys,
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
   
    const float z_i = z[i];
    const float y_0_i = ys[i].SPECIES_COMPONENT0;

    const float rho_i = rho[i];

    const float zeta_i = give_zeta(z_i, y_0_i);

    // Initialize the output
    #ifndef LOCAL_MEM_SIZE
        
        #define _GRAD_ZETA_ grad_zeta[i].XYZ
        
    #else

        #define _GRAD_ZETA_ grad_zeta_l[it]
        
        __local vec_xyz grad_zeta_l[LOCAL_MEM_SIZE];
        
        _GRAD_ZETA_ = VEC_ZERO.XYZ;        

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
            const float f_ij = kernelF(q) * CONF * m[j];
            const float z_j = z[j];
            const float y_0_j = ys[j].SPECIES_COMPONENT0;

            const float zeta_j = give_zeta(z_j, y_0_j);           
            
            //I think this is wrong
            _GRAD_ZETA_ += (z_i + z_j) / rho_j * f_ij * r_ij;         
            //that maybe   
            //_GRAD_ZETA_ += (zeta_i + zeta_j) / rho_j * f_ij * r_ij;
            // is correct. check this!         
        }
    }END_NEIGHS()

    #ifdef LOCAL_MEM_SIZE
        grad_zeta[i].XYZ = _GRAD_ZETA_;        
    #endif
}
