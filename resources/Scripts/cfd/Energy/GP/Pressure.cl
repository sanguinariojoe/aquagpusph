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
 * @brief Tool to compute the fluid pressure mechanical energy, and compressible
 * energy induced by a solid boundary in the fluid particles.
 */

#if defined(LOCAL_MEM_SIZE) && defined(NO_LOCAL_MEM)
    #error NO_LOCAL_MEM has been set.
#endif

#include "resources/Scripts/types/types.h"
#include "resources/Scripts/KernelFunctions/Kernel.h"

/** @brief Tool to compute the fluid pressure mechanical energy, and compressible
 * energy induced by a solid boundary in the fluid particles.
 *
 * See the following paper:
 *
 * M. Antuono, S. Marrone, A. Colagrossi, B. Bouscasse, "Energy balance in the
 * delta-SPH scheme". Computer methods in applied mechanincs and engineering,
 * vol 289, pp 209-226, 2015.
 *
 * @param GP_energy_degradpdt Energy induced by the ghost particles due to the
 * pressure term.
 * @param GP_energy_dedivudt Energy induced by the ghost particles due to the
 * compressibility.
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
 * @param icell Cell where each particle is located.
 * @param ihoc Head of chain for each cell (first particle found).
 * @param N Number of particles.
 * @param n_cells Number of cells in each direction
 * @param GP_energy_iset Ghost particles set to be considered.
 */
__kernel void main(__global float* GP_energy_degradpdt,
                   __global float* GP_energy_dedivudt,
                   const __global uint* iset,
                   const __global int* imove,
                   const __global vec* r,
                   const __global vec* u,
                   const __global float* rho,
                   const __global float* p,
                   const __global float* m,
                   usize N,
                   uint GP_energy_iset,
                   LINKLIST_LOCAL_PARAMS)
{
    // find position in global arrays
    const usize i = get_global_id(0);
    const usize it = get_local_id(0);
    if(i >= N)
        return;
    if(imove[i] != 1){
        GP_energy_degradpdt[i] = 0.f;
        GP_energy_dedivudt[i] = 0.f;
        return;
    }

    const vec_xyz r_i = r[i].XYZ;
    const vec_xyz u_i = u[i].XYZ;
    const float p_i = p[i];
    const float rho_i = rho[i];
    const float m_i = m[i];

    // Initialize the output
    #ifndef LOCAL_MEM_SIZE
        #define _E_GRADP_ GP_energy_degradpdt[i]
        #define _E_DIVU_ GP_energy_dedivudt[i]
    #else
        #define _E_GRADP_ degradpdt_l[i]
        #define _E_DIVU_ dedivudt_l[i]
        __local float degradpdt_l[LOCAL_MEM_SIZE];
        __local float dedivudt_l[LOCAL_MEM_SIZE];
    #endif
    _E_GRADP_ = 0.f;
    _E_DIVU_ = 0.f;

    const usize c_i = icell[i];
    BEGIN_NEIGHS(c_i, N, n_cells, icell, ihoc){
        if((iset[j] != GP_energy_iset) || (imove[j] != -1)){
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

            _E_GRADP_ -= dot(u_i, (p_i + p[j]) * f_ij * r_ij);
            _E_DIVU_ -= dot(u[j].XYZ - u_i, r_ij) * f_ij;
        }
    }END_NEIGHS()

    _E_GRADP_ *= m_i / rho_i;
    _E_DIVU_ *= p_i * m_i / rho_i;

    #ifdef LOCAL_MEM_SIZE
        GP_energy_degradpdt[i] = _E_GRADP_;
        GP_energy_dedivudt[i] = _E_DIVU_;
    #endif
}
