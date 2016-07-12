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

/** @addtogroup cfd
 * @{
 */

/** @file
 *  @brief Omega term computation
 */

#ifndef HAVE_3D
    #include "../../../../types/2D.h"
    #include "../../../../KernelFunctions/Wendland2D.hcl"
#else
    #include "../../../../types/3D.h"
    #include "../../../../KernelFunctions/Wendland3D.hcl"
#endif

/** @brief Compute the \f$ \Omega \f$ term
 *
 * \f$ \Omega_i = 1 - m_i \frac{\partial h_i}{\partial \rho_i}
 * \sum_j \frac{\partial W_{ij}}{\partial h_i} \f$
 *
 * @param imove Moving flags.
 *   - imove > 0 for regular fluid/solid particles.
 *   - imove = 0 for sensors.
 *   - imove < 0 for boundary elements/particles.
 * @param imirrored 0 if the particle has not been mirrored, 1 otherwise.
 * @param r Position \f$ \mathbf{r} \f$.
 * @param rho Density \f$ \rho \f$.
 * @param m Mass \f$ m \f$.
 * @param h_var variable kernel lenght \f$ h \f$.
 * @param Omega \f$ \Omega \f$ term.
 * @param icell Cell where each particle is located.
 * @param ihoc Head of chain for each cell (first particle found).
 * @param N Number of particles.
 * @param n_cells Number of cells in each direction
 * @see Iason Zisis, Bas van der Linden, Christina Giannopapa, Barry Koren. On
 * the derivation of SPH schemes for shocks through inhomogeneous media. Int.
 * Jnl. of Multiphysics. Vol. 9, Number 2. 2015
 */
__kernel void entry(const __global int* imove,
                    const __global int* imirrored,
                    const __global vec* r,
                    const __global float* rho,
                    const __global float* m,
                    const __global float* h_var,
                    __global float* Omega,
                    // Link-list data
                    const __global uint *icell,
                    const __global uint *ihoc,
                    // Simulation data
                    uint N,
                    uivec4 n_cells)
{
    const uint i = get_global_id(0);
    const uint it = get_local_id(0);
    if(i >= N)
        return;
    if((imove[i] != 1) || (!imirrored[i]))
        return;

    const vec_xyz r_i = r[i].XYZ;
    const float h_i = h_var[i];
    const float m_i = m[i];
    // The partial derivative of the kernel length with respect the density is
    // the same kernel length divided by the density and the number of
    // dimensions
    const float dhdrho_i = h_i / (DIMS * rho[i]);

    // The partial derivative of the kernel with respect to h multiplier is
    // 1 / h^(d + 1)
    #ifndef HAVE_3D
        const float conh = 1.f / (h_i * h_i * h_i);
    #else
        const float conh = 1.f / (h_i * h_i * h_i * h_i);
    #endif

    // Initialize the output
    #ifndef LOCAL_MEM_SIZE
        #define _OMEGA_ Omega[i]
    #else
        #define _OMEGA_ Omega_l[it]
        __local float Omega_l[LOCAL_MEM_SIZE];
        _OMEGA_ = 1.f;
    #endif

    BEGIN_LOOP_OVER_NEIGHS(){
        if((imove[j] != 1) || (imirrored[j])){
            j++;
            continue;
        }
        const vec_xyz r_ij = r[j].XYZ - r_i;
        const float q = length(r_ij) / h_i;
        if(q >= SUPPORT)
        {
            j++;
            continue;
        }
        {
            // n-scheme
            _OMEGA_ -= m_i * dhdrho_i * conh * kernelH(q);
        }
    }END_LOOP_OVER_NEIGHS()

    #ifdef LOCAL_MEM_SIZE
        Omega[i] = _OMEGA_;
    #endif
}

/*
 * @}
 */
