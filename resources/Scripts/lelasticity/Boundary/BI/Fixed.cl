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

/** @addtogroup lela
 * @{
 */

/** @file
 * @brief Boundary elements deformation gradient computation.
 */

#if defined(LOCAL_MEM_SIZE) && defined(NO_LOCAL_MEM)
    #error NO_LOCAL_MEM has been set.
#endif

#ifndef HAVE_3D
    #include "../../../types/2D.h"
    #include "../../../KernelFunctions/Wendland2D.hcl"
#else
    #include "../../../types/3D.h"
    #include "../../../KernelFunctions/Wendland3D.hcl"
#endif

/** @brief Boundary elements deformation gradient computation due to the
 * interaction with the solid particles.
 *
 * Compute the gradient of the deformation vector:
 * \f[ \nabla \mathbf{r}^{*} = \mathbf{r}^{*} \otimes \nabla \f]
 *
 * @see https://en.wikipedia.org/wiki/Matrix_calculus
 * @see https://en.wikipedia.org/wiki/Outer_product
 * 
 * @param iset Set of particles index.
 * @param imove Moving flags.
 *   - imove = 2 for regular solid particles.
 *   - imove = 0 for sensors (ignored by this preset).
 *   - imove < 0 for boundary elements/particles.
 * @param r Position \f$ \mathbf{r} \f$.
 * @param r_r0 Deformation \f$ \mathbf{r}^{*} = \mathbf{r} - \mathbf{r}_0 \f$.
 * @param rho Density \f$ \rho \f$.
 * @param m Mass \f$ m \f$.
 * @param grad_r Gradient of the deformation \f$ \nabla \mathbf{r}^{*} \f$.
 * @param icell Cell where each particle is located.
 * @param ihoc Head of chain for each cell (first particle found).
 * @param N Number of particles.
 * @param n_cells Number of cells in each direction
 * @param BIfixed_iset Set of particles affected
 */
__kernel void vol(const __global uint* iset,
                  const __global int* imove,
                  const __global vec* r,
                  const __global vec* r_r0,
                  const __global float* rho,
                  const __global float* m,
                  __global matrix* grad_r,
                  // Link-list data
                  const __global uint *icell,
                  const __global uint *ihoc,
                  // Simulation data
                  uint N,
                  uivec4 n_cells,
                  uint BIfixed_iset)
{
    const uint i = get_global_id(0);
    const uint it = get_local_id(0);
    if(i >= N)
        return;
    if((imove[i] != -3) || (iset[i] != BIfixed_iset)){
        return;
    }

    const vec_xyz r_i = r[i].XYZ;

    // Initialize the output
    #ifndef LOCAL_MEM_SIZE
        #define _GRADR_ grad_r[i]
    #else
        #define _GRADR_ grad_r_l[it]
        __local matrix grad_r_l[LOCAL_MEM_SIZE];
        _GRADR_ = MAT_ZERO;
    #endif

    BEGIN_LOOP_OVER_NEIGHS(){
        if(imove[j] != 2){
            j++;
            continue;
        }
        const vec_xyz r_ij = r[j] - r_i.XYZ;
        const float q = length(r_ij) / H;
        if(q >= SUPPORT)
        {
            j++;
            continue;
        }
        {
            const float f_ij = kernelF(q) * CONF * m[j] / rho[j];
            _GRADR_ += outer(r_r0[j].XYZ, f_ij * r_ij);
        }
    }END_LOOP_OVER_NEIGHS()

    #ifdef LOCAL_MEM_SIZE
        grad_r[i] = _GRADR_;
    #endif
}

/** @brief Boundary elements deformation gradient resulting from the interaction
 * with the boundary.
 *
 * Compute the gradient of the deformation vector:
 * \f[ \nabla \mathbf{r}^{*} = \mathbf{r}^{*} \otimes \nabla \f]
 *
 * @see https://en.wikipedia.org/wiki/Matrix_calculus
 * @see https://en.wikipedia.org/wiki/Outer_product
 * 
 * @param iset Set of particles index.
 * @param imove Moving flags.
 *   - imove = 2 for regular solid particles.
 *   - imove = 0 for sensors (ignored by this preset).
 *   - imove < 0 for boundary elements/particles.
 * @param r Position \f$ \mathbf{r} \f$.
 * @param r_r0 Deformation \f$ \mathbf{r}^{*} = \mathbf{r} - \mathbf{r}_0 \f$.
 * @param normal Normal \f$ \mathbf{n} \f$.
 * @param m Area of the boundary element \f$ s \f$.
 * @param grad_r Gradient of the deformation \f$ \nabla \mathbf{r}^{*} \f$.
 * @param icell Cell where each particle is located.
 * @param ihoc Head of chain for each cell (first particle found).
 * @param N Number of particles.
 * @param n_cells Number of cells in each direction
 * @param BIfixed_iset Set of particles affected
 */
__kernel void area(const __global uint* iset,
                   const __global int* imove,
                   const __global vec* r,
                   const __global vec* r_r0,
                   const __global vec* normal,
                   const __global float* m,
                   __global matrix* grad_r,
                   // Link-list data
                   const __global uint *icell,
                   const __global uint *ihoc,
                   // Simulation data
                   uint N,
                   uivec4 n_cells,
                   uint BIfixed_iset)
{
    const uint i = get_global_id(0);
    const uint it = get_local_id(0);
    if(i >= N)
        return;
    if((imove[i] != -3) || (iset[i] != BIfixed_iset)){
        return;
    }

    const vec_xyz r_i = r[i].XYZ;

    // Initialize the output
    #ifndef LOCAL_MEM_SIZE
        #define _GRADR_ grad_r[i]
    #else
        #define _GRADR_ grad_r_l[it]
        __local matrix grad_r_l[LOCAL_MEM_SIZE];
        _GRADR_ = grad_r[i];
    #endif

    BEGIN_LOOP_OVER_NEIGHS(){
        if(imove[j] != -3){
            j++;
            continue;
        }
        const vec_xyz r_ij = r[j] - r_i.XYZ;
        const float q = length(r_ij) / H;
        if(q >= SUPPORT)
        {
            j++;
            continue;
        }
        {
            const vec_xyz n_j = normal[j].XYZ;  // Assumed outwarding oriented
            const float area_j = m[j];
            const float w_ij = kernelW(q) * CONW * area_j;
            _GRADR_ += outer(r_r0[j].XYZ, w_ij * n_j);
        }
    }END_LOOP_OVER_NEIGHS()

    #ifdef LOCAL_MEM_SIZE
        grad_r[i] = _GRADR_;
    #endif
}

/** @brief Renormalization of the deformation gradient.
 *
 * The main drawback of the boundary integrals formulation is the requirement
 * of the renormalization of the computed differentiqal operators, which is
 * destroying several conservation properties.
 *
 * @param iset Set of particles index.
 * @param imove Moving flags.
 *   - imove = 2 for regular solid particles.
 *   - imove = 0 for sensors.
 *   - imove < 0 for boundary elements/particles.
 * @param shepard Shepard term
 * \f$ \gamma(\mathbf{x}) = \int_{\Omega}
 *     W(\mathbf{y} - \mathbf{x}) \mathrm{d}\mathbf{x} \f$.
 * @param grad_r Gradient of the deformation \f$ \nabla \mathbf{r}^{*} \f$.
 * @param N Total number of particles and boundary elements.
 * @param BIfixed_iset Set of particles affected
 */
__kernel void renormalization(const __global uint* iset,
                              const __global int* imove,
                              const __global float* shepard,
                              __global matrix* grad_r,
                              uint N,
                              uint BIfixed_iset)
{
    uint i = get_global_id(0);
    if(i >= N)
        return;
    if((imove[i] != -3) || (iset[i] != BIfixed_iset)){
        return;
    }
    grad_r[i] /= shepard[i];
}

/*
 * @}
 */