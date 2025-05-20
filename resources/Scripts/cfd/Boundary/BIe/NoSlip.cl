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
 * @brief Boundary integral friction term.
 */

#if defined(LOCAL_MEM_SIZE) && defined(NO_LOCAL_MEM)
    #error NO_LOCAL_MEM has been set.
#endif

#include "resources/Scripts/types/types.h"
#include "resources/Scripts/KernelFunctions/Kernel.h"

/** @brief Performs the boundary friction effect on the fluid particles.
 * @param iset Set of particles index.
 * @param imove Moving flags.
 *   - imove > 0 for regular fluid particles.
 *   - imove = 0 for sensors.
 *   - imove < 0 for boundary elements/particles.
 * @param r Position \f$ \mathbf{r} \f$.
 * @param normal Normal \f$ \mathbf{n} \f$.
 * @param u Velocity \f$ \mathbf{u} \f$.
 * @param rho Density \f$ \rho \f$.
 * @param m Area of the boundary element \f$ s \f$.
 * @param lap_u Velocity laplacian \f$ \frac{\Delta \mathbf{u}}{rho} \f$.
 * @param icell Cell where each particle is located.
 * @param ihoc Head of chain for each cell (first particle found).
 * @param N Number of particles.
 * @param n_cells Number of cells in each direction
 * @param noslip_iset Set of boundary elements that should be considered
 * @param dr Distance between particles \f$ \Delta r \f$.
 */
__kernel void entry(const __global uint* iset,
                    const __global int* imove,
                    const __global vec* r,
                    const __global vec* normal,
                    const __global vec* u,
                    const __global float* rho,
                    const __global float* m,
                    __global vec* lap_u,
                    usize N,
                    uint noslip_iset,
                    float dr,
                    LINKLIST_LOCAL_PARAMS)
{
    const usize i = get_global_id(0);
    const usize it = get_local_id(0);
    if(i >= N)
        return;
    if(imove[i] != 1)
        return;

    const vec_xyz r_i = r[i].XYZ;
    const vec_xyz u_i = u[i].XYZ;
    const float rho_i = rho[i];

    // Initialize the output
    #ifndef LOCAL_MEM_SIZE
        #define _LAPU_ lap_u[i].XYZ
    #else
        #define _LAPU_ lap_u_l[it]
        __local vec_xyz lap_u_l[LOCAL_MEM_SIZE];
        _LAPU_ = lap_u[i].XYZ;
    #endif

    const usize c_i = icell[i];
    BEGIN_NEIGHS(c_i, N, n_cells, icell, ihoc){
        if((imove[j] != -3) || (iset[j] != noslip_iset)){
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
            const float area_j = m[j];
            const float w_ij = kernelW(q) * CONW * area_j;

            const vec_xyz du = u[j].XYZ - u_i;
            const float dr_n = max(fabs(dot(r_ij, n_j)), dr);
            const vec_xyz du_t = du - dot(du, n_j) * n_j;
            _LAPU_ += 2.f * w_ij / (rho_i * dr_n) * du_t;
        }
    }END_NEIGHS()

    #ifdef LOCAL_MEM_SIZE
        lap_u[i].XYZ = _LAPU_;
    #endif
}

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
 * @param iset Set of particles index.
 * @param imove Moving flags.
 *   - imove > 0 for regular fluid particles.
 *   - imove = 0 for sensors.
 *   - imove < 0 for boundary elements/particles.
 * @param force_visc Viscous force on each boundary element
 * @param r Position \f$ \mathbf{r} \f$.
 * @param normal Normal \f$ \mathbf{n} \f$.
 * @param u Velocity \f$ \mathbf{u} \f$.
 * @param rho Density \f$ \rho \f$.
 * @param m Mass \f$ m \f$.
 * @param visc_dyn Dynamic viscosity \f$ \mu \f$.
 * @param icell Cell where each particle is located.
 * @param ihoc Head of chain for each cell (first particle found).
 * @param N Number of particles.
 * @param n_cells Number of cells in each direction
 * @param dr Distance between particles \f$ \Delta r \f$.
 */
__kernel void force(const __global uint* iset,
                    const __global int* imove,
                    __global vec* force_visc,
                    const __global vec* r,
                    const __global vec* normal,
                    const __global vec* u,
                    const __global float* rho,
                    const __global float* m,
                     __constant float* visc_dyn,
                    usize N,
                    float dr,
                    LINKLIST_LOCAL_PARAMS)
{
    const usize i = get_global_id(0);
    const usize it = get_local_id(0);
    if(i >= N)
        return;
    if(imove[i] != -3)
        return;
    
    const vec_xyz r_i = r[i].XYZ;
    const vec_xyz n_i = normal[i].XYZ;
    const vec_xyz u_i = u[i].XYZ;
    const float area_i = m[i];

    // Initialize the output
    #ifndef LOCAL_MEM_SIZE
        #define _F_ force_visc[i].XYZ
    #else
        #define _F_ f_l[it]
        __local vec_xyz f_l[LOCAL_MEM_SIZE];
    #endif
    _F_ = VEC_ZERO.XYZ;

    const usize c_i = icell[i];
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

        {
            const float rho_j = rho[j];
            const float m_j = m[j];
            const float visc_dyn_j = visc_dyn[iset[j]];
            const vec_xyz du = u[j].XYZ - u_i;
            const float w_ij = kernelW(q) * CONW * area_i;
            
            const float dr_n = max(fabs(dot(r_ij, n_i)), dr);
            const vec_xyz du_t = du - dot(du, n_i) * n_i;
    
            _F_ += 2.f * visc_dyn_j * m_j * w_ij / (rho_j * dr_n) * du_t;
        }
    }END_NEIGHS()

    #ifdef LOCAL_MEM_SIZE
        force_visc[i].XYZ = _F_;
    #endif
}

/** @brief Compute the force and torque at the boundary.
 *
 * @param imove Moving flags.
 *   - imove > 0 for regular fluid particles.
 *   - imove = 0 for sensors.
 *   - imove < 0 for boundary elements/particles.
 * @param r Position \f$ \mathbf{r} \f$.
 * @param force_visc Viscous force on the boundary element
 * @param moment_visc Viscous moment on the boundary element
 * @param forces_r Point with respect the moments are computed
 * \f$ \mathbf{r}_0 \f$.
 * @param N Number of particles.
 */
__kernel void moment(const __global int* imove,
                     const __global vec* r,
                     const __global vec* force_visc,
                     __global vec4* moment_visc,
                     vec forces_r,
                     usize N)
{
    const usize i = get_global_id(0);
    if(i >= N)
        return;
    if(imove[i] != -3)
        return;

    // Expand the force and the arm as 3D variables
    vec4 F = (vec4)(0.f);
    vec4 R = (vec4)(0.f);
    F.XYZ = force_visc[i].XYZ;
    R.XYZ = r[i].XYZ - forces_r.XYZ;
    moment_visc[i] = cross(R, F);    
}

/** @brief Filter out other boundaries from the viscous forces
 *
 * This is used to compute the force on an specific iset
 * @param iset Particle set
 * @param imove Moving flags.
 *   - imove > 0 for regular fluid particles.
 *   - imove = 0 for sensors.
 *   - imove < 0 for boundary elements/particles.
 * @param force_visc Viscous force on the boundary element
 * @param moment_visc Viscous moment on the boundary element
 * @param forces_iset Particles set of interest.
 * @param N Number of particles.
 */
__kernel void filter_force(const __global uint* iset,
                           const __global int* imove,
                           __global vec* force_visc,
                           __global vec4* moment_visc,
                           unsigned int forces_iset,
                           usize N)
{
    const usize i = get_global_id(0);
    if(i >= N)
        return;
    if(imove[i] != -3) {
        return;
    }

    if(iset[i] != forces_iset) {
        force_visc[i] = VEC_ZERO;
        moment_visc[i] = (vec4)(0.f);
    }
}
