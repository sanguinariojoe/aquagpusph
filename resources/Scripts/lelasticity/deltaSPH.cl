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
 * @brief delta-SPH methods, including the correction terms.
 */

#if defined(LOCAL_MEM_SIZE) && defined(NO_LOCAL_MEM)
    #error NO_LOCAL_MEM has been set.
#endif

#ifndef HAVE_3D
    #include "../types/2D.h"
    #include "../KernelFunctions/Wendland2D.hcl"
#else
    #include "../types/3D.h"
    #include "../KernelFunctions/Wendland3D.hcl"
#endif

/** @brief Simple hidrostatic based correction term.
 *
 * @param iset Set of particles index.
 * @param imove Moving flags.
 *   - imove = 2 for regular solid particles.
 *   - imove = 0 for sensors (ignored by this preset).
 *   - imove < 0 for boundary elements/particles.
 * @param lap_p_corr Correction term for the Morrison Laplacian formula.
 * @param refd Density of reference of the fluid \f$ \rho_0 \f$.
 * @param N Number of particles.
 * @param g Gravity acceleration \f$ \mathbf{g} \f$.
 */
__kernel void simple(const __global unsigned int* iset,
                     const __global int* imove,
                     __global vec* lap_p_corr,
                     __constant float* refd,
                     unsigned int N,
                     vec g)
{
    unsigned int i = get_global_id(0);
    if(i >= N)
        return;
    if(imove[i] != 2)
        return;

    lap_p_corr[i] = refd[iset[i]] * g;
}

/** @brief MLS based correction term.
 *
 * @param imove Moving flags.
 *   - imove > 0 for regular fluid particles.
 *   - imove = 0 for sensors.
 *   - imove < 0 for boundary elements/particles.
 * @param r Position \f$ \mathbf{r} \f$.
 * @param rho Density \f$ \rho \f$.
 * @param m Mass \f$ m \f$.
 * @param p Pressure \f$ p \f$.
 * @param mls_solid Kernel MLS transformation matrix \f$ L \f$.
 * @param lap_p_corr Correction term for the Morrison Laplacian formula.
 * @param icell Cell where each particle is located.
 * @param ihoc Head of chain for each cell (first particle found).
 * @param N Number of particles.
 * @param n_cells Number of cells in each direction
 */
__kernel void full(const __global int* imove,
                   const __global vec* r,
                   const __global float* rho,
                   const __global float* m,
                   const __global float* p,
                   const __global matrix* mls_solid,
                   __global vec* lap_p_corr,
                   const __global uint *icell,
                   const __global uint *ihoc,
                   uint N,
                   uivec4 n_cells)
{
    const uint i = get_global_id(0);
    const uint it = get_local_id(0);
    if(i >= N)
        return;
    if(imove[i] != 2){
        return;
    }

    const vec_xyz r_i = r[i].XYZ;
    const float p_i = p[i];

    // Initialize the output
    #ifndef LOCAL_MEM_SIZE
        #define _GRADP_ lap_p_corr[i].XYZ
    #else
        #define _GRADP_ lap_p_corr_l[it]
        __local vec_xyz lap_p_corr_l[LOCAL_MEM_SIZE];
        _GRADP_ = VEC_ZERO.XYZ;
    #endif

    BEGIN_LOOP_OVER_NEIGHS(){
        if( (i == j) || (imove[j] != 2)){
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
            _GRADP_ += (p[j] - p_i) * f_ij * r_ij;
        }
    }END_LOOP_OVER_NEIGHS()

    lap_p_corr[i] = MATRIX_DOT(mls_solid[i], _GRADP_);
}

/** @brief Laplacian of the pressure computation.
 *
 * @param imove Moving flags.
 *   - imove > 0 for regular fluid particles.
 *   - imove = 0 for sensors.
 *   - imove < 0 for boundary elements/particles.
 * @param r Position \f$ \mathbf{r} \f$.
 * @param rho Density \f$ \rho \f$.
 * @param m Mass \f$ m \f$.
 * @param p Pressure \f$ p \f$.
 * @param lap_p_corr Correction term for the Morrison Laplacian formula.
 * @param lap_p Pressure laplacian \f$ \Delta p \f$.
 * @param icell Cell where each particle is located.
 * @param ihoc Head of chain for each cell (first particle found).
 * @param N Number of particles.
 * @param n_cells Number of cells in each direction
 */
__kernel void lapp(const __global int* imove,
                   const __global vec* r,
                   const __global float* rho,
                   const __global float* m,
                   const __global float* p,
                   const __global vec* lap_p_corr,
                   __global float* lap_p,
                   const __global uint *icell,
                   const __global uint *ihoc,
                   uint N,
                   uivec4 n_cells)
{
    const uint i = get_global_id(0);
    const uint it = get_local_id(0);
    if(i >= N)
        return;
    if(imove[i] != 2){
        return;
    }

    const vec_xyz r_i = r[i].XYZ;
    const float p_i = p[i];
    const vec_xyz gradp_i = lap_p_corr[i].XYZ;

    // Initialize the output
    #ifndef LOCAL_MEM_SIZE
        #define _LAPP_ lap_p[i]
    #else
        #define _LAPP_ lap_p_l[it]
        __local float lap_p_l[LOCAL_MEM_SIZE];
        _LAPP_ = 0.f;
    #endif

    BEGIN_LOOP_OVER_NEIGHS(){
        if( (i == j) || (imove[j] != 2)){
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
            const vec_xyz gradp_ij = lap_p_corr[j].XYZ + gradp_i;
            const float f_ij = kernelF(q) * CONF * m[j] / rho[j];
            _LAPP_ += ((p[j] - p_i) - 0.5f * dot(gradp_ij, r_ij)) * f_ij;
        }
    }END_LOOP_OVER_NEIGHS()

    #ifdef LOCAL_MEM_SIZE
        lap_p[i] = _LAPP_;
    #endif
}

/** @brief Density variation rates delta-SPH term.
 *
 * @param iset Set of particles index.
 * @param imove Moving flags.
 *   - imove = 2 for regular solid particles.
 *   - imove = 0 for sensors (ignored by this preset).
 *   - imove < 0 for boundary elements/particles.
 * @param rho Density \f$ \rho_{n+1} \f$.
 * @param lap_p Pressure laplacian \f$ \Delta p \f$.
 * @param drhodt Density rate of change \f$ \frac{d \rho}{d t} \f$.
 * @param refd Density of reference of the fluid \f$ \rho_0 \f$.
 * @param delta Diffusive term \f$ \delta \f$ multiplier.
 * @param N Number of particles.
 * @param dt Time step \f$ \Delta t \f$.
 */
__kernel void deltaSPH(const __global unsigned int* iset,
                    const __global int* imove,
                    const __global float* rho,
                    const __global float* lap_p,
                    __global float* drhodt,
                    __constant float* refd,
                    __constant float* delta,
                    unsigned int N,
                    float dt)
{
    unsigned int i = get_global_id(0);
    if(i >= N)
        return;
    if(imove[i] != 2)
        return;

    const uint set_i = iset[i];
    const float delta_f = delta[set_i] * dt * rho[i] / refd[set_i];

    drhodt[i] += delta_f * lap_p[i];
}

/*
 * @}
 */