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
 * @brief MPI syncing point
 */

#if defined(LOCAL_MEM_SIZE) && defined(NO_LOCAL_MEM)
    #error NO_LOCAL_MEM has been set.
#endif

#include "resources/Scripts/types/types.h"
#include "resources/Scripts/KernelFunctions/Kernel.h"

#if __LAP_FORMULATION__ == __LAP_MONAGHAN__
    #ifndef HAVE_3D
        #define __CLEARY__ 8.f
    #else
        #define __CLEARY__ 10.f
    #endif
#endif

/** @brief Copy the data
 *
 * To make the operation as asynchronous as possible, we are copying all the
 * available data, regardless it is useful or not. That is a bit computationally
 * less efficient, but allows to start transfering data ASAP, while we still
 * operate to neglect the particles
 *
 * It is intended that the fields copied here are the same subsequently
 * exchanged later with the mpi-sync tool
 *
 * @param mpi_iset Particles set MPI copy
 * @param mpi_r Position \f$ \mathbf{r} \f$ MPI copy
 * @param mpi_u Velocity \f$ \mathbf{u} \f$ MPI copy
 * @param mpi_dudt Velocity rate of change \f$ \frac{d \mathbf{u}}{d t} \f$ MPI copy
 * @param mpi_rho Density \f$ \rho \f$ MPI copy
 * @param mpi_drhodt Density rate of change \f$ \frac{d \rho}{d t} \f$ MPI copy
 * @param mpi_m Mass \f$ m \f$ MPI copy
 * @param iset Particles set
 * @param r Position \f$ \mathbf{r} \f$
 * @param u Velocity \f$ \mathbf{u} \f$
 * @param dudt Velocity rate of change \f$ \frac{d \mathbf{u}}{d t} \f$
 * @param rho Density \f$ \rho \f$
 * @param drhodt Density rate of change \f$ \frac{d \rho}{d t} \f$
 * @param m Mass \f$ m \f$
 * @param N Number of particles
 */
__kernel void copy(__global unsigned int* mpi_iset,
                   __global vec* mpi_r,
                   __global vec* mpi_u,
                   __global vec* mpi_dudt,
                   __global float* mpi_rho,
                   __global float* mpi_drhodt,
                   __global float* mpi_m,
                   const __global unsigned int* iset,
                   const __global vec* r,
                   const __global vec* u,
                   const __global vec* dudt,
                   const __global float* rho,
                   const __global float* drhodt,
                   const __global float* m,
                   usize N)
{
    const usize i = get_global_id(0);
    if(i >= N)
        return;

    mpi_iset[i] = iset[i];
    mpi_r[i] = r[i];
    mpi_u[i] = u[i];
    mpi_dudt[i] = dudt[i];
    mpi_rho[i] = rho[i];
    mpi_drhodt[i] = drhodt[i];
    mpi_m[i] = m[i];
}

/** @brief Append the particles received from other processes
 *
 * The restored particles has always imove=0 flag
 *
 * @param imove Moving flags
 *   - imove > 0 for regular fluid particles
 *   - imove = 0 for sensors
 *   - imove < 0 for boundary elements/particles
 * @param r_in Position \f$ \mathbf{r} \f$
 * @param u_in Velocity \f$ \mathbf{u} \f$
 * @param dudt_in Velocity rate of change \f$ \frac{d \mathbf{u}}{d t} \f$
 * @param rho_in Density \f$ \rho \f$
 * @param drhodt_in Density rate of change \f$ \frac{d \rho}{d t} \f$
 * @param m Mass \f$ m \f$
 * @param mpi_local_mask Incoming processes mask
 * @param mpi_r Position \f$ \mathbf{r} \f$ MPI copy
 * @param mpi_u Velocity \f$ \mathbf{u} \f$ MPI copy
 * @param mpi_dudt Velocity rate of change \f$ \frac{d \mathbf{u}}{d t} \f$ MPI copy
 * @param mpi_rho Density \f$ \rho \f$ MPI copy
 * @param mpi_drhodt Density rate of change \f$ \frac{d \rho}{d t} \f$ MPI copy
 * @param mpi_m Mass \f$ m \f$ MPI copy
 * @param mpi_rank MPI process index
 * @param nbuffer Number of buffer particles
 * @param N Number of particles
 */
__kernel void append(__global unsigned int* iset,
                     __global vec* r,
                     __global vec* u,
                     __global vec* dudt,
                     __global float* rho,
                     __global float* drhodt,
                     __global float* m,
                     __global int* imove,
                     const __global size_t* mpi_local_mask,
                     const __global unsigned int* mpi_iset,
                     const __global vec* mpi_r,
                     const __global vec* mpi_u,
                     const __global vec* mpi_dudt,
                     const __global float* mpi_rho,
                     const __global float* mpi_drhodt,
                     const __global float* mpi_m,
                     unsigned int mpi_rank,
                     usize nbuffer,
                     usize N)
{
    const usize i = get_global_id(0);
    if(i >= N)
        return;
    if(mpi_local_mask[i] == mpi_rank)
        return;

    // The new particles are at the beginning of the mpi_* arrays. However, we
    // would consume particles from the buffer, which is at the end of the list
    // of particles. Never consume particles straight from the end, or the
    // buffer will be ruined for future executions
    const usize i_out = N - nbuffer + i;
    imove[i_out] = 1;
    iset[i_out] = mpi_iset[i];
    r[i_out] = mpi_r[i];
    u[i_out] = mpi_u[i];
    dudt[i_out] = mpi_dudt[i];
    rho[i_out] = mpi_rho[i];
    drhodt[i_out] = mpi_drhodt[i];
    m[i_out] = mpi_m[i];
}

/** @brief Remove the particles sent to other processes
 *
 * @param imove Moving flags
 *   - imove > 0 for regular fluid particles
 *   - imove = 0 for sensors
 *   - imove < 0 for boundary elements/particles
 * @param r Position \f$ \mathbf{r} \f$
 * @param mpi_local_mask Incoming processes mask
 * @param mpi_rank The current process id
 * @param N Number of particles
 */
__kernel void remove(__global int* imove,
                     __global vec* r,
                     const __global size_t* mpi_local_mask,
                     unsigned int mpi_rank,
                     usize N)
{
    const usize i = get_global_id(0);
    if(i >= N)
        return;
    if(mpi_local_mask[i] == mpi_rank)
        return;

    r[i] = VEC_INFINITY;
    // We mark the particle as -256, to be turned into a buffer one (-255)
    // later on
    imove[i] = -256;
}

/** @brief Backup the position of the MPI particles to sort that later.
 *
 * Just the particles coming from other processes are backuped
 *
 * @param mpi_neigh_mask Incoming processes mask
 * @param mpi_r Position \f$ \mathbf{r} \f$ MPI copy
 * @param mpi_r_in Position \f$ \mathbf{r} \f$ MPI copy
 * @param r_max The maximum position detected by the Link-List
 * @param mpi_rank The current process id
 * @param N Number of particles
 * @param n_radix Number of elements on @p mpi_r_in
 */
__kernel void backup_r(const __global size_t* mpi_neigh_mask,
                       const __global vec* mpi_r,
                       __global vec* mpi_r_in,
                       vec r_max,
                       unsigned int mpi_rank,
                       usize N,
                       usize n_radix)
{
    const usize i = get_global_id(0);
    if(i >= n_radix)
        return;
    if((i >= N) || (mpi_neigh_mask[i] == mpi_rank)) {
        mpi_r_in[i] = r_max;
        return;
    }
    mpi_r_in[i] = mpi_r[i];
}

/** @brief Sort all the particle variables by the cell indexes
 *
 * Due to the large number of registers consumed (21 global memory arrays are
 * loaded), it is safer carrying out the sorting process in to stages. This
 * is the first stage.
 *
 * @param mpi_r_in Unsorted position \f$ \mathbf{r} \f$.
 * @param mpi_r Sorted position \f$ \mathbf{r} \f$.
 * @param mpi_u_in Unsorted velocity \f$ \mathbf{u} \f$.
 * @param mpi_u Sorted velocity \f$ \mathbf{u} \f$.
 * @param mpi_rho_in Unsorted density \f$ \rho \f$.
 * @param mpi_rho Sorted density \f$ \rho \f$.
 * @param mpi_m_in Unsorted mass \f$ m \f$.
 * @param mpi_m Sorted mass \f$ m \f$.
 * @param mpi_id_sorted Permutations list from the unsorted space to the sorted
 * one.
 * @param N Number of particles.
 */
__kernel void sort(const __global unsigned int *mpi_iset_in, __global unsigned int *mpi_iset,
                   const __global vec *mpi_r_in, __global vec *mpi_r,
                   const __global vec *mpi_u_in, __global vec *mpi_u,
                   const __global float *mpi_rho_in, __global float *mpi_rho,
                   const __global float *mpi_m_in, __global float *mpi_m,
                   const __global usize *mpi_id_sorted,
                   usize N)
{
    const usize i = get_global_id(0);
    if(i >= N)
        return;

    const usize i_out = mpi_id_sorted[i];

    mpi_iset[i_out] = mpi_iset_in[i];
    mpi_r[i_out] = mpi_r_in[i];
    mpi_u[i_out] = mpi_u_in[i];
    mpi_rho[i_out] = mpi_rho_in[i];
    mpi_m[i_out] = mpi_m_in[i];
}

/** @brief Stiff Equation Of State (EOS) computation
 *
 * The equation of state relates the pressure and the density fields,
 * \f$ p = p_0 + c_s^2 \left(\rho - \rho_0 \right) \f$
 *
 * @param mpi_iset Set of particles index.
 * @param mpi_rho Density \f$ \rho_{n+1/2} \f$.
 * @param mpi_p Pressure \f$ p_{n+1/2} \f$.
 * @param refd Density of reference of the fluid \f$ \rho_0 \f$.
 * @param N Number of particles.
 * @param cs Speed of sound \f$ c_s \f$.
 * @param p0 Background pressure \f$ p_0 \f$.
 */
__kernel void eos(__global unsigned int* mpi_iset,
                  __global float* mpi_rho,
                  __global float* mpi_p,
                  __constant float* refd,
                  usize N,
                  float cs,
                  float p0)
{
    const usize i = get_global_id(0);
    if(i >= N)
        return;

    mpi_p[i] = p0 + cs * cs * (mpi_rho[i] - refd[mpi_iset[i]]);
}

/** @brief Fluid particles interactions with the neighbours from other
 * processes.
 *
 * Compute the differential operators involved in the numerical scheme, taking
 * into account just the fluid-fluid interactions with the MPI revceived.
 *
 * @param imove Moving flags.
 *   - imove > 0 for regular fluid particles.
 *   - imove = 0 for sensors.
 *   - imove < 0 for boundary elements/particles.
 * @param r Position \f$ \mathbf{r} \f$.
 * @param u Velocity \f$ \mathbf{u} \f$.
 * @param rho Density \f$ \rho \f$.
 * @param p Pressure \f$ p \f$.
 * @param mpi_r Neighs position \f$ \mathbf{r} \f$.
 * @param mpi_u Neighs velocity \f$ \mathbf{u} \f$.
 * @param mpi_rho Neighs density \f$ \rho \f$.
 * @param mpi_p Neighs pressure \f$ p \f$.
 * @param mpi_m Neighs mass \f$ m \f$.
 * @param grad_p Pressure gradient \f$ \frac{\nabla p}{rho} \f$.
 * @param lap_u Velocity laplacian \f$ \frac{\Delta \mathbf{u}}{rho} \f$.
 * @param div_u Velocity divergence \f$ \rho \nabla \cdot \mathbf{u} \f$.
 * @param N Number of particles.
 * @param icell Cell where each particle is located.
 * @param mpi_icell Cell where each neighbour is located.
 * @param mpi_ihoc Head of chain for each cell (first neighbour found).
 * @param n_cells Number of cells in each direction
 */
__kernel void interactions(const __global int* imove,
                           const __global vec* r,
                           const __global vec* u,
                           const __global float* rho,
                           const __global float* p,
                           const __global vec* mpi_r,
                           const __global vec* mpi_u,
                           const __global float* mpi_rho,
                           const __global float* mpi_p,
                           const __global float* mpi_m,
                           __global vec* grad_p,
                           __global vec* lap_u,
                           __global float* div_u,
                           usize N,
                           LINKLIST_REMOTE_PARAMS)
{
    const usize i = get_global_id(0);
    const usize it = get_local_id(0);
    if(i >= N)
        return;
    if(imove[i] != 1){
        return;
    }

    const vec_xyz r_i = r[i].XYZ;
    const vec_xyz u_i = u[i].XYZ;
    const float p_i = p[i];
    const float rho_i = rho[i];

    // Initialize the output
    #ifndef LOCAL_MEM_SIZE
        #define _GRADP_ grad_p[i].XYZ
        #define _LAPU_ lap_u[i].XYZ
        #define _DIVU_ div_u[i]
    #else
        #define _GRADP_ grad_p_l[it]
        #define _LAPU_ lap_u_l[it]
        #define _DIVU_ div_u_l[it]
        __local vec_xyz grad_p_l[LOCAL_MEM_SIZE];
        __local vec_xyz lap_u_l[LOCAL_MEM_SIZE];
        __local float div_u_l[LOCAL_MEM_SIZE];
        _GRADP_ = VEC_ZERO.XYZ;
        _LAPU_ = VEC_ZERO.XYZ;
        _DIVU_ = 0.f;
    #endif

    const usize c_i = icell[i];
    BEGIN_NEIGHS(c_i, N, n_cells, mpi_icell, mpi_ihoc){
        const vec_xyz r_ij = mpi_r[j].XYZ - r_i;
        const float q = length(r_ij) / H;
        if(q >= SUPPORT)
        {
            j++;
            continue;
        }
        {
            const float rho_j = mpi_rho[j];
            const float p_j = mpi_p[j];
            const float udr = dot(mpi_u[j].XYZ - u_i, r_ij);
            const float f_ij = kernelF(q) * CONF * mpi_m[j];

            _GRADP_ += (p_i + p_j) / (rho_i * rho_j) * f_ij * r_ij;

            #if __LAP_FORMULATION__ == __LAP_MONAGHAN__
                const float r2 = (q * q + 0.01f) * H * H;
                _LAPU_ += f_ij * __CLEARY__ * udr / (r2 * rho_i * rho_j) * r_ij;
            #elif __LAP_FORMULATION__ == __LAP_MORRIS__
                _LAPU_ += f_ij * 2.f / (rho_i * rho_j) * (u[j].XYZ - u_i);
            #else
                #error Unknown Laplacian formulation: __LAP_FORMULATION__
            #endif

            _DIVU_ += udr * f_ij * rho_i / rho_j;
        }
    }END_NEIGHS()

    #ifdef LOCAL_MEM_SIZE
        grad_p[i].XYZ += _GRADP_;
        lap_u[i].XYZ += _LAPU_;
        div_u[i] += _DIVU_;
    #endif
}
