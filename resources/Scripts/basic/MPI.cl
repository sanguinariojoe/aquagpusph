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

#include "resources/Scripts/types/types.h"

/** @brief Compute the particles sending mask
 *
 * The sending mask is just an array where it is stored the process each
 * particle will be sent for, just in case it shall be sent at all
 *
 * This method assumes equally spaced interface planes
 *
 * No special treatment is required at the bounding processes, since setting the
 * mask with an invalid process is just ignoring it. However, such a situation
 * would indicate problems. Thus, it is strongly recommended to use MPI sync
 * with the computational domain module, in such a way mpi_planes_orig matchs
 * the minimum computational domain point
 *
 * @param mpi_mask Output processes mask
 * @param imove Moving flags
 *   - imove > 0 for regular fluid particles
 *   - imove = 0 for sensors
 *   - imove < 0 for boundary elements/particles
 * @param r_in Position \f$ \mathbf{r} \f$
 * @param mpi_rank MPI process index
 * @param mpi_planes_orig Center of the first interface (between procs) plane
 * @param mpi_planes_dist Distance between interface planes
 * @param N Number of particles
 */
__kernel void mask_planes(__global unsigned int* mpi_mask,
                          const __global unsigned int* imove,
                          const __global vec* r_in,
                          unsigned int mpi_rank,
                          vec mpi_planes_orig,
                          vec mpi_planes_dist,
                          unsigned int N)
{
    unsigned int i = get_global_id(0);
    if(i >= N)
        return;
    if(imove[i] <= 0)
        return;

    // We need to check the distance, so we need to normalize
    const vec_xyz normal = normalize(mpi_planes_dist.XYZ);
    
#ifdef CL_VERSION_2_0
    __attribute__((opencl_unroll_hint(2)))
#endif
    for (unsigned int i=0; i < 2; i++) {
        const vec_xyz orig = mpi_planes_orig.XYZ +
                             (mpi_rank + i) * mpi_planes_dist.XYZ;
        const float normal_sign = (float)(mpi_rank) - (float)(i) - 1.f;
        const float dist = dot(r_in[i].XYZ - orig, normal_sign * normal);
        if(dist > -SUPPORT * H) {
            mpi_mask[i] = mpi_rank + i - 1;
        }
    }
}

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
 * @param mpi_r Position \f$ \mathbf{r} \f$ MPI copy
 * @param mpi_u Velocity \f$ \mathbf{u} \f$ MPI copy
 * @param mpi_dudt Velocity rate of change \f$ \frac{d \mathbf{u}}{d t} \f$ MPI copy
 * @param mpi_rho Density \f$ \rho \f$ MPI copy
 * @param mpi_drhodt Density rate of change \f$ \frac{d \rho}{d t} \f$ MPI copy
 * @param r_in Position \f$ \mathbf{r} \f$
 * @param u_in Velocity \f$ \mathbf{u} \f$
 * @param dudt_in Velocity rate of change \f$ \frac{d \mathbf{u}}{d t} \f$
 * @param rho_in Density \f$ \rho \f$
 * @param drhodt_in Density rate of change \f$ \frac{d \rho}{d t} \f$
 * @param N Number of particles
 */
__kernel void copy(__global vec* mpi_r,
                   __global vec* mpi_u,
                   __global vec* mpi_dudt,
                   __global float* mpi_rho,
                   __global float* mpi_drhodt,
                   const __global vec* r_in,
                   const __global vec* u_in,
                   const __global vec* dudt_in,
                   const __global float* rho_in,
                   const __global float* drhodt_in,
                   unsigned int N)
{
    unsigned int i = get_global_id(0);
    if(i >= N)
        return;

    mpi_r[i] = r_in[i];
    mpi_u[i] = u_in[i];
    mpi_dudt[i] = dudt_in[i];
    mpi_rho[i] = rho_in[i];
    mpi_drhodt[i] = drhodt_in[i];
}

/** @brief Drop particles belonging to different processes
 *
 * Transfer the particles belonging to different processes, i.e. outside the
 * current process bounding planes, to "the buffer". It should be always bear in
 * mind that the buffer particles are not made available until next time step
 *
 * @param imove Moving flags
 *   - imove > 0 for regular fluid particles
 *   - imove = 0 for sensors
 *   - imove < 0 for boundary elements/particles
 * @param r Position \f$ \mathbf{r} \f$
 * @param mpi_mask Output processes mask
 * @param mpi_rank MPI process index
 * @param mpi_planes_orig Center of the first interface (between procs) plane
 * @param mpi_planes_dist Distance between interface planes
 * @param domain_max Top-left-frontal corner of the computational domain
 * @param N Number of particles
 */
__kernel void drop_planes(__global unsigned int* imove,
                          __global vec* r,
                          unsigned int mpi_rank,
                          vec mpi_planes_orig,
                          vec mpi_planes_dist,
                          vec domain_max,
                          unsigned int N)
{
    unsigned int i = get_global_id(0);
    if(i >= N)
        return;
    if(imove[i] <= 0)
        return;

    // We just need the sign, so no need to normalize
    const vec_xyz normal = mpi_planes_dist.XYZ;
    
#ifdef CL_VERSION_2_0
    __attribute__((opencl_unroll_hint(2)))
#endif
    for (unsigned int i=0; i < 2; i++) {
        const vec_xyz orig = mpi_planes_orig.XYZ +
                             (mpi_rank + i) * mpi_planes_dist.XYZ;
        const float normal_sign = (float)(mpi_rank) - (float)(i) - 1.f;
        if(dot(r[i].XYZ - orig, normal_sign * normal) > 0) {
            r[i] = domain_max;
            imove[i] = -256;
        }
    }
}

/** @brief Add the particles received from other processes
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
 * @param mpi_mask Incoming processes mask
 * @param mpi_r Position \f$ \mathbf{r} \f$ MPI copy
 * @param mpi_u Velocity \f$ \mathbf{u} \f$ MPI copy
 * @param mpi_dudt Velocity rate of change \f$ \frac{d \mathbf{u}}{d t} \f$ MPI copy
 * @param mpi_rho Density \f$ \rho \f$ MPI copy
 * @param mpi_drhodt Density rate of change \f$ \frac{d \rho}{d t} \f$ MPI copy
 * @param mpi_rank MPI process index
 * @param nbuffer Number of buffer particles
 * @param N Number of particles
 */
__kernel void restore(__global vec* r_in,
                      __global vec* u_in,
                      __global vec* dudt_in,
                      __global float* rho_in,
                      __global float* drhodt_in,
                      __global unsigned int* imove,
                      const __global unsigned int* mpi_mask,
                      const __global vec* mpi_r,
                      const __global vec* mpi_u,
                      const __global vec* mpi_dudt,
                      const __global float* mpi_rho,
                      const __global float* mpi_drhodt,
                      unsigned int mpi_rank,
                      unsigned int nbuffer,
                      unsigned int N)
{
    unsigned int i = get_global_id(0);
    if(i >= N)
        return;
    if(mpi_mask[i] == mpi_rank)
        return;

    // The new particles are at the beginning of the mpi_* arrays. However, we
    // would consume particles from the buffer, which is at the end of the list
    // of particles. Never consume particles straight from the end, or the
    // buffer will be ruined for future executions
    const unsigned int i_out = N - nbuffer + i;
    imove[i_out] = 1;
    r_in[i_out] = mpi_r[i];
    u_in[i_out] = mpi_u[i];
    dudt_in[i_out] = mpi_dudt[i];
    rho_in[i_out] = mpi_rho[i];
    drhodt_in[i_out] = mpi_drhodt[i];
}
