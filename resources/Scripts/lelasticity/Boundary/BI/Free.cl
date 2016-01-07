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
 * @brief Fixed boundary elements methods.
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

/** @brief Pressure initialization.
 *
 * Since the pressure field is resulting from the density field (applying the
 * EOS), letting free the pressure at the boundary may be dangerous due to the
 * value must be unassigned. Therefore the pressure field will be initialized
 * as the background pressure, letting the user to don't sdpecifically assign a
 * pressure value without crashing the simulation. 
 *
 * @param iset Set of particles index.
 * @param imove Moving flags.
 *   - imove = 2 for regular solid particles.
 *   - imove = 0 for sensors.
 *   - imove < 0 for boundary elements/particles.
 * @param p Pressure \f$ p \f$.
 * @param p0 Background pressure \f$ p_0 \f$.
 * @param N Total number of particles and boundary elements.
 * @param BIfree_iset Set of particles affected
 */
__kernel void press(const __global uint* iset,
                    const __global int* imove,
                    __global float* p,
                    uint N,
                    float p0,
                    uint BIfree_iset)
{
    uint i = get_global_id(0);
    if(i >= N)
        return;
    if((imove[i] != -3) || (iset[i] != BIfree_iset)){
        return;
    }

    p[i] = p0;
}


/** @brief Pressure and stress deviation interpolation at the boundary elements.
 *
 * The values are computed using just the fluid information. The resulting
 * interpolated values are not renormalized yet.
 *
 * Just the elements with the flag imove = -3 are considered boundary elements.
 *
 * @param iset Set of particles index.
 * @param imove Moving flags.
 *   - imove > 0 for regular fluid particles.
 *   - imove = 0 for sensors.
 *   - imove < 0 for boundary elements/particles.
 * @param r Position \f$ \mathbf{r} \f$.
 * @param m Mass \f$ m \f$.
 * @param rho Density \f$ \rho \f$.
 * @param u Velocity \f$ \mathbf{u} \f$.
 * @param S Deviatory stress \f$ S \f$.
 * @param icell Cell where each particle is located.
 * @param ihoc Head of chain for each cell (first particle found).
 * @param N Number of particles.
 * @param n_cells Number of cells in each direction
 * @param BIfree_iset Set of particles affected
 */
__kernel void interpolation(const __global uint* iset,
                            const __global int* imove,
                            const __global vec* r,
                            const __global float* m,
                            const __global float* rho,
                            __global vec* u,
                            __global matrix* S,
                            // Link-list data
                            const __global uint *icell,
                            const __global uint *ihoc,
                            // Simulation data
                            uint N,
                            uivec4 n_cells,
                            uint BIfree_iset)
{
    const uint i = get_global_id(0);
    const uint it = get_local_id(0);
    if(i >= N)
        return;
    if((imove[i] != -3) || (iset[i] != BIfree_iset)){
        return;
    }

    const vec_xyz r_i = r[i].XYZ;

    // Initialize the output
    #ifndef LOCAL_MEM_SIZE
        #define _U_ u[i].XYZ
        #define _S_ S[i]
    #else
        #define _U_ u_l[it]
        __local vec_xyz u_l[LOCAL_MEM_SIZE];
        #define _S_ S_l[it]
        __local matrix S_l[LOCAL_MEM_SIZE];
    #endif
    _U_ = VEC_ZERO.XYZ;
    _S_ = MAT_ZERO;

    BEGIN_LOOP_OVER_NEIGHS(){
        if(imove[j] != 2){
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
            const float w_ij = kernelW(q) * CONW * m[j] / rho[j];
            _U_ += u[j].XYZ * w_ij;
            _S_ += S[j] * w_ij;
        }
    }END_LOOP_OVER_NEIGHS()

    #ifdef LOCAL_MEM_SIZE
        u[i].XYZ = _U_;
        S[i] = _S_;
    #endif
}

/** @brief Velocity and stress deviation renormalization.
 *
 * @param iset Set of particles index.
 * @param imove Moving flags.
 *   - imove = 2 for regular solid particles.
 *   - imove = 0 for sensors.
 *   - imove < 0 for boundary elements/particles.
 * @param shepard Shepard term
 * \f$ \gamma(\mathbf{x}) = \int_{\Omega}
 *     W(\mathbf{y} - \mathbf{x}) \mathrm{d}\mathbf{x} \f$.
 * @param u Velocity \f$ \mathbf{u} \f$.
 * @param S Deviatory stress \f$ S \f$.
 * @param N Total number of particles and boundary elements.
 * @param BIfree_iset Set of particles affected
 */
__kernel void shepard(const __global uint* iset,
                      const __global int* imove,
                      const __global float* shepard,
                      __global vec* u,
                      __global matrix* S,
                      uint N,
                      uint BIfree_iset)
{
    uint i = get_global_id(0);
    if(i >= N)
        return;
    if((imove[i] != -3) || (iset[i] != BIfree_iset)){
        return;
    }
    float shepard_i = shepard[i];
    if(shepard_i < 1.0E-6f){
        // It will be considered that there are not enough
        // particles to interpolate
        shepard_i = 1.f;
    }

    u[i] /= shepard_i;
    S[i] /= shepard_i;
}

/** @brief Inverse EOS to get the density from the interpolated presseure.
 *
 * @param iset Set of particles index.
 * @param imove Moving flags.
 *   - imove = 2 for regular solid particles.
 *   - imove = 0 for sensors.
 *   - imove < 0 for boundary elements/particles.
 * @param shepard Shepard term
 * \f$ \gamma(\mathbf{x}) = \int_{\Omega}
 *     W(\mathbf{y} - \mathbf{x}) \mathrm{d}\mathbf{x} \f$.
 * @param p Pressure \f$ p \f$.
 * @param S Deviatory stress \f$ S \f$.
 * @param N Total number of particles and boundary elements.
 * @param cs Speed of sound \f$ c_s \f$.
 * @param p0 Background pressure \f$ p_0 \f$.
 * @param BIfree_iset Set of particles affected
 */
__kernel void eos(const __global uint* iset,
                  const __global int* imove,
                  const __global float* p,
                  __global float* rho,
                  __constant float* refd,
                  uint N,
                  float cs,
                  float p0,
                  uint BIfree_iset)
{
    uint i = get_global_id(0);
    if(i >= N)
        return;
    if((imove[i] != -3) || (iset[i] != BIfree_iset)){
        return;
    }
    rho[i] = refd[iset[i]] + (p[i] - p0) / (cs * cs);
}

/** @brief Simple Euler time integration of the velocity
 *
 * Since the velocity is resulting from the interpolation of the solid particle,
 * it does not make any sense to integrate it with an improved Euler scheme.
 * 
 *   \f$ \mathbf{u}_{n+1} = \mathbf{u}_{n} + \Delta t
 *      \left. \frac{\mathrm{d} \mathbf{u}}{\mathrm{d}t} \right\vert_{n+1/2}
 *   \right) \f$
 *
 * @param iset Set of particles index.
 * @param imove Moving flags.
 *   - imove = 2 for regular solid particles.
 *   - imove = 0 for sensors.
 *   - imove < 0 for boundary elements/particles.
 * @param r Position \f$ \mathbf{r} \f$.
 * @param u Velocity \f$ \mathbf{u} \f$.
 * @param N Number of particles.
 * @param dt Time step \f$ \Delta t \f$.
 * @param BIfree_iset Set of particles affected
 */
__kernel void euler(const __global uint* iset,
                    const __global int* imove,
                    __global vec* r,
                    const __global vec* u,
                    unsigned int N,
                    float dt,
                    uint BIfree_iset)
{
    unsigned int i = get_global_id(0);
    if(i >= N)
        return;
    if((imove[i] != -3) || (iset[i] != BIfree_iset)){
        return;
    }

    r[i] += dt * u[i];
}

/*
 * @}
 */