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
 * @brief Particles interactions computation.
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

#if __LAP_FORMULATION__ == __LAP_MONAGHAN__
    #ifndef HAVE_3D
        #define __CLEARY__ 8.f
    #else
        #define __CLEARY__ 10.f
    #endif
#endif

/** @brief MLS based correction term, due to the particles at the other
 * portal side.
 *
 * The term computed with this function should be later renormalized by MLS.
 *
 * @param imove Moving flags.
 *   - imove > 0 for regular fluid particles.
 *   - imove = 0 for sensors.
 *   - imove < 0 for boundary elements/particles.
 * @param imirrored 0 if the particle has not been mirrored, 1 otherwise.
 * @param r Position \f$ \mathbf{r} \f$.
 * @param rmirrored Mirrored position of the particle, \a r if \a imirrored is
 * false (0).
 * @param rho Density \f$ \rho \f$.
 * @param m Mass \f$ m \f$.
 * @param p Pressure \f$ p \f$.
 * @param lap_p_corr Correction term for the Morris Laplacian formula.
 * @param icell Cell where each particle is located.
 * @param ihoc Head of chain for each cell (first particle found).
 * @param N Number of particles.
 * @param n_cells Number of cells in each direction
 */
__kernel void full(const __global int* imove,
                   const __global int* imirrored,
                   const __global vec* r,
                   const __global vec* rmirrored,
                   const __global float* rho,
                   const __global float* m,
                   const __global float* p,
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
    if((!imirrored[i]) || (imove[i] != 1)){
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
        _GRADP_ = lap_p_corr[i].XYZ;
    #endif

    BEGIN_LOOP_OVER_NEIGHS(){
        if(!imirrored[j] || (imove[j] != 1)){
            j++;
            continue;
        }
        const vec_xyz r_ij = rmirrored[j].XYZ - r_i;
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

    #ifdef LOCAL_MEM_SIZE
        lap_p_corr[i] = _GRADP_;
    #endif
}

/** @brief Laplacian of the pressure computation.
 *
 * @param imove Moving flags.
 *   - imove > 0 for regular fluid particles.
 *   - imove = 0 for sensors.
 *   - imove < 0 for boundary elements/particles.
 * @param imirrored 0 if the particle has not been mirrored, 1 otherwise.
 * @param r Position \f$ \mathbf{r} \f$.
 * @param rmirrored Mirrored position of the particle, \a r if \a imirrored is
 * false (0).
 * @param rho Density \f$ \rho \f$.
 * @param m Mass \f$ m \f$.
 * @param p Pressure \f$ p \f$.
 * @param lap_p Pressure laplacian \f$ \Delta p \f$.
 * @param icell Cell where each particle is located.
 * @param ihoc Head of chain for each cell (first particle found).
 * @param N Number of particles.
 * @param n_cells Number of cells in each direction
 */
__kernel void lapp(const __global int* imove,
                   const __global int* imirrored,
                   const __global vec* r,
                   const __global vec* rmirrored,
                   const __global float* rho,
                   const __global float* m,
                   const __global float* p,
                   __global float* lap_p,
                   // Link-list data
                   __global uint *icell,
                   __global uint *ihoc,
                   // Simulation data
                   uint N,
                   uivec4 n_cells)
{
    const uint i = get_global_id(0);
    const uint it = get_local_id(0);
    if(i >= N)
        return;
    if((!imirrored[i]) || (imove[i] != 1))
        return;

    const vec_xyz r_i = r[i].XYZ;
    const float p_i = p[i];

    // Initialize the output
    #ifndef LOCAL_MEM_SIZE
        #define _LAPP_ lap_p[i]
    #else
        #define _LAPP_ lap_p_l[it]
        __local float lap_p_l[LOCAL_MEM_SIZE];
        _LAPP_ = lap_p[i];
    #endif

	BEGIN_LOOP_OVER_NEIGHS(){
		if(!imirrored[j] || (imove[j] != 1)){
			j++;
			continue;
		}
		const vec_xyz r_ij = rmirrored[j].XYZ - r_i;
		const float q = length(r_ij) / H;
		if(q >= SUPPORT)
		{
			j++;
			continue;
		}

		{
            const float f_ij = kernelF(q) * CONF * m[j] / rho[j];
            _LAPP_ += (p[j] - p_i) * f_ij;
		}
	}END_LOOP_OVER_NEIGHS()

    #ifdef LOCAL_MEM_SIZE
        lap_p[i] = _LAPP_;
    #endif
}

/** @brief Compute the mirrored position of the fluid particles.
 *
 * The mirrored particles (the ones close enough to the symmetry plane) will be
 * marked with \a imirrored = 1.
 * 
 * @param imirrored 0 if the particle has not been mirrored, 1 otherwise.
 * @param lap_p_corr Correction term for the Morris Laplacian formula.
 * @param lap_p_corr_mirrored Mirrored correction term for the Morris
 * Laplacian formula.
 * @param N Number of particles.
 * @param symmetry_n Normal of the symmetry plane. It is assumed as normalized.
 */
__kernel void mirror(const __global int* imirrored,
                     const __global vec* lap_p_corr,
                     __global vec* lap_p_corr_mirrored,
                     unsigned int N,
                     vec symmetry_n)
{
    // find position in global arrays
    unsigned int i = get_global_id(0);
    if(i >= N)
        return;

    if(!imirrored[i]){
        lap_p_corr_mirrored[i] = lap_p_corr[i];
        return;
    }

    const float v_n = dot(-lap_p_corr[i], symmetry_n);
    lap_p_corr_mirrored[i] = lap_p_corr[i] + 2.f * v_n * symmetry_n;
}

/** @brief Laplacian of the pressure correction.
 *
 * @param imove Moving flags.
 *   - imove > 0 for regular fluid particles.
 *   - imove = 0 for sensors.
 *   - imove < 0 for boundary elements/particles.
 * @param imirrored 0 if the particle has not been mirrored, 1 otherwise.
 * @param r Position \f$ \mathbf{r} \f$.
 * @param rmirrored Mirrored position of the particle, \a r if \a imirrored is
 * false (0).
 * @param rho Density \f$ \rho \f$.
 * @param m Mass \f$ m \f$.
 * @param p Pressure \f$ p \f$.
 * @param lap_p_corr Correction term for the Morris Laplacian formula.
 * @param lap_p_corr_mirrored Mirrored correction term for the Morris
 * Laplacian formula.
 * @param lap_p Pressure laplacian \f$ \Delta p \f$.
 * @param icell Cell where each particle is located.
 * @param ihoc Head of chain for each cell (first particle found).
 * @param N Number of particles.
 * @param n_cells Number of cells in each direction
 */
__kernel void lapp_corr(const __global int* imove,
                        const __global int* imirrored,
                        const __global vec* r,
                        const __global vec* rmirrored,
                        const __global float* rho,
                        const __global float* m,
                        const __global vec* lap_p_corr,
                        const __global vec* lap_p_corr_mirrored,
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
    if((!imirrored[i]) || (imove[i] != 1))
        return;

    const vec_xyz r_i = r[i].XYZ;
    const vec_xyz gradp_i = lap_p_corr[i].XYZ;

    // Initialize the output
    #ifndef LOCAL_MEM_SIZE
        #define _LAPP_ lap_p[i]
    #else
        #define _LAPP_ lap_p_l[it]
        __local float lap_p_l[LOCAL_MEM_SIZE];
        _LAPP_ = lap_p[i];
    #endif

    BEGIN_LOOP_OVER_NEIGHS(){
		if(!imirrored[j] || (imove[j] != 1)){
            j++;
            continue;
        }
        const vec_xyz r_ij = rmirrored[j].XYZ - r_i;
        const float q = length(r_ij) / H;
        if(q >= SUPPORT)
        {
            j++;
            continue;
        }
        {
            const vec_xyz gradp_ij = lap_p_corr_mirrored[j].XYZ + gradp_i;
            const float f_ij = kernelF(q) * CONF * m[j] / rho[j];
            _LAPP_ -= 0.5f * dot(gradp_ij, r_ij) * f_ij;
        }
    }END_LOOP_OVER_NEIGHS()

    #ifdef LOCAL_MEM_SIZE
        lap_p[i] = _LAPP_;
    #endif
}
