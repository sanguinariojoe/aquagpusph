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
    #include "../../../../types/2D.h"
    #include "../../../../KernelFunctions/Wendland2D.hcl"
#else
    #include "../../../../types/3D.h"
    #include "../../../../KernelFunctions/Wendland3D.hcl"
#endif

#if __LAP_FORMULATION__ == __LAP_MONAGHAN__
    #ifndef HAVE_3D
        #define __CLEARY__ 8.f
    #else
        #define __CLEARY__ 10.f
    #endif
#endif

/** @brief Particles interactions computation.
 *
 * Compute the rates of variation due to the fluid (fixed particles will be
 * included here).
 *
 * During this stage some other operations are performed as well, like the
 * values interpolation in the boundaries (for DeLeffe boundary conditions),
 * the sensors meassurement, or the Shepard factor computation.
 * @param imove Moving flags.
 *   - imove > 0 for regular fluid particles.
 *   - imove = 0 for sensors.
 *   - imove < 0 for boundary elements/particles.
 * @param imirrored 0 if the particle has not been mirrored, 1 otherwise.
 * @param r Position \f$ \mathbf{r} \f$.
 * @param rmirrored Mirrored position of the particle, \a r if \a imirrored is
 * false (0).
 * @param u Velocity \f$ \mathbf{u} \f$.
 * @param umirrored Mirrored velocity of the particle, \a u if \a imirrored is
 * false (0).
 * @param rho Density \f$ \rho \f$.
 * @param m Mass \f$ m \f$.
 * @param p Pressure \f$ p \f$.
 * @param h_var variable kernel lenght \f$ h \f$.
 * @param Omega \f$ \Omega \f$ term.
 * @param grad_p Pressure gradient \f$ \frac{\nabla p}{rho} \f$.
 * @param lap_u Velocity laplacian \f$ \frac{\Delta \mathbf{u}}{rho} \f$.
 * @param div_u Velocity divergence \f$ \rho \nabla \cdot \mathbf{u} \f$.
 * @param shepard Shepard term
 * \f$ \gamma(\mathbf{x}) = \int_{\Omega}
 *     W(\mathbf{y} - \mathbf{x}) \mathrm{d}\mathbf{x} \f$.
 * @param icell Cell where each particle is located.
 * @param ihoc Head of chain for each cell (first particle found).
 * @param N Number of particles.
 * @param n_cells Number of cells in each direction
 */
__kernel void entry(const __global int* imove,
                    const __global int* imirrored,
                    const __global vec* r,
                    const __global vec* rmirrored,
                    const __global vec* u,
                    const __global vec* umirrored,
                     const __global float* rho,
                    const __global float* m,
                    const __global float* p,
                    const __global float* h_var,
                    const __global float* Omega,
                    __global vec* grad_p,
                    __global vec* lap_u,
                    __global float* div_u,
                    __global float* shepard,
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
    const vec_xyz u_i = u[i].XYZ;
    const float rho_i = rho[i];
    const float m_i = m[i];
    const float V_i = m_i / rho_i;
    const float h_i = h_var[i];
    #ifndef HAVE_3D
        const float conw_i = 1.f / (h_i * h_i);
        const float conf_i = 1.f / (h_i * h_i * h_i * h_i);
    #else
        const float conw_i = 1.f / (h_i * h_i * h_i);
        const float conf_i = 1.f / (h_i * h_i * h_i * h_i * h_i);
    #endif
    const float Omega_i = Omega[i];
    const float p_i = p[i] * m_i * m_i / (rho_i * rho_i) / Omega_i;

    // Initialize the output
    #ifndef LOCAL_MEM_SIZE
        #define _GRADP_ grad_p[i].XYZ
        #define _LAPU_ lap_u[i].XYZ
        #define _DIVU_ div_u[i]
        #define _SHEPARD_ shepard[i]
    #else
        #define _GRADP_ grad_p_l[it]
        #define _LAPU_ lap_u_l[it]
        #define _DIVU_ div_u_l[it]
        #define _SHEPARD_ shepard_l[it]
        __local vec_xyz grad_p_l[LOCAL_MEM_SIZE];
        __local vec_xyz lap_u_l[LOCAL_MEM_SIZE];
        __local float div_u_l[LOCAL_MEM_SIZE];
        __local float shepard_l[LOCAL_MEM_SIZE];
        _GRADP_ = grad_p[i].XYZ;
        _LAPU_ = lap_u[i].XYZ;
        _DIVU_ = div_u[i];
        _SHEPARD_ = shepard[i];
    #endif

	BEGIN_LOOP_OVER_NEIGHS(){
		if((!imirrored[j]) || ((imove[j] != 1) && (imove[j] != -1))){
			j++;
			continue;
		}
        const vec_xyz r_ij = rmirrored[j].XYZ - r_i;
        const float h_j = h_var[j];
        const float l_ij = length(r_ij);
        const float q_i = min(l_ij / h_i, SUPPORT);
        const float q_j = min(l_ij / h_j, SUPPORT);
        if((q_i == SUPPORT) && (q_j == SUPPORT))
        {
            j++;
            continue;
        }
		{
            const float rho_j = rho[j];
            const float m_j = m[j];
            const float udr = dot(umirrored[j].XYZ - u_i, r_ij);
            #ifndef HAVE_3D
                const float conf_j = 1.f / (h_j * h_j * h_j * h_j);
            #else
                const float conf_j = 1.f / (h_j * h_j * h_j * h_j * h_j);
            #endif
            const float p_j = p[j] * m_j * m_j / (rho_j * rho_j) / Omega[j];

			const float wi_ij = conw_i * kernelW(q_i);
            const float fi_ij = conf_i * kernelF(q_i);
            const float fj_ij = conf_j * kernelF(q_j);

            _GRADP_ += (p_i * fi_ij + p_j * fj_ij) / m_i * r_ij;

            #if __LAP_FORMULATION__ == __LAP_MONAGHAN__
                const float r2 = (q * q + 0.01f) * H * H;
                _LAPU_ += 0.5f * (fi_ij + fj_ij) * m_j * __CLEARY__ *
                          udr / (r2 * rho_i * rho_j) * r_ij;
            #elif __LAP_FORMULATION__ == __LAP_MORRIS__
                _LAPU_ += (fi_ij + fj_ij) * m_j / (rho_i * rho_j) *
                          (umirrored[j].XYZ - u_i);
            #else
                #error Unknown Laplacian formulation: __LAP_FORMULATION__
            #endif

            _DIVU_ += udr * fi_ij * m_j * rho_i / (rho_j * Omega_i);
			_SHEPARD_ += wi_ij * m_j / rho_j;
		}
	}END_LOOP_OVER_NEIGHS()

    #ifdef LOCAL_MEM_SIZE
        grad_p[i].XYZ = _GRADP_;
        lap_u[i].XYZ = _LAPU_;
        div_u[i] = _DIVU_;
        shepard[i] = _SHEPARD_;
    #endif
}
