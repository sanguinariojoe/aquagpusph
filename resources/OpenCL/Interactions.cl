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

#ifndef HAVE_3D
    #include "types/2D.h"
    #include "KernelFunctions/Wendland2D.hcl"
#else
    #include "types/3D.h"
    #include "KernelFunctions/Wendland3D.hcl"
#endif

#ifndef HAVE_3D
    #ifndef NEIGH_CELLS
        /** @def NEIGH_CELLS
         * @brief Number of neigh cells.
         *
         * In 2D cases 9 cells must be computed, while in 3D simulations 27
         * cells must be computed.
         */ 
        #define NEIGH_CELLS 9
    #endif
#else
    #ifndef NEIGH_CELLS
        /** @def NEIGH_CELLS
         * @brief Number of neigh cells.
         *
         * In 2D cases 9 cells must be computed, while in 3D simulations 27
         * cells must be computed.
         */ 
        #define NEIGH_CELLS 27
    #endif
#endif

#ifndef uint
    /** @def uint
     * @brief Short alias for unsigned integer type.
     */ 
    #define uint unsigned int
#endif

/** @brief Particles interactions computation.
 *
 * Compute the rates of variation due to the fluid (fixed particles will be
 * included here).
 *
 * During this stage some other operations are performed as well, like the
 * values interpolation in the boundaries (for DeLeffe boundary conditions),
 * the sensors meassurement, or the Shepard factor computation.
 * @param iset Set of particles index.
 * @param imove Moving flags.
 *   - imove > 0 for regular fluid particles.
 *   - imove = 0 for sensors.
 *   - imove < 0 for boundary elements/particles.
 * @param pos Position \f$ \mathbf{r} \f$.
 * @param v Velocity \f$ \mathbf{u} \f$.
 * @param rho Density \f$ \rho \f$.
 * @param m Mass \f$ m \f$.
 * @param p Pressure \f$ p \f$.
 * @param refd Density of reference \f$ \rho_0 \f$ (one per set of particles)
 * @param grad_p Pressure gradient \f$ \frac{\nabla p}{rho} \f$.
 * @param lap_u Velocity laplacian \f$ \frac{\Delta \mathbf{u}}{rho} \f$.
 * @param div_u Velocity divergence \f$ \rho \nabla \cdot \mathbf{u} \f$.
 * @param lap_p Pressure laplacian \f$ \Delta p \f$.
 * @param shepard Shepard term
 * \f$ \gamma(\mathbf{x}) = \int_{\Omega}
 *     W(\mathbf{y} - \mathbf{x}) \mathrm{d}\mathbf{x} \f$.
 * @param icell Cell where each particle is located.
 * @param ihoc Head of chain for each cell (first particle found).
 * @param N Number of particles.
 * @param n_cells Number of cells in each direction
 * @param h Kernel characteristic length \f$ h \f$.
 * @param support Kernel support \f$ s \f$, such that
 * \f$ W \left( s \cdot h \right) = 0 \f$.
 * @param g Gravity acceleration \f$ \mathbf{g} \f$.
 */
__kernel void main(const __global uint* iset,
                   const __global int* imove,
                   const __global vec* pos,
                   const __global vec* v,
                   const __global float* rho,
                   const __global float* m,
                   const __global float* p,
                   __constant float* refd,
                   __global vec* grad_p,
                   __global vec* lap_u,
                   __global float* div_u,
                   __global float* lap_p,
                   __global float* shepard,
                   // Link-list data
                   __global uint *icell,
                   __global uint *ihoc,
                   // Simulation data
                   uint N,
                   uivec4 n_cells,
                   float h,
                   float support,
                   vec g)
{
    const uint i = get_global_id(0);
    const uint it = get_local_id(0);
    if(i >= N)
        return;

    const uint c_i = icell[i];
    const int move_i = imove[i];
    const vec_xyz pos_i = pos[i].XYZ;
    const vec_xyz v_i = v[i].XYZ;
    const float p_i = p[i];
    const float rho_i = rho[i];
    const float refd_i = refd[iset[i]];

    const float prfac_i = p_i / (rho_i * rho_i);

    #ifndef HAVE_3D
        const float conw = 1.f/(h*h);
        const float conf = 1.f/(h*h*h*h);
    #else
        const float conw = 1.f/(h*h*h);
        const float conf = 1.f/(h*h*h*h*h);
    #endif

    // Initialize the output
    #ifndef LOCAL_MEM_SIZE
        #define _GRADP_ grad_p[i].XYZ
        #define _LAPU_ lap_u[i].XYZ
        #define _DIVU_ div_u[i]
        #define _LAPP_ lap_p[i]
        #define _SHEPARD_ shepard[i]
    #else
        #define _GRADP_ grad_p_l[it]
        #define _LAPU_ lap_u_l[it]
        #define _DIVU_ div_u_l[it]
        #define _LAPP_ lap_p_l[it]
        #define _SHEPARD_ shepard_l[it]
        __local vec_xyz grad_p_l[LOCAL_MEM_SIZE];
        __local vec_xyz lap_u_l[LOCAL_MEM_SIZE];
        __local float div_u_l[LOCAL_MEM_SIZE];
        __local float lap_p_l[LOCAL_MEM_SIZE];
        __local float shepard_l[LOCAL_MEM_SIZE];
        _GRADP_ = VEC_ZERO.XYZ;
        _LAPU_ = VEC_ZERO.XYZ;
        _DIVU_ = 0.f;
        _LAPP_ = 0.f;
        _SHEPARD_ = 0.f;
    #endif

    // Loop over neighbour particles
    // =============================
    {
        uint j;
        // Home cell, starting from the next particle
        // ==========================================
        j = i + 1;
        while((j < N) && (icell[j] == c_i) ) {
            if(move_i == 0){
                #include "InteractionsSensors.hcl"
            }
            else if((move_i == 1) || (move_i == -1)){
                #include "Interactions.hcl"
            }
            else{
                #include "InteractionsBounds.hcl"
            }
            j++;
        }

        // Neighbour cells
        // ===============
        for(uint cell = 1; cell < NEIGH_CELLS; cell++) {
            uint c_j;
            switch(cell) {
                case 0: c_j = c_i + 0; break;
                case 1: c_j = c_i + 1; break;
                case 2: c_j = c_i - 1; break;
                case 3: c_j = c_i + n_cells.x; break;
                case 4: c_j = c_i + n_cells.x + 1; break;
                case 5: c_j = c_i + n_cells.x - 1; break;
                case 6: c_j = c_i - n_cells.x; break;
                case 7: c_j = c_i - n_cells.x + 1; break;
                case 8: c_j = c_i - n_cells.x - 1; break;
                #ifdef HAVE_3D
                    case 9 : c_j = c_i + 0             - n_cells.x*n_cells.y; break;
                    case 10: c_j = c_i + 1             - n_cells.x*n_cells.y; break;
                    case 11: c_j = c_i - 1             - n_cells.x*n_cells.y; break;
                    case 12: c_j = c_i + n_cells.x     - n_cells.x*n_cells.y; break;
                    case 13: c_j = c_i + n_cells.x + 1 - n_cells.x*n_cells.y; break;
                    case 14: c_j = c_i + n_cells.x - 1 - n_cells.x*n_cells.y; break;
                    case 15: c_j = c_i - n_cells.x     - n_cells.x*n_cells.y; break;
                    case 16: c_j = c_i - n_cells.x + 1 - n_cells.x*n_cells.y; break;
                    case 17: c_j = c_i - n_cells.x - 1 - n_cells.x*n_cells.y; break;

                    case 18: c_j = c_i + 0             + n_cells.x*n_cells.y; break;
                    case 19: c_j = c_i + 1             + n_cells.x*n_cells.y; break;
                    case 20: c_j = c_i - 1             + n_cells.x*n_cells.y; break;
                    case 21: c_j = c_i + n_cells.x     + n_cells.x*n_cells.y; break;
                    case 22: c_j = c_i + n_cells.x + 1 + n_cells.x*n_cells.y; break;
                    case 23: c_j = c_i + n_cells.x - 1 + n_cells.x*n_cells.y; break;
                    case 24: c_j = c_i - n_cells.x     + n_cells.x*n_cells.y; break;
                    case 25: c_j = c_i - n_cells.x + 1 + n_cells.x*n_cells.y; break;
                    case 26: c_j = c_i - n_cells.x - 1 + n_cells.x*n_cells.y; break;
                #endif
            }

            j = ihoc[c_j];
            while((j < N) && (icell[j] == c_j)) {
                if(move_i == 0){
                    #include "InteractionsSensors.hcl"
                }
                else if((move_i == 1) || (move_i == -1)){
                    #include "Interactions.hcl"
                }
                else{
                    #include "InteractionsBounds.hcl"
                }
                j++;
            }
        }

        // Home cell, starting from the head of chain
        // ==========================================
        j = ihoc[c_i];
        while(j < i) {
            if(move_i == 0){
                #include "InteractionsSensors.hcl"
            }
            else if((move_i == 1) || (move_i == -1)){
                #include "Interactions.hcl"
            }
            else{
                #include "InteractionsBounds.hcl"
            }
            j++;
        }
    }

    // Self particle effect
    // ====================
    if((move_i == 1) || (move_i == -1)){
        const float wab = kernelW(0.f) * conw * m[i];
        _SHEPARD_ += wab / rho_i;
    }

    #ifdef LOCAL_MEM_SIZE
        grad_p[i].XYZ = _GRADP_;
        lap_u[i].XYZ = _LAPU_;
        div_u[i] = _DIVU_;
        lap_p[i] = _LAPP_;
        shepard[i] = _SHEPARD_;
    #endif
}
