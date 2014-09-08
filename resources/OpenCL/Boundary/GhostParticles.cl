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
 * @brief Ghost particles fluid extension based boundary condition.
 * (See Aqua::CalcServer::Boundary::GhostParticles for details)
 */

#ifndef HAVE_3D
    #include "../types/2D.h"
    #include "../KernelFunctions/Wendland2D.hcl"
#else
    #include "../types/3D.h"
    #include "../KernelFunctions/Wendland3D.hcl"
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

#include "Wall.hcl"

/** @brief Compute the effect of the boundary based on the virtual ghost
 * particles.
 *
 * In the virtual ghost particles model the particles are mirrored respect to
 * the wall "on the fly", i.e. During the interactions computation.
 * @param ifluid Fluid index.
 * @param imove Moving flags.
 *   - imove > 0 for regular fluid particles.
 *   - imove = 0 for sensors.
 *   - imove < 0 for boundary elements/particles.
 * @param pos Position \f$ \mathbf{r} \f$.
 * @param v Velocity \f$ \mathbf{u} \f$.
 * @param dens Density \f$ \rho \f$.
 * @param mass Mass \f$ m \f$.
 * @param press Pressure \f$ p \f$.
 * @param viscdyn Dynamic viscosity \f$ \nu \f$ (one per fluid)
 * @param refd Density of reference of the fluid \f$ \rho_0 \f$.
 * @param dvdt Velocity rate of change \f$ \frac{d \mathbf{u}}{d t} \f$.
 * @param drdt Density rate of change \f$ \frac{d \rho}{d t} \f$.
 * @param drdt_F Density rate of change restricted to the diffusive term
 * \f$ \left. \frac{d \rho}{d t} \right\vert_F \f$.
 * @param shepard Shepard term
 * \f$ \gamma(\mathbf{x}) = \int_{\Omega}
 *     W(\mathbf{y} - \mathbf{x}) \mathrm{d}\mathbf{x} \f$.
 * @param icell Cell where each particle is located.
 * @param ihoc Head of chain for each cell (first particle found).
 * @param N Number of particles.
 * @param lvec Number of cells in each direction
 * @param grav Gravity acceleration \f$ \mathbf{g} \f$.
 */
__kernel void Boundary(__global int* ifluid, __global int* imove,
                       __global vec* pos, __global vec* v,
                       __global float* dens, __global float* mass,
                       __global float* press, __constant float* viscdyn,
                       __constant float* refd, __global vec* dvdt,
                       __global float* drdt, __global float* drdt_F,
                       __global float* shepard,
                       // Link-list data
                       __global uint *icell, __global uint *ihoc,
                       // Simulation data
                       uint N, uivec lvec, vec grav,
                       // Wall specific data
                       #ifdef HAVE_3D
                           vec p1_w, vec p2_w, vec p3_w, vec p4_w, vec n_w,
                           vec v1_w, vec v2_w, vec v3_w, vec v4_w
                       #else
                           vec p1_w, vec p2_w, vec n_w, vec v1_w, vec v2_w
                       #endif
                       // Continuity equation diffusive term data
                       #ifdef __DELTA_SPH__
                           , __constant float* delta
                           , float dt, float cs
                       #endif
                       )
{
    const uint i = get_global_id(0);
    const uint it = get_local_id(0);
    if(i >= N)
        return;
    if(imove[i] <= 0){
        return;
    }

    // ---- | ------------------------ | ----
    // ---- V ---- Your code here ---- V ----

    // Ghost particles must be computed only for moving particles
    // In order to compute ghost particles wall must be nearer than kernel total height
    const vec pos_i = pos[i];
    vec pos_w = wallProjection(pos_i);
    vec dir = pos_w - pos_i;
    if(dot(dir, dir) / (h * h) >= sep * sep){
        return;
    }

    const uint c_i = icell[i];
    const vec v_i = v[i];
    const float press_i = press[i];
    const float dens_i = dens[i];
    const float refd_i = refd[ifluid[i]];
    #ifndef __FREE_SLIP__
        const float viscdyn_i = viscdyn[ifluid[i]];
    #endif
    #ifdef __DELTA_SPH__
        const float delta_i = delta[ifluid[i]];
    #endif

    const float prfac_i = press_i / (dens_i * dens_i);

    #ifndef HAVE_3D
        const float conw = 1.f/(h*h);
        const float conf = 1.f/(h*h*h*h);
    #else
        const float conw = 1.f/(h*h*h);
        const float conf = 1.f/(h*h*h*h*h);
    #endif

    // Initialize the output
    #ifndef LOCAL_MEM_SIZE
        #define _F_ dvdt[i]
        #define _DRDT_ drdt[i]
        #define _DRDT_F_ drdt_F[i]
        #define _SHEPARD_ shepard[i]
    #else
        #define _F_ dvdt_l[it]
        #define _DRDT_ drdt_l[it]
        #define _DRDT_F_ drdt_F_l[it]
        #define _SHEPARD_ shepard_l[it]
        __local vec dvdt_l[LOCAL_MEM_SIZE];
        __local float drdt_l[LOCAL_MEM_SIZE];
        __local float drdt_F_l[LOCAL_MEM_SIZE];
        __local float shepard_l[LOCAL_MEM_SIZE];
        dvdt_l[it] = dvdt[i];
        drdt_l[it] = drdt[i];
        drdt_F_l[it] = drdt_F[i];
        shepard_l[it] = shepard[i];
    #endif

    // Loop over neighbour particles
    // =============================
    {
        uint j;
        // Home cell, starting from the self particle
        // ==========================================
        j = i;
        while((j < N) && (icell[j] == c_i) ) {
            // Sensor specific computation
            #include "GhostParticles.hcl"
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
                case 3: c_j = c_i + lvec.x; break;
                case 4: c_j = c_i + lvec.x + 1; break;
                case 5: c_j = c_i + lvec.x - 1; break;
                case 6: c_j = c_i - lvec.x; break;
                case 7: c_j = c_i - lvec.x + 1; break;
                case 8: c_j = c_i - lvec.x - 1; break;
                #ifdef HAVE_3D
                    case 9 : c_j = c_i + 0          - lvec.x*lvec.y; break;
                    case 10: c_j = c_i + 1          - lvec.x*lvec.y; break;
                    case 11: c_j = c_i - 1          - lvec.x*lvec.y; break;
                    case 12: c_j = c_i + lvec.x     - lvec.x*lvec.y; break;
                    case 13: c_j = c_i + lvec.x + 1 - lvec.x*lvec.y; break;
                    case 14: c_j = c_i + lvec.x - 1 - lvec.x*lvec.y; break;
                    case 15: c_j = c_i - lvec.x     - lvec.x*lvec.y; break;
                    case 16: c_j = c_i - lvec.x + 1 - lvec.x*lvec.y; break;
                    case 17: c_j = c_i - lvec.x - 1 - lvec.x*lvec.y; break;

                    case 18: c_j = c_i + 0          + lvec.x*lvec.y; break;
                    case 19: c_j = c_i + 1          + lvec.x*lvec.y; break;
                    case 20: c_j = c_i - 1          + lvec.x*lvec.y; break;
                    case 21: c_j = c_i + lvec.x     + lvec.x*lvec.y; break;
                    case 22: c_j = c_i + lvec.x + 1 + lvec.x*lvec.y; break;
                    case 23: c_j = c_i + lvec.x - 1 + lvec.x*lvec.y; break;
                    case 24: c_j = c_i - lvec.x     + lvec.x*lvec.y; break;
                    case 25: c_j = c_i - lvec.x + 1 + lvec.x*lvec.y; break;
                    case 26: c_j = c_i - lvec.x - 1 + lvec.x*lvec.y; break;
                #endif
            }

            j = ihoc[c_j];
            while((j < N) && (icell[j] == c_j)) {
                #include "GhostParticles.hcl"
                j++;
            }            
        }
        // Home cell, starting from the head of chain
        // ==========================================
        j = ihoc[c_i];
        while(j < i) {
            #include "GhostParticles.hcl"
            j++;
        }
    }

    #ifdef LOCAL_MEM_SIZE
        dvdt[i] = _F_;
        drdt[i] = _DRDT_;
        drdt_F[i] = _DRDT_F_;
        shepard[i] = _SHEPARD_;
    #endif

    // ---- A ---- Your code here ---- A ----
    // ---- | ------------------------ | ----

}
