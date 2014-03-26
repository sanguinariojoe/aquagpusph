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

#ifndef HAVE_3D
    #include "../types/2D.h"
#else
    #include "../types/3D.h"
#endif

#ifndef HAVE_3D
	#ifndef NEIGH_CELLS
		/** @def NEIGH_CELLS Number of neighbour cells. In 2D case 8,
		 * and the main cells must be computed, but in 3D 27 cells,
		 * must be computed.
		 */ 
		#define NEIGH_CELLS 9
	#endif
#else
	#ifndef NEIGH_CELLS
		#define NEIGH_CELLS 27
	#endif
#endif

/* Mathematica
 */
#ifndef M_PI
	#define M_PI 3.14159265359f
#endif
#ifndef iM_PI
	#define iM_PI 0.318309886f
#endif

#ifndef uint
	#define uint unsigned int
#endif

/** Compute the boundary effect, based on an elastic bounce when the particles
 * will tresspass the wall.
 * @param imove Moving flags.
 * @param pos Positions.
 * @param normal Normals.
 * @param v Velocities.
 * @param f Accelerations.
 * @param fin Accelerations (from the previous time step).
 * @param icell Cell where each particle is located.
 * @param ihoc Head particle of chain for each cell.
 * @param N Number of particles.
 * @param dt Time step.
 * @param lvec Number of cells
 * @param grav Gravity acceleration.
 * @param r_element The considered elements radius.
 */
__kernel void Boundary(__global int* imove,
                       __global vec* pos, __global vec* normal,
                       __global vec* v, __global vec* f,
                       // Link-list data
                       __global uint *icell, __global uint *ihoc,
                       // Simulation data
                       uint N, float dt, uivec lvec, vec grav,
                       float r_element )
{
	const uint i = get_global_id(0);
	const uint it = get_local_id(0);
	if(i >= N)
		return;
	if(imove[i] <= 0)
		return;

	// ---- | ------------------------ | ----
	// ---- V ---- Your code here ---- V ----

    const uint c_i = icell[i];
    const vec pos_i = pos[i];
    vec v_i = v[i];
    const vec f_i = f[i];

    // Loop over neighbour particles
    // =============================
    {
        uint j;
        // Home cell, starting from the next particle
        // ==========================================
        j = i + 1;
        while((j < N) && (icell[j] == c_i) ) {
            #include "ElasticBounce.hcl"
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
                #include "ElasticBounce.hcl"
                j++;
            }            
        }
        // Home cell, starting from the head of chain
        // ==========================================
        j = ihoc[c_i];
        while(j < i) {
            #include "ElasticBounce.hcl"
            j++;
        }
    }

	// ---- A ---- Your code here ---- A ----
	// ---- | ------------------------ | ----

}
