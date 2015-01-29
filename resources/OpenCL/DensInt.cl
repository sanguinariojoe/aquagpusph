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
 * @brief Density field geometrical interpolation.
 * (See Aqua::CalcServer::DensityInterpolation for details)
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
     * @brief Unsigned integer short alias.
	 */ 
	#define uint unsigned int
#endif

/** @brief Interpolate the density as
 * \f$ \rho_i = \frac{1}{\gamma_i} \sum_j W \left(
        \mathbf{r}_j - \mathbf{r}_i
   \right) m_j \f$.
 * @param dens Density \f$ \rho \f$ to be computed.
 * @param imove Moving flags.
 *   - imove > 0 for regular fluid particles.
 *   - imove = 0 for sensors.
 *   - imove < 0 for boundary elements/particles.
 * @param r Position \f$ \mathbf{r} \f$.
 * @param mass Mass \f$ m \f$.
 * @param shepard Shepard term
 * \f$ \gamma(\mathbf{x}) = \int_{\Omega}
 *     W(\mathbf{y} - \mathbf{x}) \mathrm{d}\mathbf{x} \f$.
 * @param icell Cell where each particle is located.
 * @param ihoc Head of chain for each cell (first particle found).
 * @param N Number of particles.
 * @param lvec Number of cells in each direction
 * @see Aqua::CalcServer::DensityInterpolation
 */
__kernel void DensityInterpolation(__global float* dens,
                                   __global int* imove,
                                   __global vec* r,
                                   __global float* mass,
                                   __global float* shepard,
                                   __global uint *icell,
                                   __global uint *ihoc,
                                   uint N,
                                   uivec lvec)
{
	const uint i = get_global_id(0);
	const uint it = get_local_id(0);
	if(i >= N)
		return;
	#if __BOUNDARY__==0 || __BOUNDARY__==2
		if(imove[i] <= 0)
			return;
	#endif

	// ---- | ------------------------ | ----
	// ---- V ---- Your code here ---- V ----

    const uint c_i = icell[i];
    const vec r_i = r[i];
    const float shepard_i = shepard[i];

	#ifndef HAVE_3D
		const float conw = 1.f/(h*h);
	#else
		const float conw = 1.f/(h*h*h);
	#endif

	// Initialize the output
    #ifndef LOCAL_MEM_SIZE
	    #define _DENS_ dens[i]
    #else
	    #define _DENS_ dens_l[it]
        __local float dens_l[LOCAL_MEM_SIZE];
    #endif
	_DENS_ = 0.f;


	// Loop over neighbour particles
    // =============================
	{
        uint j;
        // Home cell, starting from the self particle
        // ==========================================
        j = i;
		while((j < N) && (icell[j] == c_i) ) {
			#include "DensInt.hcl"
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
    			#include "DensInt.hcl"
				j++;
			}			
		}
		// Home cell, starting from the head of chain
        // ==========================================
		j = ihoc[c_i];
		while(j < i) {
			#include "DensInt.hcl"
			j++;
		}
    }

	// Write resulting density output
	dens[i] = _DENS_ / shepard[i];

	// ---- A ---- Your code here ---- A ----
	// ---- | ------------------------ | ----

}
