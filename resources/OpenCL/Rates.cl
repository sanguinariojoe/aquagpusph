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
 * (See Aqua::CalcServer::Rates for details)
 */

// To use Gaussian kernel please compile AQUAgpusph with Gauss kernel option
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
 * @param refd Density of reference of the fluid \f$ \rho_0 \f$.
 * @param delta Delta factor \f$ \delta \f$.
 * @param dt Time step \f$ \Delta t \f$.
 * @param cs Speed of sound \f$ c_s \f$.
 * @see Aqua::CalcServer::Rates
 */
__kernel void Rates(__global int* ifluid,
                    __global int* imove,
                    __global vec* pos,
                    __global vec* v,
                    __global float* dens,
                    __global float* mass,
                    __global float* press,
                    __constant float* viscdyn,
                    __global vec* dvdt,
                    __global float* drdt,
                    __global float* drdt_F,
                    __global float* shepard,
                    // Link-list data
                    __global uint *icell,
                    __global uint *ihoc,
                    // Simulation data
                    uint N,
                    uivec lvec,
                    vec grav
                    #ifdef __DELTA_SPH__
                        // Continuity equation diffusive term data
                        , __constant float* refd
                        , __constant float* delta
                        , float dt
                        , float cs
                    #endif
                    )
{
	const uint i = get_global_id(0);
	const uint it = get_local_id(0);
	if(i >= N)
		return;

	// ---- | ------------------------ | ----
	// ---- V ---- Your code here ---- V ----

    const uint c_i = icell[i];
    const int move_i = imove[i];
    const vec pos_i = pos[i];
    const vec v_i = v[i];
    const float press_i = press[i];
    const float dens_i = dens[i];
    const float viscdyn_i = viscdyn[ifluid[i]];
	#ifdef __DELTA_SPH__
        const float delta_i = delta[ifluid[i]];
        const float refd_i = refd[ifluid[i]];
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
    #endif
	_F_       = VEC_ZERO;
	_DRDT_    = 0.f;
	_DRDT_F_  = 0.f;
	_SHEPARD_ = 0.f;

	// Loop over neighbour particles
    // =============================
	{
        uint j;
        // Home cell, starting from the next particle
        // ==========================================
        j = i + 1;
		while((j < N) && (icell[j] == c_i) ) {
			// Sensor specific computation
			if(!move_i){
				#include "RatesSensors.hcl"
			}
			else{
				#if __BOUNDARY__==0
					// ElasticBounce boundary condition
					if(move_i<0){
						#include "RatesBounds.hcl"
					}
					else{
						#include "Rates.hcl"
					}
				#elif __BOUNDARY__==1
					// Fix particles
					#include "Rates.hcl"
				#elif __BOUNDARY__==2
					// DeLeffe boundary condition
					if(move_i<0){
						#include "RatesBounds.hcl"
					}
					else{
						#include "Rates.hcl"
					}
				#else
					#error Unknow boundary condition
				#endif
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
				if(!move_i){
					#include "RatesSensors.hcl"
				}
				else{
					#if __BOUNDARY__==0
						// ElasticBounce boundary condition
						if(move_i < 0){
							#include "RatesBounds.hcl"
						}
						else{
							#include "Rates.hcl"
						}
					#elif __BOUNDARY__==1
						// Fix particles
						#include "Rates.hcl"
					#elif __BOUNDARY__==2
						// DeLeffe boundary condition
						if(move_i < 0){
							#include "RatesBounds.hcl"
						}
						else{
							#include "Rates.hcl"
						}
					#else
						#error Unknow boundary condition
					#endif
		        }
				j++;
			}			
		}
		// Home cell, starting from the head of chain
        // ==========================================
		j = ihoc[c_i];
		while(j < i) {
			if(!move_i){
				#include "RatesSensors.hcl"
			}
			else{
				#if __BOUNDARY__==0
					// ElasticBounce boundary condition
					if(move_i < 0){
						#include "RatesBounds.hcl"
					}
					else{
						#include "Rates.hcl"
					}
				#elif __BOUNDARY__==1
					// Fix particles
					#include "Rates.hcl"
				#elif __BOUNDARY__==2
					// DeLeffe boundary condition
					if(move_i < 0){
						#include "RatesBounds.hcl"
					}
					else{
						#include "Rates.hcl"
					}
				#else
					#error Unknow boundary condition
				#endif
			}
			j++;
		}
    }
	// Self particle effect
	// ====================
	if(move_i){
		#if __BOUNDARY__==0 || __BOUNDARY__==2
			// Contour not included
			if(move_i > 0){
				const float wab = kernelW(0.f) * conw * mass[i];
				_SHEPARD_ += wab / dens_i;
			}
		#else
			const float wab = kernelW(0.f) * conw * mass[i];
			_SHEPARD_ += wab / dens_i;
		#endif
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
