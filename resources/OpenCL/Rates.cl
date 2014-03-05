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

// To use Gaussian kernel please compile AQUAgpusph with Gauss kernel option
#ifndef HAVE_3D
	#include "KernelFunctions/Wendland2D.hcl"
#else
	#include "KernelFunctions/Wendland3D.hcl"
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

#ifndef M_PI
	#define M_PI 3.14159265359f
#endif
#ifndef iM_PI
	#define iM_PI 0.318309886f
#endif

#ifndef uint
	#define uint unsigned int
#endif

/** Compute the rates of variation due to the fluid (fixed particles will be
 * included here). During this stage some other operations are performed as
 * well, like the values interpolation in the boundaries (for DeLeffe boundary
 * conditions), the sensors meassurement, or the Shepard factor computation.
 * @param iFluid Fluid identifier.
 * @param pos Position.
 * @param v Velocity.
 * @param dens Density.
 * @param hp Kernel height.
 * @param pmass Mass
 * @param press Pressure.
 * @param Visckin Kinetic viscosity (one per fluid)
 * @param Viscdyn Dynamic viscosity (one per fluid)
 * @param shepard Shepard term (0th correction).
 * @param f Acceleration.
 * @param drdt Rate of change of the density.
 * @param shepard Shepard factor.
 * @param icell Cell where the particle is situated.
 * @param ihoc Head particle of cell chain.
 * @param N Number of particles & sensors.
 * @param lvec Number of cells at each direction.
 * @param grav Gravity acceleration.
 */
__kernel void Rates(__global int* ifluid, __global int* imove,
                    __global vec* pos, __global vec* v, __global float* dens,
                    __global float* mass, __global float* press,
                    __constant float* Visckin, __constant float* viscdyn,
                    __global vec* f, __global float* drdt,
                    __global float* drdt_F, __global float* shepard,
                    // Link-list data
                    __global uint *icell, __global uint *ihoc,
                    // Simulation data
                    uint N, uivec lvec, vec grav
                    #ifdef __DELTA_SPH__
                        // Continuity equation diffusive term data
                        , __constant float* refd, __constant float* delta
                        , float dt, float cs
                    #endif
                    )
{
	// find position in global arrays
	uint i = get_global_id(0);			// Particle at sorted space
	uint it = get_local_id(0);			// Particle at local memory (temporal storage)
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
	    #define _F_ f[i]
	    #define _DRDT_ drdt[i]
    	#define _DRDT_F_ drdt_F[i]
    	#define _SHEPARD_ shepard[i]
    #else
	    #define _F_ f_l[it]
	    #define _DRDT_ drdt_l[it]
	    #define _DRDT_F_ drdt_F_l[it]
	    #define _SHEPARD_ shepard_l[it]
        __local vec f_l[LOCAL_MEM_SIZE];
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
		f[i] = _F_;
		drdt[i] = _DRDT_ + _DRDT_F_;
		drdt_F[i] = _DRDT_F_;
		shepard[i] = _SHEPARD_;
	#else
		drdt[i] += _DRDT_F_;
	#endif

	// ---- A ---- Your code here ---- A ----
	// ---- | ------------------------ | ----

}
