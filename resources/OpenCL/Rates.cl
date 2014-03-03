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

#ifdef _g
	#error '_g' is already defined.
#endif
#define _g __global

#ifdef _c
	#error '_c' is already defined.
#endif
#define _c __constant

#ifdef _l
	#error '_l' is already defined.
#endif
#define _l __local

/** Sort all input data by cells (the output data will remain unsorted).
 * @param iFluidin Unsorted fluid identifier.
 * @param iFluid Sorted fluid identifier.
 * @param posin Unsorted positions.
 * @param vin Unsorted velocity.
 * @param hpin Unsorted kernels height.
 * @param densin Unsorted density.
 * @param pressin Unsorted pressure.
 * @param pmassin Unsorted mass.
 * @param pos Sorted positions.
 * @param v Sorted velocity.
 * @param hp Sorted kernels height.
 * @param dens Sorted density.
 * @param press Sorted pressure.
 * @param pmass Sorted mass.
 * @param dPermut Sorted space -> unsorted space permutations.
 * @param iPermut Unsorted space -> sorted space permutations.
 * @param N Number of particles.
 */
__kernel void SortData( _g int* iFluidin, _g int* iFluid,
                        _g int* iMovein, _g int* iMove,
                        _g vec* posin, _g vec* vin,
                        _g float* densin, _g float* pressin, _g float* pmassin,
                        _g vec* pos, _g vec* v,
                        _g float* dens, _g float* press, _g float* pmass,
                        _g uint *dPermut, _g uint *iPermut,
                        uint N)
{
	uint i = get_global_id(0);
	if(i >= N)
		return;

	// ---- | ------------------------ | ----
	// ---- V ---- Your code here ---- V ----

	// We assume i in the unsorted space, and labp in the sorted one 
	uint labp = iPermut[i];

	iFluid[labp] = iFluidin[i];
	iMove[labp]  = iMovein[i];
	v[labp]      = vin[i];
	dens[labp]   = densin[i];
	press[labp]  = pressin[i];
	pmass[labp]  = pmassin[i];
	pos[labp]    = posin[i];

	// ---- A ---- Your code here ---- A ----
	// ---- | ------------------------ | ----
}

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
 * @param lcell Cell where the particle is situated.
 * @param ihoc Head particle of cell chain.
 * @param dPermut Sorted space -> unsorted space permutations.
 * @param iPermut Unsorted space -> sorted space permutations.
 * @param n Number of particles.
 * @param N Number of particles & sensors.
 * @param lvec Number of cells at each direction.
 * @param grav Gravity acceleration.
 */
__kernel void Rates( _g int* iFluid, _g int* iMove,
                     _g vec* pos, _g vec* v, _g float* dens,
                     _g float* pmass, _g float* press, _c float* Visckin,
                     _c float* Viscdyn, _g vec* f, _g float* drdt,
                     _g float* drdt_F, _g float* shepard,
                     // Link-list data
                     _g uint *lcell, _g uint *ihoc,
                     _g uint *dPermut, _g uint *iPermut,
                     // Simulation data
                     uint n, uint N, uivec lvec, vec grav
                     #ifdef __DELTA_SPH__
                         // Continuity equation diffusive term data
                         , _c float* refd, _c float* delta
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

	/* All the data has previously been sorted, so two spaces must be
	 * considereed:
	 * # Sorted space, where i is the index of the particle.
	 * # Unsorted space, where labp is the index of the particle.
	 * In the sorted space are stored the input data in order to read
	 * the variables in a convenient coalescing way, while in the unsorted
	 * space are stored the output data to can easily manage the time
	 * integration and output files writting.
	 */

	// Kernel variables
	float dist, conw, conf, wab, fab;
	// Particle data
	int iIFluid, iIMove;
	uint j,labp, lc;
	vec iPos, iV;
	float iDens, iPress, iVisckin, iViscdyn;
	// Neighbours data
	vec r,dv;
	float r1, vdr, pDens, pMass, prfac, viscg;

	j = i;              // Backup of the variable, in order to compare later
	labp = dPermut[i];  // Particle index at unsorted space
	lc = lcell[i];      // Cell of the particle

    #ifndef LOCAL_MEM_SIZE
	    #define _F_ f[labp]
	    #define _DRDT_ drdt[labp]
    	#define _DRDT_F_ drdt_F[labp]
    	#define _SHEPARD_ shepard[labp]
    #else
	    #define _F_ f_l[it]
	    #define _DRDT_ drdt_l[it]
	    #define _DRDT_F_ drdt_F_l[it]
	    #define _SHEPARD_ shepard_l[it]
        _l vec f_l[LOCAL_MEM_SIZE];
        _l float drdt_l[LOCAL_MEM_SIZE];
        _l float drdt_F_l[LOCAL_MEM_SIZE];
        _l float shepard_l[LOCAL_MEM_SIZE];
    #endif
	// Output initialization
	_F_       = VEC_ZERO;
	_DRDT_    = 0.f;
	_DRDT_F_  = 0.f;
	_SHEPARD_ = 0.f;

	// Particle data
	iPos = pos[i];
	iV = v[i];
	iDens = dens[i];
	iPress = press[i];
	iIFluid = iFluid[i];
	iVisckin = Visckin[iIFluid];
	iViscdyn = Viscdyn[iIFluid];
	iIMove = iMove[i];
	#ifdef __DELTA_SPH__
		float rDens, iDelta, drfac, rdr;
		rDens  = refd[iIFluid];
		iDelta = delta[iIFluid];
	#endif
	//! 3th.- Loop over all neightbour particles
	{
		//! 3.a.- Home cell, starting by next particle
		i++;
		while( (i<N) && (lcell[i]==lc) ) {
			// Sensor specific computation
			if(!iIMove){
				#include "RatesSensors.hcl"
			}
			else{
				#if __BOUNDARY__==0
					// ElasticBounce boundary condition
					if(iIMove<0){
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
					if(iIMove<0){
						#include "RatesBounds.hcl"
					}
					else{
						#include "Rates.hcl"
					}
				#else
					#error Unknow boundary condition
				#endif
			}
			i++;
		}
		//! 3.b.- Neighbour cells
		for(uint cell = 1; cell < NEIGH_CELLS; cell++) {
			// Get the neighbour cell index
            uint c_j;
			switch(cell) {
				// Cells at the same Z than main cell
				case 0: c_j = lc + 0; break;
				case 1: c_j = lc + 1; break;
				case 2: c_j = lc - 1; break;
				case 3: c_j = lc + lvec.x; break;
				case 4: c_j = lc + lvec.x + 1; break;
				case 5: c_j = lc + lvec.x - 1; break;
				case 6: c_j = lc - lvec.x; break;
				case 7: c_j = lc - lvec.x + 1; break;
				case 8: c_j = lc - lvec.x - 1; break;
				#ifdef HAVE_3D
					// Cells bellow main cell
					case 9 : c_j = lc + 0          - lvec.x*lvec.y; break;
					case 10: c_j = lc + 1          - lvec.x*lvec.y; break;
					case 11: c_j = lc - 1          - lvec.x*lvec.y; break;
					case 12: c_j = lc + lvec.x     - lvec.x*lvec.y; break;
					case 13: c_j = lc + lvec.x + 1 - lvec.x*lvec.y; break;
					case 14: c_j = lc + lvec.x - 1 - lvec.x*lvec.y; break;
					case 15: c_j = lc - lvec.x     - lvec.x*lvec.y; break;
					case 16: c_j = lc - lvec.x + 1 - lvec.x*lvec.y; break;
					case 17: c_j = lc - lvec.x - 1 - lvec.x*lvec.y; break;
					// Cells over main cell
					case 18: c_j = lc + 0          + lvec.x*lvec.y; break;
					case 19: c_j = lc + 1          + lvec.x*lvec.y; break;
					case 20: c_j = lc - 1          + lvec.x*lvec.y; break;
					case 21: c_j = lc + lvec.x     + lvec.x*lvec.y; break;
					case 22: c_j = lc + lvec.x + 1 + lvec.x*lvec.y; break;
					case 23: c_j = lc + lvec.x - 1 + lvec.x*lvec.y; break;
					case 24: c_j = lc - lvec.x     + lvec.x*lvec.y; break;
					case 25: c_j = lc - lvec.x + 1 + lvec.x*lvec.y; break;
					case 26: c_j = lc - lvec.x - 1 + lvec.x*lvec.y; break;
				#endif
			}

			i = ihoc[c_j];
			while((i < N) && (lcell[i] == c_j)) {
				if(!iIMove){
					#include "RatesSensors.hcl"
				}
				else{
					#if __BOUNDARY__==0
						// ElasticBounce boundary condition
						if(iIMove<0){
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
						if(iIMove<0){
							#include "RatesBounds.hcl"
						}
						else{
							#include "Rates.hcl"
						}
					#else
						#error Unknow boundary condition
					#endif
		                }
				i++;
			}			
		}
		//! 3.c.- Home cell, starting from head of chain.
		i = ihoc[lc];
		while(i < j) {
			if(!iIMove){
				#include "RatesSensors.hcl"
			}
			else{
				#if __BOUNDARY__==0
					// ElasticBounce boundary condition
					if(iIMove<0){
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
					if(iIMove<0){
						#include "RatesBounds.hcl"
					}
					else{
						#include "Rates.hcl"
					}
				#else
					#error Unknow boundary condition
				#endif
			}
			i++;
		}
	}
	//! 4th.- Append own particle as part of shepard term
	// Sensors not included
	if(iIMove){
		#if __BOUNDARY__==0 || __BOUNDARY__==2
			// Contour not included
			if(iIMove>0){
				#ifndef HAVE_3D
					conw = 1.f/(h*h);				// Different for 1d and 3d
				#else
					conw = 1.f/(h*h*h);			// Different for 1d and 3d
				#endif
				wab = kernelW(0.f)*conw*pmass[j];
				_SHEPARD_ += wab/iDens;
			}
		#else
			#ifndef HAVE_3D
				conw = 1.f/(h*h);				// Different for 1d and 3d
			#else
				conw = 1.f/(h*h*h);			// Different for 1d and 3d
			#endif
			wab = kernelW(0.f)*conw*pmass[j];
			_SHEPARD_ += wab/iDens;
		#endif
	}
	//! 5th.- Write output into global memory (at unsorted space)
	#ifdef LOCAL_MEM_SIZE
		f[labp] = _F_;
		drdt[labp] = _DRDT_ + _DRDT_F_;
		drdt_F[labp] = _DRDT_F_;
		shepard[labp] = _SHEPARD_;
	#else
		drdt[labp] += _DRDT_F_;
	#endif

	// ---- A ---- Your code here ---- A ----
	// ---- | ------------------------ | ----

}
