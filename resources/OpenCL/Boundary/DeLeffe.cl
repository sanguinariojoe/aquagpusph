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

/* Kernel to use. Take care with support distance (sep), that may vary.
 */
#ifndef HAVE_3D
	#include "../KernelFunctions/Wendland2D.hcl"
#else
	#include "../KernelFunctions/Wendland3D.hcl"
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

#ifdef __NO_LOCAL_MEM__
	#define _F_ f[labp]
	#define _DRDT_ drdt[labp]
	#define _GRADW_ gradW[labp]
#else
	#define _F_ lF[it]
	#define _DRDT_ lDrdt[it]
	#define _GRADW_ lGradW[it]
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

/** Set new values of vertices.
 * @param f Forces over particles.
 * @param iMove Movement flags.
 * @param drdt Density evolution of particles.
 * @param press Pressure at unsorted space.
 * @param pressin Pressure at sorted space.
 * @param dens Density of particles at unsorted space.
 * @param densin Density of particles at sorted space.
 * @param refd Reference density.
 * @param ifluid Fluid identifier.
 * @param gamma Gamma.
 * @param normal Particle normal at unsorted space.
 * @param normalin Particle normal at sorted space.
 * @param shepard Shepard term (0th correction).
 * @param shepardin Shepard term backup.
 * @param iPermut Transform each unsorted space index into their sorted space index.
 * @param N Total number of vertices.
 * @param cs Speed of sound.
 */
__kernel void Vertices( _g int* iMove,
                        _g vec* f, _g float* drdt, _g float* press, _g float* pressin, _g float* dens, _g float* densin,
                        _g float* refd, _g uint* ifluid, _g float* gamma, 
                        _g vec* normal, _g vec* normalin, 
                        _g float* shepard, _g uint *iPermut, uint N, float cs )
{
	uint i = get_global_id(0);
	if(i >= N)
		return;

	/** To solve vertexes pressure and density values are to alternatives (define only
	 * one of them): \n
	 * _ALL_INTERPOLATED_: Purposed by Ferrand, density and pressure are interpolated
	 * by nearest particles. \n
	 * _DENSITY_BASED_: Density is interpolated, but pressure is computed as the
	 * maximum value of interpolated value and state equation application. \n
	 * _PRESSURE_BASED_: Pressure field is interpolated, and density is computed
	 * using state equation. \n
	 * Pressure field is corrected (in all cases) by hydrostatic field. Since density
	 * is not corrected, if you are using weakly compressible (where fluid pressure and
	 * density are related) consideer using _PRESSURE_BASED_ algorithm. At the other
	 * hand _ALL_INTERPOLATED_ is not recommended.
	 */
	#define _PRESSURE_BASED_

	// ---- | ------------------------ | ----
	// ---- V ---- Your code here ---- V ----

	// Sort normals and kernel gradient (for every particle)
	normalin[iPermut[i]] = normal[i];

	// Now compute only vertexes
	if(iMove[i]>=0)
		return;
	float iShepard = shepard[i];
	if(iShepard < 0.001f)		// Alone particle will consider
		iShepard=1.f;
	#ifdef _ALL_INTERPOLATED_
		dens[i]     = drdt[i]/iShepard;
		press[i]    = (f[i].x + f[i].y)/iShepard;
	#elif defined _DENSITY_BASED_
		dens[i]  = drdt[i]/iShepard;
		// Batchelor 1967
		float rdenf = refd[ifluid[i]];
		float gammf = gamma[ifluid[i]];
		float ddenf = dens[i]/rdenf;
		float prb   = cs*cs*rdenf/gammf;
		press[i]    = max(prb*(pow(ddenf,gammf)-1.f), f[i].x/iShepard) + f[i].y/iShepard;
	#elif defined _PRESSURE_BASED_
		press[i]     = (f[i].x + f[i].y)/iShepard;
		press[i]     = max(press[i], 0.f);
		// Reversed Batchelor 1967
		float rdenf  = refd[ifluid[i]];
		float gammf  = gamma[ifluid[i]];
		float iprb   = gammf/(cs*cs*rdenf);
		dens[i]      = rdenf*pow( 1.f + iprb*press[i], 1.f/gammf );
	#else
		#error "Unknow vertexes field computation algorithm"
	#endif
	densin[iPermut[i]]  = dens[i];
	pressin[iPermut[i]] = press[i];
	// Restore zero values
	f[i]     = VEC_ZERO;
	drdt[i]  = 0.f;

	// ---- A ---- Your code here ---- A ----
	// ---- | ------------------------ | ----

}

/** Performs the boundary effect over particles.
 * @param iFluid Fluid identifier.
 * @param iMove Movement flags.
 * @param pos Position of particles.
 * @param normal Particle normal (fluid particles can have null normal).
 * @param v Velocity of particles.
 * @param dens Density of particles.
 * @param hp Kernel height of particles.
 * @param press Pressure.
 * @param mass Mass.
 * @param Viscdyn Dynamic viscosity (one per fluid)
 * @param f Forces over particles.
 * @param drdt Density evolution of particles.
 * @param gradW Shepard term gradient.
 * @param lcell Cell where the particle is situated.
 * @param ihoc Head particle of cell chain.
 * @param dPermut Transform each sorted space index into their unsorted space index.
 * @param iPermut Transform each unsorted space index into their sorted space index.
 * @param N Number of particles.
 * @param hfac Kernel height factor
 * @param lvec Number of cells
 */
__kernel void Boundary(	_g int* iFluid, _g int* iMove,
                        _g vec* pos, _g vec* normal, _g vec* v, _g float* dens, _g float* hp,
                        _g float* press, _g float* mass, _g float* Viscdyn,
                        _g vec* f, _g float* drdt, _g vec* gradW,
                        // Link-list data
                        _g uint *lcell, _g uint *ihoc, _g uint *dPermut, _g uint *iPermut,
                        // Simulation data
                        uint N, float hfac, uivec lvec
                        #ifdef __NO_LOCAL_MEM__
                        	)
                        #else
                        	// Local memory to accelerate writes
                        	, _l vec* lF, _l float* lDrdt, _l vec* lGradW)
                        #endif
{
	// find position in global arrays
	uint i = get_global_id(0);			// Particle at sorted space
	uint it = get_local_id(0);			// Particle at local memory (temporal storage)
	if(i >= N)
		return;
	// Don't compute vertexes or sensors
	if(iMove[i]<=0)
		return;

	// ---- | ------------------------ | ----
	// ---- V ---- Your code here ---- V ----

	/* All data has been sorted previously, so two spaces must be considereed:
	 * Sorted space, where i is the index of the particle.
	 * Unsorted space, where labp is the index of the particle.
	 *
	 * Sorted space is used usually, in order to read variables coalescing, and unsorted space is used
	 * eventually to write data at the original arrays. Aiming to avoid write several times into the unsorted global
	 * memory address, local memory is employed to the output.
	 */

	// Kernel variables
	float dist, conw, wab;
	// Main particle data
	int iIFluid;
	uint j,labp, lc;
	vec iPos, iV;
	float iHp, iDens, iPress, iViscdyn;
	// Neighbour particles data
	uint cellCount, lcc;
	vec r,n,dv, prfac, viscg;
	float r1,r0, vdr, sDens, sPress, sArea, q;
	//! 1nd.- Particle data (reads it now in order to avoid read it from constant memory several times)
	j = i;                              // Backup of the variable, in order to compare later
	labp = dPermut[i];                  // Particle index at unsorted space
	lc = lcell[i];                      // Cell of the particle
	iHp = hp[i];                        // Kernel height of the particle [m]
	iPos = pos[i];                      // Position of the particle
	iV = v[i];                          // Velocity of the particle
	iDens = dens[i];                    // Density of the particle
	iPress = press[i];                  // Pressure of the particle
	iIFluid = iFluid[i];                // Fluid index of the particle
	iViscdyn = Viscdyn[iIFluid];        // Dynamic viscosity (clamped with artificial viscosity)
	//! 2nd.- Initialize output
	#ifndef __NO_LOCAL_MEM__
		_F_       = f[labp];
		_DRDT_    = drdt[labp];
		_GRADW_   = gradW[labp];
	#endif
	//! 3th.- Loop over all neightbour particles
	{
		//! 3.a.- Home cell, starting by next particle
		i++;
		while( (i<N) && (lcell[i]==lc) ) {
			if(iMove[i] < 0){
				#include "DeLeffe.hcl"
			}
			i++;
		}
		//! 3.b.- Neighbour cells
		for(cellCount=1;cellCount<NEIGH_CELLS;cellCount++) {
			// Loop over 8 neighbour cells, taking cell index (lcc)
			switch(cellCount) {
				// Cells at the same Z than main cell
				case 0: lcc = lc + 0; break;
				case 1: lcc = lc + 1; break;
				case 2: lcc = lc - 1; break;
				case 3: lcc = lc + lvec.x; break;
				case 4: lcc = lc + lvec.x + 1; break;
				case 5: lcc = lc + lvec.x - 1; break;
				case 6: lcc = lc - lvec.x; break;
				case 7: lcc = lc - lvec.x + 1; break;
				case 8: lcc = lc - lvec.x - 1; break;
				#ifdef HAVE_3D
					// Cells bellow main cell
					case 9 : lcc = lc + 0          - lvec.x*lvec.y; break;
					case 10: lcc = lc + 1          - lvec.x*lvec.y; break;
					case 11: lcc = lc - 1          - lvec.x*lvec.y; break;
					case 12: lcc = lc + lvec.x     - lvec.x*lvec.y; break;
					case 13: lcc = lc + lvec.x + 1 - lvec.x*lvec.y; break;
					case 14: lcc = lc + lvec.x - 1 - lvec.x*lvec.y; break;
					case 15: lcc = lc - lvec.x     - lvec.x*lvec.y; break;
					case 16: lcc = lc - lvec.x + 1 - lvec.x*lvec.y; break;
					case 17: lcc = lc - lvec.x - 1 - lvec.x*lvec.y; break;
					// Cells over main cell
					case 18: lcc = lc + 0          + lvec.x*lvec.y; break;
					case 19: lcc = lc + 1          + lvec.x*lvec.y; break;
					case 20: lcc = lc - 1          + lvec.x*lvec.y; break;
					case 21: lcc = lc + lvec.x     + lvec.x*lvec.y; break;
					case 22: lcc = lc + lvec.x + 1 + lvec.x*lvec.y; break;
					case 23: lcc = lc + lvec.x - 1 + lvec.x*lvec.y; break;
					case 24: lcc = lc - lvec.x     + lvec.x*lvec.y; break;
					case 25: lcc = lc - lvec.x + 1 + lvec.x*lvec.y; break;
					case 26: lcc = lc - lvec.x - 1 + lvec.x*lvec.y; break;
				#endif
			}
			// Sub-loop over particles into neighbour cells
			i = ihoc[lcc];
			while( (i<N) && (lcell[i]==lcc) ) {
				if(iMove[i] < 0){
					#include "DeLeffe.hcl"
				}
				i++;
			}			
		}
		//! 3.c.- Home cell, starting from head of chain.
		i = ihoc[lc];
		while(i < j) {
			if(iMove[i] < 0){
				#include "DeLeffe.hcl"
			}
			i++;
		}
	}
	//! 4th.- Write output into global memory (at unsorted space)
	#ifndef __NO_LOCAL_MEM__
		f[labp]       = _F_;
		drdt[labp]    = _DRDT_;
		gradW[labp]   = _GRADW_;
	#endif

	// ---- A ---- Your code here ---- A ----
	// ---- | ------------------------ | ----

}

