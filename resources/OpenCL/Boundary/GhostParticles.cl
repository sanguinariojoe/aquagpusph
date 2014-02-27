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
	#define _DRDT_F_ drdt_F[labp]
	#define _SHEPARD_ shepard[labp]
#else
	#define _F_ lF[it]
	#define _DRDT_ lDrdt[it]
	#define _DRDT_F_ lDrdt_F[it]
	#define _SHEPARD_ lShepard[it]
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

#include "Wall.hcl"

/** Method called outside to compute forces and density rate due to ghost particles.
 * @param iFluid Fluid identifier.
 * @param pos Position of particles.
 * @param v Velocity of particles.
 * @param f Forces over particles.
 * @param dens Density of particles.
 * @param drdt Density evolution of particles.
 * @param hp Kernel height of particles.
 * @param pmass mass of particles.
 * @param press pressure of particles.
 * @param spsound Speed of sound.
 * @param sigma Viscosity time step term.
 * @param shepard Shepard term (0th correction).
 * @param gradW Shepard term gradient.
 * @param Visckin Kinetic viscosity (one per fluid)
 * @param Viscdyn Dynamic viscosity (one per fluid)
 * @param refd Reference density (one per fluid)
 * @param lcell Cell where the particle is situated.
 * @param ihoc Head particle of cell chain.
 * @param validCell Mark cells that have at least one fluid particle.
 * @param dPermut Transform each sorted space index into their unsorted space index.
 * @param iPermut Transform each unsorted space index into their sorted space index.
 * @param n Number of particles.
 * @param N Number of particles & sensors.
 * @param hfac Kernel height factor.
 * @param lvec Number of cells at each direction.
 * @param grav Gravity force.
 */
__kernel void Boundary(_g int* iFluid, _g int* iMove, _g vec* pos, _g vec* v,
                       _g float* dens, _g float* pmass,
                       _g float* press, _c float* Visckin, _c float* Viscdyn,
                       _c float* refd, _g vec* f, _g float* drdt,
                       _g float* drdt_F, _g float* shepard,
                       // Link-list data
                       _g uint *lcell, _g uint *ihoc, _g short* validCell,
                       _g uint *dPermut, _g uint *iPermut,
                       // Simulation data
                       uint n, uint N, float hfac, uivec lvec, vec grav,
                       // Wall specific data
                       _c struct Wall* wall
                       #ifdef __DELTA_SPH__
                           // Continuity equation diffusive term data
                           , _c float* delta
                           , float dt, float cs
                       #endif
                       #ifdef __NO_LOCAL_MEM__
                           )
                       #else
                           // Local memory to accelerate writes
                           , _l vec* lF, _l float* lDrdt
                           , _l float* lDrdt_F, _l float* lShepard)
                       #endif
{
	// find position in global arrays
	uint i = get_global_id(0);			// Particle at sorted space
	uint it = get_local_id(0);			// Particle at local memory (temporal storage)
	if(i >= N)
		return;

	// ---- | ------------------------ | ----
	// ---- V ---- Your code here ---- V ----

	// Ghost particles must be computed only for moving particles
	int iIMove = iMove[i];
	if(iIMove <= 0){
		return;
	}
	// In order to compute ghost particles wall must be nearer than kernel total height
	vec iPos   = pos[i];
	vec wPos   = wallProjection(iPos, wall);
	vec wDir   = wPos-iPos;
	float dist = sep*h;
	if(dot(wDir, wDir) >= dist*dist){
		return;
	}

	/* All data has been sorted previously, so two spaces must be considereed:
	 * Sorted space, where i is the index of the particle.
	 * Unsorted space, where labp is the index of the particle.
	 *
	 * Sorted space is used usually, in order to read variables coalescing, and unsorted space is used
	 * eventually to write data at the original arrays. Aiming to avoid write several times into the unsorted global
	 * memory address, local memory is employed to the output.
	 */
	// Kernel variables
	float conw, conf, wab, fab;
	// Particle data
	int iIFluid;
	uint j,labp, lc;
	vec iV;
	float iDens, iPress, iVisckin, iViscdyn, rDens;
	// Neighbours data
	uint cellCount, lcc;
	vec pPos, pV,pVn,pVt, r,dv;
	float r1, vdr, pDens, pPress, pMass, prfac, viscg;
	j = i;							// Backup of the variable, in order to compare later
	labp = dPermut[i];					// Particle index at unsorted space
	lc = lcell[i];						// Cell of the particle
	if(!validCell[lc]){					// Don't waste time computing boundaries without fluid.
			shepard[labp] = 0.f;
			return;
	}
	//! 1st.- Initialize output
	#ifndef __NO_LOCAL_MEM__
		_F_       = f[labp];
		_DRDT_    = drdt[labp];
		_DRDT_F_  = drdt_F[labp];
		_SHEPARD_ = shepard[labp];
	#endif

	//! 2nd.- Particle data (reads it now in order to avoid read it from constant memory several times)
	iV       = v[i];                // Velocity of the particle
	iDens    = dens[i];             // Density of the particle
	iPress   = press[i];            // Pressure of the particle
	iIFluid  = iFluid[i];           // Fluid index of the particle
	iVisckin = Visckin[iIFluid];    // Kinematic viscosity of the particle
	iViscdyn = Viscdyn[iIFluid];    // Dynamic viscosity (clamped with artificial viscosity)
	rDens    = refd[iIFluid];       // Density of reference
	#ifdef __DELTA_SPH__
		float iDelta, drfac, rdr;
		iDelta = delta[iIFluid];
	#endif
	//! 3th.- Loop over all neightbour particles
	{
		//! 3.a.- Home cell, starting by own particle
		while( (i<N) && (lcell[i]==lc) ) {
			#include "GhostParticles.hcl"
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
				#include "GhostParticles.hcl"
				i++;
			}			
		}
		//! 3.c.- Home cell, starting from head of chain.
		i = ihoc[lc];
		while(i < j) {
			#include "GhostParticles.hcl"
			i++;
		}

	}
	//! 5th.- Write output into global memory (at unsorted space)
	#ifndef __NO_LOCAL_MEM__
		f[labp]       =  _F_;
		drdt[labp]    =  _DRDT_ + _DRDT_F_;
		drdt_F[labp]  =  _DRDT_F_;
		shepard[labp] =  _SHEPARD_;
	#else
		drdt[labp]   +=  _DRDT_F_;
	#endif

	// ---- A ---- Your code here ---- A ----
	// ---- | ------------------------ | ----

}
