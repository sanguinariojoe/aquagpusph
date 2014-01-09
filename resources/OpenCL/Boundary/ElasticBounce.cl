/*
 *  This file is part of AQUA-gpusph, a free CFD program based on SPH.
 *  Copyright (C) 2012  Jose Luis Cercos Pita <jl.cercos@upm.es>
 *
 *  AQUA-gpusph is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  (at your option) any later version.
 *
 *  AQUA-gpusph is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with AQUA-gpusph.  If not, see <http://www.gnu.org/licenses/>.
 */

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

/* OpenCL related
 */
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

/** Performs the boundary effect over particles.
 * @param iMove Movement flags.
 * @param pos Position of particles.
 * @param v Velocity of particles.
 * @param f Forces over particles.
 * @param fin Previous time step forces over particles.
 * @param hp Kernel height of the particles.
 * @param outPos Position of the particles (at unsorted space, so can be used as output variable).
 * @param lcell Cell where the particle is situated.
 * @param ihoc Head particle of cell chain.
 * @param dPermut Transform each sorted space index into their unsorted space index.
 * @param iPermut Transform each unsorted space index into their sorted space index.
 * @param N Number of particles.
 * @param lvec Number of cells
 * @param grav Gravity vector
 * @param r_element Tangential distance to the element in order to considerate that the particle is passing through him
 */
__kernel void Boundary( _g int* iMove,
                        _g vec* pos, _g vec* v, _g vec* f, _g vec* fin, _g vec* normal,
                        _g float* hp, _g vec* outPos,
                        // Link-list data
                        _g uint *lcell, _g uint *ihoc, _g uint *dPermut, _g uint *iPermut,
                        // Simulation data
                        uint N, float dt, uivec lvec, vec grav, float r_element )
{
	// find position in global arrays
	uint i = get_global_id(0);			// Particle at sorted space
	uint it = get_local_id(0);			// Particle at local memory (temporal storage)
	if(i >= N)
		return;
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

	// Particle data
	uint j,labp, lc;
	vec iPos, iV, iF, iFin;
	float iHp, nV, nF, nFin, nG, dist;
	// Neighbours data
	uint cellCount, lcc;
	vec n,p,pV,r,rt;
	float r0;
	//! 1nd.- Particle data (reads it now in order to avoid read it from constant memory several times)
	j = i;							// Backup of the variable, in order to compare later
	labp = dPermut[i];					// Particle index at unsorted space
	lc   = lcell[i];					// Cell of the particle
	iPos = pos[i];						// Position of the particle
	iV   = v[labp];						// Velocity of the particle
	iF   = f[labp];						// Acceleration of the particle
	iFin = fin[labp];					// Acceleration of the particle (previous time step)
	iHp  = hp[i];						// Kernel height of the particle
	//! Loop over all neightbour particles
	{
		//! a.- Home cell, starting by next particle
		i++;
		while( (i<N) && (lcell[i]==lc) ) {
			if(iMove[i] < 0){
				#include "ElasticBounce.hcl"
			}
			i++;
		}
		//! b.- Neighbour cells
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
					#include "ElasticBounce.hcl"
				}
				i++;
			}
		}
		//! c.- Home cell, starting from head of chain.
		i = ihoc[lc];
		while(i < j) {
			if(iMove[i] < 0){
				#include "ElasticBounce.hcl"
			}
			i++;
		}
	}

	// ---- A ---- Your code here ---- A ----
	// ---- | ------------------------ | ----

}
