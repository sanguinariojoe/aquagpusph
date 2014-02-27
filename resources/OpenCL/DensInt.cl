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
	#include "KernelFunctions/Wendland2D.hcl"
#else
	#include "KernelFunctions/CubicSpline3D.hcl"
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

/** Method called outside to interpolate density.
 * @param dens Density of particles (output).
 * @param iMove Moving particle flag.
 * @param pos Position of particles.
 * @param hp Kernel height of particles.
 * @param pmass particle mass
 * @param shepard Shepard term (0th correction).
 * @param lcell Cell where the particle is situated.
 * @param ihoc Head particle of cell chain.
 * @param dPermut Transform each sorted space index into their unsorted space index.
 * @param iPermut Transform each unsorted space index into their sorted space index.
 * @param N Number of particles.
 * @param hfac Kernel height factor
 * @param lvec Number of cells in each direction
 */
__kernel void DensityInterpolation( _g float* dens, _g int* iMove, _g vec* pos,
                                    _g float* pmass, _g float* shepard,
                                    // Link-list data
                                    _g uint *lcell, _g uint *ihoc, _g uint *dPermut, _g uint *iPermut,
                                    // Simulation data
                                    uint N, float hfac, uivec lvec
                                	)
{
	// find position in global arrays
	uint i = get_global_id(0);			// Particle at sorted space
	uint it = get_local_id(0);			// Particle at local memory (temporal storage)
	if(i >= N)
		return;
	#if __BOUNDARY__==0 || __BOUNDARY__==2
		if(iMove[i]<=0)
			return;
	#endif

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
	// Particle data
	uint j,labp, lc;
	float iShepard;
	vec iPos;
	// Neighbours data
	uint cellCount, lcc;
	vec r;
	float r1, pMass;

	//! 1nd.- Particle data (reads it now in order to avoid read it from constant memory several times)
	j = i;							// Backup of the variable, in order to compare later
	labp = dPermut[i];					// Particle index at unsorted space
	lc = lcell[i];						// Cell of the particle
	iShepard = shepard[labp];				// Shepard value
	iPos = pos[i];						// Position of the particle
	//! 2nd.- Initialize output
    #ifndef LOCAL_MEM_SIZE
	    #define _DENS_ dens[labp]
    #else
	    #define _DENS_ dens_l[it]
        _l float dens_l[LOCAL_MEM_SIZE];
    #endif


	_DENS_ = 0.f;
	//! 3th.- Loop over all neightbour particles
	{
		//! 3.a.- Home cell, starting by next particle
		i++;
		while( (i<N) && (lcell[i]==lc) ) {
			#include "DensInt.hcl"
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
				#include "DensInt.hcl"
				i++;
			}			
		}
		//! 3.c.- Home cell, starting from head of chain.
		i = ihoc[lc];
		while(i < j) {                  // Own particle must be included
			#include "DensInt.hcl"
			i++;
		}

	}
	//! 4th.- Append the self particle density interpolation
	if(iMove[j]){
		#ifndef HAVE_3D
			conw = 1.f/(h*h);
		#else
			conw = 1.f/(h*h*h);
		#endif
		wab = kernelW(0.f)*conw*pmass[j];
		_DENS_ += wab;
	}

	//! 5th.- Write output into global memory (at unsorted space)
	dens[labp] = _DENS_/iShepard;

	// ---- A ---- Your code here ---- A ----
	// ---- | ------------------------ | ----

}
