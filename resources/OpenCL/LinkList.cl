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
    #include "types/2D.h"
#else
    #include "types/3D.h"
#endif

#ifndef M_PI
	#define M_PI 3,14159265359
#endif

#ifdef _g
	#error '_g' is already defined.
#endif
#define _g __global

#ifdef _c
	#error '_c' is already defined.
#endif
#define _c __constant

/** Method called to init ihoc array as all cells empty.
 * @param ihoc I Head Of Cell array.
 * @param n Number of cells
 * @param N Number of particles
 */
__kernel void InitIhoc(_g unsigned int* ihoc, unsigned int n, unsigned int N)
{
	// find position in global arrays
	unsigned int i = get_global_id(0);
	if(i >= n)
		return;

	// ---- | ------------------------ | ----
	// ---- V ---- Your code here ---- V ----

	// Set as empty cell (send to a particle out of bounds)
	ihoc[i] = N;

	// ---- A ---- Your code here ---- A ----
	// ---- | ------------------------ | ----

}

/** Method called outside to locate every particle in their cell.
 * @param lcell Cell where the particle is situated.
 * @param pos Particle position
 * @param n Number of particles
 * @param N Next power of two of n, that marks the dimension of lcell. All
 * particles next to N must be assigned to a cell > maximum cell.
 * @param l Number of cells in each direction
 * @param lxy Total number of cells (lx·ly)
 * @param posmin Minimum position of a particle
 * @param rdist Cell dimension
 */
__kernel void LCell(_g unsigned int* lcell, _g vec* pos, unsigned int n, unsigned int N, vec posmin, float rdist, unsigned int lxy, uivec l)
{
	// find position in global arrays
	unsigned int i = get_global_id(0);
	if(i >= N)
		return;

	// ---- | ------------------------ | ----
	// ---- V ---- Your code here ---- V ----

	uivec cell;
	unsigned int cellID;
	// Normal particles
	if(i < n) {
		cell.x = (unsigned int)((pos[i].x - posmin.x) * rdist) + (unsigned int)3;
		cell.y = (unsigned int)((pos[i].y - posmin.y) * rdist) + (unsigned int)3;
		#ifndef HAVE_3D
			cellID  =  cell.x-(unsigned int)1
			        + (cell.y-(unsigned int)1) * l.x;
			lcell[i] = cellID;
		#else
			cell.z   = (unsigned int)((pos[i].z - posmin.z) * rdist) + (unsigned int)3;
			cellID   = cell.x-(unsigned int)1
			        + (cell.y-(unsigned int)1) * l.x
			        + (cell.z-(unsigned int)1) * l.x * l.y;
			lcell[i] = cellID;
		#endif
	}
	// Particles out of bounds (extra particles)
	else {
		lcell[i] = lxy;
	}

	// ---- A ---- Your code here ---- A ----
	// ---- | ------------------------ | ----

}

/** Method called to set ihoc in the sorted space.
 * @param lcell Cell where the particle is situated.
 * @param ihoc I Head Of Cell array.
 * @param N Number of particles (minus one)
 */
__kernel void LinkList(_g unsigned int* lcell, _g unsigned int* ihoc, unsigned int N)
{
    // find position in global arrays
    unsigned int i = get_global_id(0);
	if(i >= N)
		return;

	// ---- | ------------------------ | ----
	// ---- V ---- Your code here ---- V ----

	// We are looking the first particle on each cell.
	// A particle is the first of the cell if the
	// previous particle is allocated on other cell.
	// As special case, first particle is ever a head
	// of chain.
	unsigned int lc,lc2;
	lc = lcell[i];
	lc2 = lcell[i+1];
	if(i==0){
		ihoc[lc] = 0;
	}
	if(lc2 != lc){
		ihoc[lc2] = i+1;
	}

	// ---- A ---- Your code here ---- A ----
	// ---- | ------------------------ | ----

}
