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
 * @brief Link-list computation tool.
 * (See Aqua::CalcServer::LinkList for details)
 */

#ifndef HAVE_3D
    #include "types/2D.h"
#else
    #include "types/3D.h"
#endif

/** @brief "Head of Chain" array initialization.
 *
 * The head of chain of the cells will be set as @paramname{N}, which means
 * that no particles have been found.
 *
 * @param ihoc List of first particle in each cell.
 * @param n Number of cells.
 * @param N Number of particles.
 */
__kernel void InitIhoc(__global unsigned int* ihoc,
                       unsigned int n,
                       unsigned int N)
{
	// find position in global arrays
	unsigned int i = get_global_id(0);
	if(i >= n)
		return;

	// ---- | ------------------------ | ----
	// ---- V ---- Your code here ---- V ----

	// Set as an empty cell (send to a particle out of bounds)
	ihoc[i] = N;

	// ---- A ---- Your code here ---- A ----
	// ---- | ------------------------ | ----

}

/** @brief Cell where each particle is located computation.
 * @param icell Cell where the particle is situated.
 * @param pos Position \f$ \mathbf{r} \f$.
 * @param n Number of particles.
 * @param N Next power of two of @paramname{n}.
 * It is the dimension of @paramname{icell} array (requirement for
 * Aqua::CalcServer::RadixSort). All particles next to @paramname{N} must be
 * assigned to a cell greater than @paramname{lxy}.
 * @param posmin Minimum coordinates of the bounds box.
 * @param rdist Cell length
 * @param lxy Total number of cells (l.x \f$ \cdot \f$ l.y)
 * @param l Number of cells in each direction
 */
__kernel void LCell(__global unsigned int* icell,
                    __global vec* pos,
                    unsigned int n,
                    unsigned int N,
                    vec posmin,
                    float rdist,
                    unsigned int lxy,
                    uivec l)
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
			icell[i] = cellID;
		#else
			cell.z   = (unsigned int)((pos[i].z - posmin.z) * rdist) + (unsigned int)3;
			cellID   = cell.x-(unsigned int)1
			        + (cell.y-(unsigned int)1) * l.x
			        + (cell.z-(unsigned int)1) * l.x * l.y;
			icell[i] = cellID;
		#endif
	}
	// Particles out of bounds (extra particles)
	else {
		icell[i] = lxy;
	}

	// ---- A ---- Your code here ---- A ----
	// ---- | ------------------------ | ----

}

/** @brief "Head of Chain" computation.
 *
 * With the particles sorted it is easy to find the "Heads of chains" ones,
 * just looking for the particles such that the previous one is not located in
 * the same cell.
 *
 * @param icell Cell where the particle is situated.
 * @param ihoc List of first particle in each cell.
 * @param N Number of particles (minus one).
 */
__kernel void LinkList(__global unsigned int* icell,
                       __global unsigned int* ihoc,
                       unsigned int N)
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
	lc = icell[i];
	lc2 = icell[i+1];
	if(i==0){
		ihoc[lc] = 0;
	}
	if(lc2 != lc){
		ihoc[lc2] = i+1;
	}

	// ---- A ---- Your code here ---- A ----
	// ---- | ------------------------ | ----

}
