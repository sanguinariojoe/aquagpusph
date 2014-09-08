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
 * @brief Outdated data, just ignore it.
 */

#ifndef HAVE_3D
    #include "../types/2D.h"
#else
    #include "../types/3D.h"
#endif

/** @struct Portal
* @brief Specific portal storage, that get inlet or outlet.
* @note Interior normals.
*/
struct Portal
{
	/// Start corner
	vec corner;
	/// Up vector
	vec up;
	/// Side vector
	vec side;
	/// Normal
	vec normal;
};

/** @struct PortalPair
 * @brief Full portal data structure, that contains a pair of portals as
 * inlet/outlet.
 */
struct PortalPair
{
	struct Portal in,out;
};

/** @brief Teleport particles that pass trought a portal.
 * @param iMove Movement flags. Fixed particles will not considered.
 * @param pos Position of particles.
 * @param N Number of particles.
 * @param in First portal plane.
 * @param out Second portal plane.
 */
__kernel void Portal( __global int* iMove, __global vec* pos, uint N,
                      struct PortalPair portal)
{
	// find position in global arrays
	uint i = get_global_id(0);
	if(i >= N)
		return;
	if(iMove[i]<=0)
		return;

	// ---- | ------------------------ | ----
	// ---- V ---- Your code here ---- V ----

	struct Portal in  = portal.in;
	struct Portal out = portal.out;
	// Get some usefull data
	float inX  = fast_length(in.side);
	float inY  = fast_length(in.up);
	float outX = fast_length(out.side);
	float outY = fast_length(out.up);
	// Test if particle pass throught out portal
	vec relPos = pos[i] - out.corner;
	float relY = dot(relPos, out.up) / (outY*outY);
	float relX = 0.f;
	#ifdef HAVE_3D
		relX = dot(relPos, out.side) / (outX*outX);
	#endif
	float n    = dot(relPos, out.normal);
	if( (n < 0.f) && (relX >= 0.f) && (relX <= 1.f) && (relY >= 0.f) && (relY <= 1.f) ){
		// Teleport the particle to in portal
		pos[i] = in.corner + relX*in.side + relY*in.up + n*in.normal;
		return;
	}

	// ---- A ---- Your code here ---- A ----
	// ---- | ------------------------ | ----

}
