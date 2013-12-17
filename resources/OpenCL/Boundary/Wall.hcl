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

#ifdef HAVE_3D
	typedef struct Wall
	{
		/// 1st Corner
		vec p1;
		/// 2nd Corner
		vec p2;
		/// 3rd Corner
		vec p3;
		/// 4th Corner
		vec p4;
		/// Normal to face
		vec n;
		/// 1st Corner velocity
		vec v1;
		/// 2nd Corner velocity
		vec v2;
		/// 3rd Corner velocity
		vec v3;
		/// 4th Corner velocity
		vec v4;
	} Wall;

	/** Compute point projected over wall
	 * (useful to mirror the particle or
	 * get if particle is in bounds).
	 * @param pos Particle position.
	 * @param wall Wall to reflect.
	 * @return Reflection point.
	 */
	vec wallProjection(vec pos, __constant struct Wall* wall){
		vec p = pos - wall->p1;
		vec n = dot(p, wall->n)*wall->n;
		return pos - n;
	}

	/** Get if point is into wall bounds.
	 * @param pos Particle position.
	 * @param wall Wall to reflect.
	 * @return true if point is on wall bounds,
	 * false otherwise.
	 * @warning Point out of wall plane will
	 * not analized.
	 */
	vec isOnWallBounds(vec pos, __constant struct Wall* wall){
		// Test if point is in 1->2->3 triangle
		float u = dot( (pos      - wall->p2), (wall->p1 - wall->p2) ) /
		          dot( (wall->p1 - wall->p2), (wall->p1 - wall->p2) );
		float v = dot( (pos      - wall->p2), (wall->p3 - wall->p2) ) /
		          dot( (wall->p3 - wall->p2), (wall->p3 - wall->p2) );
		if( (u+v >= 0.f) && (u+v <= 1.f) )
			return true;
		// Before test 1->3->4 triangle we
		// must ensure that it is not a triangular
		// wall.
		float lv =  dot( (wall->p1 - wall->p4), (wall->p1 - wall->p4) );
		if(!lv)
			return false;
		// Test if the point is in the 1->3->4 triangle
		float u = dot( (pos      - wall->p4), (wall->p3 - wall->p4) ) /
		          dot( (wall->p1 - wall->p4), (wall->p3 - wall->p4) );
		float v = dot( (pos      - wall->p4), (wall->p1 - wall->p4) ) /
		          lv;
		if( (u+v >= 0.f) && (u+v <= 1.f) )
			return true;
		// Tes failed
		return false;
	}

	/** Compute wall point velocity
	 * (useful to mirror the particle as ASM).
	 * @param pos Wall point.
	 * @param wall Wall to reflect.
	 * @return Velocity at wall point.
	 * @remarks Point into wall will not tested.
	 * @note Phong approach will selected. To do it
	 * we launch a ray from p1 to pos, getting the
	 * point on opposite segment, interpolating value
	 * into them. Then linear interpolation is
	 * performed on initial ray launched.
	 */
	vec wallVelocity(vec pos, __constant struct Wall* wall){
		// Find opposite segment point
		vec d      = pos - wall->p1;
		float ld1  = length(d);
		vec d     /= ld;
		float ld2  = dot(wall->p3 - pos, d);
		vec p2     = pos + ld2*d;
		// We may try to find the opposite point
		// on segment 3->2, due to triangular
		// face will accept it ever.
		vec td     = wall->p2 - wall->p3;
		vec vv     = wall->v2;
		if(dot( td, d ) < 0.f)
			// We fails, the point is in the
			// segment 4->3.
			td = wall->p4 - wall->p3;
			vv = wall->v4;
		// Compute distances over segment
		float ltd  = length(td);
		float lt   = dot( d, td )*(ld1 + ld2)/ltd;
		// Interpolate value on segment
		float f    = lt /ltd;
		vv         = f*vv + (1.f - f)*wall->v3;
		// Interpolate value on launched ray
		f          = ld1 / (ld1+ld2);
		return f*vv + (1.f - f)*wall->v1;
	}
#else
	typedef struct Wall
	{
		/// 1st Corner
		vec p1;
		/// 2nd Corner
		vec p2;
		/// Normal to face
		vec n;
		/// 1st Corner velocity
		vec v1;
		/// 2nd Corner velocity
		vec v2;
	} Wall;

	/** Compute point projected over wall
	 * (useful to mirror the particle or
	 * get if particle is in bounds).
	 * @param pos Particle position.
	 * @param wall Wall to reflect.
	 * @return Reflection point.
	 */
	vec wallProjection(vec pos, __constant struct Wall* wall){
		vec p = pos - wall->p1;
		vec n = dot(p, wall->n)*wall->n;
		return pos - n;
	}

	/** Get if point is into wall bounds.
	 * @param pos Particle position.
	 * @param wall Wall to reflect.
	 * @return true if point is on wall bounds,
	 * false otherwise.
	 * @warning Point out of wall plane will
	 * not analized.
	 */
	bool isOnWallBounds(vec pos, __constant struct Wall* wall){
		vec p    = pos - wall->p1;
		vec d    = wall->p2 - wall->p1;
		float l2 = dot( p, d );
		if( (l2 < 0.f) || (l2 > dot(d,d)) )
			return false;
		return true;
	}

	/** Compute wall point velocity
	 * (useful to mirror the particle as ASM).
	 * @param pos Wall point.
	 * @param wall Wall to reflect.
	 * @return Velocity at wall point.
	 * @warning Point into wall will not tested.
	 */
	vec wallVelocity(vec pos, __constant struct Wall* wall){
		vec p   = pos - wall->p1;
		vec d   = wall->p2 - wall->p1;
		float f = length(p)/length(d);
		return f*wall->v2 + (1.f - f)*wall->v1;
	}
#endif
