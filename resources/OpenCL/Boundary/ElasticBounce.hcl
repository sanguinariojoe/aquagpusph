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

// ------------------------------------------------------------------
// face properties
// ------------------------------------------------------------------
p  = pos[i];
n  = normal[dPermut[i]];
pV = v[dPermut[i]];
// ------------------------------------------------------------------
// Vertex relation
// ------------------------------------------------------------------
r  = iPos - p;                                                                   // Vector from vertex to particle               [m]
r0 = dot(r,n);                                                                   // Normal distance between particle and vertex  [m]
// Test for swap normal (that must be internal oriented)
if(r0 < 0.f){
	n  = -n;
	r0 = -r0;                                                                // Positive at normal direction                 [m]
}
// ------------------------------------------------------------------
// Movement data
// ------------------------------------------------------------------
nV   = dot(iV-pV,n);
nF   = dot(iF,n);
nFin = dot(iFin,n);
nG   = dot(grav,n);
dist = dt*nV + 0.5f*dt*dt*(nF + nG);                                  // Distance moved, positive at normal direction [m]
// ------------------------------------------------------------------
// Since normal has been internally oriented, if dist < 0, the
// particle is running against the wall, and then two cases can be
// discriminated:
// * The particle is placed in the effect zone of the wall.
// * The partilce is placed outside the effect zone, but will enter
//   inside it.
// ------------------------------------------------------------------
if( (dist < 0.f) && (r0 + dist <= __MIN_BOUND_DIST__*iHp) ){
	// ------------------------------------------------------------------
	// Reflect particle velocity (using elastic factor)
	// ------------------------------------------------------------------
	// As first approach, particle can be completely fliped, but in order to
	// don't perturbate the moments and forces computation is better choice
	// modify only the velocity.
	v[labp]   = iV - (1.f+__ELASTIC_FACTOR__)*(
	            nV
	          + 0.5f*dt*(nF + nG)
                    )*n;
	// Modify value to next wall tests.
	iV = v[labp];
	// As second approach, an acceleration can be imposed in order to ensure
	// that the particle can't move against the wall.
	/*
	f[labp]   = iF - (1.f+__ELASTIC_FACTOR__)*(
		    nV/(0.5f*dt)
		  + nF
		  + nG
	            )*n;
	// Modify value to next wall tests.
	iF = f[labp];
	*/
	// As thord approach, the acceleration can be set as the gravity (no net
	// effect over the forces), and the velocity corrected in oreder to avoid
	// the particle's wall tresspasing.
	/*
	f[labp]   = iF - (nF + nG)*n;
	v[labp]   = iV - (1.f+__ELASTIC_FACTOR__)*(
	            nV
		  + 0.5f*dt*(nF + nG)
	            )*n
	          + 0.5f*dt*(nF + nG)*n;
	// Modify the values to next wall tests.
	iF = f[labp];
	iV = v[labp];
	*/
}
// ------------------------------------------------------------------
// In some cases user can want to grant that the particles are enough
// far from the boundaries.
// ------------------------------------------------------------------
#ifdef __FORCE_MIN_BOUND_DIST__
if(r0 < __MIN_BOUND_DIST__*iHp){
	// Move ever the same distance in order to avoid several particles
	// placed in the same point near to the corners.
	r             = __MIN_BOUND_DIST__*iHp*n;
	iPos         += r;
	outPos[labp] += r;
}
#endif

