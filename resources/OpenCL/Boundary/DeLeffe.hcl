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

// ------------------------------------------------------------------
// face properties
// ------------------------------------------------------------------
n = normal[i];                                                // Face normal                                    [m]
sDens  = dens[i];                                             // Density field at face                          [kg/m3]
if(sDens <= 0.0001f){ // No fluid present, avoid this vertex
	i++;
	continue;
}
sPress = press[i];                                            // Pressure field at face.                        [Pa]
sArea  = mass[i];                                             // Face area.                                     [m2]
// ------------------------------------------------------------------
// other properties
// ------------------------------------------------------------------
dist = sep*iHp;                                                // Maximum interaction distance                   [m]
// ------------------------------------------------------------------
// Vertex-particle interaction parameters
// ------------------------------------------------------------------
r  = iPos - pos[i];                                            // Vector from vertex to particle                 [m]
r1 = fast_length(r);                                           // Distance between particle and vertex           [m]
if(r1 >= dist){
	i++;
	continue;
}
r0 = dot(r,n);                                                 // Normal distance between particle and vertex     [m]
// Test for swap normal (that must be internal oriented)
if(r0 < 0.f){
	n  = -n;
	r0 = -r0;
}
q = r1/iHp;

// ------------------------------------------------------------------
// Vertex computation
// ------------------------------------------------------------------
#ifndef HAVE_3D
	conw = 1.f/(iHp*iHp);
#else
	conw = 1.f/(iHp*iHp*iHp);
#endif
{
	dv   = iV - v[i];                                       // Delta of velocity                                [m/s]
	vdr  = sDens * dot(dv, n);                              // Projected velocity over normal                   [m/s]
	//---------------------------------------------------------------
	//       calculate the kernel wab
	//---------------------------------------------------------------
	wab = kernelW(q)*conw*sArea;
	//---------------------------------------------------------------
	//       calculate the pressure factor
	//---------------------------------------------------------------
	prfac = sDens * ( iPress/(iDens*iDens) + sPress/(sDens*sDens) ) * n;
	//---------------------------------------------------------------
	//       calculate viscosity terms (Macia etal PTP-2012)
	//---------------------------------------------------------------
	viscg = VEC_ZERO;
	#ifdef __NO_SLIP__
	if(r0 > 0.01f*iHp){
		viscg = 2.f*iViscdyn/iDens * dv * r0/dot(r,r);
	}
	#endif
	//---------------------------------------------------------------
	//       force computation
	//---------------------------------------------------------------
	_F_ += wab*(prfac - viscg);                             // Particle force                                   [m/s2]
	//---------------------------------------------------------------
	//       rate of change of density
	//---------------------------------------------------------------
	_DRDT_ -= wab*vdr;                                      // Density varaition rate                           [kg/m3/s]
	//---------------------------------------------------------------
	//       Shepard term gradient
	//---------------------------------------------------------------
	_GRADW_ += n*wab/sDens;                                 // Kernel gradient                                  [1/m]
}
