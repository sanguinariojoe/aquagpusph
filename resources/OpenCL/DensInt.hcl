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
// Study if the particles can interact
// ------------------------------------------------------------------
#if __BOUNDARY__==0 || __BOUNDARY__==2
	if(iMove[i] < 0){
		i++;
		continue;
	}
#endif
if(!iMove[i]){
	i++;
	continue;
}
hav = 0.5f*(iHp + hp[i]);                      // Average kernel heights for interactions [m]
dist = sep*hav;                                // Maximum interaction distance            [m]
r = iPos - pos[i];
r1  = fast_length(r);                          // Distance between particles              [m]
if( r1 <= dist )
{
	#ifndef HAVE_3D
		conw = 1.f/(hav*hav);                  // Different for 1d and 3d
	#else
		conw = 1.f/(hav*hav*hav);              // Different for 1d and 3d
	#endif
	//---------------------------------------------------------------
	//       calculate the kernel wab and the function fab
	//---------------------------------------------------------------
	pMass = pmass[i];                          // Mass of neighbour particle              [kg]
	wab = kernelW(r1/hav)*conw*pMass;
	//---------------------------------------------------------------
	// 	density computation
	//---------------------------------------------------------------
	_DENS_ += wab;                             //                                         [kg/m3]
}
