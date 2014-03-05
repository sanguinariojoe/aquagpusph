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
// Study if two particles can interact
// ------------------------------------------------------------------
if(!imove[j]){
	j++;
	continue;
}
#if __BOUNDARY__==0 || __BOUNDARY__==2
	// ElasticBounce or DeLeffe boundary condition
	if(imove[j]<0){
		j++;
		continue;
	}
#endif

const vec r = pos_i - pos[j];
const float q = fast_length(r) / h;
if(q <= sep)
{
	//---------------------------------------------------------------
	//       calculate the kernel wab and the function fab
	//---------------------------------------------------------------
	const float dens_j = dens[j];
	const float mass_j = mass[j];
    const float press_j = press[j];
	const float wab = kernelW(q) * conw * mass_j;
	//---------------------------------------------------------------
	//       pressure computation (stored on f)
	//---------------------------------------------------------------
	_F_.x += press_j / dens_j * wab;
	_F_.y += dot(grav, r) * wab;
	//---------------------------------------------------------------
	// 	density computation (stored on drdt)
	//---------------------------------------------------------------
	_DRDT_ += wab;
	//---------------------------------------------------------------
	// 	Shepard term
	//---------------------------------------------------------------
	_SHEPARD_ += wab / dens_j;
}
