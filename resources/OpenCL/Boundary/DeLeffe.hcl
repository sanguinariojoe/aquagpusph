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
 * @brief Boundary element - Fluid particle interaction.
 * (See DeLeffe.cl for details)
 *
 * It is prefearable to use a header to be included instead of generating a
 * function for thye particles interaction, which imply more registries
 * consumption.
 */

if(imove[j] >= 0){
	j++;
	continue;
}

// ------------------------------------------------------------------
// face properties
// ------------------------------------------------------------------
vec n_j = normal[j];
const float dens_j = dens[j];
if(dens_j <= 0.0001f){
	j++;
	continue;
}
const float area_j = mass[j];
// ------------------------------------------------------------------
// Vertex-particle interaction parameters
// ------------------------------------------------------------------
const vec r = pos_i - pos[j];
const float q = fast_length(r) / h;
if(q >= sep){
    j++;
    continue;
}
// Test for swap normal (which must be internally oriented)
float r0 = dot(r, n_j);
if(r0 < 0.f){
	n_j  = -n_j;
	r0 = -r0;
}

// ------------------------------------------------------------------
// Boundary element computation
// ------------------------------------------------------------------
{
	const vec dv = v_i - v[j];
	const float vdr = dens_j * dot(dv, n_j);
	//---------------------------------------------------------------
	//       calculate the kernel wab
	//---------------------------------------------------------------
	const float wab = kernelW(q) * conw * area_j;
	//---------------------------------------------------------------
	//       calculate the pressure factor
	//---------------------------------------------------------------
	const vec prfac = dens_j * (prfac_i + press[j] / (dens_j * dens_j)) * n_j;
	//---------------------------------------------------------------
	//       calculate viscosity terms (Macia etal PTP-2012)
	//---------------------------------------------------------------
	#ifdef __NO_SLIP__
        const float r2 = (q * q + 0.01f) * h * h;
	    const vec viscg = 2.f * viscdyn_i * dv * r0 / (dens_i * r2);
    #else
    	const vec viscg = VEC_ZERO;
	#endif
	//---------------------------------------------------------------
	//       force computation
	//---------------------------------------------------------------
	_F_ += wab * (prfac - viscg);
	//---------------------------------------------------------------
	//       rate of change of density
	//---------------------------------------------------------------
	_DRDT_ -= wab * vdr;
}
