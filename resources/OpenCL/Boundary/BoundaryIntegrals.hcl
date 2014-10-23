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

// Artificial viscosity factor
#ifndef __CLEARY__
    #ifndef HAVE_3D
        #define __CLEARY__ 8.f
    #else
        #define __CLEARY__ 15.f
    #endif
#endif

if(imove[j] != -3){
	j++;
	continue;
}

// ------------------------------------------------------------------
// face properties
// ------------------------------------------------------------------
const vec n_j = normal[j];  // Assumed outwarding oriented
const float rho_j = rho[j];
if(rho_j <= 0.01f * refd_i){
	j++;
	continue;
}
const float area_j = m[j];

const vec r = pos[j] - pos_i;
const float q = fast_length(r) / h;
if(q >= support){
    j++;
    continue;
}

// ------------------------------------------------------------------
// Boundary element computation
// ------------------------------------------------------------------
{
    const float p_j = p[j];
	const vec dv = v[j] - v_i;
	const float vdr = rho_j * dot(dv, n_j);
	//---------------------------------------------------------------
	//       calculate the kernel wab
	//---------------------------------------------------------------
	const float wab = kernelW(q) * conw * area_j;
	//---------------------------------------------------------------
	//       calculate the pressure factor
	//---------------------------------------------------------------
	const vec prfac = rho_j * (prfac_i + p_j / (rho_j * rho_j)) * n_j;
	//---------------------------------------------------------------
    //       calculate viscosity terms
	//---------------------------------------------------------------
    const float r2 = (q * q + 0.01f) * h * h;
    const vec lapufac = __CLEARY__ * vdr / (r2 * rho_i * rho_j) * n_j;
    //---------------------------------------------------------------
    //     Momentum equation (grad(p)/rho and lap(u)/rho)
    //---------------------------------------------------------------
    _GRADP_ += prfac * wab;
    // _LAPU_ += lapufac * wab;
    //---------------------------------------------------------------
    //     Continuity equation (rho*div(u))
    //---------------------------------------------------------------
    _DIVU_ += vdr * wab;
    //---------------------------------------------------------------
    //     Density diffusion term (lap(p))
    //---------------------------------------------------------------
    // const float ndr = - rho_j * dot(r, n_j) / r2;
    // const float drfac = (p_j - p_i) - refd_i * dot(g, r);
    // _LAPP_ -= drfac * ndr * wab;
}
