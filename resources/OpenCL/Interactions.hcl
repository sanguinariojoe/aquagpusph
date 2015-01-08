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
 * @brief Fluid particles interaction computation.
 * (See Rates.cl)
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
        #define __CLEARY__ 10.f
    #endif
#endif

//---------------------------------------------------------------
//       calculate the kernel wab and the function fab
//---------------------------------------------------------------
const float rho_j = rho[j];
const float m_j = m[j];
const float p_j = p[j];
const float wab = kernelW(q) * conw * m_j;
const float fab = kernelF(q) * conf * m_j;
//---------------------------------------------------------------
//       calculate the pressure factor
//---------------------------------------------------------------
const float prfac = prfac_i + p_j / (rho_j * rho_j);
//---------------------------------------------------------------
//       calculate viscosity terms
//---------------------------------------------------------------
const float vdr = dot(v[j].XYZ - v_i, r);
float lapufac = 0.f;
if(move_j > 0){
    const float r2 = (q * q + 0.01f) * h * h;
    lapufac = __CLEARY__ * vdr / (r2 * rho_i * rho_j);
}
//---------------------------------------------------------------
//     Momentum equation (grad(p)/rho and lap(u)/rho)
//---------------------------------------------------------------
_GRADP_ += r * fab * prfac;
_LAPU_ += r * fab * lapufac;
//---------------------------------------------------------------
//     Conserving mass equation (rho*div(u))
//---------------------------------------------------------------
_DIVU_ += vdr * fab;
//---------------------------------------------------------------
//     Density diffusion term (lap(p))
//---------------------------------------------------------------
const float drfac = (p_j - p_i) - refd_i * dot(g.XYZ, r);
_LAPP_ += drfac * fab / rho_j;
//---------------------------------------------------------------
//     Shepard term
//---------------------------------------------------------------
_SHEPARD_ += wab / rho_j;
