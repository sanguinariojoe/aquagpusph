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
        #define __CLEARY__ 15.f
    #endif
#endif

if((imove[j] != -1) &&
   (imove[j] != 1)){
    j++;
    continue;
}

const vec r = pos_i - pos[j];
const float q = fast_length(r) / h;
if(q < support)
{
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
    //       calculate viscosity terms (Cleary's viscosity)
    //---------------------------------------------------------------
    const float vdr = dot(v_i - v[j], r);
    float viscg = 0.f;
    if(imove[j] > 0){
        const float r2 = (q * q + 0.01f) * h * h;
        viscg = -__CLEARY__ * visc_dyn_i * vdr / (r2 * rho_i * rho_j);
    }
    //---------------------------------------------------------------
    //       force computation
    //---------------------------------------------------------------
    _DVDT_ -= r * fab * (prfac + viscg);
    //---------------------------------------------------------------
    //     rate of change of density
    //---------------------------------------------------------------
    _DRDT_ += vdr * fab;
    //---------------------------------------------------------------
    //     Shepard term
    //---------------------------------------------------------------
    _SHEPARD_ += wab / rho_j;
}
