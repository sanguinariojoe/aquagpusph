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

// Artificial viscosity factor
#ifndef __CLEARY__
    #ifndef HAVE_3D
        #define __CLEARY__ 8.f
    #else
        #define __CLEARY__ 15.f
    #endif
#endif

if(!imove[j]){
    j++;
    continue;
}
#if __BOUNDARY__==0 || __BOUNDARY__==2
    if(imove[j] < 0){
        j++;
        continue;
    }
#endif

const vec r = pos_i - pos[j];
const float q = fast_length(r) / h;
if(q < sep)
{
    //---------------------------------------------------------------
    //       calculate the kernel wab and the function fab
    //---------------------------------------------------------------
    const float dens_j = dens[j];
    const float mass_j = mass[j];
    const float press_j = press[j];
    const float wab = kernelW(q) * conw * mass_j;
    const float fab = kernelF(q) * conf * mass_j;
    //---------------------------------------------------------------
    //       calculate the pressure factor
    //---------------------------------------------------------------
    const float prfac = prfac_i + press_j / (dens_j * dens_j);
    /*
    if(imove[j] < 0){
        prfac = max(prfac, 0.f);
    }
    */
    //---------------------------------------------------------------
    //       calculate viscosity terms (Cleary's viscosity)
    //---------------------------------------------------------------
    const float vdr = dot(v_i - v[j], r);
    float viscg = 0.f;
    if(imove[j] > 0){
        const float r2 = (q * q + 0.01f) * h * h;
        viscg = -__CLEARY__ * viscdyn_i * vdr / (r2 * dens_i * dens_j);
    }
    //---------------------------------------------------------------
    //       force computation
    //---------------------------------------------------------------
    _F_ -= r * fab * (prfac + viscg);
    //---------------------------------------------------------------
    //     rate of change of density
    //---------------------------------------------------------------
    _DRDT_ += vdr * fab;
    //---------------------------------------------------------------
    //       Density diffusion term
    //---------------------------------------------------------------
    #ifdef __DELTA_SPH__
        const float drfac = (press_i - press_j) - refd_i * dot(grav, r);
        // Ferrari
        // _DRDT_F_ += delta_i * drfac * r1 * fab / (cs * dens_j);
        // Molteni
        // _DRDT_F_ += delta_i * h * drfac * fab / (cs * dens_j);
        // Cercos
        _DRDT_F_ += delta_i * dt * dens_i * drfac * fab / (refd_i * dens_j);
    #endif
    //---------------------------------------------------------------
    //     Shepard term
    //---------------------------------------------------------------
    _SHEPARD_ += wab / dens_j;
}
