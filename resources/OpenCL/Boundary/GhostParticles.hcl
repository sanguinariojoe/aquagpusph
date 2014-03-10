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

if(imove[j] <= 0){
    j++;
    continue;
}

pos_w = wallProjection(pos[j]);
if(!isOnWallBounds(pos_w)){
    j++;
    continue;
}
const vec pos_j = 2.f * pos_w - pos[j];
const vec r = pos_i - pos_j;
const float q = fast_length(r) / h;
if(q < sep)
{
    //---------------------------------------------------------------
    //       calculate mirrored values of p,v
    //---------------------------------------------------------------
    #if __PRESS_MODEL__ == 0
        // ASM (antisymmetric model)
        const float press_j = -press[j];
    #elif __PRESS_MODEL__ == 1
        // SSM (symmetric model)
        dir = pos_j - pos[j];
        const float press_j = press[j] + dot(dir, grav) * refd_i;
    #elif __PRESS_MODEL__ == 2
        // Takeda
        #error Takeda not implemented yet!
    #else
        #error Unknow pressure extension model
    #endif

    vec v_n = dot(v[j], n) * n;
    vec v_t = v[j] - v_n;
    #if __NORMAL_U_MODEL__ == 0
        // ASM (antisymmetric model)
        v_n = 2.f * dot(wallVelocity(pos_w), n) * n - v_n;
    #elif __NORMAL_U_MODEL__ == 1
        // SSM (symmetric model)
    #elif __NORMAL_U_MODEL__ == 2
        // Takeda
        #error Takeda not implemented yet!
    #elif __NORMAL_U_MODEL__ == 3
        // U0M (No velocity)
        v_n = dot(wallVelocity(pos_w), n) * n;
    #else
        #error Unknow normal velocity extension model
    #endif
    #if __TANGENT_U_MODEL__ == 0
        // ASM (antisymmetric model)
        const vec v_w = wallVelocity(pos_w);
        v_t = 2.f * (v_w - dot(v_w, n) * n) - v_t;
    #elif __TANGENT_U_MODEL__ == 1
        // SSM (symmetric model)
    #elif __TANGENT_U_MODEL__ == 2
        // Takeda
        #error Takeda not implemented yet!
    #elif __TANGENT_U_MODEL__ == 3
        // U0M (No velocity)
        const vec v_w = wallVelocity(pos_w);
        v_t = v_w - dot(v_w, n) * n;
    #else
        #error Unknow tangent velocity extension model
    #endif
    const vec v_j = v_n + v_t;

    //---------------------------------------------------------------
    //       calculate the kernel wab and the function fab
    //---------------------------------------------------------------
    const float dens_j = dens[j];
    const float mass_j = mass[j];
    const float wab = kernelW(q) * conw * mass_j;
    const float fab = kernelF(q) * conf * mass_j;
    //---------------------------------------------------------------
    //       calculate the pressure factor
    //---------------------------------------------------------------
    const float prfac = prfac_i + press_j / (dens_j * dens_j);
    //---------------------------------------------------------------
    //       calculate viscosity terms (Cleary's viscosity)
    //---------------------------------------------------------------
    const float vdr = dot(v_i - v_j, r);
    #ifdef __FREE_SLIP__
        const float viscg = 0.f;
    #else
        const float r2 = (q * q + 0.01f) * h * h;
        const float viscg = -__CLEARY__ * viscdyn_i * vdr / (r2 * dens_i * dens_j);
    #endif
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
        // _DRDT_F_ += delta_i * drfac * fab / (cs * dens_j);
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
