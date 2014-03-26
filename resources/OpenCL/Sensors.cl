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

#ifndef HAVE_3D
    #include "types/2D.h"
#else
    #include "types/3D.h"
#endif

#ifndef M_PI
    #define M_PI 3,14159265359
#endif
#ifndef iM_PI
    #define iM_PI 0.318309886f
#endif

#ifndef uint
    #define uint unsigned int
#endif

#ifdef _g
    #error '_g' is already defined.
#endif
#define _g __global

#ifdef _c
    #error '_c' is already defined.
#endif
#define _c __constant

#ifdef _l
    #error '_l' is already defined.
#endif
#define _l __local

/** Set the sensors interpolated values.
 * @param f Pressure components, it means the interpolated value and
 * the hydrostatic compensation.
 * @param drdt Interpolated density value.
 * @param drdt_F Density variance (2nd order momentum).
 * @param press Pressure.
 * @param pressin Pressure backup.
 * @param dens Density.
 * @param dens_var Density variance (Just for sensors).
 * @param densin Density backup.
 * @param refd Density of reference.
 * @param ifluid Fluid identifier.
 * @param gamma Gamma EOS exponent.
 * @param shepard Shepard term (0th correction).
 * @param cs Speed of sound.
 * @param i0 First particle which is a sensor.
 * @param N Number of sensors.
 */
__kernel void Sensors(_g vec* f, _g float* drdt, _g float* drdt_F,
                      _g float* press, _g float* pressin, _g float* dens,
                      _g float* densin, _g float* dens_var, _g float* refd,
                      _g uint* ifluid, _c float* gamma, _g float* shepard,
                      float cs, uint i0, uint N )
{
    uint i = get_global_id(0);
    if(i >= N)
        return;
    i += i0;    

    // ---- | ------------------------ | ----
    // ---- V ---- Your code here ---- V ----

    float shepard_i = shepard[i];
    if(shepard_i < 0.01f){
        // It will be considered that there are not enough
        // particles to interpolate
        shepard_i = 1.f;
        f[i] = VEC_ZERO;
        drdt[i] = 0.f;
        drdt_F[i] = 0.f;
    }

    // Meassured density
    // =================
    dens[i] = drdt[i] / shepard_i;
    densin[i] = dens[i];

    // Meassured pressure
    // ==================
    const float refd_i = refd[ifluid[i]];
    const float gamma_i = gamma[ifluid[i]];
    const float ddenf_i = dens[i] / refd_i;
    const float prb_i = cs * cs * refd_i / gamma_i;
    // Based on the density value (Applying the EOS)
    // press[i] = prb_i * (pow(ddenf_i, gamma_i) - 1.f) + f[i].y / shepard_i;
    // All interpolated
    // press[i] = (f[i].x + f[i].y) / shepard_i;
    // Maximum value (EOS vs interpolation)
    press[i] = f[i].y / shepard_i + max(prb_i * (pow(ddenf_i, gamma_i) - 1.f),
                                        f[i].x / shepard_i);
    pressin[i] = press[i];

    // Meassured pressure variance
    // ===========================
    dens_var[i - i0] = drdt_F[i] / shepard_i - dens[i] * dens[i];

    // Restore zero values
    f[i] = VEC_ZERO;
    drdt[i] = 0.f;
    drdt_F[i] = 0.f;

    // ---- A ---- Your code here ---- A ----
    // ---- | ------------------------ | ----
}

