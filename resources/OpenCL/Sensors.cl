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
 * @brief Sensors fields computation.
 * (See Aqua::CalcServer::Sensors for details)
 */

#ifndef HAVE_3D
    #include "types/2D.h"
#else
    #include "types/3D.h"
#endif

/** @brief Set the sensors interpolated values.
 * @param dvdt Pressure components, it means the interpolated value and
 * the hydrostatic compensation.
 * @param drdt Interpolated density value.
 * @param drdt_F Density variance (2nd order momentum).
 * @param press Pressure \f$ p \f$.
 * @param pressin Pressure backup \f$ p \f$.
 * @param dens Density \f$ \rho \f$.
 * @param dens_var Density variance (Just for sensors) \f$ \sigma_\rho \f$.
 * @param densin Density backup \f$ \rho \f$.
 * @param refd Density of reference \f$ \rho_0 \f$.
 * @param ifluid Fluid identifier.
 * @param gamma Gamma EOS exponent \f$ \gamma \f$.
 * @param shepard Shepard term (0th correction) \f$ \gamma_a \f$.
 * @param cs Speed of sound \f$ c_s \f$.
 * @param i0 First particle which is a sensor.
 * @param N Number of sensors.
 * @see Aqua::CalcServer::Sensors
 * @see Aqua::CalcServer::Rates
 */
__kernel void Sensors(__global vec* dvdt,
                      __global float* drdt,
                      __global float* drdt_F,
                      __global float* press,
                      __global float* pressin,
                      __global float* dens,
                      __global float* densin,
                      __global float* dens_var,
                      __global float* refd,
                      __global unsigned int* ifluid,
                      __constant float* gamma,
                      __global float* shepard,
                      float cs,
                      unsigned int i0,
                      unsigned int N )
{
    unsigned int i = get_global_id(0);
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
        dvdt[i] = VEC_ZERO;
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
    // press[i] = prb_i * (pow(ddenf_i, gamma_i) - 1.f) + dvdt[i].y / shepard_i;
    // All interpolated
    // press[i] = (dvdt[i].x + dvdt[i].y) / shepard_i;
    // Maximum value (EOS vs interpolation)
    press[i] = dvdt[i].y / shepard_i + max(prb_i * (pow(ddenf_i, gamma_i) - 1.f),
                                        dvdt[i].x / shepard_i);
    pressin[i] = press[i];

    // Meassured pressure variance
    // ===========================
    dens_var[i - i0] = drdt_F[i] / shepard_i - dens[i] * dens[i];

    // Restore zero values
    dvdt[i] = VEC_ZERO;
    drdt[i] = 0.f;
    drdt_F[i] = 0.f;

    // ---- A ---- Your code here ---- A ----
    // ---- | ------------------------ | ----
}

