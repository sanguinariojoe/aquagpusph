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
 * @brief Boundary integral term computation.
 * (See Aqua::CalcServer::Boundary::DeLeffe for details)
 */

#ifndef HAVE_3D
    #include "../../../types/2D.h"
#else
    #include "../../../types/3D.h"
#endif

/** @brief Set the fields for the boundary elements using the pressure
 * interpolated by Boundary/BI/Interpolation.cl.
 *
 * Therefore the pressure value will be renormalizaed, and the density computed
 * using the inverse of the EOS.
 *
 * @param iset Set of particles index.
 * @param imove Moving flags.
 *   - imove > 0 for regular fluid particles.
 *   - imove = 0 for sensors.
 *   - imove < 0 for boundary elements/particles.
 * @param shepard Shepard term
 * \f$ \gamma(\mathbf{x}) = \int_{\Omega}
 *     W(\mathbf{y} - \mathbf{x}) \mathrm{d}\mathbf{x} \f$.
 * @param rho Density \f$ \rho \f$.
 * @param p Pressure \f$ p \f$.
 * @param refd Density of reference of the fluid \f$ \rho_0 \f$.
 * @param N Total number of particles and boundary elements.
 * @param cs Speed of sound \f$ c_s \f$.
 * @param p0 Background pressure \f$ p_0 \f$.
 * @see Boundary/BI/Interpolation.cl
 */
__kernel void entry(const __global uint* iset,
                    const __global int* imove,
                    const __global float* shepard,
                    __global float* rho,
                    __global float* p,
                    __constant float* refd,
                    uint N,
                    float cs,
                    float p0)
{
    uint i = get_global_id(0);
    if(i >= N)
        return;

    if(imove[i] != -3){
        return;
    }

    float shepard_i = shepard[i];
    if(shepard_i < 1.0E-6f){
        // It will be considered that there are not enough
        // particles to interpolate
        shepard_i = 1.f;
    }

    p[i] /= shepard_i;
    // Reversed EOS
    rho[i] = refd[iset[i]] + (p[i] - p0) / (cs * cs);
}