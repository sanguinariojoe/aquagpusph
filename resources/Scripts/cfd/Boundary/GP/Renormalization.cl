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

#include "resources/Scripts/types/types.h"

/** @brief Set the fields for the boundary elements using the pressure
 * interpolated by Boundary/BI/Interpolation.cl.
 *
 * Therefore the pressure value will be renormalizaed, and the density computed
 * using the inverse of the EOS.
 *
 * @param imove Moving flags.
 *   - imove > 0 for regular fluid particles.
 *   - imove = 0 for sensors.
 *   - imove < 0 for boundary elements/particles.
 * @param shepard Shepard term
 * \f$ \gamma(\mathbf{x}) = \int_{\Omega}
 *     W(\mathbf{y} - \mathbf{x}) \mathrm{d}\mathbf{x} \f$.
 * @param gp_rho Interpolated density in the mirrored position \f$ \rho \f$.
 * @param gp_p Interpolated pressure in the mirrored position \f$ p \f$.
 * @param gp_u Interpolated velocity in the mirrored position \f$ \mathbf{u} \f$.
 * @param N Total number of particles and boundary elements.
 */
__kernel void entry(const __global int* imove,
                    const __global float* shepard,
                    __global float* gp_rho,
                    __global float* gp_p,
                    __global vec* gp_u,
                    usize N)
{
    const usize i = get_global_id(0);
    if(i >= N)
        return;

    if(imove[i] != -1){
        return;
    }

    float shepard_i = shepard[i];
    if(shepard_i < 1.0E-6f){
        // It will be considered that there are not enough
        // particles to interpolate
        shepard_i = 1.f;
    }

    gp_p[i] /= shepard_i;
    gp_rho[i] /= shepard_i;
    gp_u[i] /= shepard_i;
}
