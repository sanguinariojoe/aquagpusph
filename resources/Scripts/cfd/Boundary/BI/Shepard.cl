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

/** @brief Renormalize the differential operators.
 *
 * The main drawback of the boundary integrals formulation is the requirement
 * of the renormalization of the computed differentiqal operators, which is
 * destroying several conservation properties.
 *
 * @param imove Moving flags.
 *   - imove > 0 for regular fluid particles.
 *   - imove = 0 for sensors.
 *   - imove < 0 for boundary elements/particles.
 * @param shepard Shepard term
 * \f$ \gamma(\mathbf{x}) = \int_{\Omega}
 *     W(\mathbf{y} - \mathbf{x}) \mathrm{d}\mathbf{x} \f$.
 * @param grad_p Pressure gradient \f$ \frac{\nabla p}{rho} \f$.
 * @param lap_u Velocity laplacian \f$ \frac{\Delta \mathbf{u}}{rho} \f$.
 * @param div_u Velocity divergence \f$ \rho \nabla \cdot \mathbf{u} \f$.
 * @param lap_p Pressure laplacian \f$ \Delta p \f$.
 * @param N Total number of particles and boundary elements.
 * @param cs Speed of sound \f$ c_s \f$.
 * @see Boundary/BI/Interactions.cl
 */
__kernel void entry(const __global int* imove,
                    const __global float* shepard,
                    __global vec* grad_p,
                    __global vec* lap_u,
                    __global float* div_u,
                    __global float* lap_p,
                    uint N,
                    float cs)
{
    uint i = get_global_id(0);
    if(i >= N)
        return;

    if(imove[i] != 1){
        return;
    }

    float shepard_i = shepard[i];

    grad_p[i] /= shepard_i;
    lap_u[i] /= shepard_i;
    // Following line is causing problems at the free surface
    // div_u[i] /= shepard_i;
}
