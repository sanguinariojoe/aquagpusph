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
    #include "types/2D.h"
#else
    #include "types/3D.h"
#endif

/** @def _PRESSURE_BASED_
 * @brief Set the pressure and density values from the interpolated (and
 * hydrostatically corrected) pressure value in the point.
 *
 * To set the boundary elements pressure and density values two approaches
 * may be considred (just one of them must be activated):
 *   - _ALL_INTERPOLATED_: Proposed by Ferrand, the density and pressure
 *     are interpolated from the neighbour particles.
 *   - _DENSITY_BASED_: Just the density is interpolated, computing the
 *     pressure using the EOS.
 *   - _PRESSURE_BASED_: The pressure value is interpolated, and the
 *     density is computed the EOS.
 *
 * Pressure field is corrected (in all the cases) by the hydrostatic field.
 * Since density is not corrected, if you are using weakly compressible
 * (where the fluid pressure and density are related) consideer using
 * the _PRESSURE_BASED_ algorithm.
 */
#define _PRESSURE_BASED_

/** @brief Set the fields for the boundary elements using the data interpolated
 * by Aqua::CalcServer::Rates.
 * @param iset Set of particles index.
 * @param imove Moving flags.
 *   - imove > 0 for regular fluid particles.
 *   - imove = 0 for sensors.
 *   - imove < 0 for boundary elements/particles.
 * @param grad_p Pressure gradient \f$ \nabla p \f$.
 * @param div_u Velocity divergence \f$ \nabla \cdot \mathbf{u} \f$.
 * @param rho Density \f$ \rho \f$.
 * @param p Pressure \f$ p \f$.
 * @param shepard Shepard term
 * \f$ \gamma(\mathbf{x}) = \int_{\Omega}
 *     W(\mathbf{y} - \mathbf{x}) \mathrm{d}\mathbf{x} \f$.
 * @param refd Density of reference of the fluid \f$ \rho_0 \f$.
 * @param gamma Gamma EOS exponent \f$ \gamma \f$.
 * @param N Total number of particles and boundary elements.
 * @param cs Speed of sound \f$ c_s \f$.
 * @see BoundaryIntegrals.cl
 */
__kernel void main(const __global uint* iset,
                   const __global int* imove,
                   __global vec* grad_p,
                   __global float* div_u,
                   __global float* rho,
                   __global float* p,
                   const __global float* shepard,
                   __constant float* refd,
                   __constant float* gamma,
                   uint N,
                   float cs)
{
    uint i = get_global_id(0);
    if(i >= N)
        return;

    if((imove[i] == -1) || (imove[i] > 0)){
        return;
    }

    float shepard_i = shepard[i];
    if(shepard_i < 1.0E-6f){
        // It will be considered that there are not enough
        // particles to interpolate
        shepard_i = 1.f;
        grad_p[i] = VEC_ZERO;
        div_u[i] = 0.f;
    }

    #ifdef _ALL_INTERPOLATED_
        rho[i] = div_u[i] / shepard_i;
        p[i] = (grad_p[i].x + grad_p[i].y) / shepard_i;
    #elif defined _DENSITY_BASED_
        rho[i] = div_u[i] / shepard_i;
        // Batchelor 1967
        const float rdenf = refd[iset[i]];
        const float gammf = gamma[iset[i]];
        const float ddenf = rho[i] / rdenf;
        const float prb = cs * cs * rdenf / gammf;
        p[i] = max(prb * (pow(ddenf, gammf) - 1.f),
                   grad_p[i].x / shepard_i) + grad_p[i].y / shepard_i;
    #elif defined _PRESSURE_BASED_
        // p[i] = max(0.f, (grad_p[i].x + grad_p[i].y) / shepard_i);
        p[i] = (grad_p[i].x + grad_p[i].y) / shepard_i;
        // Reversed Batchelor 1967
        const float rdenf = refd[iset[i]];
        const float gammf = gamma[iset[i]];
        const float iprb = gammf / (cs * cs * rdenf);
        rho[i] = rdenf * pow(1.f + iprb * p[i], 1.f / gammf);
    #else
        #error "Unknow boundary elements field computation algorithm"
    #endif
    // Rinitializate output variables
    grad_p[i] = VEC_ZERO;
    div_u[i] = 0.f;
}
