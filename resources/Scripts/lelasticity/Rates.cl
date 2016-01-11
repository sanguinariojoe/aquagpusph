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

/** @addtogroup lela
 * @{
 */

/** @file
 * @brief Velocity and density variation rates computation.
 */

#ifndef HAVE_3D
    #include "../types/2D.h"
#else
    #include "../types/3D.h"
#endif

/** @brief Velocity and density variation rates computation.
 *
 * The mass conservation and momentum equations are applied from the already
 * computed differential operators:
 *
 *   - \f$ \frac{\mathrm{d} \mathbf{u}}{\mathrm{d} t} =
 *     - \frac{\nabla \cdot \sigma}{rho}
 *     + \mathbf{g}\f$
 *   - \f$ \frac{\mathrm{d} \rho}{\mathrm{d} t} =
 *     - \rho \nabla \cdot \mathbf{u}\f$
 *
 * @param iset Set of particles index.
 * @param imove Moving flags.
 *   - imove = 2 for regular solid particles.
 *   - imove = 0 for sensors (ignored by this preset).
 *   - imove < 0 for boundary elements/particles.
 * @param div_sigma Divergence of the stress tensor
 * 	   \f$ \frac{\nabla \cdot \sigma}{rho} \f$.
 * @param div_u Velocity divergence \f$ \rho \nabla \cdot \mathbf{u} \f$.
 * @param grad_u Velocity gradient \f$ \nabla \mathbf{u} \f$.
 * @param S Deviatory stress \f$ S \f$.
 * @param dudt Velocity rate of change \f$ \frac{d \mathbf{u}}{d t} \f$.
 * @param drhodt Density rate of change \f$ \frac{d \rho}{d t} \f$.
 * @param dSdt Deviatory stress rate of change \f$ \frac{d S}{d t} \f$.
 * @param shear_mod Shear modulus \f$ \mu \f$.
 * @param N Number of particles.
 * @param g Gravity acceleration \f$ \mathbf{g} \f$.
 */
__kernel void entry(const __global unsigned int* iset,
                    const __global int* imove,
                    const __global vec* div_sigma,
                    const __global float* div_u,
                    const __global matrix* grad_u,
                    const __global matrix* S,
                    __global vec* dudt,
                    __global float* drhodt,
                    __global matrix* dSdt,
                    __constant float* shear_mod,
                    unsigned int N,
                    vec g)
{
    unsigned int i = get_global_id(0);
    if(i >= N)
        return;
    if(imove[i] != 2)
        return;

    const float mu = shear_mod[iset[i]];

    // Momentum equation
    dudt[i] = -div_sigma[i] + g;
    // Conservation of mass equation
    drhodt[i] = -div_u[i];
    // Deviatory stress rate of change
    const matrix epsilon = 0.5f * (grad_u[i] + grad_u[i].TRANSPOSE);
    const matrix omega = 0.5f * (grad_u[i] - grad_u[i].TRANSPOSE);
    dSdt[i] = 2.f * mu * (epsilon - MATRIX_TRACE(epsilon) / DIMS * MAT_EYE)
              + MATRIX_MUL(S[i], omega.TRANSPOSE) + MATRIX_MUL(omega, S[i]);
}

/*
 * @}
 */