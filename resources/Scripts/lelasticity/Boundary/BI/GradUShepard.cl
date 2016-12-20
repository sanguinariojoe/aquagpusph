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
 * @brief Renormalization of the deformation gradient
 */

#include "resources/Scripts/types/types.h"

/** @brief Renormalization of the deformation gradient.
 *
 * The main drawback of the boundary integrals formulation is the requirement
 * of the renormalization of the computed differentiqal operators, which is
 * destroying several conservation properties.
 *
 * @param imove Moving flags.
 *   - imove = 2 for regular solid particles.
 *   - imove = 0 for sensors.
 *   - imove < 0 for boundary elements/particles.
 * @param shepard Shepard term
 * \f$ \gamma(\mathbf{x}) = \int_{\Omega}
 *     W(\mathbf{y} - \mathbf{x}) \mathrm{d}\mathbf{x} \f$.
 * @param grad_u Gradient of the velocity \f$ \nabla \mathbf{u} \f$.
 * @param N Total number of particles and boundary elements.
 */
__kernel void entry(const __global int* imove,
                    const __global float* shepard,
                    __global matrix* grad_u,
                    uint N)
{
    uint i = get_global_id(0);
    if(i >= N)
        return;
    if(imove[i] != 2){
        return;
    }
    grad_u[i] /= shepard[i];
}

/*
 * @}
 */