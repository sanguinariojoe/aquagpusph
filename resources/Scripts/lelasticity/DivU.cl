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
 * @brief Velocity divergence computation
 */

#ifndef HAVE_3D
    #include "../types/2D.h"
#else
    #include "../types/3D.h"
#endif

/** @brief Velocity divergence computation
 *
 * Since we already computed the gradient of the velocity, we can get the
 * divergence as the sum of the diagonal components:
 * \f[ \nabla \cdot \mathbf{u} = \sum_i^d \left(\nabla \mathbf{u}\right)_ii \f]
 *
 * @param imove Moving flags.
 *   - imove = 2 for regular solid particles.
 *   - imove = 0 for sensors (ignored by this preset).
 *   - imove < 0 for boundary elements/particles.
 * @param rho Density \f$ \rho \f$.
 * @param grad_u Gradient of the deformation \f$ \nabla \mathbf{r}^{*} \f$.
 * @param div_u Velocity divergence \f$ \rho \nabla \cdot \mathbf{u} \f$.
 * @param N Number of particles.
 */
__kernel void entry(const __global int* imove,
                    const __global float* rho,
                    const __global matrix* grad_u,
                    __global float* div_u,
                    unsigned int N)
{
    unsigned int i = get_global_id(0);
    if(i >= N)
        return;
    if(imove[i] != 2)
        return;

    div_u[i] = rho[i] * dot(VEC_ONE.XYZ, grad_u[i].DIAG);
}

/*
 * @}
 */