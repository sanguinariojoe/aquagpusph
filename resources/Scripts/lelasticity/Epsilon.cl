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
 * @brief Strain-displacement computation
 */

#ifndef HAVE_3D
    #include "../types/2D.h"
#else
    #include "../types/3D.h"
#endif

/** @brief Strain-displacement computation.
 *
 * \f$ \varepsilon = \frac{1}{2} \left(
 *   \nabla \mathbf{r}^{*} + \left( \nabla \mathbf{r}^{*} \right)^{T}
 * \right) \f$
 *
 * @see https://en.wikipedia.org/wiki/Linear_elasticity
 * @see https://en.wikipedia.org/wiki/Hooke's_law
 *
 * @param grad_r Gradient of the deformation \f$ \nabla \mathbf{r}^{*} \f$.
 * @param epsilon Strain-displacement \f$ \varepsilon \f$.
 * @param N Number of particles.
 */
__kernel void entry(const __global matrix* grad_r,
                    __global matrix* epsilon,
                    unsigned int N)
{
    unsigned int i = get_global_id(0);
    if(i >= N)
        return;

    epsilon[i] = 0.5f * (grad_r[i] + grad_r[i].TRANSPOSE);
}

/*
 * @}
 */