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

/** @defgroup lela Linear elasticity preset
 *
 * @brief Preset to perform Linear elasticity simulations
 * 
 * @{
 */

/** @file
 * @brief Equation Of State (EOS) computation
 */

#ifndef HAVE_3D
    #include "../types/2D.h"
#else
    #include "../types/3D.h"
#endif

/** @brief Deformation computation
 *
 * The Strain-displacement and stress are computed on top of the deformations,
 * i.e. the position of the particles relative to these equilibrium one:
 * \f$ \mathbf{r}^{*} = \mathbf{r} - \mathbf{r}_0 \f$
 *
 * @param r_in Position \f$ \mathbf{r}_{n+1/2} \f$.
 * @param r0 Equilibrium position \f$ \mathbf{r}_{0} \f$.
 * @param r_r0 Deformation \f$ \mathbf{r}^{*} = \mathbf{r} - \mathbf{r}_0 \f$.
 * @param N Number of particles.
 */
__kernel void entry(const __global vec* r_in,
                    const __global vec* r0,
                    __global vec* r_r0,
                    unsigned int N)
{
    unsigned int i = get_global_id(0);
    if(i >= N)
        return;

    r_r0[i] = r_in[i] - r0[i];
}

/*
 * @}
 */