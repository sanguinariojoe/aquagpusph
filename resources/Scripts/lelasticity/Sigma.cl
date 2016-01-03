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
 * @brief Stress tensor computation
 */

#ifndef HAVE_3D
    #include "../types/2D.h"
#else
    #include "../types/3D.h"
#endif

/** @brief Stress tensor computation.
 *
 * \f[ \lambda = \frac{E \nu}{(1 + \nu) (1 - 2\nu)} \f]
 * \f[ \mu = \frac{E}{2 (1 + \nu)} \f]
 * \f[ \sigma_{ii} = 2 \mu \varepsilon_{ii} +
 *                   \lambda \sum_{j=1}^{DIMS} \varepsilon_{jj}\f]
 * \f[ \sigma_{ij} = 2 \mu \varepsilon_{ij} \f]
 *
 * @see https://en.wikipedia.org/wiki/Linear_elasticity
 * @see https://en.wikipedia.org/wiki/Hooke's_law
 *
 * @param iset Set of particles index.
 * @param epsilon Strain-displacement \f$ \varepsilon \f$.
 * @param sigma Stress tensor \f$ \sigma \f$.
 * @param young Young's modulus \f$ E \f$.
 * @param poisson Poisson's ratio \f$ \nu \f$.
 * @param N Number of particles.
 */
__kernel void entry(const __global uint* iset,
                    const __global matrix* epsilon,
                    __global matrix* sigma,
                    __constant float* young,
                    __constant float* poisson,
                    unsigned int N)
{
    unsigned int i = get_global_id(0);
    if(i >= N)
        return;

    // Compute the Lame constants
    const float E = young[iset[i]];
    const float nu = poisson[iset[i]];
    const float mu = E / (1.f + nu); 
    const float lambda = mu * nu / (1.f - 2.f * nu); 

    sigma[i] = mu * epsilon[i] + 
    		   lambda * dot(VEC_ONE.XYZ, epsilon[i].DIAG) * MAT_EYE;
}

/*
 * @}
 */