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

/**
 * \addtogroup ideal_gas
 * @{
 */

/** @file
 * @brief Internal energy variation rates computation.
 *
 * This is an extension of resources/Scripts/cfd/Rates.cl
 */

#include "resources/Scripts/types/types.h"

/** @brief Internal energy variation rates computation.
 *
 * The internal energy equation is applied from the already computed
 * differential operators:
 *
 *   - \f$ \frac{\mathrm{d} \e}{\mathrm{d} t} =
 *     - \frac{p}{\rho} \nabla \cdot \mathbf{u} \f$
 *
 * @param iset Set of particles index.
 * @param imove Moving flags.
 *   - imove > 0 for regular fluid particles.
 *   - imove = 0 for sensors.
 *   - imove < 0 for boundary elements/particles.
 * @param rho Density \f$ \rho_{n+1} \f$.
 * @param p Pressure \f$ p_{n+1/2} \f$.
 * @param div_u Velocity divergence \f$ \rho \nabla \cdot \mathbf{u} \f$.
 * @param deintdt Internal energy rate of change
 * \f$ \left. \frac{d e}{d t} \right\vert_{n+1} \f$.
 * @param N Number of particles.
 */
__kernel void entry(const __global uint* iset,
                    const __global int* imove,
                    const __global float* rho,
                    const __global float* p,
                    const __global float* div_u,
                    __global float* deintdt,
                    const usize N)
{
    const usize i = get_global_id(0);
    if(i >= N)
        return;
    if(imove[i] != 1)
        return;

    deintdt[i] = -p[i] / (rho[i] * rho[i]) * div_u[i];
}

/*
 * @}
 */
