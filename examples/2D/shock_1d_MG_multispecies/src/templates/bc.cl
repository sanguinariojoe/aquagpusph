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

/** @defgroup basic Basic preset
 *
 * @brief Basic preset of tools to build more complex sets of tools later
 * 
 * @{
 */

/** @file
 *  @brief 1st order Euler time integration scheme
 *
 * Time integration is based in the following 1st order integration scheme:
 *   - \f$ \mathbf{u}_{n+1} = \mathbf{u}_{n} + \Delta t
        \left. \frac{\mathrm{d}\mathbf{u}}{\mathrm{d}t} \right\vert_{n}
     \f$
 *   - \f$ \mathbf{r}_{n+1} = \mathbf{r}_{n} + \Delta t \, \mathbf{u}_{n}
     + \frac{\Delta t^2}{2}
        \left. \frac{\mathrm{d}\mathbf{u}}{\mathrm{d}t} \right\vert_{n}
     \f$
 *   - \f$ \rho_{n+1} = \rho_{n} + \Delta t
        \left. \frac{\mathrm{d}\rho}{\mathrm{d}t} \right\vert_{n}
     \f$
 */

#include "resources/Scripts/types/types.h"

__kernel void set_fixed(__global int* imove,
                        const __global vec* r,
                        const usize N,
                        const float L)
{
    const usize i = get_global_id(0);
    if(i >= N)
        return;
    if(imove[i] < 0) {
        // Forget about the buffer particles
        return;
    }

    if(fabs(r[i].x) > 0.5 * L - SUPPORT * H)
        imove[i] = 0;
}

__kernel void unset_fixed(__global int* imove,
                          const usize N)
{
    const usize i = get_global_id(0);
    if(i >= N)
        return;
    if(imove[i] < 0) {
        // Forget about the buffer particles
        return;
    }

    imove[i] = 1;
}

__kernel void set_1d(__global vec* dudt,
                     const usize N)
{
    const usize i = get_global_id(0);
    if(i >= N)
        return;
    dudt[i].y = 0.f;
}

