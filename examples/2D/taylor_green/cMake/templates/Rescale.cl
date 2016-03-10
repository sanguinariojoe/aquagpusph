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
 * @brief Set the velocity according to the initialization stage.
 */

#ifndef HAVE_3D
    #include "@RESOURCES_DIR@/Scripts/types/2D.h"
#else
    #include "@RESOURCES_DIR@/Scripts/types/3D.h"
#endif

/** @brief Rescale the pressure and the velocity of the particles to restore
 * the energy dissipated.
 *
 * @param iset Set of particles index.
 * @param u Velocity \f$ \mathbf{u} \f$.
 * @param p Pressure \f$ p \f$.
 * @param visc_dyn Dynamic viscosity \f$ \mu \f$.
 * @param refd Density of reference of the fluid \f$ \rho_0 \f$.
 * @param N Number of particles.
 * @param dt Time step \f$ \Delta t \f$.
 */
__kernel void entry(const __global uint* iset,
                    __global vec* u,
                    __global float* p,
                    __constant float* visc_dyn,
                    __constant float* refd,
                    unsigned int N,
                    float dt)
{
    unsigned int i = get_global_id(0);
    if(i >= N)
        return;

    const float visc_kin = visc_dyn[iset[i]] / refd[iset[i]];
    
    u[i].XYZ += u[i].XYZ * expm1(2.f * visc_kin * dt);
    p[i] += p[i] * expm1(4.f * visc_kin * dt);
}
