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

/** @brief Set the velocity field for the channel boundaries, and the
 * non-inertial terms for the fluid particles.
 *
 * @param iset Set of particles index.
 * @param imove Moving flags.
 *   - imove > 0 for regular fluid particles.
 *   - imove = 0 for sensors.
 *   - imove < 0 for boundary elements/particles.
 * @param u Velocity \f$ \mathbf{u} \f$.
 * @param dudt Velocity rate of change \f$ \frac{d \mathbf{u}}{d t} \f$.
 * @param U Cylinder velocity \f$ \mathbf{U} \f$.
 * @param dUdt Cylinder acceleration \f$ \frac{d \mathbf{U}}{d t} \f$.
 * @param N Number of particles.
 */
__kernel void entry(const __global uint* iset,
                    const __global int* imove,
                    __global vec* u,
                    __global vec* dudt,
                    float U,
                    float dUdt,
                    unsigned int N)
{
    unsigned int i = get_global_id(0);
    if(i >= N)
        return;
    if(imove[i] == 1){
        // It is a fluid particle, affected just by the non-inertial term
        dudt[i].x += dUdt;
    }
    else if((imove[i] == -3) && (iset[i] == 0)){
        // It is a channel boundary, set the velocity
        u[i].x = U;
    }
}
