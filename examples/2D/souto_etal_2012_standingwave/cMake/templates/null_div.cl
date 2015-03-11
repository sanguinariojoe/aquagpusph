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

#ifndef HAVE_3D
    #include "@RESOURCES_DIR@/OpenCL/types/2D.h"
#else
    #include "@RESOURCES_DIR@/OpenCL/types/3D.h"
#endif

/** @brief Set a helper array where the z coordinate of the particles far away
 * from the sensor is vanished.
 * @param imove Moving flags.
 *   - imove > 0 for regular fluid particles.
 *   - imove = 0 for sensors.
 *   - imove < 0 for boundary elements/particles.
 * @param div_u Velocity divergence \f$ \rho \nabla \cdot \mathbf{u} \f$.
 * @param N Number of particles.
 */
__kernel void main(const __global int* imove,
                   __global float *div_u,
                   // Simulation data
                   uint N)
{
    const uint i = get_global_id(0);
    if(i >= N)
        return;
    if(imove[i] <= 0){
        return;
    }

    div_u[i] = 0.f;
}
 
