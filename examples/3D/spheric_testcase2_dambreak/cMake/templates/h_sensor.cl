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
    #include "@RESOURCES_DIR@/Scripts/cfd/types/2D.h"
#else
    #include "@RESOURCES_DIR@/Scripts/cfd/types/3D.h"
#endif

/** @brief Set a helper array where the z coordinate of the particles far away
 * from the sensor is vanished.
 * @param imove Moving flags.
 *   - imove > 0 for regular fluid particles.
 *   - imove = 0 for sensors.
 *   - imove < 0 for boundary elements/particles.
 * @param r Position \f$ \mathbf{r} \f$.
 * @param h_sensor_z Helper array of z coordinates, where just the particles
 * close enough to the sensor are preserved.
 * @param N Number of particles.
 * @param h_sensor_x Position of the sensor (h_sensor_y = 0.0).
 * @param dr Distance between particle \f$ \Delta r \f$.
 */
__kernel void entry(const __global int* imove,
                    const __global vec* r,
                    __global float *h_sensor_z,
                    // Simulation data
                    uint N,
                    float h_sensor_x,
                    float dr)
{
    const uint i = get_global_id(0);
    if(i >= N)
        return;

    if(imove[i] <= 0){
        h_sensor_z[i] = 0.f;
        return;
    }

    const float x = r[i].x - h_sensor_x;
    const float y = r[i].y;

    if((fabs(x) > 2.f * dr) || (fabs(y) > 2.f * dr)){
        h_sensor_z[i] = 0.f;
        return;        
    }

    h_sensor_z[i] = r[i].z + 0.5f * dr;
}
 
