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

#include "resources/Scripts/types/types.h"

/** @brief Set a helper array where the z coordinate of the particles far away
 * from the sensor is vanished.
 * @param imove Moving flags.
 *   - imove > 0 for regular fluid particles.
 *   - imove = 0 for sensors.
 *   - imove < 0 for boundary elements/particles.
 * @param r Position \f$ \mathbf{r} \f$.
 * @param h_sensory Helper array of y coordinates, where just the particles
 * close enough to the sensor are preserved.
 * @param N Number of particles.
 * @param h_sensorx Position of the sensor (h_sensor_y = 0.0).
 * @param dr Distance between particle \f$ \Delta r \f$.
 * @param depth The maximum depth of the tank. No value below -depth will be
 * accepted
 */
__kernel void entry(const __global int* imove,
                    const __global vec* r,
                    __global float *h_sensory,
                    // Simulation data
                    uint N,
                    float h_sensorx,
                    float dr,
                    float depth)
{
    const uint i = get_global_id(0);
    if(i >= N)
        return;

    if(imove[i] <= 0){
        h_sensory[i] = -depth;
        return;
    }

    const float x = r[i].x - h_sensorx;

    if(fabs(x) > 2.f * dr){
        h_sensory[i] = -depth;
        return;        
    }

    h_sensory[i] = r[i].y + 0.5f * dr;
}
 
