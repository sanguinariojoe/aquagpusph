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
 * @brief Tool to compute the fluid global energy components.
 */

#ifndef HAVE_3D
    #include "../../types/2D.h"
#else
    #include "../../types/3D.h"
#endif

/** @brief Tool to compute the fluid kinetic energy.
 *
 * \f$ E^{kin} = \sum_{a \in Fluid} \frac{1}{2} m_a
 * \mathbf{u}_a \cdot \mathbf{u}_a\f$
 *
 * @param energy_kin Particle kinetic energy:
 * \f$ E^{kin}_a = \frac{1}{2} m_a \mathbf{u}_a \cdot \mathbf{u}_a\f$
 * @param imove Moving flags.
 *   - imove = 1 for regular fluid particles.
 *   - imove = 0 for sensors.
 *   - imove < 0 for boundary elements/particles.
 * @param u Velocity \f$ \mathbf{u} \f$.
 * @param m Mass \f$ m \f$.
 * @param N Number of particles.
 */
__kernel void entry(__global float* energy_kin,
                    const __global int* imove,
                    const __global vec* u,
                    const __global float* m,
                    unsigned int N)
{
    // find position in global arrays
    unsigned int i = get_global_id(0);
    if(i >= N)
        return;
    if(imove[i] != 1){
        energy_kin[i] = 0.f;
        return;
    }

    energy_kin[i] = 0.5f * m[i] * dot(u[i], u[i]);
}
