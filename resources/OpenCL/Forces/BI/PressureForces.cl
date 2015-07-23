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
 * @brief Tool to compute the fluid pressure force and moment.
 */

#ifndef HAVE_3D
    #include "../../types/2D.h"
#else
    #include "../../types/3D.h"
#endif

/** @brief Tool to compute the pressure force and moment for an especific body.
 *
 * In this approach the following operation is performed for the boundary
 * elements:
 * \f$$ f_a = p_a \, \mathbf{n}_a \, s_a \f$$
 * where \f$ s_a \f$ is the area of the element, stored in the masses array.
 * The moment is computed therefore as:
 * \f$$ m_a = p_a \, \mathbf{n}_a \, s_a \times \left(
 *   \mathbf{r}_a - \mathbf{r}_0 \right) \f$$
 * becoming \f$ \mathbf{r}_0 \f$ the reference point where the moment should be
 * computed.
 *
 * @param pressureForces_f Force of each boundary element to be computed [N].
 * @param pressureForces_m Moment of each boundary element to be computed
 * [N \f$ \cdot \f$ m].
 * @param iset Set of particles index.
 * @param imove Moving flags.
 *   - imove > 0 for regular fluid particles.
 *   - imove = 0 for sensors.
 *   - imove < 0 for boundary elements/particles.
 * @param r Position \f$ \mathbf{r} \f$.
 * @param normal Normal \f$ \mathbf{n} \f$.
 * @param p Pressure \f$ p \f$.
 * @param m Mass \f$ m \f$.
 * @param pressureForces_iset Particles set to be computed.
 * @param pressureForces_r Point with respect the moments are computed
 * \f$ \mathbf{r}_0 \f$.
 */
__kernel void main(__global vec* pressureForces_f,
                   __global vec4* pressureForces_m,
                   const __global uint* iset,
                   const __global int* imove,
                   const __global vec* r,
                   const __global vec* normal,
                   const __global float* p,
                   const __global float* m,
                   unsigned int N,
                   unsigned int pressureForces_iset,
                   vec pressureForces_r)
{
    // find position in global arrays
    unsigned int i = get_global_id(0);
    if(i >= N)
        return;
    if((iset[i] != pressureForces_iset) || (imove[i] != -3)){
        pressureForces_f[i] = VEC_ZERO;
        pressureForces_m[i] = (vec4)(0.f, 0.f, 0.f, 0.f);
        return;
    }

    const vec arm = r[i] - pressureForces_r;
    const vec f = p[i] * m[i] * normal[i];
    pressureForces_f[i] = f;
    pressureForces_m[i].z = arm.x * f.y - arm.y * f.x;
    pressureForces_m[i].w = 0.f;
    #ifdef HAVE_3D
        pressureForces_m[i].x = arm.y * f.z - arm.z * f.y;
        pressureForces_m[i].y = arm.z * f.x - arm.x * f.z;
    #else
        pressureForces_m[i].x = 0.f;
        pressureForces_m[i].y = 0.f;
    #endif
}
