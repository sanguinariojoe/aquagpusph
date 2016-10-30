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

/** @addtogroup basic
 * @{
 */

/** @file
 *  @brief Splitting particles methods
 */

#ifndef HAVE_3D
    #include "../../types/2D.h"
#else
    #include "../../types/3D.h"
#endif

/** @brief Set the refinement level of the particles inside.
 *
 * @param imove Moving flags.
 *   - imove > 0 for regular fluid/solid particles.
 *   - imove = 0 for sensors.
 *   - imove < 0 for boundary elements/particles.
 * @param ilevel Current refinement level of the particle.
 * @param r Position \f$ \mathbf{r} \f$.
 * @param level Target refinement level of the particle.
 * @param N Number of particles.
 * @param multiresolution_sphere_center Center of the refinement area.
 * @param multiresolution_sphere_radius Radius of the refinement area.
 * @param multiresolution_sphere_level Refinement level inside the area.
 */
__kernel void entry(__global const int* imove,
                    __global const unsigned int* ilevel,
                    __global const vec* r,
                    __global unsigned int* level,
                    unsigned int N,
                    vec multiresolution_sphere_center,
                    float multiresolution_sphere_radius,
                    unsigned int multiresolution_sphere_level)
{
    unsigned int i = get_global_id(0);
    if(i >= N)
        return;

    if(imove[i] <= 0) {
        // Neglect boundary elements/particles
        return;
    }

    vec_xyz r_c = r[i].XYZ - multiresolution_sphere_center.XYZ;
    float radius = multiresolution_sphere_radius;
    // We must be careful at the time of selecting the target level of daughter
    // particles, because a mother particle right at the border of the
    // refinement area will generate daughters out of it (which will be removed
    // inmediately after this).
    if(multiresolution_sphere_level == ilevel[i]){
        // The particle seems to be a daughter particle, let's keep it for a
        // while
        radius += H;
    }

    if(dot(r_c, r_c) > radius * radius){
        // The particle is just simply outside the refinement area
        return;
    }

    if(level[i] > multiresolution_sphere_level){
        // A more dominant refinement area is already managing the particle
        return;
    }

    // We can safely say that the refinement target should be set
    level[i] = multiresolution_sphere_level;
}

/*
 * @}
 */