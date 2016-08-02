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
 * @param ilevel0 Level of refinement of the particle, by construction.
 * @param r Position \f$ \mathbf{r} \f$.
 * @param level Target refinement level of the particle.
 * @param gamma_m Mass multiplier \f$ \gamma_m \f$.
 * @param N Number of particles.
 * @param multiresolution_cube_min Minimum point of the refinement area.
 * @param multiresolution_cube_max Maximum point of the refinement area.
 * @param multiresolution_cube_level Refinement level inside the area.
 * @see P.N. Sun, A. Colagrossi, S. Marrone, A.M. Zhang. Multi-resolution
 * delta-SPH model. 2016
 */
__kernel void entry(__global const int* imove,
                    __global const unsigned int* ilevel0,
                    __global const vec* r,
                    __global unsigned int* level,
                    __global float* gamma_m,
                    unsigned int N,
                    vec multiresolution_cube_min,
                    vec multiresolution_cube_max,
                    unsigned int multiresolution_cube_level)
{
    unsigned int i = get_global_id(0);
    if(i >= N)
        return;

    if(imove[i] <= 0) {
        // Neglect boundary elements/particles
        return;
    }

    // We must be careful at the time of selecting the target level of daughter
    // particles, because a mother particle right at the border of the
    // refinement area will generate daughters out of it (which will be removed
    // later).
    vec_xyz r_min = multiresolution_cube_min.XYZ;
    vec_xyz r_max = multiresolution_cube_max.XYZ;
    if(multiresolution_cube_level == ilevel0[i]){
        // The particle seems to be a daughter particle, let's keep it for a
        // while
        r_min -= H * VEC_ONE.XYZ;
        r_max += H * VEC_ONE.XYZ;
    }

    if(    (r[i].x < r_min.x)
        || (r[i].y < r_min.y)
        || (r[i].x > r_max.x)
        || (r[i].y > r_max.y)
        #ifdef HAVE_3D
        || (r[i].z < r_min.z)
        || (r[i].z > r_max.z)
        #endif
      )
    {
        // The particles is just simply not inside the refinement area, forgive
        // it
        return;
    }

    if(level[i] > multiresolution_cube_level){
        // A more dominant refinement area is already managing the particle
        return;
    }

    // We can safely say that the refinement target should be set
    level[i] = multiresolution_cube_level;

    // To select gamma_m, we must take care about eventually overlapped
    // refinement areas. To do that, we should operate in a different way the
    // mother and the daughter particles.
    const vec r_diff = min(r[i].XYZ - multiresolution_cube_min.XYZ,
                           multiresolution_cube_max.XYZ - r[i].XYZ);
    #ifndef HAVE_3D
        const float gamma = min(1.f, max(0.f,
            min(r_diff.x, r_diff.y) / (SUPPORT * H)));
    #else
        const float gamma = min(1.f, max(0.f,
            min(min(r_diff.x, r_diff.y), r_diff.z) / (SUPPORT * H)));
    #endif
    
    if(ilevel0[i] == multiresolution_cube_level - 1){
        // A mother particle, which should vanish during the entrance in this
        // area
        gamma_m[i] = min(gamma_m[i], 1.f - gamma);
    }
    else if(ilevel0[i] < multiresolution_cube_level){
        // An outdated mother particle
        gamma_m[i] = 0.f;
    }
    else if(ilevel0[i] >= multiresolution_cube_level){
        // A daughter particle, generated in this refinement level. We should
        // try to grow the particle as much as possible
        gamma_m[i] = max(gamma_m[i], gamma);
    }
}

/*
 * @}
 */