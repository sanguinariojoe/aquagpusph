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
 * @brief Euler-XYZ based velocity computation.
 */

#include "resources/Scripts/types/types.h"

/** @brief Compute the boundary elements velocity applying Euler-XYZ motion.
 *
 * In Euler-XYZ the following transformation is applied to a particle \f$ a \f$:
 * \f[ R_z \cdot R_y \cdot R_x \cdot \mathbf{x_a} + \mathbf{cor}, \f]
 * where \f$ \mathbf{cor} \f$ is the position of the center of rotation (global
 * translations), and \f$ \mathbf{x_a} \f$ is the constant position of the
 * boundary element with respect to \f$ \mathbf{cor} \f$, and
 * \f$ R_x, R_y, R_z \f$ are the rotation matrices:
 * \f[ R_x = \left[ \begin{matrix}
       1 &  0                  &  0                  \\
       0 &  \mathrm{cos}(\phi) & -\mathrm{sin}(\phi) \\
       0 &  \mathrm{sin}(\phi) &  \mathrm{cos}(\phi) \\
   \end{matrix} \right], \f]
 * \f[ R_y = \left[ \begin{matrix}
        \mathrm{cos}(\theta) &  0 &  \mathrm{sin}(\theta) \\
        0                    &  1 &  0                    \\
       -\mathrm{sin}(\theta) &  0 &  \mathrm{cos}(\theta) \\
   \end{matrix} \right], \f]
 * \f[ R_z = \left[ \begin{matrix}
        \mathrm{cos}(\psi) & -\mathrm{sin}(\psi) & 0 \\
        \mathrm{sin}(\psi) &  \mathrm{cos}(\psi) & 0 \\
        0                  &  0                  & 1 \\
   \end{matrix} \right]. \f]
 *
 * To compute the velocity the following process can be followed:
 *   -# The velocity due to the rotations is computed in the local coordinates:
 *      \f$ \omega \times \mathbf{x_a} \f$, with \f$ \omega =
 *      \left[ \dot \phi, \dot \theta, \dot \psi \rigth] \f$
 *   -# Then the vector is rotated using the rotation matrix.
 *   -# Finally the linear velocity, \f$ \dot \mathbf{cor} \f$ is added.
 *
 * @param iset Set of particles index.
 * @param imove Moving flags.
 *   - imove > 0 for regular fluid particles.
 *   - imove = 0 for sensors.
 *   - imove < 0 for boundary elements/particles.
 * @param r Position \f$ \mathbf{r} \f$.
 * @param u Velocity \f$ \mathbf{u} \f$.
 * @param N Number of particles.
 * @param motion_iset Set of particles affected.
 * @param motion_drdt Center of rotation velocity.
 * @param motion_a Rotation angles \f$ \phi, \theta, \psi \f$.
 * @param motion_dadt Angular velocities.
 * @see MotionTransform.cl
 */
__kernel void entry(const __global uint* iset,
                    const __global int* imove,
                    __global vec* r,
                    __global vec* u,
                    usize N,
                    unsigned int motion_iset,
                    vec motion_drdt,
                    vec4 motion_a,
                    vec4 motion_dadt)
{
    // find position in global arrays
    const usize i = get_global_id(0);
    if(i >= N)
        return;
    if((iset[i] != motion_iset) || (imove[i] == 1)){
        return;
    }

    // Compute the velocity due to the rotation in the local frame of reference
    #ifndef HAVE_3D
        vec u_i = (vec)(-motion_dadt.z * r[i].y, motion_dadt.z * r[i].x);
    #else
        vec u_i = cross(motion_dadt, r[i]);
    #endif
    vec uu;

    // Transform it to the global coordinates
    const float cphi = cos(motion_a.x);
    const float sphi = sin(motion_a.x);
    const float ctheta = cos(motion_a.y);
    const float stheta = sin(motion_a.y);
    const float cpsi = cos(motion_a.z);
    const float spsi = sin(motion_a.z);

    #ifdef HAVE_3D
        // Rotate along x
        uu = u_i;
        u_i.y = cphi * uu.y - sphi * uu.z;
        u_i.z = sphi * uu.y + cphi * uu.z;
        // Rotate along y
        uu = u_i;
        u_i.x = ctheta * uu.x + stheta * uu.z;
        u_i.z = -stheta * uu.x + ctheta * uu.z;
    #endif
    // Rotate along z
    uu = u_i;
    u_i.x = cpsi * uu.x - spsi * uu.y;
    u_i.y = spsi * uu.x + cpsi * uu.y;

    // Add the linear velocity
    u[i] = u_i + motion_drdt;
}

