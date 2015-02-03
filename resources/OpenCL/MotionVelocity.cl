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
 * @brief Euler XYZ based untransformation script.
 */

#ifndef HAVE_3D
    #include "types/2D.h"
#else
    #include "types/3D.h"
#endif

/** @brief Compute the boundary points velocity applying EulerXYZ motion.
 *
 * In EulerXYZ the following transformation is applied to a particle \f$ a \f$:
 * \f[ R_z \cdot R_y \cdot R_x \cdot \mathbf{x_a} + \mathbf{cor}, \f]
 * where \f$ \mathbf{cor} \f$ is the position of the center of rotation (global
 * translations), and \f$ R_x, R_y, R_z \f$ are the rotation matrices:
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
 * Therefore the velocity could be computed as the combination of four terms:
 * \f[ \mathbf[u] =
       \mathbf[u]_1 + \mathbf[u]_2 + \mathbf[u]_3 + \mathbf[u]_4, \f]
 * with:
 * \f[ \mathbf[u]_1 = \dot{R_z} \cdot R_y \cdot R_x \cdot \mathbf{x_a}, \f]
 * \f[ \mathbf[u]_2 = R_z \cdot \dot{R_y} \cdot R_x \cdot \mathbf{x_a}, \f]
 * \f[ \mathbf[u]_3 = R_z \cdot R_y \cdot \dot{R_x} \cdot \mathbf{x_a}, \f]
 * \f[ \mathbf[u]_4 = \dot{\mathbf{cor}}. \f]
 *
 * @param imove Moving flags.
 *   - imove > 0 for regular fluid particles.
 *   - imove = 0 for sensors.
 *   - imove < 0 for boundary elements/particles.
 * @param r Position \f$ \mathbf{r} \f$.
 * @param u Velocity \f$ \mathbf{u} \f$.
 * @param N Number of particles.
 * @param motion_r Center of rotation.
 * @param motion_drdt Center of rotation velocity.
 * @param motion_a Rotation angles \f$ \phi, \theta, \psi \f$.
 * @param motion_dadt Angular velocities.
 * @see MotionUnTransform.cl
 * @see MotionVelocity.cl
 */
__kernel void main(const __global int* imove,
                   __global vec* r,
                   __global vec* u,
                   unsigned int N,
                   vec motion_r,
                   vec motion_drdt,
                   vec4 motion_a,
                   vec4 motion_dadt)
{
    // find position in global arrays
    int i = get_global_id(0);
    if(i >= N)
        return;
    if(imove[i]>0)
        return;

    const vec r_i = r[i];
    vec u1, u2, u3, uu;

    const float cphi = cos(motion_a.x);
    const float sphi = sin(motion_a.x);
    const float ctheta = cos(motion_a.y);
    const float stheta = sin(motion_a.y);
    const float cpsi = cos(motion_a.z);
    const float spsi = sin(motion_a.z);
    const float dphidt = motion_dadt.x;
    const float dthetadt = motion_dadt.y;
    const float dpsidt = motion_dadt.z;

    //---------------------------------------------
    // Compute u1 (acceleration along z)
    //---------------------------------------------
    u1 = r_i;
    #ifdef HAVE_3D
        u1.z = 0.f;
        uu = u1;
        // Rotate along x
        u1.y = cphi * uu.y - sphi * uu.z;
        // Rotate along y
        uu = u1;
        u1.x = ctheta * uu.x + stheta * uu.z;
    #endif
    // Rotate along z
    uu = u1;
    u1.x = dpsidt * (-spsi * uu.x - cpsi * uu.y);
    u1.y = dpsidt * (cpsi * uu.x - spsi * uu.y);

    //---------------------------------------------
    // Compute u2 (acceleration along y)
    //---------------------------------------------
    #ifdef HAVE_3D
        u2 = r_i;
        u2.y = 0.f;
        uu = u2;
        // Rotate along x
        u2.z = sphi * uu.y + cphi * uu.z;
        // Rotate along y
        uu = u2;
        u2.x = dthetadt * (-stheta * uu.x + ctheta * uu.z);
        u2.z = dthetadt * (-ctheta * uu.x - stheta * uu.z);
        // Rotate along z
        uu = u2;
        u2.x = cpsi * uu.x - spsi * uu.y;
    #else
        u2 = VEC_ZERO;
    #endif

    //---------------------------------------------
    // Compute u3 (acceleration along x)
    //---------------------------------------------
    #ifdef HAVE_3D
        u3 = r_i;
        u3.x = 0.f;
        uu = u3;
        // Rotate along x
        u3.y = dphidt * (-sphi * uu.y - cphi * uu.z);
        u3.z = dphidt * (cphi * uu.y - sphi * uu.z);
        // Rotate along y
        uu = u3;
        u3.z = -stheta * uu.x + ctheta * uu.z;
        // Rotate along z
        uu = u3;
        u3.y = spsi * uu.x + cpsi * uu.y;
    #else
        u3 = VEC_ZERO;
    #endif

    // COR velocity
    u[i] = u1 + u2 + u3 + motion_drdt;
}

