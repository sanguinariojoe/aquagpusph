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

/** @brief Transform the boundaries applying EulerXYZ motion.
 *
 * Just the position and the normal of the particle are modified, but not the
 * velocity which is changed by MotionVelocity.cl.
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
 * @param imove Moving flags.
 *   - imove > 0 for regular fluid particles.
 *   - imove = 0 for sensors.
 *   - imove < 0 for boundary elements/particles.
 * @param r Position \f$ \mathbf{r} \f$.
 * @param normal Normal \f$ \mathbf{n} \f$.
 * @param N Number of particles.
 * @param motion_r Center of rotation.
 * @param motion_a Rotation angles \f$ \phi, \theta, \psi \f$.
 * @see MotionUnTransform.cl
 * @see MotionVelocity.cl
 */
__kernel void main(const __global int* imove,
                   __global vec* r,
                   __global vec* normal,
                   unsigned int N,
                   vec motion_r,
                   vec4 motion_a)
{
    // find position in global arrays
    int i = get_global_id(0);
    if(i >= N)
        return;
    if(imove[i]>0)
        return;

    vec r_i, rr, n_i, nn;

    const float cphi = cos(motion_a.x);
    const float sphi = sin(motion_a.x);
    const float ctheta = cos(motion_a.y);
    const float stheta = sin(motion_a.y);
    const float cpsi = cos(motion_a.z);
    const float spsi = sin(motion_a.z);

    //---------------------------------------------
    // Transform the point
    //---------------------------------------------
    n_i = normal[i];
    r_i = r[i];
    #ifdef HAVE_3D
        // Rotate along x
        rr = r_i;
        nn = n_i;
        r_i.y = cphi * rr.y - sphi * rr.z;
        r_i.z = sphi * rr.y + cphi * rr.z;
        n_i.y = cphi * nn.y - sphi * nn.z;
        n_i.z = sphi * nn.y + cphi * nn.z;
        // Rotate along y
        rr = r_i;
        nn = n_i;
        r_i.x = ctheta * rr.x + stheta * rr.z;
        r_i.z = -stheta * rr.x + ctheta * rr.z;
        n_i.x = ctheta * nn.x + stheta * nn.z;
        n_i.z = -stheta * nn.x + ctheta * nn.z;
    #endif
    // Rotate along z
    rr = r_i;
    nn = n_i;
    r_i.x = cpsi * rr.x - spsi * rr.y;
    r_i.y = spsi * rr.x + cpsi * rr.y;
    n_i.x = cpsi * nn.x - spsi * nn.y;
    n_i.y = spsi * nn.x + cpsi * nn.y;

    normal[i] = n_i;
    // Displace the point
    r[i] = r_i + motion_r;
}

