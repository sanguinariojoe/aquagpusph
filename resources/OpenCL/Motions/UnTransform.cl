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
 * @brief Euler-XYZ based untransformation script.
 */

#ifndef HAVE_3D
    #include "../types/2D.h"
#else
    #include "../types/3D.h"
#endif

/** @brief Untransform the previously applied EulerXYZ motion.
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
 * Therefore to invert the transformation the following expression can be
 * applied:
 * \f[ R_x^{-1} \cdot R_y^{-1} \cdot R_z^{-1} \cdot
   \left(\mathbf{x_a} - \mathbf{cor}\right), \f]
 * where the inverse rotation matrices are obtained just using
 * \f$ -\phi, -\theta, -\psi \f$ angles.
 *
 * @param iset Set of particles index.
 * @param imove Moving flags.
 *   - imove > 0 for regular fluid particles.
 *   - imove = 0 for sensors.
 *   - imove < 0 for boundary elements/particles.
 * @param r Position \f$ \mathbf{r} \f$.
 * @param normal Normal \f$ \mathbf{n} \f$.
 * @param N Number of particles.
 * @param motion_iset Set of particles affected.
 * @param motion_r_in Center of rotation.
 * @param motion_a_in Rotation angles \f$ \phi, \theta, \psi \f$.
 * @see MotionTransform.cl
 */
__kernel void entry(const __global uint* iset,
                    const __global int* imove,
                    __global vec* r,
                    __global vec* normal,
                    unsigned int N,
                    unsigned int motion_iset,
                    vec motion_r_in,
                    vec4 motion_a_in)
{
    // find position in global arrays
    int i = get_global_id(0);
    if(i >= N)
        return;
    if((iset[i] != motion_iset) || (imove[i] > 0)){
        return;
    }

    vec r_i, rr, n_i, nn;

    const float cphi = cos(motion_a_in.x);
    const float sphi = -sin(motion_a_in.x);
    const float ctheta = cos(motion_a_in.y);
    const float stheta = -sin(motion_a_in.y);
    const float cpsi = cos(motion_a_in.z);
    const float spsi = -sin(motion_a_in.z);

    //---------------------------------------------
    // Untransform the point
    //---------------------------------------------
    n_i = normal[i];
    // Undisplace the point
    r_i = r[i] - motion_r_in;
    // Unrotate along z
    rr = r_i;
    nn = n_i;
    r_i.x = cpsi * rr.x - spsi * rr.y;
    r_i.y = spsi * rr.x + cpsi * rr.y;
    n_i.x = cpsi * nn.x - spsi * nn.y;
    n_i.y = spsi * nn.x + cpsi * nn.y;
    #ifdef HAVE_3D
        // Unrotate along y
        rr = r_i;
        nn = n_i;
        r_i.x = ctheta * rr.x + stheta * rr.z;
        r_i.z = -stheta * rr.x + ctheta * rr.z;
        n_i.x = ctheta * nn.x + stheta * nn.z;
        n_i.z = -stheta * nn.x + ctheta * nn.z;
        // Unrotate along x
        rr = r_i;
        nn = n_i;
        r_i.y = cphi * rr.y - sphi * rr.z;
        r_i.z = sphi * rr.y + cphi * rr.z;
        n_i.y = cphi * nn.y - sphi * nn.z;
        n_i.z = sphi * nn.y + cphi * nn.z;
    #endif

    normal[i] = n_i;
    r[i] = r_i;
}