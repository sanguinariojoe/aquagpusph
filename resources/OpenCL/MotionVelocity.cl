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
 * @param v Velocity \f$ \mathbf{u} \f$.
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
                   __global vec* v,
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
    vec v1, v2, v3, vv;

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
    // Compute v1 (acceleration along z)
    //---------------------------------------------
    v1 = r_i;
    #ifdef HAVE_3D
        v1.z = 0.f;
        vv = v1;
        // Rotate along x
        v1.y = cphi * vv.y - sphi * vv.z;
        // Rotate along y
        vv = v1;
        v1.x = ctheta * vv.x + stheta * vv.z;
    #endif
    // Rotate along z
    vv = v1;
    v1.x = dpsidt * (-spsi * vv.x - cpsi * vv.y);
    v1.y = dpsidt * (cpsi * vv.x - spsi * vv.y);

    //---------------------------------------------
    // Compute v2 (acceleration along y)
    //---------------------------------------------
    #ifdef HAVE_3D
        v2 = r_i;
        v2.y = 0.f;
        vv = v2;
        // Rotate along x
        v2.z = sphi * vv.y + cphi * vv.z;
        // Rotate along y
        vv = v2;
        v2.x = dthetadt * (-stheta * vv.x + ctheta * vv.z);
        v2.z = dthetadt * (-ctheta * vv.x - stheta * vv.z);
        // Rotate along z
        vv = v2;
        v2.x = cpsi * vv.x - spsi * vv.y;
    #else
        v2 = VEC_ZERO;
    #endif

    //---------------------------------------------
    // Compute v3 (acceleration along x)
    //---------------------------------------------
    #ifdef HAVE_3D
        v3 = r_i;
        v3.x = 0.f;
        vv = v3;
        // Rotate along x
        v3.y = dphidt * (-sphi * vv.y - cphi * vv.z);
        v3.z = dphidt * (cphi * vv.y - sphi * vv.z);
        // Rotate along y
        vv = v3;
        v3.z = -stheta * vv.x + ctheta * vv.z;
        // Rotate along z
        vv = v3;
        v3.y = spsi * vv.x + cpsi * vv.y;
    #else
        v3 = VEC_ZERO;
    #endif

    // COR velocity
    v[i] = v1 + v2 + v3 + motion_drdt;
}

