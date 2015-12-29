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
 * @brief Euler XYZ based acceleration computation.
 */

#ifndef HAVE_3D
    #include "../types/2D.h"
#else
    #include "../types/3D.h"
#endif

/** @brief Compute the boundary elements acceleration applying Euler-XYZ motion.
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
 * To compute the acceleration the following process can be followed:
 *   -# The acceleration due to the rotations is computed in the local
 *      coordinates: \f$ \dot \omega \times \mathbf{x_a} \f$, with
 *      \f$ \dot \omega = \left[ \ddot \phi, \ddot \theta, \ddot \psi \rigth]\f$
 *   -# Then the vector is rotated using the rotation matrix.
 *   -# Finally the linear acceleration, \f$ \ddot \mathbf{cor} \f$ is added.
 *
 * @param iset Set of particles index.
 * @param imove Moving flags.
 *   - imove > 0 for regular fluid particles.
 *   - imove = 0 for sensors.
 *   - imove < 0 for boundary elements/particles.
 * @param r Position \f$ \mathbf{r} \f$.
 * @param dudt Velocity rate of change \f$ \frac{d \mathbf{u}}{d t} \f$.
 * @param N Number of particles.
 * @param motion_iset Set of particles affected.
 * @param motion_r Center of rotation.
 * @param motion_ddrddt Center of rotation aceleration.
 * @param motion_a Rotation angles \f$ \phi, \theta, \psi \f$.
 * @param motion_ddaddt Angular accelerations.
 * @see MotionVelocity.cl
 */
__kernel void entry(const __global uint* iset,
                    const __global int* imove,
                    __global vec* r,
                    __global vec* dudt,
                     unsigned int N,
                    unsigned int motion_iset,
                    vec motion_r,
                    vec motion_ddrddt,
                    vec4 motion_a,
                    vec4 motion_ddaddt)
{
    // find position in global arrays
    int i = get_global_id(0);
    if(i >= N)
        return;
    if((iset[i] != motion_iset) || (imove[i] > 0)){
        return;
    }

    // Compute the velocity due to the rotation in the local frame of reference
    #ifndef HAVE_3D
        vec dudt_i = (vec)(-motion_ddaddt.z * r[i].y, motion_ddaddt.z * r[i].x);
    #else
        vec dudt_i = cross(motion_ddaddt, r[i]);
    #endif
    vec duudt;

    // Transform it to the global coordinates
    const float cphi = cos(motion_a.x);
    const float sphi = sin(motion_a.x);
    const float ctheta = cos(motion_a.y);
    const float stheta = sin(motion_a.y);
    const float cpsi = cos(motion_a.z);
    const float spsi = sin(motion_a.z);

    #ifdef HAVE_3D
        // Rotate along x
        duudt = dudt_i;
        dudt_i.y = cphi * duudt.y - sphi * duudt.z;
        dudt_i.z = sphi * duudt.y + cphi * duudt.z;
        // Rotate along y
        duudt = dudt_i;
        dudt_i.x = ctheta * duudt.x + stheta * duudt.z;
        dudt_i.z = -stheta * duudt.x + ctheta * duudt.z;
    #endif
    // Rotate along z
    duudt = dudt_i;
    dudt_i.x = cpsi * duudt.x - spsi * duudt.y;
    dudt_i.y = spsi * duudt.x + cpsi * duudt.y;

    // Add the linear velocity
    dudt[i] = dudt_i + motion_ddrddt;
}

