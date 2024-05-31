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
 * @brief Spiky kernel definition (3D version). The spiky kernel is a common
 * Cubic Spline kernel, where the middle know is displaced to grant that the
 * maximum of the gradient lies on 1 / hfac
 */

#ifndef _KERNEL_H_INCLUDED_
#define _KERNEL_H_INCLUDED_

#ifndef M_PI
    /** @def M_PI
     * \f$ \pi \f$ value.
     */
    #define M_PI 3.14159265359f
#endif

/// Spline middle knot
#define A 2.f / (4.f * HFAC - 1.f)
static float a = A;
/// Renormalization factors
#define WCON 15.f * (A - 2.f) * A * A * (A + 2.f) /                            \
    (M_PI * (10.f * A * A * A - 28.f * A * A + 15.f * A - 4.f))
#define FCON 3.f / 2.f * WCON
static float wcon = WCON;
static float fcon = FCON;

/** @brief The kernel value
 * \f$ W \left(\mathbf{r_j} - \mathbf{r_i}; h\right) \f$.
 * @param q Normalized distance \f$ \frac{\mathbf{r_j} - \mathbf{r_i}}{h} \f$.
 * @return Kernel value.
 */
inline float kernelW(float q)
{
    if (q < a)
        return wcon * ((a + 2.f) * q * q * q - 6.f * a * q * q + 4.f * a * a) /
            (2.f * a * a * (a + 2.f));
    return wcon * (2.f - q) * (2.f - q) * (2.f - q);
}

/** @brief The kernel gradient factor
 * \f$ F \left(\mathbf{r_j} - \mathbf{r_i}; h\right) \f$
 *
 * The factor \f$ F \f$ is defined such that
 * \f$ \nabla W \left(\mathbf{r_j} - \mathbf{r_i}; h\right) =
       \frac{\mathbf{r_j} - \mathbf{r_i}}{h^d} \cdot
       F \left(\mathbf{r_j} - \mathbf{r_i}; h\right)
   \f$.
 *
 * @param q Normalized distance \f$ \frac{\mathbf{r_j} - \mathbf{r_i}}{h} \f$.
 * @return Kernel gradient factor value
 */
inline float kernelF(float q)
{
    if (q < a)
        return -fcon * ((a + 2.f) * q - 4.f * a) / (a * a * (a + 2.f));
    return -fcon * ((2.f - q) * (2.f - q)) / ((a - 2.f) * (a + 2.f) * q);
}

#endif    // _KERNEL_H_INCLUDED_
