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
 * @brief Spiky kernel definition (2D version). The main feature of the spiky
 * kernel is the uniform sign of the gradient, which should naturally avoid the
 * particles clamping.
 */

#ifndef _KERNEL_H_INCLUDED_
#define _KERNEL_H_INCLUDED_

#define WFAC 0.10807913094111594f
#define FFAC 0.05403956547055797f
#define __SQRT2 1.4142135623730951f
#define __SQRT2POW5 5.656854249492382f

/** @brief The kernel value
 * \f$ W \left(\mathbf{r_j} - \mathbf{r_i}; h\right) \f$.
 * @param q Normalized distance \f$ \frac{\mathbf{r_j} - \mathbf{r_i}}{h} \f$.
 * @return Kernel value.
 */
inline float kernelW(float q)
{
    return WFAC * (__SQRT2 - sqrt(q)) * (2.f - q) * (2.f - q);
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
    if (q == 0.f)
        return 0.f;
    const float sqrt_q = sqrt(q);
    return FFAC * (-5.f * q + __SQRT2POW5 * sqrt_q + 2.f) * (2.f - q)
        / (q * sqrt_q);
}

#endif    // _KERNEL_H_INCLUDED_
