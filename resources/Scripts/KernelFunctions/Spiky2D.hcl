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

#ifndef M_PI
    /** @def M_PI
     * \f$ \pi \f$ value.
     */
    #define M_PI 3.14159265359f
#endif
#ifndef iM_PI
    /** @def iM_PI
     * \f$ \frac{1}{\pi} \f$ value.
     */
    #define iM_PI 0.318309886f
#endif

/** @brief The kernel value
 * \f$ W \left(\mathbf{r_j} - \mathbf{r_i}; h\right) \f$.
 * @param q Normalized distance \f$ \frac{\mathbf{r_j} - \mathbf{r_i}}{h} \f$.
 * @return Kernel value.
 */
float kernelW(float q)
{
    const float wcon = 0.357143f * iM_PI;  // 0.357143f = 5/14
    return wcon * (2.f - q) * (2.f - q) * (2.f - q);
}

/** @brief The kernel gradient factor
 * \f$ F \left(\mathbf{r_j} - \mathbf{r_i}; h\right) \f$
 *
 * The factor \f$ F \f$ is defined such that
 * \f$ \nabla W \left(\mathbf{r_j} - \mathbf{r_i}; h\right) =
       \frac{\mathbf{r_j} - \mathbf{r_i}}{h} \cdot
       F \left(\mathbf{r_j} - \mathbf{r_i}; h\right)
   \f$.
 *
 * @param q Normalized distance \f$ \frac{\mathbf{r_j} - \mathbf{r_i}}{h} \f$.
 * @return Kernel gradient factor value
 */
float kernelF(float q)
{
    const float wcon = 1.071429f * iM_PI;  // 1.071429f = 3*5/14
    return wcon * (2.f/q - 1.f) * (2.f/q - 1.f);
}

/** @brief An equivalent kernel function to compute the Shepard factor using the
 * boundary integral elements instead of the fluid volume.
 *
 * The kernel is defined as follows:
 * \f$ \hat{W} \left(\rho; h\right) =
 * \frac{1}{\rho^d} \int \rho^{d - 1} W \left(\rho; h\right) d\rho \f$
 *
 * @param q Normalized distance \f$ \frac{\mathbf{r_j} - \mathbf{r_i}}{h} \f$.
 * @return Equivalent kernel
 */
float kernelS(float q)
{
    float wcon = 4.f * iM_PI / 112.f;
    float q2 = q * q;
    float q3 = q2 * q;
    return wcon * (-2.f * q3 - 15.f * q2 + 40.f * q - 40.f - 16.f / q2 );
}

/** @brief The kernel partial derivative with respect to the characteristic
 * height \f$ \frac{\partial W}{\partial h} \f$
 *
 * The result returned by this function should be multiplied by
 * \f$ \frac{1}{h^{d + 1}} \f$, where d is 2,3 for 2D and 3D respectively.
 *
 * @param q Normalized distance \f$ \frac{\mathbf{r_j} - \mathbf{r_i}}{h} \f$.
 * @return Kernel partial derivative factor
 * @note This function is useful for variable resolution (non-constant kernel
 * height)
 * @see Iason Zisis, Bas van der Linden, Christina Giannopapa and Barry Koren,
 * On the derivation of SPH schemes for shocks through inhomogeneous media. Int.
 * Jnl. of Multiphysics (2015).
 */
float kernelH(float q)
{
    const float wcon = 2.142857f * iM_PI;  // 2.142857f = 15/7
    return wcon * (2.f - q) * (2.f - q) * (1.f - q);
}

#endif    // _KERNEL_H_INCLUDED_
