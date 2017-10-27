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
 * @brief Wendland kernel definition (2D version).
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
    float wcon = 0.109375f*iM_PI;  // 0.109375f = 7/64
    return wcon*(1.f+2.f*q) * (2.f-q)*(2.f-q)*(2.f-q)*(2.f-q);
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
    float wcon = 1.09375f*iM_PI;  // 1.09375f = 2*5*7/64
    return wcon*(2.f-q)*(2.f-q)*(2.f-q);
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
    float wcon = 0.109375f*iM_PI;  // 0.109375f = 7/64
    return wcon * (2.f-q)*(2.f-q)*(2.f-q) * (
        4.f * q * (1.f + 2.f * q) -
        DIMS * (2.f-q) * (1.f + 2.f * q) -
        2.f * q * (2.f-q));
}

/** @brief An equivalent kernel function to compute the Shepard factor using the
 * boundary integral elements instead of the fluid volume.
 *
 * For practical purposes, the kernel computation is split in 2 parts: The
 * polynomial part, where trucation errors are acceptable, and the divergent
 * part, which requires an analytical solution.
 * This function computes the polynomial part.
 *
 * The kernel is defined as follows:
 * \f$ \hat{W} \left(\rho; h\right) =
 * \frac{1}{\rho^d} \int \rho^{d - 1} W \left(\rho; h\right) d\rho \f$
 *
 * @param q Normalized distance \f$ \frac{\mathbf{r_j} - \mathbf{r_i}}{h} \f$.
 * @return Equivalent kernel polynomial part
 * @see kernelS_D
 */
float kernelS_P(float q)
{
    float wcon = 0.109375f*iM_PI;  // 0.109375f = 7/64
    float q2 = q * q;
    float q3 = q2 * q;
    return wcon * (  0.285714f * q3 * q2  // 0.285714f = 2/7
                   - 2.5f * q2 * q2       // 2.5f = 5/2
                   + 8.f * q3
                   - 10.f * q2
                   + 8.f
                  );
}

/** @brief An equivalent kernel function to compute the Shepard factor using the
 * boundary integral elements instead of the fluid volume.
 *
 * For practical purposes, the kernel computation is split in 2 parts: The
 * polynomial part, where trucation errors are acceptable, and the divergent
 * part, which requires an analytical solution.
 * This function computes the divergent part.
 *
 * The kernel is defined as follows:
 * \f$ \hat{W} \left(\rho; h\right) =
 * \frac{1}{\rho^d} \int \rho^{d - 1} W \left(\rho; h\right) d\rho \f$
 *
 * @param d Normal distance to the wall,
 * \f$ (\mathbf{r_j} - \mathbf{r_i}) \cdot \mathbf{n_j} \f$.
 * @param t Tangential distance to the boundary element,
 * \f$ \vert (\mathbf{r_j} - \mathbf{r_i}) - \mathbf{n_j} \left(
 *      \left( \mathbf{r_j} - \mathbf{r_i} \right) \cdot \mathbf{n_j} 
 * \right) \vert \f$.
 * @param s Area of the boundary element, \f$ 2 * \Delta r \f$.
 * @return Equivalent kernel divergent part
 * @see kernelS_P
 * @warning Due to the analytical nature of the solution, this kernel should not
 * be multiplied by the element area, nor divided by \f$ h^2 \f$
 */
float kernelS_D(d, t, s)
{
    const float wcon = 0.5f * iM_PI;
    const float dr = 0.5f * s;
    return wcon / d * \
           (atan((t + dr) / d) - atan((t - dr) / d));
}

#endif    // _KERNEL_H_INCLUDED_
