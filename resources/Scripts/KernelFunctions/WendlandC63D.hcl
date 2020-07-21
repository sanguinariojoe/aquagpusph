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
 * @brief Wendland kernel definition (3D version).
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
const float kernelW(float q)
{
    const float wcon  = 0.01041412353515625f*iM_PI;  // 0.01041412353515625f = 1365 / 64 / 2**11
    const float tmq = 2.f - q;
    const float facq8 = tm2 * tmq * tmq * tmq * tmq * tmq * tmq * tmq;
    return wcon * (1.f + 4.f*q + 6.25f*q*q + 4.f*q*q*q) * facq8;
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
 * @return Kernel amount
 */
const float kernelF(float q)
{
    const float wcon  = 0.01041412353515625f*iM_PI;  // 0.01041412353515625f = 1365 / 64 / 2**11
    const float tmq = 2.f - q;
    const float facq7 = tmq * tmq * tmq * tmq * tmq * tmq * tmq;
    return wcon * facq7 * (-11.f - 38.5f*q - 44.f*q*q);
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
const float kernelS_P(float q)
{
    const float wcon = 0.005440848214285714f*iM_PI;  // 0.00544085f = 7/14336 = 7/(78 * 2**11)  
    const float q2   = q * q;
    const float q5   = q2 * q2 * q;
    return wcon * ( 0.2857142857142857f * q5 * q5 * q   - // 0.2857142857142857 = 2/7
		    4.4423076923076925f * q5 * q5       + // 4.4423076923076925 = 231/52
		    29.333333333333333.f * q5 * q2 * q2 - // 29.333333333333333 = 88/3
		    105.f * q5 * q2 * q                 + // 
		    211.2f * q5 * q2                    - // 211.2 = 2112/10
		    205.33333333333334.f * q5 * q       + // 205.33333333333333 = 616/3
		    150.85714285714286f * q2 * q2       - // 150.85714285714286 = 1056/7
		    140.8f * q2                         + // 140.8f = 704/5
		    85.33333333333333f                    // 85.33333333333333f = 256/3
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
 * \f$ (\mathbf{r_j} - \mathbf{r_i}) \cdot \mathbf{t_j} \f$.
 * @param b 0 in 2D simulations, distance along the normal direction in 3D
 * simulations, \f$ (\mathbf{r_j} - \mathbf{r_i}) \cdot \mathbf{b_j} \f$.
 * @param s Area of the boundary element, \f$ 2 * \Delta r \f$.
 * @return Equivalent kernel divergent part
 * @see kernelS_P
 * @warning Due to the analytical nature of the solution, this kernel should not
 * be multiplied by the element area, nor divided by \f$ h^2 \f$
 */
const float kernelS_D(float d, float t, float b, float s)
{
    const float wcon = 0.5f * iM_PI;
    const float dr = 0.5f * s;
    return -wcon * (atan((t + dr) / d) - atan((t - dr) / d));
}

/** @brief Helper function to compute the solid angle of a rectangular patch,
 * with a corner placed in the projection of the origin into the patch plane.
 *
 * b
 *  A
 *  |
 *  |XXXXXXX
 *  |XXXXXXX
 * -+-----------> t
 *
 * @param d Normal distance of the patch to the origin.
 * @param a Width of the rectangular patch.
 * @param b Height of the rectangular patch.
 * @return Equivalent kernel divergent part
 * @see kernelS_D
 */
const float _Omega(const float a, const float b)
{
    const float a2 = a * a;
    const float b2 = b * b;
    return acospi(min(sqrt(
        (1.f + a2 + b2) / ((1.f + a2) * (1.f + b2))
    ), 1.f));
}

/* @brief An equivalent kernel function to compute the Shepard factor using the
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
 * @param d Unsigned normal distance to the wall,
 * \f$ \vert (\mathbf{r_j} - \mathbf{r_i}) \cdot \mathbf{n_j} \vert \f$.
 * @param t Unsigned Tangential distance to the boundary element,
 * \f$ \vert (\mathbf{r_j} - \mathbf{r_i}) \cdot \mathbf{t_j} \vert \f$.
 * @param b 0 in 2D simulations, unsigned distance along the binormal direction
 * in 3D simulations,
 * \f$ \vert (\mathbf{r_j} - \mathbf{r_i}) \cdot \mathbf{b_j} \vert  \f$.
 * @param s Area of the boundary element, \f$ \Delta r^2 \f$.
 * @return Equivalent kernel divergent part
 * @see kernelS_P
 * @warning Due to the analytical nature of the solution, this kernel should not
 * be multiplied by the element area, nor divided by \f$ h^2 \f$
 */
const float kernelS_D(float d, float t, float b, float s)
{
    const float wcon = 0.25f;
    const float dr = 0.5f * sqrt(s);
    const float t1 = (t - dr) / d, t2 = (t + dr) / d;
    const float b1 = (b - dr) / d, b2 = (b + dr) / d;
    const float st = sign(t1), sb = sign(b1);
    return -wcon * (_Omega(t2, b2)
                    - st * _Omega(t1, b2)
                    - sb * _Omega(t2, b1)
                    + st * sb * _Omega(t1, b1));
}

#endif    // _KERNEL_H_INCLUDED_
