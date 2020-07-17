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
 * @brief Wendland kernel C6 definition (2D version).
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
	float wcon  = 0.010881696428571428f*iM_PI;  // 0.0108817f = 7/7168 = 7/(78 * 2**10)
	float facq8 = (2.f-q)*(2.f-q)*(2.f-q)*(2.f-q)*(2.f-q)*(2.f-q)*(2.f-q)*(2.f-q);
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
 * @return Kernel gradient factor value
 */
float kernelF(float q)
{
	float wcon  = 0.010881696428571428f*iM_PI;  // 0.0108817f = 7/7168 = 7/(78 * 2**10)
	float facq7 = (2.f-q)*(2.f-q)*(2.f-q)*(2.f-q)*(2.f-q)*(2.f-q)*(2.f-q);
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
float kernelS_P(float q)
{
    float wcon = 0.010881696428571428f*iM_PI;  // 0.0108817f = 78/7168 = 78/(7 * 2**10)
    float q2 = q * q;
    float q5 = q2 * q2 * q;
    return wcon * ( 0.3076923076923077f * q5 * q5 * q  - // 0.3076923076923077f = 4/13
		    4.8125f * q5 * q5                  + // 4.8125f = 77/16
		    32.f * q5 * q2 * q2                - //
		    115.5f * q5 * q2 * q               + // 115.5f = 1155/10
		    234.66666666666667f * q5 * q2      - // 234.66666667f = 704/3
		    231.f * q5 * q                     + //
		    176.f * q2 * q2                    - //
		    176.f * q2                         + //
		    128.f
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
float kernelS_D(float d, float t, float b, float s)
{
    const float wcon = 0.5f * iM_PI;
    const float dr = 0.5f * s;
    return -wcon * (atan((t + dr) / d) - atan((t - dr) / d));
}

#endif	// _KERNEL_H_INCLUDED_
