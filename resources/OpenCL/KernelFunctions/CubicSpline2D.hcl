/*
 *  This file is part of AQUA-gpusph, a free CFD program based on SPH.
 *  Copyright (C) 2012  Jose Luis Cercos Pita <jl.cercos@upm.es>
 *
 *  AQUA-gpusph is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  (at your option) any later version.
 *
 *  AQUA-gpusph is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with AQUA-gpusph.  If not, see <http://www.gnu.org/licenses/>.
 */

#ifndef _KERNEL_H_INCLUDED_
#define _KERNEL_H_INCLUDED_

#ifndef M_PI
	#define M_PI 3.14159265359f
#endif
#ifndef iM_PI
	#define iM_PI 0.318309886f
#endif

/// @def sep Amount of kernel heights h into the maximum effect distance. 
#ifndef sep
	#define sep 2.f
#endif

/** Method that returns kernel amount with given distance / kernel height rate.
 * @param q distance over kernel height
 * @return Kernel amount
 */
float kernelW(float q)
{
	float wcona = 15.f*iM_PI/7.f;
	float wconb = 5.f*iM_PI/14.f;
	if(q <= 1.f)
		return wcona*(2.f/3.f - q*q + 0.5f*q*q*q);
	else if (q < 2.f)
		return wconb*(2.f-q)*(2.f-q)*(2.f-q);
	return 0.f;
}

/** Method that returns kernel derivative with given distance / kernel height rate.
 * @param q distance over kernel height
 * @return Kernel amount
 */
float kernelF(float q)
{
	float wcona = 15.f*iM_PI/7.f;
	float wconb = 5.f*iM_PI/14.f;
	if(q <= 1.f)
		return wcona*(-2.f + 1.5f*q);
	else if (q < 2.f)
		return -3.f*wconb*(2.f-q)*(2.f-q)/q;
	return 0.f;
}

#endif	// _KERNEL_H_INCLUDED_
