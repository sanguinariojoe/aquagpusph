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
	float wcon = 0.08203125f*iM_PI;                                             // 0.08203125f = 21/256
	return wcon*(1.f+2.f*q) * (2.f-q)*(2.f-q)*(2.f-q)*(2.f-q);
}

/** Method that returns kernel derivative with given distance / kernel height rate.
 * @param q distance over kernel height
 * @return Kernel amount
 */
float kernelF(float q)
{
	float wcon = 0.8203125f*iM_PI;                                              // 0.8203125f = 10*21/256
	return -wcon*(2.f-q)*(2.f-q)*(2.f-q);                                        // Take care, one q will be added later as r/h
}

#endif	// _KERNEL_H_INCLUDED_
