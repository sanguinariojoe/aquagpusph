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
 * @brief Density field interpolation interaction computation.
 * (See DensInt.cl)
 *
 * It is prefearable to use a header to be included instead of generating a
 * function for thye particles interaction, which imply more registries
 * consumption.
 */

if(!imove[j]){
    j++;
    continue;
}
#if __BOUNDARY__==0 || __BOUNDARY__==2
    if(imove[j] < 0){
        j++;
        continue;
    }
#endif

const vec r = pos_i - pos[j];
const float q = fast_length(r) / h;
if(q < sep)
{
	//---------------------------------------------------------------
	//       calculate the kernel wab and the function fab
	//---------------------------------------------------------------
	const float wab = kernelW(q) * conw * mass[j];
	//---------------------------------------------------------------
	// 	density computation
	//---------------------------------------------------------------
	_DENS_ += wab;
}
