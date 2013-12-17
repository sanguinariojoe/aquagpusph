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

#ifndef DENSITYINTERPOLATION_H_INCLUDED
#define DENSITYINTERPOLATION_H_INCLUDED

// ----------------------------------------------------------------------------
// Include Generic kernel
// ----------------------------------------------------------------------------
#include <CalcServer/Kernel.h>

namespace Aqua{ namespace CalcServer{

/** @class DensityInterpolation DensityInterpolation.h CalcServer/DensityInterpolation.h
 * @brief Density interpolation is a smoothing density field technique. SPH perform density
 * as evolution process, but in order to solve some instabilities a field interpolation
 * can be performed.
 * @warning Don't use Density interpolation if you don't know what are you doing.
 */
struct DensityInterpolation : public Aqua::CalcServer::Kernel
{
	/** Constructor.
	 */
	DensityInterpolation();

	/** Destructor.
	 */
	~DensityInterpolation();

	/** Performs density interpolation.
	 * @return false if all gone right. \n true otherwise.
	 */
	bool execute();

private:
	/** Setup OpenCL kernel
	 * @return false if all gone right. \n true otherwise.
	 */
	bool setupOpenCL();

	/// OpenCL script path
	char* mPath;
	/// OpenCL program
	cl_program clProgram;
	/// OpenCL kernel
	cl_kernel clKernel;
	/// Global work size
	size_t clGlobalWorkSize;
	/// Local work size
	size_t clLocalWorkSize;
	/// true if local memory can be used on kernel.
	bool isLocalMemory;
};

}}  // namespace

#endif // DENSITYINTERPOLATION_H_INCLUDED
