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

#ifndef ELASTICBOUNCE_H_INCLUDED
#define ELASTICBOUNCE_H_INCLUDED

// ----------------------------------------------------------------------------
// Include Generic kernel
// ----------------------------------------------------------------------------
#include <CalcServer/Kernel.h>

namespace Aqua{ namespace CalcServer{ namespace Boundary{

/** @class ElasticBounce ElasticBounce.h CalcServer/Boundary/ElasticBounce.h
 * @brief Simplest boundary condition.
 * This boundary condition only test if a particle
 * will pass trought a wall, and reflect it as elastic
 * bounce, affected by a factor that set the amount of
 * energy lost in the process. \n
 * This boundary condition is ussually used clompementary
 * to others boundary conditions because guarantee that
 * particles can't pass trought walls.
 */
class ElasticBounce : public Aqua::CalcServer::Kernel
{
public:
	/** Constructor.
	 */
	ElasticBounce();

	/** Destructor.
	 */
	~ElasticBounce();

	/** Boundary effect computation.
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
	/// OpenCL vertices set kernel
	cl_kernel clKernel;
	/// Global work size (calculated with clLocalWorkSize).
	size_t clGlobalWorkSize;
	/// Local work size (default value = 256)
	size_t clLocalWorkSize;
};

}}}  // namespace

#endif // ELASTICBOUNCE_H_INCLUDED
