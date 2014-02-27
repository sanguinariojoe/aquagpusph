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

#ifndef GHOSTPARTICLES_H_INCLUDED
#define GHOSTPARTICLES_H_INCLUDED

// ----------------------------------------------------------------------------
// Include Generic kernel
// ----------------------------------------------------------------------------
#include <CalcServer/Kernel.h>

// ----------------------------------------------------------------------------
// Include standar libraries
// ----------------------------------------------------------------------------
#include <vector>
#include <deque>

namespace Aqua{ namespace CalcServer{ namespace Boundary{

/** @class GhostParticles GhostParticles.h CalcServer/Boundary/GhostParticles.h
 * @brief Boundary condition based on extending the fluid with a mirroring
 * process of the particles with respect to the wall, conveniently setting
 * modified fields. In order to improve the consistency the model used for the
 * pressure, normal velocity, and tangent velocity can be modified, but the
 * following model is set by default (Free slip):
 *   - pressModel = "SSM"
 *   - normalUModel = "ASM"
 *   - tangentUModel = "SSM"
 */
class GhostParticles : public Aqua::CalcServer::Kernel
{
public:
	/** Constructor.
	 */
	GhostParticles();

	/** Destructor.
	 */
	~GhostParticles();

	/** Compute the boundary condition.
	 * @return false if all gone right, true otherwise.
	 */
	bool execute();

private:
	/** Setup the OpenCL stuff
	 * @return false if all gone right, true otherwise.
	 */
	bool setupOpenCL();

	/** Create walls OpenCL instances.
	 * @return false if all gone right, true otherwise.
	 */
	bool createWalls();

	/// OpenCL script path
	char* _path;

	/// OpenCL program
	cl_program _program;
	/// OpenCL vertices set kernel
	cl_kernel _kernel;
	/// Global work size (calculated with local_work_size).
	size_t _global_work_size;
	/// Local work size (default value = 256)
	size_t _local_work_size;
	/// Array of walls with all parameters
	std::deque<cl_mem> _walls;
	/// true if \f$delta\f$-SPH (cont. eq. diffusive term) must be applied.
	bool _is_delta;
};

}}}  // namespace

#endif // GHOSTPARTICLES_H_INCLUDED
