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

#ifndef TIMESTEP_H_INCLUDED
#define TIMESTEP_H_INCLUDED

// ----------------------------------------------------------------------------
// Include Generic kernel
// ----------------------------------------------------------------------------
#include <CalcServer/Kernel.h>

// ----------------------------------------------------------------------------
// Include Reduction tool
// ----------------------------------------------------------------------------
#include <CalcServer/Reduction.h>

namespace Aqua{ namespace CalcServer{

/** @class TimeStep TimeStep.h CalcServer/TimeStep.h
 * @brief Time step computation. Courant condition is
 * stablished forcing maximum distance that each particle
 * can be translated in a time step to 0.1 <i>h</i>.
 */
class TimeStep : public Aqua::CalcServer::Kernel
{
public:
	/** Constructor.
	 */
	TimeStep();

	/** Destructor.
	 */
	~TimeStep();

	/** Time step computation.
	 * @return false if all gone right. \n true otherwise.
	 */
	bool execute();

private:
	/** Setup OpenCL kernel
	 * @return 0 if any error happens. \n
	 * Error code otherwise.
	 */
	bool setupOpenCL();

	/// OpenCL script path
	char* _path;

	/// OpenCL program
	cl_program program;
	/// OpenCL kernel
	cl_kernel kernel;
	/// Global work size (calculated with local_work_size).
	size_t _global_work_size;
	/// Local work size (default value = 256)
	size_t _local_work_size;

	/// Convective time step reduction tool
	Reduction *reduction;

	/// Main time step
	float MainDt;

	/// Time step clamped flag (in order to avoid annoying clamping error)
	int dtClamp;
};

}}  // namespace

#endif // TIMESTEP_H_INCLUDED
