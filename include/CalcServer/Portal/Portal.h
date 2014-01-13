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

#ifndef PORTAL_H_INCLUDED
#define PORTAL_H_INCLUDED

// ----------------------------------------------------------------------------
// Include the Problem setup manager header
// ----------------------------------------------------------------------------
#include <ProblemSetup.h>

// ----------------------------------------------------------------------------
// Include Generic kernel
// ----------------------------------------------------------------------------
#include <CalcServer/Kernel.h>

namespace Aqua{ namespace CalcServer{ namespace Portal{

/** @class Portal Portal.h CalcServer/Portal/Portal.h
 * @brief Base portal class, that only teleport particles that
 * pass throught outlet portal to the inlet portal.
 */
class Portal : public Aqua::CalcServer::Kernel
{
public:
	/** Constructor.
	 */
	Portal(InputOutput::ProblemSetup::sphPortal *portal);

	/** Destructor.
	 */
	~Portal();

	/** Teleport particles from outlet to inlet.
	 * @return false if all gone right, true otherwise.
	 */
	bool execute();

private:
	/** Setup the OpenCL stuff
	 * @return false if all gone right, true otherwise.
	 */
	bool setupOpenCL();

	/// OpenCL script path
	char* _path;
	/// Main portal
	InputOutput::ProblemSetup::sphPortal *mPortal;

	/// OpenCL program
	cl_program _program;
	/// OpenCL vertices set kernel
	cl_kernel _kernel;
	/// Global work size (calculated with local_work_size).
	size_t _global_work_size;
	/// Local work size (default value = 256)
	size_t _local_work_size;
};

}}}  // namespace

#endif // PORTAL_H_INCLUDED
