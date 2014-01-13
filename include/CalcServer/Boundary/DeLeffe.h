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

#ifndef DELEFFE_H_INCLUDED
#define DELEFFE_H_INCLUDED

#include <CalcServer/Kernel.h>

namespace Aqua{ namespace CalcServer{ namespace Boundary{

/** @class DeLeffe DeLeffe.h CalcServer/Boundary/DeLeffe.h
 * @brief Boundary integrals based boundary condition.
 */
class DeLeffe : public Aqua::CalcServer::Kernel
{
public:
	/** Constructor.
	 */
	DeLeffe();

	/** Destructor.
	 */
	~DeLeffe();

	/** Compute the boundary condition.
	 * @return false if all gone right, true otherwise.
	 */
	bool execute();

private:
	/** Setup the OpenCL stuff
	 * @return false if all gone right, true otherwise.
	 */
	bool setupOpenCL();

	/** Set the boundary elements properties.
	 * @return false if all gone right, true otherwise.
	 */
	bool elements();

	/** Compute the boundary condition.
	 * @return false if all gone right, true otherwise.
	 */
	bool boundary();

	/// OpenCL script path
	char* _path;

	/// OpenCL program
	cl_program _program;
	/// OpenCL boundary elements seter kernel
	cl_kernel _setup_kernel;
	/// OpenCL boundary condition kernel
	cl_kernel _boundary_kernel;
	/// Global work size (calculated with local_work_size).
	size_t _global_work_size;
	/// Local work size (default value = 256)
	size_t _local_work_size;
	/// true if local memory can be used on kernel.
	bool _use_local_mem;
};

}}}  // namespace

#endif // DELEFFE_H_INCLUDED
