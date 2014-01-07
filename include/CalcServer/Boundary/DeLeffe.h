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

// ----------------------------------------------------------------------------
// Include Generic kernel
// ----------------------------------------------------------------------------
#include <CalcServer/Kernel.h>

namespace Aqua{ namespace CalcServer{ namespace Boundary{

/** @class DeLeffe DeLeffe.h CalcServer/Boundary/DeLeffe.h
 * @brief Boundary condition that not consider particles
 * at the boundary. \n
 * Boundary condition requires several OpenCL kernels. \n
 * 1st the effect of walls over particles is computed. \n
 * 2nd the shepard term is applied. \n
 * 3rd the boundary vertices properties are set.
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

	/** DeLeffe boundary computation.
	 * @return false if all gone right. \n true otherwise.
	 */
	bool execute();

private:
	/** Setup OpenCL kernel
	 * @return false if all gone right. \n true otherwise.
	 */
	bool setupOpenCL();

	/** Set vertices.
	 * @return false if all gone right. \n true otherwise.
	 */
	bool vertices();

	/** Perform boundary effect.
	 * @return false if all gone right. \n true otherwise.
	 */
	bool boundary();

	/// OpenCL script path
	char* _path;

	/// OpenCL program
	cl_program _program;
	/// OpenCL vertices set kernel
	cl_kernel clVerticesKernel;
	/// OpenCL boundary effect kernel
	cl_kernel clBoundaryKernel;
	/// Global work size (calculated with local_work_size).
	size_t _global_work_size;
	/// Local work size (default value = 256)
	size_t _local_work_size;
	/// true if local memory can be used on kernel.
	bool _use_local_mem;
};

}}}  // namespace

#endif // DELEFFE_H_INCLUDED
