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

#ifndef RATES_H_INCLUDED
#define RATES_H_INCLUDED

// ----------------------------------------------------------------------------
// Include Generic kernel
// ----------------------------------------------------------------------------
#include <CalcServer/Kernel.h>

namespace Aqua{ namespace CalcServer{

/** @class Rates Rates.h CalcServer/Rates.h
 * @brief With the Link-List data the interaction between
 * particles can be performed. As a result, density rate
 * variation, forces, and some time step data will computed.
 */
class Rates : public Aqua::CalcServer::Kernel
{
public:
	/** Constructor.
	 */
	Rates();

	/** Destructor.
	 */
	~Rates();

	/** Rates of variation calculation.
	 * @return false if all gone right. \n true otherwise.
	 */
	bool execute();

private:
	/** Setup OpenCL kernel
	 * @return false if all gone right. \n true otherwise.
	 */
	bool setupOpenCL();
	/** Auxiliar method to alloc memory for link-list (when needed).
	 * @return false if all gone right. \n true otherwise.
	 */
	bool allocLinkList();

	/// OpenCL script path
	char *_path;

	/// OpenCL program
	cl_program program;
	/// OpenCL icell kernel
	cl_kernel kernel;
	/// OpenCL ihoc init kernel
	cl_kernel clSortKernel;
	/// Global work size.
	size_t _global_work_size;
	/// Local work size
	size_t _local_work_size;
	/// true if \f$delta\f$-SPH (cont. eq. diffusive term) must be applied.
	bool isDelta;
	/// true if local memory can be used on kernel.
	bool isLocalMemory;
};

}}  // namespace

#endif // RATES_H_INCLUDED
