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

#ifndef SENSORS_H_INCLUDED
#define SENSORS_H_INCLUDED

// ----------------------------------------------------------------------------
// Include Generic kernel
// ----------------------------------------------------------------------------
#include <CalcServer/Kernel.h>

namespace Aqua{ namespace CalcServer{

/** @class Sensors Sensors.h CalcServer/Sensors.h
 * @brief Perform sensors calculation.
 * @todo Fix sensors technique
 */
class Sensors : public Aqua::CalcServer::Kernel
{
public:
	/** Constructor.
	 */
	Sensors();

	/** Destructor.
	 */
	~Sensors();

	/** Sensors calculation.
	 * @return false if all gone right. \n true otherwise.
	 */
	bool execute();

protected:
	/** Retrieve data form server, and print output.
	 * @return false if all gone right. \n true otherwise.
	 */
	bool printOutput();

private:
	/** Setup OpenCL kernel
	 * @return false if all gone right. \n true otherwise.
	 */
	bool setupOpenCL();

	/** Start output
	 * @return false if all gone right. \n true otherwise.
	 */
	bool initOutput();

	/// Number of sensors
	int n;

	/// OpenCL script path
	char* _path;

	/// Output file
	FILE* Output;
	/// Last time when a file was printed
	float OutputTime;

	/// Host storage for positions
	vec *hPos;
	/// Host storage for pressures
	cl_float *hPress;
	/// Host storage for density
	cl_float *hDens;
	/// Host storage for Shepard term
	cl_float *hSumW;
	/// Host storage for Shepard term gradient
	vec *hGradW;

	/// OpenCL program
	cl_program _program;
	/// OpenCL kernel
	cl_kernel _kernel;
	/// Global work size
	size_t _global_work_size;
	/// Local work size
	size_t _local_work_size;
};

}}  // namespace

#endif // SENSORS_H_INCLUDED
