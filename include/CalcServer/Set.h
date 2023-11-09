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
 * @brief Set all the components of an array with the desired value.
 * (See Aqua::CalcServer::Set for details)
 * @note Hardcoded versions of the files CalcServer/Set.cl.in and
 * CalcServer/Set.hcl.in are internally included as a text array.
 */

#ifndef SET_H_INCLUDED
#define SET_H_INCLUDED

#include <CalcServer.h>
#include <CalcServer/SetScalar.h>
#include <CalcServer/Kernel.h>

namespace Aqua {
namespace CalcServer {

/** @class Set Set.h CalcServer/Set.h
 * @brief Set all the components of an array with the desired value.
 */
class Set final : public Aqua::CalcServer::ScalarExpression
{
  public:
	/** Constructor.
	 * @param name Tool name.
	 * @param var_name Variable to set.
	 * @param value Value to set.
	 * @param once Run this tool just once. Useful to make initializations.
	 * @note Some helpers are available for value:
	 *   - VEC_ZERO: Zeroes vector.
	 *   - VEC_ONE: Ones vector, in 3D cases the last component will be zero.
	 *   - VEC_ALL_ONE: Equal to VEC_ONE, but in 3D cases the last component
	 * will be one as well.
	 *   - VEC_INFINITY: INFINITY components vector, in 3D cases the last
	 * component will be zero.
	 *   - VEC_ALL_INFINITY: Equal to VEC_INFINITY, but in 3D cases the last
	 * component will be INFINITY as well.
	 *   - VEC_NEG_INFINITY: -VEC_INFINITY
	 *   - VEC_ALL_NEG_INFINITY: -VEC_ALL_INFINITY.
	 */
	Set(const std::string name,
	    const std::string var_name,
	    const std::string value,
	    bool once = false);

	/** Destructor.
	 */
	~Set();

	/** Initialize the tool.
	 */
	void setup();

  protected:
	/** @brief Evaluate the expression if possible
	 */
	void _solve();

	/** Execute the tool
	 * @param events List of events that shall be waited before safe execution
	 * @return OpenCL event to be waited before accessing the dependencies
	 */
	cl_event _execute(const std::vector<cl_event> events);

  private:
	/** Get the input variable
	 */
	void variable();

	/** Setup the OpenCL stuff
	 */
	void setupOpenCL();

	/** Update the input looking for changed value.
	 * @return false if all gone right, true otherwise.
	 */
	void setVariables();

	/// Input variable name
	std::string _var_name;
	/// Value to set
	std::string _value;

	/// Input variable
	InputOutput::ArrayVariable* _var;

	/// Memory object sent
	cl_mem _input;

	/// OpenCL kernel
	cl_kernel _kernel;

	/// Global work sizes in each step
	size_t _global_work_size;
	/// Local work sizes in each step
	size_t _work_group_size;
	/// Number of elements
	unsigned int _n;

	/// Memory storage for the value
	void* _data;

	/// Kernel enqueuer
	KernelEnqueuer* _kernel_enqueuer;
};

}
} // namespace

#endif // SET_H_INCLUDED
