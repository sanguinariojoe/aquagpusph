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
 * @brief UnSort Recover the original id of each particle.
 * (See Aqua::CalcServer::UnSort for details)
 * @note Hardcoded versions of the files CalcServer/UnSort.cl.in and
 * CalcServer/UnSort.hcl.in are internally included as a text array.
 */

#ifndef UNSORT_H_INCLUDED
#define UNSORT_H_INCLUDED

#include <CalcServer.h>
#include <CalcServer/Kernel.h>

namespace Aqua {
namespace CalcServer {

/** @class UnSort UnSort.h CalcServer/UnSort.h
 * @brief Recover the original id of each particle.
 *
 * This tool is not designed for the common usage but as an auxiliar tool for
 * the savers or the MPI syncing, therefore it will not be selectable for the
 * users.
 */
class UnSort : public Aqua::CalcServer::Tool
{
  public:
	/** Constructor
	 * @param name Tool name
	 * @param var_name Variable to unsort
	 * @param permutations_name Permutations to be considered
	 * @param once Run this tool just once. Useful to make initializations
	 */
	UnSort(const std::string name,
	       const std::string var_name,
	       const std::string permutations_name = "id",
	       bool once = false);

	/** Destructor.
	 */
	~UnSort();

	/** Initialize the tool.
	 */
	void setup();

	/** Get the input variable.
	 * @return The variable to become sorted.
	 */
	InputOutput::ArrayVariable* input() { return _var; }

	/** Get the memory object where the unsorted data is stored.
	 * @return The memory object where the unsorted data is stored.
	 */
	cl_mem output() { return _output; }

  protected:
	/** Execute the tool
	 * @param events List of events that shall be waited before safe execution
	 * @return OpenCL event to be waited before accessing the dependencies
	 */
	cl_event _execute(const std::vector<cl_event> events);

  private:
	/** Get the input variable
	 */
	void variables();

	/** Create the output memory object
	 */
	void setupMem();

	/** Setup the OpenCL stuff
	 */
	void setupOpenCL();

	/// Input variable name
	std::string _var_name;

	/// Input permutations name
	std::string _perms_name;

	/// Input variable
	InputOutput::ArrayVariable* _var;

	/// ID variable
	InputOutput::ArrayVariable* _id_var;

	/// ID Memory object sent
	cl_mem _id_input;

	/// Memory object sent
	cl_mem _input;

	/// Output memory object
	cl_mem _output;

	/// OpenCL kernel
	cl_kernel _kernel;

	/// Global work sizes in each step
	size_t _global_work_size;
	/// Local work sizes in each step
	size_t _work_group_size;
	/// Number of elements
	unsigned int _n;

	/// Arguments setter for the kernel
	ArgSetter* _args_setter;
};

}
} // namespace

#endif // UNSORT_H_INCLUDED
