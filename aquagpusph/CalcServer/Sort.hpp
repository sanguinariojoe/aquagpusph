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
 * @brief Methods to perform a radix sort using the GPU (or any device
 * supported by OpenCL).
 * (See Aqua::CalcServer::Sort for details)
 * @note Hardcoded versions of the files CalcServer/Sort.cl.in and
 * CalcServer/Sort.hcl.in are internally included as a text array.
 */

#ifndef SORT_H_INCLUDED
#define SORT_H_INCLUDED

#include "aquagpusph/sphPrerequisites.hpp"
#include "CalcServer.hpp"
#include "Tool.hpp"
#include "Kernel.hpp"

namespace Aqua {
namespace CalcServer {

/** @class Sort Sort.h CalcServer/Sort.h
 * @brief Methods to perform a vector sorting using the GPU (or any device
 * supported by OpenCL).
 *
 * The sorting is carried out considering a bitonic sorter:
 * https://en.wikipedia.org/wiki/Bitonic_sorter
 * @note Hardcoded versions of the files CalcServer/Sort.cl.in and
 * CalcServer/Sort.hcl.in are internally included as a text array.
 */
class Sort : public Aqua::CalcServer::Tool
{
  public:
	/** Constructor.
	 * @param tool_name Tool name.
	 * @param variable Variable to sort.
	 * @param permutations Variable where the permutations will be stored.
	 * @param inv_permutations Variable where the inverse permutations will be
	 * stored.
	 * @param once Run this tool just once. Useful to make initializations.
	 */
	Sort(const std::string tool_name,
	     const std::string variable = "icell",
	     const std::string permutations = "id_unsorted",
	     const std::string inv_permutations = "id_sorted",
	     bool once = false);

	/** Destructor
	 */
	~Sort();

	/** Initialize the tool.
	 */
	void setup();

  protected:
	/** Execute the tool
	 * @param events List of events that shall be waited before safe execution
	 * @return OpenCL event to be waited before accesing the dependencies
	 */
	cl_event _execute(const std::vector<cl_event> events);

  private:
	/** Get the variables to compute.
	 */
	void variables();

	/** Setup the OpenCL stuff
	 */
	void setupOpenCL();

	/** Compile the kernels
	 *
	 * This function is useful to give several shots with different values of
	 * Aqua::CalcServer::Sort::_local_work_size
	 * @return The list of kernels
	 */
	std::vector<cl_kernel> compileOpenCL();

	/** Setup the memory objects.
	 */
	void setupMems();

	/** Send the fixed arguments to the kernels.
	 */
	void setupArgs();

	/// Name of the variable to sort
	std::string _var_name;

	/// Permutations array name
	std::string _perms_name;

	/// Inverse permutations array name
	std::string _inv_perms_name;

	/// Variable to sort
	InputOutput::ArrayVariable* _var;

	/// Permutations array
	InputOutput::ArrayVariable* _perms;

	/// Inverse permutations array
	InputOutput::ArrayVariable* _inv_perms;

	/// Number of keys to sort
	size_t _n;

	/// Padded number of keys to sort (power of 2)
	size_t _n_padded;

	/// Type of the variable
	std::string _var_type;

	/// Maximum value the variable can take (according to the type)
	std::string _var_max;

	/// Permutations initialization kernel
	cl_kernel _init_kernel;
	/// Bitonic merge start
	cl_kernel _start_kernel;
	/// Bitonic local merge
	cl_kernel _local_kernel;
	/// Bitonic global merge
	cl_kernel _global_kernel;
	/// Bitonic global merge
	cl_kernel _inv_perms_kernel;

	/// Values, already padded
	cl_mem _vals;
	/// Permutations, already padded
	cl_mem _permut;

	/** @brief Maximum local work size to run the task
	 *
	 * It is limited by the maximum number of threads reported by
	 * clGetKernelWorkGroupInfo() as well as the local memory available by the
	 * device.
	 */
	size_t _local_work_size;
	/// Global work size (assuming the maximum local work size) to compute _n
	/// threads.
	size_t _global_work_size;
};

}
} // namespace

#endif // SORT_H_INCLUDED
