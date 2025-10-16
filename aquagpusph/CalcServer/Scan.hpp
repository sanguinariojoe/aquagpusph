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
 * @brief Prefix sum
 * (See Aqua::CalcServer::Scan for details)
 * @note Hardcoded versions of the files CalcServer/Scan.cl.in and
 * CalcServer/Scan.hcl.in are internally included as a text array.
 */

#ifndef SCAN_H_INCLUDED
#define SCAN_H_INCLUDED

#include <vector>
#include "CalcServer.hpp"
#include "Kernel.hpp"

/** @def _WGS Work group size
 * @note Must be power of 2, and in some devices greather than 32.
 */
#ifndef _WGS
#define _WGS 64U
#endif

namespace Aqua {
namespace CalcServer {

/** @class Scan Scan.h CalcServer/Scan.h
 * @brief A prefix sum, also known as cumsum
 * @see Scan.cl
 * @note Hardcoded versions of the files CalcServer/Scan.cl.in and
 * CalcServer/Scan.hcl.in are internally included as a text array.
 * @warning This tool will not work with ihoc array, since it might be
 * reallocated
 * @todo The Scan algorithm allows a segmented prefix sum, however such option
 * is not still available for the user.
 */
class Scan : public Aqua::CalcServer::Tool
{
  public:
	/** @brief Scan definition.
	 * @param name Tool name.
	 * @param var_name Variable to be scanned.
	 * @param once Run this tool just once. Useful to make initializations.
	 */
	Scan(const std::string name,
	     const std::string var_name,
	     bool once = false);

	/// Destructor.
	~Scan();

	/** @brief Initialize the tool.
	 *
	 * This method should be called after the constructor, such that it could
	 * report errors that the application may handle quitting in a safe way.
	 */
	void setup();

  protected:
	/** Execute the tool
	 * @param events List of events that shall be waited before safe execution
	 * @return OpenCL event to be waited before accesing the dependencies
	 */
	cl_event _execute(const std::vector<cl_event> events);

  private:
	/** @brief Extract the input and output variables from the provided data in
	 * Scan().
	 * @see Aqua::InputOutput::Variables
	 */
	void variables();

	/** @brief Setup the OpenCL stuff
	 */
	void setupOpenCL();

	/** @brief Carry out a recursive step of the scan
	 * @param stage The reduction step
	 * @param events List of events to wait for
	 * @return The event marking when this stage is already finished
	 */
	cl_event scan(const unsigned int stage,
	              const std::vector<cl_event> events);

	/** @brief Create the compilation flags
	 *
	 * The compilation flags depends on the intended work group size
	 *
	 * @param local_size Work group size
	 * @return Flags string
	 */
	const std::string flags();

	/// Variable name
	std::string _var_name;

	/// Variable
	InputOutput::ArrayVariable* _var;

	/// Length of the input array
	size_t _n;

	/// Number of work groups at each substage
	std::vector<size_t> _k;
	/// Inputs at each stage
	std::vector<cl_mem> _input;
	/// Partials at each stage (segmented scan is not available yet)
	std::vector<cl_mem> _partial;
	/// Flags at each stage (segmented scan is not available yet)
	std::vector<cl_mem> _flags;

	/// Init kernel
	cl_kernel _init;
	/// Scan kernel
	cl_kernel _scan;
	/// Sweep-up kernel
	cl_kernel _sweepup;
	/// Sweep-down kernel
	cl_kernel _sweepdown;
};

}
} // namespace

#endif // SCAN_H_INCLUDED
