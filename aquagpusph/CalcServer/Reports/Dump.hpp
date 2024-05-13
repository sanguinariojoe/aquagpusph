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
 * @brief Runtime particles set output file.
 * (See Aqua::CalcServer::Reports::SetSetTabFile for details)
 */

#pragma once

#include "aquagpusph/sphPrerequisites.hpp"
#include "Report.hpp"

namespace Aqua {
namespace CalcServer {
namespace Reports {

/** @class Dump Dump.h CalcServer/Dump.h
 * @brief Dump an array at any time.
 *
 * Useful mainly for debugging. Several arrays can be dumped simultaneously,
 * provided that all of them have the same length
 */
class Dump : public Aqua::CalcServer::Reports::Report
{
  public:
	/** @brief Constructor
	 * @param tool_name Tool name
	 * @param fields Fields to be printed
	 * @param output_file File to be written. Several scape strings can be used,
	 * as described in Aqua::newFilePath()
	 * @remarks The output file will be cleared
	 */
	Dump(const std::string tool_name,
	     const std::string fields,
	     const std::string output_file,
	     bool binary=false);

	/** @brief Destructor
	 */
	~Dump();

	/** @brief Initialize the tool.
	 */
	void setup();

	/** Print the data on the tabulated file
	 * @note This function is public to let the OpenCL callback call it, but it
	 * is not meant to be used by regular users
	 */
	void print();

  protected:
	/** Execute the tool
	 * @param events List of events that shall be waited before safe execution
	 * @return OpenCL event to be waited before accessing the dependencies
	 */
	cl_event _execute(const std::vector<cl_event> events);

	/** Download the data from the device, and store it.
	 * @param vars Fields to download.
	 * @return The events to be waited before the data is available.
	 * @note The returned data must be manually freed.
	 */
	std::vector<cl_event> download(std::vector<InputOutput::Variable*> vars);

  private:
	/// Output file name
	std::string _output_file;

	/// File opening mode
	bool _binary;

	/// List of host allocated memory objects
	std::vector<void*> _data;

	/// Length of the arrays to be saved
	size_t _n;
};

}
}
} // namespace
