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

/** @class Save Save.h CalcServer/Save.h
 * @brief Calls for a full files saving
 *
 * Useful mainly for debugging
 */
class Save : public Aqua::CalcServer::Tool
{
  public:
	/** @brief Constructor
	 * @param tool_name Tool name
	 */
	Save(const std::string tool_name);

	/** @brief Destructor
	 */
	~Save();

	/** @brief Initialize the tool.
	 */
	void setup();

  protected:
	/** Execute the tool
	 * @param events List of events that shall be waited before safe execution
	 * @return OpenCL event to be waited before accessing the dependencies
	 */
	cl_event _execute(const std::vector<cl_event> events);
};

}
}
} // namespace
