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
 * @brief Base class for all the report file managers.
 * (See Aqua::InputOutput::Report for details)
 */

#include <sstream>
#include <iostream>
#include "Logger.hpp"
#include "Report.hpp"
#include "aquagpusph/AuxiliarMethods.hpp"

namespace Aqua {
namespace InputOutput {

Report::Report() {}

Report::~Report() {}

void
Report::file(std::string filename)
{
	_output_file = filename;
}

void
Report::file(std::string basename, unsigned int start_index)
{
	try {
		file(newFilePath(basename, start_index));
	} catch (std::invalid_argument& e) {
		file(setStrConstantsCopy(basename));
		std::ostringstream msg;
		msg << "Overwriting file '" << file() << "'" << std::endl;
		LOG(L_WARNING, msg.str());
	}
}

}
} // namespace
