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
 * @brief Runtime array dumping to an output file.
 * (See Aqua::CalcServer::Reports::Save for details)
 */

#include <sstream>
#include <fstream>
#include <iomanip>
#include "aquagpusph/AuxiliarMethods.hpp"
#include "aquagpusph/InputOutput/Logger.hpp"
#include "aquagpusph/CalcServer/CalcServer.hpp"
#include "aquagpusph/FileManager.hpp"
#include "Save.hpp"

namespace Aqua {
namespace CalcServer {
namespace Reports {

Save::Save(const std::string tool_name)
  : Tool(tool_name)
{
}

Save::~Save()
{
}

void
Save::setup()
{
	std::ostringstream msg;
	msg << "Loading the save report \"" << name() << "\"..." << std::endl;
	LOG(L_INFO, msg.str());

	Tool::setup();
}

cl_event
Save::_execute(const std::vector<cl_event> UNUSED_PARAM events)
{
	InputOutput::FileManager *file_manager =
		InputOutput::FileManager::singleton();
	file_manager->save(0.f);
	file_manager->waitForSavers();

	return NULL;
}

}
}
} // namespace
