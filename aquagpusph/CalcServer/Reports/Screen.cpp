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
 * @brief On screen runtime output.
 * (See Aqua::CalcServer::Reports::Screen for details)
 */

#include <sstream>
#include <algorithm>
#include "aquagpusph/InputOutput/Logger.hpp"
#include "aquagpusph/CalcServer/CalcServer.hpp"
#include "Screen.hpp"

namespace Aqua {
namespace CalcServer {
namespace Reports {

Screen::Screen(const std::string tool_name,
               const std::string fields,
               const std::string color,
               bool bold)
  : Report(tool_name, fields)
  , _color(color)
  , _bold(bold)
{
}

Screen::~Screen() {}

void
Screen::setup()
{
	unsigned int i;

	std::ostringstream msg;
	msg << "Loading the report \"" << name() << "\"..." << std::endl;
	LOG(L_INFO, msg.str());

	// Set the color in lowercase
	std::transform(_color.begin(), _color.end(), _color.begin(), ::tolower);

	Report::setup();
}

/** @brief Callback called when all the variables can be read by
 * Aqua::CalcServer::Screen.
 *
 * This function is just redirecting the work to
 * InputOutput::Logger::writeReport() using the
 * Aqua::CalcServer::Report::data() function to extract the information
 * @param event The triggering event
 * @param event_command_status CL_COMPLETE upon all dependencies successfully
 * fulfilled. A negative integer if one or mor dependencies failed.
 * @param user_data A casted pointer to the Aqua::CalcServer::Screen
 * tool (or the inherited one)
 */
void CL_CALLBACK
screenreport_cb(cl_event event, cl_int event_command_status, void* user_data)
{
	clReleaseEvent(event);
	auto tool = (Screen*)user_data;
	if (event_command_status != CL_COMPLETE) {
		std::stringstream msg;
		msg << "Skipping \"" << tool->name() << "\" due to dependency errors."
		    << std::endl;
		LOG(L_WARNING, msg.str());
		clSetUserEventStatus(tool->getUserEvent(), event_command_status);
		clReleaseEvent(tool->getUserEvent());
		return;
	}

	InputOutput::Logger::singleton()->writeReport(
	    tool->data(), tool->color(), tool->bold());
	cl_int err_code;
	err_code = clSetUserEventStatus(tool->getUserEvent(), CL_COMPLETE);
	CHECK_OCL_OR_THROW(err_code,
	                   std::string("Failure setting as complete the tool \"") +
	                       tool->name() + "\".");
	err_code = clReleaseEvent(tool->getUserEvent());
	CHECK_OCL_OR_THROW(
	    err_code,
	    std::string("Failure releasing the user event in tool \"") +
	        tool->name() + "\".");
}

cl_event
Screen::_execute(const std::vector<cl_event> events)
{
	std::vector<cl_event> no_events;
	return setCallback(no_events, screenreport_cb);
}

}
}
} // namespace
