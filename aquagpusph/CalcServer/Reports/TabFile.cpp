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
 * @brief Runtime output file.
 * (See Aqua::CalcServer::Reports::TabFile for details)
 */

#include "aquagpusph/AuxiliarMethods.hpp"
#include "aquagpusph/InputOutput/Logger.hpp"
#include "TabFile.hpp"

namespace Aqua {
namespace CalcServer {
namespace Reports {

TabFile::TabFile(const std::string tool_name,
                 const std::string fields,
                 const std::string output_file)
  : Report(tool_name, fields)
  , _output_file("")
{
	try {
		unsigned int i = 0;
		_output_file = newFilePath(output_file, i, 1);
	} catch (std::invalid_argument e) {
		std::ostringstream msg;
		_output_file = setStrConstantsCopy(output_file);
		msg << "Overwriting '" << _output_file << "'" << std::endl;
		LOG(L_WARNING, msg.str());
	}
}

TabFile::~TabFile()
{
	if (_f.is_open())
		_f.close();
}

void
TabFile::setup()
{
	std::ostringstream msg;
	msg << "Loading the report \"" << name() << "\"..." << std::endl;
	LOG(L_INFO, msg.str());

	_f.open(_output_file.c_str(), std::ios::out);

	Report::setup();

	// Write the header
	_f << "# ";
	std::vector<InputOutput::Variable*> vars = variables();
	for (auto var : vars) {
		_f << var->name() << " ";
	}
	_f << std::endl;
	_f.flush();
}

void
TabFile::print()
{
	std::string out = replaceAllCopy(data(false, false), "\n", " ");
	_f << out << std::endl;
	_f.flush();

	cl_int err_code;
	err_code = clSetUserEventStatus(getUserEvent(), CL_COMPLETE);
	if (err_code != CL_SUCCESS) {
		std::stringstream msg;
		msg << "Failure setting as complete the tool \"" << name() << "\"."
		    << std::endl;
		LOG(L_ERROR, msg.str());
		InputOutput::Logger::singleton()->printOpenCLError(err_code);
		return;
	}
	err_code = clReleaseEvent(getUserEvent());
	if (err_code != CL_SUCCESS) {
		std::stringstream msg;
		msg << "Failure releasing the user event at tool \"" << name() << "\"."
		    << std::endl;
		LOG(L_ERROR, msg.str());
		InputOutput::Logger::singleton()->printOpenCLError(err_code);
		return;
	}
}

/** @brief Callback called when all the variables can be read by
 * Aqua::CalcServer::TabFile.
 *
 * This function is just redirecting the work to
 * Aqua::CalcServer::TabFile::print()
 * @param event The triggering event
 * @param event_command_status CL_COMPLETE upon all dependencies successfully
 * fulfilled. A negative integer if one or mor dependencies failed.
 * @param user_data A casted pointer to the Aqua::CalcServer::Screen
 * tool (or the inherited one)
 */
void CL_CALLBACK
tabfilereport_cb(cl_event event, cl_int event_command_status, void* user_data)
{
	clReleaseEvent(event);
	auto tool = (TabFile*)user_data;
	if (event_command_status != CL_COMPLETE) {
		std::stringstream msg;
		msg << "Skipping \"" << tool->name() << "\" due to dependency errors."
		    << std::endl;
		LOG(L_WARNING, msg.str());
		clSetUserEventStatus(tool->getUserEvent(), event_command_status);
		clReleaseEvent(tool->getUserEvent());
		return;
	}

	tool->print();
}

cl_event
TabFile::_execute(const std::vector<cl_event> events)
{
	return setCallback(events, tabfilereport_cb);
}

}
}
} // namespace
