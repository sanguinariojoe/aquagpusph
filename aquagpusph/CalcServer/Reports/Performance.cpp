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
 * @brief Performance report.
 * (See Aqua::CalcServer::Reports::Performance for details)
 */

#include <algorithm>
#include <iomanip>
#include "aquagpusph/AuxiliarMethods.hpp"
#include "aquagpusph/InputOutput/Logger.hpp"
#include "aquagpusph/CalcServer/CalcServer.hpp"
#include "Performance.hpp"
#include "aquagpusph/ext/json.hpp"

using json = nlohmann::json;

namespace Aqua {
namespace CalcServer {
namespace Reports {

Performance::Performance(const std::string tool_name,
                         const std::string color,
                         bool bold,
                         const std::string output_file)
  : Report(tool_name, "dummy_fields_string")
  , _color(color)
  , _bold(bold)
  , _output_file(output_file)
  , _timer(0)
  , _user_event(NULL)
{
	if (output_file != "") {
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
}

Performance::~Performance() {}

void
Performance::setup()
{
	unsigned int i;
	auto C = CalcServer::singleton();

	std::ostringstream msg;
	msg << "Loading the report \"" << name() << "\"..." << std::endl;
	LOG(L_INFO, msg.str());

	// Set the color in lowercase
	std::transform(_color.begin(), _color.end(), _color.begin(), ::tolower);

	Tool::setup();

	setInputDependencies(
		{ "t", "end_t", "iter", "end_iter", "frame", "end_frame" });
}

void
Performance::print()
{
	std::stringstream data;
	auto C = CalcServer::singleton();

	size_t allocated_MB = computeAllocatedMemory() / (1024 * 1024);
	data << "Performance:" << std::endl
	     << "Memory=" << std::setw(18) << allocated_MB << "MB" << std::endl;

	// Add the tools time elapsed
	auto snapshots = ((ProfilingInfo*)C)->get();
	auto snapshot = snapshots.front();
	cl_long tools_elapsed = 0;
	for (auto sample : snapshot.samples) {
		tools_elapsed += ProfilingInfo::delta(sample.end, sample.start);
	}
	auto timer = CalcServer::host_timer();
	auto elapsed = (!_timer) ? tools_elapsed : ProfilingInfo::delta(
		timer, _timer);
	_timer = timer;
	auto overhead = elapsed - tools_elapsed;
	data << "Elapsed=" << std::setw(17) << elapsed * 1.e-9 << "s"
	     << std::endl;
	data << "Overhead=" << std::setw(17) << overhead * 1.e-9 << "s"
	     << std::endl;

	// Compute the progress
	InputOutput::Variables* vars = C->variables();
	float progress = 0.f;
	float t = *(float*)vars->get("t")->get_async();
	float end_t = *(float*)vars->get("end_t")->get_async();
	progress = max(progress, t / end_t);
	unsigned int iter = *(unsigned int*)vars->get("iter")->get_async();
	unsigned int end_iter = *(unsigned int*)vars->get("end_iter")->get_async();
	progress = max(progress, (float)iter / end_iter);
	unsigned int frame = *(unsigned int*)vars->get("frame")->get_async();
	unsigned int end_frame = *(unsigned int*)vars->get("end_frame")->get_async();
	progress = max(progress, (float)frame / end_frame);

	// And the estimated time to arrive
	float total_elapsed = elapsed * used_times();
	float ETA = total_elapsed * (1.f / progress - 1.f);

	data << "Percentage=" << std::setw(14) << std::setprecision(2)
	     << progress * 100.f;
	data << "   ETA=" << ETA << std::endl;

	// Replace the trailing space by a line break
	if (data.str().back() == ' ') {
		data.seekp(-1, data.cur);
		data << std::endl;
	}

	InputOutput::Logger::singleton()->writeReport(data.str(), _color, _bold);

	// Write the output file
	if (_output_file.compare("")) {
		auto f = writer();
		json jdata = { { "progress", progress },
		               { "elapsed", elapsed } };
		auto jsnapshot = json::object();
		jsnapshot["step"] = snapshot.step;
		jsnapshot["samples"] = json::array();
		for (auto sample : snapshot.samples) {
			auto jsample = json::object();
			jsample["name"] = sample.name;
			jsample["start"] = sample.start;
			jsample["end"] = sample.end;
			auto [indeps, outdeps] = sample.tool->getDependencies();
			auto jindeps = json::array();
			for (auto dep : indeps)
				jindeps.push_back(dep->name());
			jsample["dependencies"]["in"] = jindeps;
			auto joutdeps = json::array();
			for (auto dep : outdeps)
				joutdeps.push_back(dep->name());
			jsample["dependencies"]["out"] = joutdeps;
			jsnapshot["samples"].push_back(jsample);
		}
		jdata["snapshot"] = jsnapshot;
		f << jdata.dump(4) << std::endl;
		f.close();
	}

	cl_int err_code;
	err_code = clSetUserEventStatus(_user_event, CL_COMPLETE);
	if (err_code != CL_SUCCESS) {
		std::stringstream msg;
		msg << "Failure setting as complete the tool \"" << name() << "\"."
		    << std::endl;
		LOG(L_ERROR, msg.str());
		InputOutput::Logger::singleton()->printOpenCLError(err_code);
		return;
	}
	err_code = clReleaseEvent(_user_event);
	if (err_code != CL_SUCCESS) {
		std::stringstream msg;
		msg << "Failure releasing the user event at tool \"" << name() << "\"."
		    << std::endl;
		LOG(L_ERROR, msg.str());
		InputOutput::Logger::singleton()->printOpenCLError(err_code);
		return;
	}
}

size_t
Performance::computeAllocatedMemory()
{
	size_t allocated_mem = 0;
	auto C = CalcServer::singleton();

	// Get the allocated memory in the variables
	InputOutput::Variables* vars = C->variables();
	allocated_mem += vars->allocatedMemory();

	// Gwet the additionally allocated memory in the tools
	auto tools = C->tools();
	for (auto tool : tools) {
		allocated_mem += tool->allocatedMemory();
	}

	return allocated_mem;
}

/** @brief Callback called when all the command prior to
 * Aqua::CalcServer::Performance are completed.
 *
 * This function is just redirecting the work to
 * Aqua::CalcServer::Performance::print()
 * @param event The triggering event
 * @param event_command_status CL_COMPLETE upon all dependencies successfully
 * fulfilled. A negative integer if one or mor dependencies failed.
 * @param user_data A casted pointer to the Aqua::CalcServer::Performance
 * tool (or the inherited one)
 */
void CL_CALLBACK
performance_cb(cl_event event, cl_int event_command_status, void* user_data)
{
	clReleaseEvent(event);
	auto tool = (Performance*)user_data;
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
Performance::_execute(const std::vector<cl_event> events)
{
	cl_int err_code;
	auto C = CalcServer::singleton();
	cl_event event;
	// We are running everything on a parallel thread which is launched when
	// all the previous enqueued commands are completed
	err_code = clEnqueueMarkerWithWaitList(
		C->command_queue(), 0, NULL, &event);
	if (err_code != CL_SUCCESS) {
		std::stringstream msg;
		msg << "Failure setting the marker for tool \"" << name() << "\"."
		    << std::endl;
		LOG(L_ERROR, msg.str());
		InputOutput::Logger::singleton()->printOpenCLError(err_code);
		throw std::runtime_error("OpenCL execution error");
	}

	// Now we create a user event that we will set as completed when we already
	// readed the input dependencies
	_user_event = clCreateUserEvent(C->context(), &err_code);
	if (err_code != CL_SUCCESS) {
		std::stringstream msg;
		msg << "Failure creating the event for tool \"" << name() << "\"."
		    << std::endl;
		LOG(L_ERROR, msg.str());
		InputOutput::Logger::singleton()->printOpenCLError(err_code);
		throw std::runtime_error("OpenCL execution error");
	}

	// So it is time to register our callback on our trigger
	err_code = clRetainEvent(_user_event);
	if (err_code != CL_SUCCESS) {
		std::stringstream msg;
		msg << "Failure retaining the event for tool \"" << name() << "\"."
		    << std::endl;
		LOG(L_ERROR, msg.str());
		InputOutput::Logger::singleton()->printOpenCLError(err_code);
		throw std::runtime_error("OpenCL execution error");
	}
	err_code = clSetEventCallback(event, CL_COMPLETE, performance_cb, this);
	if (err_code != CL_SUCCESS) {
		std::stringstream msg;
		msg << "Failure registering the solver callback in tool \"" << name()
		    << "\"." << std::endl;
		LOG(L_ERROR, msg.str());
		InputOutput::Logger::singleton()->printOpenCLError(err_code);
		throw std::runtime_error("OpenCL execution error");
	}

	return _user_event;
}

std::ofstream
Performance::writer()
{
	std::ofstream f(_output_file, std::ios::out | std::ios::trunc);
	if (!f.is_open()) {
		std::ostringstream msg;
		msg << "Failure writing on '" << _output_file << "'" << std::endl;
		LOG(L_ERROR, msg.str());
	}
	return f;
}

}
}
} // namespace
