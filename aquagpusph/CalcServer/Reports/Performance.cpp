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
#ifdef VTK_NLOHMANN_JSON
#include <vtk_nlohmannjson.h>
#include VTK_NLOHMANN_JSON(json.hpp)
#else
#include <nlohmann/json.hpp>
#endif // VTK_NLOHMANN JSON


using json = nlohmann::json;

namespace Aqua {
namespace CalcServer {
namespace Reports {

Performance::Performance(const std::string tool_name,
                         const std::string output_file)
  : Report(tool_name, "dummy_fields_string")
  , _output_file(output_file)
  , _output_file_index(0)
  , _timer(0)
{
}

Performance::~Performance() {}

void
Performance::setup()
{
	std::ostringstream msg;
	msg << "Loading the report \"" << name() << "\"..." << std::endl;
	LOG(L_INFO, msg.str());

	Tool::setup();

	setInputDependencies(
	    { "t", "end_t", "iter", "end_iter", "frame", "end_frame" });
}

/// List of time format multipliers
static const std::vector<float> time_mults(
	{ 1.0, 1.e-3, 1.e-6, 1.e-9, 1.e-9 / 60, 1.e-9 / 3600 });
/// List of time units
static const std::vector<std::string> time_units(
	{ "ns", "μs", "ms", "s", "m", "h" });

/** @brief Reformat the time in nanoseconds as the best possible alternative
 * between nanoseconds, microseconds, milliseconds, seconds minutes or hours
 * @param t Time in nanoseconds
 * @return The reformated time as well as the resulting units, i.e. ns, μs, ms,
 * s, m or h
 */
std::tuple<float, std::string> reformat_nanosec(float t)
{
	float t_val;
	std::string t_units;
	for (unsigned int i = 0; i < time_mults.size(); i++) {
		t_val = t * time_mults[i];
		t_units = time_units[i];
		if(t_val < 1.0)
			break;
	}
	return { t_val, t_units };
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
	auto elapsed =
	    (!_timer) ? tools_elapsed : ProfilingInfo::delta(timer, _timer);
	_timer = timer;
	auto [elapsed_v, elapsed_u] = reformat_nanosec(elapsed);
	data << "Elapsed=" << std::setw(17) << elapsed_v << elapsed_u << std::endl;
#ifdef HAVE_GPUPROFILE
	auto overhead = elapsed - tools_elapsed;
	auto [overhead_v, overhead_u] = reformat_nanosec(overhead);
	data << "Overhead=" << std::setw(17) << overhead_v << overhead_u
	     << std::endl;
#endif

	// Compute the progress
	InputOutput::Variables* vars = C->variables();
	float progress = 0.f;
	float t = *(float*)vars->get("t")->get_async();
	float end_t = *(float*)vars->get("end_t")->get_async();
	progress = (std::max)(progress, t / end_t);
	unsigned int iter = *(unsigned int*)vars->get("iter")->get_async();
	unsigned int end_iter = *(unsigned int*)vars->get("end_iter")->get_async();
	progress = (std::max)(progress, (float)iter / end_iter);
	unsigned int frame = *(unsigned int*)vars->get("frame")->get_async();
	unsigned int end_frame =
	    *(unsigned int*)vars->get("end_frame")->get_async();
	progress = (std::max)(progress, (float)frame / end_frame);

	// And the estimated time to arrive
	float total_elapsed = elapsed * used_times();
	float ETA = total_elapsed * (1.f / progress - 1.f);

	data << "Percentage=" << std::setw(14) << std::setprecision(2)
	     << progress * 100.f;
	auto [ETA_v, ETA_u] = reformat_nanosec(ETA);
	data << "   ETA=" << ETA_v << ETA_u << std::endl;

	// Replace the trailing space by a line break
	if (data.str().back() == ' ') {
		data.seekp(-1, data.cur);
		data << std::endl;
	}

	InputOutput::Logger::singleton()->writeReport(data.str());

	// Write the output file
	if (_output_file.compare("")) {
		auto f = writer();
		json jdata = { { "progress", progress }, { "elapsed", elapsed } };
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

	setUserEventStatus(CL_COMPLETE);
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

/** @brief Callback called when all the commands prior to
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
		tool->setUserEventStatus(event_command_status);
		return;
	}

	tool->print();
}

cl_event
Performance::_execute(const std::vector<cl_event> UNUSED_PARAM events)
{
	auto C = CalcServer::singleton();
	return setCallback({ C->marker() }, performance_cb);
}

std::ofstream
Performance::writer()
{
	std::string output_file;
	try {
		output_file = newFilePath(_output_file, _output_file_index, 6);
	} catch (std::invalid_argument& e) {
		output_file = setStrConstantsCopy(_output_file);
	}

	std::ofstream f(output_file, std::ios::out | std::ios::trunc);
	if (!f.is_open()) {
		std::ostringstream msg;
		msg << "Failure writing on '" << output_file << "'" << std::endl;
		LOG(L_ERROR, msg.str());
	}
	return f;
}

}
}
} // namespace
