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
  , _output_file("")
  , _first_execution(true)
{
	gettimeofday(&_tic, NULL);
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

Performance::~Performance()
{
	if (_f.is_open())
		_f.close();
}

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

	// Open the output file
	if (_output_file.compare("")) {
		_f.open(_output_file.c_str(), std::ios::out);
		// Write the header
		_f << "\"t\",\"elapsed average (s)\",\"elapsed variance (s)\"";
		for (auto tool : C->tools()) {
			_f << ",\"" << tool->name() << " begin (ns)\""
			   << ",\"" << tool->name() << " begin variance (ns)\""
			   << ",\"" << tool->name() << " end (ns)\""
			   << ",\"" << tool->name() << " end variance (ns)\"";
		}
		_f << std::endl;
	}

	Tool::setup();
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

/** @brief Get the starting timers
 * @param tools The tools to analyze
 * @return The reference timers
 */
std::deque<cl_ulong>
ref_timer(std::vector<Tool*> tools)
{
	for (auto tool : tools)
		if (tool->begin().size())
			return tool->begin();
	LOG(L_WARNING, "No profiled tools were found");
	return tools.front()->begin();
}

/** @brief Get the average and standard deviation of the beginning timers
 * @param tools The tools to analyze
 * @return The reference timers
 */
std::tuple<std::vector<cl_ulong>, std::vector<cl_ulong>>
begin_timers(std::vector<Tool*> tools)
{
	std::vector<cl_ulong> v, e;
	auto t0 = ref_timer(tools);
	for (auto tool : tools) {
		if (!tool->begin().size()) {
			v.push_back(0.0);
			e.push_back(0.0);
			continue;
		}
		auto [t, t_std] = stats(Profiler::delta(tool->begin(), t0));
		v.push_back(t);
		e.push_back(t_std);
	}
	return { v, e };
}

/** @brief Get the average and standard deviation of the ending timers
 * @param tools The tools to analyze
 * @return The reference timers
 */
std::tuple<std::vector<cl_ulong>, std::vector<cl_ulong>>
end_timers(std::vector<Tool*> tools)
{
	std::vector<cl_ulong> v, e;
	auto t0 = ref_timer(tools);
	for (auto tool : tools) {
		if (!tool->end().size()) {
			v.push_back(0.0);
			e.push_back(0.0);
			continue;
		}
		auto [t, t_std] = stats(Profiler::delta(tool->end(), t0));
		v.push_back(t);
		e.push_back(t_std);
	}
	return { v, e };
}

cl_event
Performance::_execute(const std::vector<cl_event> events)
{
	auto C = CalcServer::singleton();
	clFinish(C->command_queue());
	std::stringstream data;

	size_t allocated_MB = computeAllocatedMemory() / (1024 * 1024);
	data << "Performance:" << std::endl
	     << "Memory=" << std::setw(18) << allocated_MB << "MB" << std::endl;

	// Add the tools time elapsed
	auto tools = C->tools();
	float elapsed = 0.f, elapsed_std = 0.f;
	for (auto tool : tools) {
		// Exclude the tool itself
		if (this == tool) {
			continue;
		}
		auto [t_ave, t_std] = tool->elapsed();
		elapsed += t_ave * 1.e-9;
		elapsed_std += t_std * 1.e-9;
	}

	data << "Elapsed=" << std::setw(17) << elapsed << "s  (+-" << std::setw(16)
	     << elapsed_std << "s)" << std::endl;

	// Compute the progress
	InputOutput::Variables* vars = C->variables();
	float progress = 0.f;
	float t = *(float*)vars->get("t")->get();
	float end_t = *(float*)vars->get("end_t")->get();
	progress = max(progress, t / end_t);
	unsigned int iter = *(unsigned int*)vars->get("iter")->get();
	unsigned int end_iter = *(unsigned int*)vars->get("end_iter")->get();
	progress = max(progress, (float)iter / end_iter);
	unsigned int frame = *(unsigned int*)vars->get("frame")->get();
	unsigned int end_frame = *(unsigned int*)vars->get("end_frame")->get();
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
	if (_f.is_open()) {
		_f << t << "," << elapsed << "," << elapsed_std;
		auto [t0, t0_std] = begin_timers(tools);
		auto [t1, t1_std] = end_timers(tools);
		for (unsigned int i = 0; i < tools.size(); i++) {
			_f << "," << t0[i] << "," << t0_std[i] << "," << t1[i] << ","
			   << t1_std[i];
		}
		_f << std::endl;
	}

	return NULL;
}

}
}
} // namespace
