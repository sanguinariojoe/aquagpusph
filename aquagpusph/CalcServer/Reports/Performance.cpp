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

Performance::~Performance()
{
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

/** @brief Get the reference timers
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


/** @brief Timer data holder
 */
typedef struct _timer {
	/// Average value
	cl_long v;
	/// Standard deviation
	cl_long e;
} timer;

/** @brief Get the timers of the tools
 * @param tools The tools to analyze
 * @param t0 Reference timer
 * @return The beginning and ending timers
 */
std::tuple<std::vector<timer>, std::vector<timer>>
tools_timers(std::vector<Tool*> tools, std::deque<cl_ulong> t0)
{
	std::vector<timer> begin, end;
	for (auto tool : tools) {
		if (!tool->begin().size()) {
			begin.push_back(timer({0, 0}));
			end.push_back(timer({0, 0}));
			continue;
		}
		{
			auto [t, t_std] = stats(Profiler::delta(tool->begin(), t0));
			begin.push_back(timer({t, t_std}));
		}
		{
			auto [t, t_std] = stats(Profiler::delta(tool->end(), t0));
			end.push_back(timer({t, t_std}));
		}
	}
	return { begin, end };
}

/** @brief Get the substage timers
 * @param substage The substage to analyze
 * @param t0 Reference timer
 * @return The beginning and ending timers
 */
std::tuple<timer, timer>
substage_timers(Profile* substage, std::deque<cl_ulong> t0)
{
	timer begin, end;
	{
		auto [t, t_std] = stats(Profiler::delta(substage->begin(), t0));
		begin.v = t;
		begin.e = t_std;
	}
	{
		auto [t, t_std] = stats(Profiler::delta(substage->end(), t0));
		end.v = t;
		end.e = t_std;
	}
	return { begin, end };
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
	if (_output_file.compare("")) {
		auto f = writer();
		json jdata = {
			{"progress", progress},
			{"elapsed", {
				{"avg", elapsed},
				{"std", elapsed_std}
			}}
		};
		jdata["tools"] = json::array();
		auto tref = ref_timer(tools);
		auto [t0, t1] = tools_timers(tools, tref);
		for (unsigned int i = 0; i < tools.size(); i++) {
			auto tool =tools[i];
			auto jtool = json::object();
			jtool["name"] = tool->name();
			auto [indeps, outdeps] = tool->getDependencies();
			auto jindeps = json::array();
			for (auto dep : indeps)
				jindeps.push_back(dep->name());
			jtool["dependencies"]["in"] = jindeps;
			auto joutdeps = json::array();
			for (auto dep : outdeps)
				joutdeps.push_back(dep->name());
			jtool["dependencies"]["out"] = joutdeps;
			jtool["start"]["avg"] = t0[i].v;
			jtool["start"]["std"] = t0[i].e;
			jtool["end"]["avg"] = t1[i].v;
			jtool["end"]["std"] = t1[i].e;
			jtool["substages"] = json::array();
			auto substages = tool->subinstances();
			for (auto substage : substages) {
				auto [tt0, tt1] = substage_timers(substage, tref);
				auto jsubstage = json::object();
				jsubstage["name"] = substage->name();
				jsubstage["start"]["avg"] = tt0.v;
				jsubstage["start"]["std"] = tt0.e;
				jsubstage["end"]["avg"] = tt1.v;
				jsubstage["end"]["std"] = tt1.e;
				jtool["substages"].push_back(jsubstage);
			}
			jdata["tools"].push_back(jtool);
		}
		f << jdata.dump(4) << std::endl;
		f.close();
	}

	return NULL;
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
