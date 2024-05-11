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
 * @brief Simulation configuration data structures.
 * (See Aqua::InputOutput::ProblemSetup for details)
 */

#include <limits>
#include <sstream>

#include "ProblemSetup.hpp"
#include "InputOutput/Logger.hpp"
#include "AuxiliarMethods.hpp"

namespace Aqua {
namespace InputOutput {

ProblemSetup::ProblemSetup()
  : _n_cmd_queues(15)
  , _copy(false)
{
	time_opts.sim_end_mode = __NO_OUTPUT_MODE__;
	time_opts.sim_end_time = 0.f;
	time_opts.sim_end_step = 0;
	time_opts.sim_end_frame = 0;
	time_opts.output_mode = __NO_OUTPUT_MODE__;
	time_opts.output_fps = 0.f;
	time_opts.output_ipf = 0;
}

ProblemSetup::ProblemSetup(const ProblemSetup& p)
  : settings(p.settings)
  , variables(p.variables)
  , definitions(p.definitions)
  , tools(p.tools)
  , reports(p.reports)
  , time_opts(p.time_opts)
  , sets(p.sets)
  , _copy(true)
{
}

ProblemSetup::~ProblemSetup()
{
	if (_copy)
		return;

	unsigned int i;
	for (i = tools.size(); i > 0; i--) {
		if (toolInstances(tools.at(i - 1)) == 1) {
			// This is the last remaining instance
			delete tools.at(i - 1);
		}
		tools.erase(tools.begin() + i - 1);
	}
	tools.clear();
	for (auto report : reports) {
		delete report;
	}
	reports.clear();
	for (auto set : sets) {
		delete set;
	}
	sets.clear();
}

ProblemSetup::sphSettings::sphSettings()
  : save_on_fail(true)
  , debug_tools(false)
  , base_path("")
{
}

void
ProblemSetup::sphVariables::registerVariable(std::string name,
                                             std::string type,
                                             std::string length,
                                             std::string value)
{
	names.push_back(name);
	types.push_back(type);
	lengths.push_back(length);
	values.push_back(value);
}

void
ProblemSetup::sphDefinitions::define(const std::string name,
                                     const std::string value,
                                     const bool evaluate)
{
	if (isDefined(name)) {
		undefine(name);
	}

	names.push_back(name);
	values.push_back(value);
	evaluations.push_back(evaluate);
}

bool
ProblemSetup::sphDefinitions::isDefined(const std::string name)
{
	for (auto str : names) {
		if (!name.compare(str)) {
			return true;
		}
	}

	return false;
}

void
ProblemSetup::sphDefinitions::undefine(const std::string name)
{
	unsigned int i;
	for (i = 0; i < names.size(); i++) {
		if (!name.compare(names.at(i))) {
			names.erase(names.begin() + i);
			values.erase(values.begin() + i);
			evaluations.erase(evaluations.begin() + i);
			return;
		}
	}
}

void
ProblemSetup::sphTool::set(const std::string name, const std::string value)
{
	if (has(name)) {
		_data[name] = value;
		return;
	}
	_data.insert(std::make_pair(name, value));
}

const std::string
ProblemSetup::sphTool::get(const std::string name)
{
	return _data[name].c_str();
}

const std::string
ProblemSetup::sphTool::get(unsigned int index)
{
	unsigned int i = 0;
	for (auto& d : _data) {
		if (i == index) {
			return d.second;
		}
		i++;
	}
	return NULL;
}

const std::string
ProblemSetup::sphTool::getName(unsigned int index)
{
	unsigned int i = 0;
	for (auto& d : _data) {
		if (i == index) {
			return d.first;
		}
		i++;
	}
	return NULL;
}

bool
ProblemSetup::sphTool::has(const std::string name)
{
	std::map<std::string, std::string>::iterator var = _data.find(name);
	if (var != _data.end())
		return true;
	return false;
}

unsigned int
ProblemSetup::toolInstances(ProblemSetup::sphTool* tool)
{
	unsigned int i, n = 0;
	for (auto t : tools) {
		if (tool == t) {
			n++;
		}
	}
	return n;
}

ProblemSetup::sphParticlesSet::sphParticlesSet()
  : _n(0)
{
}

ProblemSetup::sphParticlesSet::~sphParticlesSet() {}

void
ProblemSetup::sphParticlesSet::addScalar(std::string name, std::string value)
{
	_snames.push_back(name);
	_svalues.push_back(value);
}

void
ProblemSetup::sphParticlesSet::input(std::string path,
                                     std::string format,
                                     std::string fields)
{
	_in_path = setStrConstantsCopy(path);
	_in_format = format;
	// Split the fields
	std::istringstream f(replaceAllCopy(fields, " ", ""));
	std::string s;
	while (getline(f, s, ',')) {
		_in_fields.push_back(s);
	}
}

void
ProblemSetup::sphParticlesSet::output(std::string path,
                                      std::string format,
                                      std::string fields)
{
	_out_path = setStrConstantsCopy(replaceAllCopy(path, "%d", "{index}"));
	_out_format = format;
	// Split the fields
	std::istringstream f(replaceAllCopy(fields, " ", ""));
	std::string s;
	while (getline(f, s, ',')) {
		_out_fields.push_back(s);
	}
}

}
} // namespace
