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
 * @brief Resize an array.
 * (See Aqua::CalcServer::Resize for details)
 */

#include <sstream>
#include "aquagpusph/AuxiliarMethods.hpp"
#include "aquagpusph/InputOutput/Logger.hpp"
#include "Resize.hpp"
#include "CalcServer.hpp"

namespace Aqua {
namespace CalcServer {

Resize::Resize(const std::string name,
         const std::string var_name,
         const std::string length,
         bool shrink,
         bool once)
  : ScalarExpression(name, length, "float", once)
  , _var_name(var_name)
  , _value(length)
  , _shrink(shrink)
  , _var(NULL)
{
	auto profilers = Profiler::substages();
	profilers.push_back(new EventProfile("Kernel", this));
	Profiler::substages(profilers);
}

Resize::~Resize()
{
}

void
Resize::setup()
{
	std::ostringstream msg;
	msg << "Loading the tool \"" << name() << "\"..." << std::endl;
	LOG(L_INFO, msg.str());

	variable();
	_var->reallocatable(true);

	auto C = CalcServer::singleton();
	setOutputType(
		C->device_addr_bits() == 64 ? "unsigned long" : "unsigned int");
	ScalarExpression::setup();
	ScalarExpression::_solve();
}

void
Resize::_solve()
{
	cl_int err_code;
	auto C = CalcServer::singleton();
	ScalarExpression::_solve();

	// Get the new size
	size_t n;
	if (C->device_addr_bits() == 64) {
		ulcl n_typed;
		memcpy(&n_typed, ScalarExpression::getValue(), sizeof(ulcl));
		n = (size_t)n_typed;
	} else {
		uicl n_typed;
		memcpy(&n_typed, ScalarExpression::getValue(), sizeof(uicl));
		n = (size_t)n_typed;
	}

	// Check if we can skip reallocating
	const size_t typesize = InputOutput::Variables::typeToBytes(_var->type());
	const size_t n_current = _var->size() / typesize;
	if ((n == n_current) ||
		(!_shrink && (n < n_current)))
	{
		return;
	}

	// Reallocate
	cl_mem mem = *(cl_mem*)_var->get();
	if (mem) {
		err_code = clReleaseMemObject(mem);
		CHECK_OCL_OR_THROW(err_code,
			std::string("Failure releasing ") +
			std::to_string(_var->size()) +
			" bytes on the device memory for tool \"" + name() + "\".");
	}
	mem = NULL;

	mem = clCreateBuffer(C->context(),
	                     CL_MEM_READ_WRITE,
	                     n * typesize,
	                     NULL,
	                     &err_code);
	CHECK_OCL_OR_THROW(err_code,
		std::string("Failure allocating ") +
		std::to_string(n * typesize) +
		" bytes on the device memory for tool \"" + name() + "\".");
	_var->set_async(&mem);
}

cl_event
Resize::_execute(const std::vector<cl_event> events)
{
	cl_int err_code;
	auto C = CalcServer::singleton();

	// Ask ScalarExpression subtask to do the job, which will call our ::_solve
	// function
	cl_event out_event = ScalarExpression::_execute(events);
	err_code = clFlush(C->command_queue());
	CHECK_OCL_OR_THROW(
	    err_code,
	    std::string("Failure flushing the command queue at \"") + name() +
	        "\".");

	return out_event;
}

void
Resize::variable()
{
	InputOutput::Variables* vars = CalcServer::singleton()->variables();
	if (!vars->get(_var_name)) {
		std::stringstream msg;
		msg << "The tool \"" << name()
		    << "\" is asking the undeclared variable \"" << _var_name << "\"."
		    << std::endl;
		LOG(L_ERROR, msg.str());
		throw std::runtime_error("Invalid variable");
	}
	if (!vars->get(_var_name)->isArray()) {
		std::stringstream msg;
		msg << "The tool \"" << name() << "\" is asking the variable \""
		    << _var_name << "\", which is a scalar." << std::endl;
		LOG(L_ERROR, msg.str());
		throw std::runtime_error("Invalid variable type");
	}
	_var = (InputOutput::ArrayVariable*)vars->get(_var_name);

	std::vector<InputOutput::Variable*> deps = { _var };
	setOutputDependencies(deps);
}

}
} // namespaces
