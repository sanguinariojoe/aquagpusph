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
 * @brief Set a scalar variable.
 * (See Aqua::CalcServer::SetScalar for details)
 */

#include <InputOutput/Logger.h>
#include <CalcServer/SetScalar.h>
#include <CalcServer.h>

namespace Aqua {
namespace CalcServer {

SetScalar::SetScalar(const std::string name,
                     const std::string var_name,
                     const std::string value,
                     bool once)
  : Tool(name, once)
  , _var_name(var_name)
  , _value(value)
  , _var(NULL)
  , _event(NULL)
{
}

SetScalar::~SetScalar() {}

void
SetScalar::setup()
{
	std::ostringstream msg;
	msg << "Loading the tool \"" << name() << "\"..." << std::endl;
	LOG(L_INFO, msg.str());

	Tool::setup();
	variables();
}

void CL_CALLBACK
solver(cl_event event, cl_int event_command_status, void* user_data)
{
	auto tool = (SetScalar*)user_data;

	clReleaseEvent(event);
	InputOutput::Variables* vars = CalcServer::singleton()->variables();

	auto var = tool->getOutputVariable();
	cl_event user_event = tool->getEvent();

	void* data = malloc(var->typesize());
	if (!data) {
		std::stringstream msg;
		msg << "Failure allocating " << var->typesize()
		    << " bytes for the variable \"" << var->name() << "\"."
		    << std::endl;
		LOG(L_ERROR, msg.str());
		clSetUserEventStatus(user_event, -1);
		return;
	}

	try {
		vars->solve(var->type(), tool->getExpression(), data, "__NONE");
	} catch (...) {
		free(data);
		clSetUserEventStatus(user_event, -2);
		return;
	}

	var->set_async(data);
	free(data);
	// Ensure that the variable is populated
	vars->populate(var);

	clSetUserEventStatus(user_event, CL_COMPLETE);
	clReleaseEvent(user_event);
}

cl_event
SetScalar::_execute(const std::vector<cl_event> events)
{
	cl_int err_code;
	cl_event trigger;
	CalcServer* C = CalcServer::singleton();

	// We create a trigger event to be marked as completed when all the
	// dependencies are fullfiled
	cl_uint num_events_in_wait_list = events.size();
	const cl_event* event_wait_list = events.size() ? events.data() : NULL;

	err_code = clEnqueueMarkerWithWaitList(
	    C->command_queue(), num_events_in_wait_list, event_wait_list, &trigger);
	if (err_code != CL_SUCCESS) {
		std::stringstream msg;
		msg << "Failure setting the trigger for tool \"" << name() << "\"."
		    << std::endl;
		LOG(L_ERROR, msg.str());
		InputOutput::Logger::singleton()->printOpenCLError(err_code);
		throw std::runtime_error("OpenCL execution error");
	}

	// Now we create a user event that we will set as completed when we already
	// solved the equation and set the varaible value
	auto event = clCreateUserEvent(C->context(), &err_code);
	if (err_code != CL_SUCCESS) {
		std::stringstream msg;
		msg << "Failure creating the event for tool \"" << name() << "\"."
		    << std::endl;
		LOG(L_ERROR, msg.str());
		InputOutput::Logger::singleton()->printOpenCLError(err_code);
		throw std::runtime_error("OpenCL execution error");
	}
	err_code = clRetainEvent(event);
	if (err_code != CL_SUCCESS) {
		std::stringstream msg;
		msg << "Failure retaining the event for tool \"" << name() << "\"."
		    << std::endl;
		LOG(L_ERROR, msg.str());
		InputOutput::Logger::singleton()->printOpenCLError(err_code);
		throw std::runtime_error("OpenCL execution error");
	}
	_event = event;

	// So it is time to register our callback on our trigger
	err_code = clSetEventCallback(trigger, CL_COMPLETE, solver, this);
	if (err_code != CL_SUCCESS) {
		std::stringstream msg;
		msg << "Failure registering the solver callback in tool \"" << name()
		    << "\"." << std::endl;
		LOG(L_ERROR, msg.str());
		InputOutput::Logger::singleton()->printOpenCLError(err_code);
		throw std::runtime_error("OpenCL execution error");
	}

	return _event;
}

InputOutput::Variable*
SetScalar::variable(const std::string& var_name) const
{
	InputOutput::Variables* vars = CalcServer::singleton()->variables();
	auto var = vars->get(var_name);
	if (!var) {
		std::stringstream msg;
		msg << "The tool \"" << name()
		    << "\" is querying the undeclared variable \"" << var_name << "\"."
		    << std::endl;
		LOG(L_ERROR, msg.str());
		throw std::runtime_error("Invalid variable");
	}
	if (var->type().find('*') != std::string::npos) {
		std::stringstream msg;
		msg << "The tool \"" << name() << "\" is considering the variable \""
		    << _var_name << "\", which is an array." << std::endl;
		LOG(L_ERROR, msg.str());
		throw std::runtime_error("Invalid variable type");
	}
	return var;
}

void
SetScalar::variables()
{
	InputOutput::Variables* vars = CalcServer::singleton()->variables();
	_in_vars = vars->exprVariables(_value);
	_var = variable(_var_name);
	std::vector<InputOutput::Variable*> out_vars{ _var };
	setDependencies(_in_vars, out_vars);
}

}
} // namespaces
