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

#include "aquagpusph/InputOutput/Logger.hpp"
#include "SetScalar.hpp"
#include "CalcServer.hpp"
#include <sys/time.h>

namespace Aqua {
namespace CalcServer {

ScalarExpression::ScalarExpression(const std::string name,
                                   const std::string expr,
                                   const std::string type,
                                   bool once)
  : Tool(name, once)
  , _value(expr)
  , _output(NULL)
  , _output_type(type)
{
	setOutputType(type);
	Profiler::subinstances( { new ScalarProfile("solver") } );
}

ScalarExpression::~ScalarExpression()
{
	free(_output);
}

void
ScalarExpression::setup()
{
	std::ostringstream msg;
	msg << "Loading the tool \"" << name() << "\"..." << std::endl;
	LOG(L_INFO, msg.str());

	Tool::setup();
	variables();
}

void
ScalarExpression::setOutputType(const std::string type)
{
	_output_type = type;
	size_t typesize = InputOutput::Variables::typeToBytes(type);
	free(_output);
	_output = malloc(typesize);
	if (!_output) {
		std::ostringstream msg;
		msg << "Failure allocating " << typesize << " bytes (" << type
		    << ")for the output." << std::endl;
		LOG(L_ERROR, msg.str());
		throw std::bad_alloc();
	}
}

void
ScalarExpression::solve()
{
	cl_int err_code;

	dynamic_cast<ScalarProfile*>(Profiler::subinstances().front())->start();

	cl_event user_event = getUserEvent();

	try {
		_solve();
	} catch (...) {
		std::stringstream msg;
		LOG(L_ERROR, msg.str());
		err_code = clSetUserEventStatus(user_event, -1);
		if (err_code != CL_SUCCESS) {
			std::stringstream msg;
			msg << "Failure setting the error on the tool \"" << name()
			    << "\" user event." << std::endl;
			LOG(L_ERROR, msg.str());
		}
		return;
	}

	err_code = clSetUserEventStatus(user_event, CL_COMPLETE);
	if (err_code != CL_SUCCESS) {
		std::stringstream msg;
		msg << "Failure setting as complete the tool \"" << name() << "\""
		    << std::endl;
		LOG(L_ERROR, msg.str());
		InputOutput::Logger::singleton()->printOpenCLError(err_code);
		return;
	}

	err_code = clReleaseEvent(user_event);
	if (err_code != CL_SUCCESS) {
		std::stringstream msg;
		msg << "Failure releasing the user even of tool \"" << name() << "\""
		    << std::endl;
		LOG(L_ERROR, msg.str());
		InputOutput::Logger::singleton()->printOpenCLError(err_code);
		return;
	}

	dynamic_cast<ScalarProfile*>(Profiler::subinstances().front())->end();
}

/** @brief Callback called when all the dependencies of the
 * Aqua::CalcServer::ScalarExpression tool are fulfilled.
 *
 * This function is just redirecting the work to
 * Aqua::CalcServer::ScalarExpression::solve()
 * @param event The triggering event
 * @param event_command_status CL_COMPLETE upon all dependencies successfully
 * fulfilled. A negative integer if one or mor dependencies failed.
 * @param user_data A casted pointer to the Aqua::CalcServer::ScalarExpression
 * tool (or the inherited one)
 */
void CL_CALLBACK
solver(cl_event event, cl_int event_command_status, void* user_data)
{
	clReleaseEvent(event);
	auto tool = (ScalarExpression*)user_data;
	if (event_command_status != CL_COMPLETE) {
		std::stringstream msg;
		msg << "Skipping \"" << tool->name() << "\" due to dependency errors."
		    << std::endl;
		LOG(L_WARNING, msg.str());
		clSetUserEventStatus(tool->getUserEvent(), event_command_status);
		clReleaseEvent(tool->getUserEvent());
		return;
	}

	tool->solve();
}

void
ScalarExpression::_solve()
{
	InputOutput::Variables* vars = CalcServer::singleton()->variables();
	try {
		vars->solve(getOutputType(), getExpression(), _output, "__NONE");
	} catch (...) {
		throw;
	}
}

cl_event
ScalarExpression::_execute(const std::vector<cl_event> events)
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
	_user_event = event;

	// So it is time to register our callback on our trigger
	err_code = clRetainEvent(event);
	if (err_code != CL_SUCCESS) {
		std::stringstream msg;
		msg << "Failure retaining the event for tool \"" << name() << "\"."
		    << std::endl;
		LOG(L_ERROR, msg.str());
		InputOutput::Logger::singleton()->printOpenCLError(err_code);
		throw std::runtime_error("OpenCL execution error");
	}
	err_code = clSetEventCallback(trigger, CL_COMPLETE, solver, this);
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

void
ScalarExpression::variables()
{
	InputOutput::Variables* vars = CalcServer::singleton()->variables();
	_in_vars = vars->exprVariables(_value);
	setInputDependencies(_in_vars);
}

// We start setting the outpu variable type as float, and modify it later
SetScalar::SetScalar(const std::string name,
                     const std::string var_name,
                     const std::string value,
                     bool once)
  : ScalarExpression(name, value, "float", once)
  , _var_name(var_name)
  , _var(NULL)
{
}

SetScalar::~SetScalar() {}

void
SetScalar::setup()
{
	ScalarExpression::setup();
	_var = variable(_var_name);
	if (!_var) {
		std::stringstream msg;
		msg << "Invalid output variable \"" << _var_name << " for tool \""
		    << name() << "\"." << std::endl;
		LOG(L_ERROR, msg.str());
		throw std::invalid_argument("Invalid variable");
	}
	setOutputType(_var->type());
	std::vector<InputOutput::Variable*> out_vars{ _var };
	setOutputDependencies(out_vars);
}

void
SetScalar::_solve()
{
	ScalarExpression::_solve();
	InputOutput::Variables* vars = CalcServer::singleton()->variables();
	_var->set_async((void*)getValue());
	vars->populate(_var);
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

}
} // namespaces
