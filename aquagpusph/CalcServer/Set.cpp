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
 * @brief Set all the components of an array with the desired value.
 * (See Aqua::CalcServer::Set for details)
 * @note Hardcoded versions of the files CalcServer/Set.cl.in and
 * CalcServer/Set.hcl.in are internally included as a text array.
 */

#include <sstream>
#include "aquagpusph/AuxiliarMethods.hpp"
#include "aquagpusph/InputOutput/Logger.hpp"
#include "Set.hpp"
#include "CalcServer.hpp"

namespace Aqua {
namespace CalcServer {

#ifndef DOXYGEN_SHOULD_SKIP_THIS
#include "aquagpusph/CalcServer/Set.hcl"
#include "aquagpusph/CalcServer/Set.cl"
#endif
std::string SET_INC = xxd2string(Set_hcl_in, Set_hcl_in_len);
std::string SET_SRC = xxd2string(Set_cl_in, Set_cl_in_len);

Set::Set(const std::string name,
         const std::string var_name,
         const std::string value,
         bool once)
  : ScalarExpression(name, value, "float", once)
  , _var_name(var_name)
  , _value(value)
  , _var(NULL)
  , _input(NULL)
  , _kernel(NULL)
  , _global_work_size(0)
  , _work_group_size(0)
  , _n(0)
  , _data(NULL)
{
	auto profilers = Profiler::substages();
	profilers.push_back(new EventProfile("Kernel", this));
	Profiler::substages(profilers);
}

Set::~Set()
{
	if (_kernel)
		clReleaseKernel(_kernel);
	if (_data)
		free(_data);
}

void
Set::setup()
{
	std::ostringstream msg;
	msg << "Loading the tool \"" << name() << "\"..." << std::endl;
	LOG(L_INFO, msg.str());

	variable();

	size_t typesize = InputOutput::Variables::typeToBytes(_var->type());
	_data = malloc(typesize);
	if (!_data) {
		std::stringstream msg;
		msg << "Failure allocating " << typesize << " bytes for the variable \""
		    << _var->name() << "\" value." << std::endl;
		LOG(L_ERROR, msg.str());
		throw std::bad_alloc();
	}
	try {
		ScalarExpression::setup();
		setOutputType(_var->type());
		ScalarExpression::_solve();
	} catch (...) {
		free(_data);
		_data = NULL;
		LOG(L_INFO, "Falling back to definition mode...\n");
	}

	_input = *(cl_mem*)_var->get();
	_n = _var->size() / typesize;
	setupOpenCL();
}

void
Set::_solve()
{
	if (_data) {
		ScalarExpression::_solve();
		memcpy(_data,
		       getValue(),
		       InputOutput::Variables::typeToBytes(_var->type()));
	}

	setVariables();
}

cl_event
Set::_execute(const std::vector<cl_event> events)
{
	cl_int err_code;
	auto C = CalcServer::singleton();

	// Eventually evaluate the expression and set the variables to the kernel
	// To do that we do not need to wait until the output variable is ready
	auto in_events = getEvents(dep_events::in);
	cl_event args_event = ScalarExpression::_execute(in_events);
	err_code = clFlush(C->command_queue());
	CHECK_OCL_OR_THROW(
	    err_code,
	    std::string("Failure flushing the command queue at \"") + name() +
	        "\".");
	for (auto event : in_events) {
		// Retain the event until we work with it
		err_code = clReleaseEvent(event);
		CHECK_OCL_OR_THROW(err_code,
		                   std::string("Failure releasing an event in \"") +
		                       name() + "\" tool.");
	}
	err_code = clWaitForEvents(1, &args_event);
	CHECK_OCL_OR_THROW(
	    err_code,
	    std::string("Failure waiting for the args setter in tool \"") +
	        name() + "\".");
	err_code = clReleaseEvent(args_event);
	CHECK_OCL_OR_THROW(
	    err_code,
	    std::string("Failure releasing the args setter event in tool \"") +
	        name() + "\".");

	// And enqueue the kernel execution
	cl_event event;
	const cl_event* wait_events = events.size() ? events.data() : NULL;
	err_code = clEnqueueNDRangeKernel(C->command_queue(),
	                                  _kernel,
	                                  1,
	                                  NULL,
	                                  &_global_work_size,
	                                  &_work_group_size,
	                                  events.size(),
	                                  wait_events,
	                                  &event);
	CHECK_OCL_OR_THROW(err_code,
	                   std::string("Failure executing the tool \"") + name() +
	                       "\".");
	auto profiler = dynamic_cast<EventProfile*>(Profiler::substages().back());
	profiler->start(event);
	profiler->end(event);

	return event;
}

void
Set::variable()
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

void
Set::setupOpenCL()
{
	cl_int err_code;
	cl_kernel kernel;
	CalcServer* C = CalcServer::singleton();

	// Create a header for the source code where the operation will be placed
	std::ostringstream source;
	source << SET_INC;
	if (!_data) {
		// No valid equation is available, so let's try to use _value as a
		// straight definition
		source << " #define VALUE " << _value;
	}
	source << SET_SRC;

	// Compile the kernel
	std::ostringstream flags;
	if (!_var->type().compare("unsigned int*")) {
		// Spaces are not a good business into definitions passed as args
		flags << " -DT=uint";
	} else if (!_var->type().compare("unsigned long*")) {
		// Spaces are not a good business into definitions passed as args
		flags << " -DT=ulong";
	} else {
		std::string t = trimCopy(_var->type());
		t.pop_back(); // Remove the asterisk
		flags << " -DT=" << t;
	}
	kernel = compile_kernel(source.str(), "set", flags.str());
	err_code = clGetKernelWorkGroupInfo(kernel,
	                                    C->device(),
	                                    CL_KERNEL_WORK_GROUP_SIZE,
	                                    sizeof(size_t),
	                                    &_work_group_size,
	                                    NULL);
	CHECK_OCL_OR_THROW(
	    err_code,
	    std::string("Failure querying CL_KERNEL_WORK_GROUP_SIZE in tool \"") +
	        name() + "\".");
	if (_work_group_size < __CL_MIN_LOCALSIZE__) {
		LOG(L_ERROR, "insufficient local memory.\n");
		std::stringstream msg;
		msg << "\t" << _work_group_size
		    << " local work group size with __CL_MIN_LOCALSIZE__="
		    << __CL_MIN_LOCALSIZE__ << std::endl;
		LOG0(L_DEBUG, msg.str());
		throw std::runtime_error("OpenCL error");
	}

	_global_work_size = roundUp<size_t>(_n, _work_group_size);
	_kernel = kernel;
	err_code = clSetKernelArg(kernel, 0, _var->typesize(), _var->get());
	CHECK_OCL_OR_THROW(
	    err_code,
	    std::string("Failure sending the array argument to tool \"") + name() +
	        "\".");
	err_code = C->setKernelSizeArg(kernel, 1, _n);
	CHECK_OCL_OR_THROW(
	    err_code,
	    std::string("Failure sending the array size argument to tool \"") +
	        name() + "\".");
	if (!_data)
		return;
	err_code = clSetKernelArg(
	    kernel, 2, InputOutput::Variables::typeToBytes(_var->type()), _data);
	CHECK_OCL_OR_THROW(
	    err_code,
	    std::string("Failure sending the value argument to tool \"") + name() +
	        "\".");
}

void
Set::setVariables()
{
	cl_int err_code;

	if (_data) {
		err_code =
		    clSetKernelArg(_kernel,
		                   2,
		                   InputOutput::Variables::typeToBytes(_var->type()),
		                   _data);
		CHECK_OCL_OR_THROW(
		    err_code,
		    std::string("Failure setting the value to the tool \"") + name() +
		        "\".");
	}

	if (_input != *(cl_mem*)_var->get(false)) {
		// For some reason the input variable has changed...
		err_code = clSetKernelArg(_kernel, 0, _var->typesize(), _var->get());
		CHECK_OCL_OR_THROW(err_code,
		                   std::string("Failure setting the variable \"") +
		                       _var->name() + "\" to the tool \"" + name() +
		                       "\".");
		_input = *(cl_mem*)_var->get(false);
	}
}

}
} // namespaces
