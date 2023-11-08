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

#include <AuxiliarMethods.h>
#include <InputOutput/Logger.h>
#include <CalcServer/Set.h>
#include <CalcServer.h>

namespace Aqua {
namespace CalcServer {

#ifndef DOXYGEN_SHOULD_SKIP_THIS
#include "CalcServer/Set.hcl"
#include "CalcServer/Set.cl"
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
  , _local_work_size(0)
  , _n(0)
  , _data(NULL)
{
}

Set::~Set()
{
	if (_kernel)
		clReleaseKernel(_kernel);
	_kernel = NULL;
	if (_data)
		free(_data);
	_data = NULL;
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
	if(_data) {
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
	unsigned int i;
	cl_int err_code;
	cl_event event;
	CalcServer* C = CalcServer::singleton();

	// Eventually evaluate the expression and set the variables to the kernel
	cl_event args_event = ScalarExpression::_execute(events);

	std::vector<cl_event> wait_list = events;
	wait_list.push_back(args_event);

	err_code = clEnqueueNDRangeKernel(C->command_queue(),
	                                  _kernel,
	                                  1,
	                                  NULL,
	                                  &_global_work_size,
	                                  &_local_work_size,
	                                  wait_list.size(),
	                                  wait_list.data(),
	                                  &event);
	if (err_code != CL_SUCCESS) {
		std::stringstream msg;
		msg << "Failure executing the tool \"" << name() << "\"." << std::endl;
		LOG(L_ERROR, msg.str());
		InputOutput::Logger::singleton()->printOpenCLError(err_code);
		throw std::runtime_error("OpenCL execution error");
	}
	err_code = clReleaseEvent(args_event);
	if (err_code != CL_SUCCESS) {
		std::stringstream msg;
		msg << "Failure releasing the args setter event in \"" << name()
			<< "\" tool." << std::endl;
		LOG(L_ERROR, msg.str());
		InputOutput::Logger::singleton()->printOpenCLError(err_code);
		throw std::runtime_error("OpenCL execution error");
	}

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
	if (vars->get(_var_name)->type().find('*') == std::string::npos) {
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
		flags << "-DT=uint";
	} else {
		std::string t = trimCopy(_var->type());
		t.pop_back(); // Remove the asterisk
		flags << "-DT=" << t;
	}
	kernel = compile_kernel(source.str(), "set", flags.str());
	err_code = clGetKernelWorkGroupInfo(kernel,
	                                    C->device(),
	                                    CL_KERNEL_WORK_GROUP_SIZE,
	                                    sizeof(size_t),
	                                    &_local_work_size,
	                                    NULL);
	if (err_code != CL_SUCCESS) {
		LOG(L_ERROR, "Failure querying the work group size.\n");
		InputOutput::Logger::singleton()->printOpenCLError(err_code);
		clReleaseKernel(kernel);
		throw std::runtime_error("OpenCL error");
	}
	if (_local_work_size < __CL_MIN_LOCALSIZE__) {
		LOG(L_ERROR, "insufficient local memory.\n");
		std::stringstream msg;
		msg << "\t" << _local_work_size
		    << " local work group size with __CL_MIN_LOCALSIZE__="
		    << __CL_MIN_LOCALSIZE__ << std::endl;
		LOG0(L_DEBUG, msg.str());
		throw std::runtime_error("OpenCL error");
	}

	_global_work_size = roundUp(_n, _local_work_size);
	_kernel = kernel;
	err_code = clSetKernelArg(kernel, 0, _var->typesize(), _var->get());
	if (err_code != CL_SUCCESS) {
		LOG(L_ERROR, "Failure sending the array argument\n");
		InputOutput::Logger::singleton()->printOpenCLError(err_code);
		throw std::runtime_error("OpenCL error");
	}
	err_code = clSetKernelArg(kernel, 1, sizeof(unsigned int), (void*)&_n);
	if (err_code != CL_SUCCESS) {
		LOG(L_ERROR, "Failure sending the array size argument\n");
		InputOutput::Logger::singleton()->printOpenCLError(err_code);
		throw std::runtime_error("OpenCL error");
	}
	if (!_data)
		return;
	err_code = clSetKernelArg(
	    kernel, 2, InputOutput::Variables::typeToBytes(_var->type()), _data);
	if (err_code != CL_SUCCESS) {
		LOG(L_ERROR, "Failure sending the value argument\n");
		InputOutput::Logger::singleton()->printOpenCLError(err_code);
		throw std::runtime_error("OpenCL error");
	}
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
		if (err_code != CL_SUCCESS) {
			std::stringstream msg;
			msg << "Failure setting the value to the tool \"" << name()
			    << "\"." << std::endl;
			LOG(L_ERROR, msg.str());
			InputOutput::Logger::singleton()->printOpenCLError(err_code);
			throw std::runtime_error("OpenCL error");
		}
	}

	if (_input != *(cl_mem*)_var->get()) {
		// For some reason the input variable has changed...
		err_code = clSetKernelArg(_kernel, 0, _var->typesize(), _var->get());
		if (err_code != CL_SUCCESS) {
			std::stringstream msg;
			msg << "Failure setting the variable \"" << _var->name()
			    << "\" to the tool \"" << name() << "\"." << std::endl;
			LOG(L_ERROR, msg.str());
			InputOutput::Logger::singleton()->printOpenCLError(err_code);
			throw std::runtime_error("OpenCL error");
		}

		_input = *(cl_mem*)_var->get();
	}
}

}
} // namespaces
