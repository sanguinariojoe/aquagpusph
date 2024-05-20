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
 * @brief Reductions, like scans, prefix sums, maximum or minimum, etc...
 * (See Aqua::CalcServer::Reduction for details)
 * @note Hardcoded versions of the files CalcServer/Reduction.cl.in and
 * CalcServer/Reduction.hcl.in are internally included as a text array.
 */

#include <sstream>
#include "aquagpusph/AuxiliarMethods.hpp"
#include "aquagpusph/InputOutput/Logger.hpp"
#include "Reduction.hpp"
#include "Kernel.hpp"
#include "CalcServer.hpp"

namespace Aqua {
namespace CalcServer {

#ifndef DOXYGEN_SHOULD_SKIP_THIS
#include "aquagpusph/CalcServer/Reduction.hcl"
#include "aquagpusph/CalcServer/Reduction.cl"
#endif
std::string REDUCTION_INC = xxd2string(Reduction_hcl_in, Reduction_hcl_in_len);
std::string REDUCTION_SRC = xxd2string(Reduction_cl_in, Reduction_cl_in_len);

Reduction::Reduction(const std::string name,
                     const std::string input_name,
                     const std::string output_name,
                     const std::string operation,
                     const std::string null_val,
                     bool once)
  : Tool(name, once)
  , _input_name(input_name)
  , _output_name(output_name)
  , _operation(operation)
  , _null_val(null_val)
  , _input_var(NULL)
  , _output_var(NULL)
  , _input(NULL)
{
}

Reduction::~Reduction()
{
	for (auto mem : _mems) {
		if (mem && (mem != _mems.front())) {
			// The first element can't be removed
			clReleaseMemObject(mem);
		}
	}
	_mems.clear();
	for (auto kernel : _kernels) {
		if (kernel)
			clReleaseKernel(kernel);
	}
	_kernels.clear();
	_global_work_sizes.clear();
	_local_work_sizes.clear();
}

void
Reduction::setup()
{
	std::ostringstream msg;
	msg << "Loading the tool \"" << name() << "\"..." << std::endl;
	LOG(L_INFO, msg.str());

	Tool::setup();
	variables();

	_mems.push_back(*(cl_mem*)_input_var->get());
	_input = *(cl_mem*)_input_var->get();
	size_t n = _input_var->size() /
	           InputOutput::Variables::typeToBytes(_input_var->type());
	_n.push_back(n);
	setupOpenCL();
	std::vector<Profile*> profilers;
	for (size_t i = 0; i < _kernels.size(); i++) {
		std::stringstream reduction_stage;
		reduction_stage << "step " << i + 1 << "/" << _kernels.size();
		profilers.push_back(new EventProfile(reduction_stage.str(), this));
	}
	profilers.push_back(new EventProfile("read", this));
	Profiler::substages(profilers);
}

/** @brief Callback called when Aqua::CalcServer::Reduction already set the new
 * value on the output value.
 *
 * This function is just populating the variable
 * @param event The triggering event
 * @param event_command_status CL_COMPLETE upon all dependencies successfully
 * fulfilled. A negative integer if one or mor dependencies failed.
 * @param user_data A casted pointer to the Aqua::CalcServer::ScalarExpression
 * tool (or the inherited one)
 */
void CL_CALLBACK
populator(cl_event event, cl_int event_command_status, void* user_data)
{
	cl_int err_code;
	clReleaseEvent(event);
	auto tool = (Reduction*)user_data;
	if (event_command_status != CL_COMPLETE) {
		std::stringstream msg;
		msg << "Skipping \"" << tool->name()
		    << "\" variable population due to dependency errors." << std::endl;
		LOG(L_WARNING, msg.str());
		clSetUserEventStatus(tool->getUserEvent(), event_command_status);
		clReleaseEvent(tool->getUserEvent());
		return;
	}

	auto var = tool->getOutputDependencies()[0];
	auto vars = CalcServer::singleton()->variables();
	vars->populate(var);
	err_code = clSetUserEventStatus(tool->getUserEvent(), CL_COMPLETE);
	CHECK_OCL_OR_THROW(err_code,
	                   "Failure setting the variable population event");
	clReleaseEvent(tool->getUserEvent());
}

cl_event
Reduction::_execute(const std::vector<cl_event> events)
{
	cl_event event, out_event;
	cl_int err_code;
	auto C = CalcServer::singleton();
	auto vars = C->variables();

	// We must execute several kernels in a sequential way. The process is
	// though seeded by a reduction which just needs to read the input array
	event = _input_var->getWritingEvent();
	err_code = clRetainEvent(event);
	CHECK_OCL_OR_THROW(
	    err_code,
	    std::string("Failure retaining the input event of tool \"") + name() +
	        "\".");
	for (size_t i = 0; i < _kernels.size(); i++) {
		const size_t global_work_size = _global_work_sizes.at(i);
		const size_t local_work_size = _local_work_sizes.at(i);
		err_code = clEnqueueNDRangeKernel(C->command_queue(),
		                                  _kernels.at(i),
		                                  1,
		                                  NULL,
		                                  &global_work_size,
		                                  &local_work_size,
		                                  1,
		                                  &event,
		                                  &out_event);
		CHECK_OCL_OR_THROW(err_code,
		                   std::string("Failure executing the step ") +
		                       std::to_string(i) + " in tool \"" + name() +
		                       "\".");

		auto profiler =
		    dynamic_cast<EventProfile*>(Profiler::substages().at(i));
		profiler->start(out_event);
		profiler->end(out_event);

		// Replace the event by the new one
		err_code = clReleaseEvent(event);
		CHECK_OCL_OR_THROW(
		    err_code,
		    std::string("Failure releasing the input event for the step ") +
		        std::to_string(i) + " in tool \"" + name() + "\".");
		event = out_event;
	}

	// Get back the result
	std::vector<cl_event> read_events = events;
	read_events.push_back(event);
	err_code = clEnqueueReadBuffer(C->command_queue(),
	                               _mems.at(_mems.size() - 1),
	                               CL_FALSE,
	                               0,
	                               _output_var->typesize(),
	                               _output_var->get(),
	                               read_events.size(),
	                               read_events.data(),
	                               &out_event);
	CHECK_OCL_OR_THROW(
	    err_code,
	    std::string("Failure reading back the result in tool \"") + name() +
	        "\".");

	auto profiler = dynamic_cast<EventProfile*>(Profiler::substages().back());
	profiler->start(out_event);

	err_code = clReleaseEvent(event);
	CHECK_OCL_OR_THROW(
	    err_code,
	    std::string("Failure releasing the transactional event in tool \"") +
	        name() + "\".");
	event = out_event;

	// Although the variable will be correctly set when clEnqueueReadBuffer()
	// finish its job, we want to populate it so other parts of the code are
	// aware of the change, like the math solver Aqua::Tokenizer
	_user_event = clCreateUserEvent(C->context(), &err_code);
	CHECK_OCL_OR_THROW(
	    err_code,
	    std::string("Failure creating the user event for tool \"") + name() +
	        "\".");

	err_code = clRetainEvent(_user_event);
	CHECK_OCL_OR_THROW(
	    err_code,
	    std::string("Failure retaining the user event for tool \"") + name() +
	        "\".");
	err_code = clSetEventCallback(event, CL_COMPLETE, populator, this);
	CHECK_OCL_OR_THROW(
	    err_code,
	    std::string("Failure registering the solver callback in tool \"") +
	        name() + "\".");

	// It is not really a good business to return an user event, which is kind
	// of special. e.g. it cannot be profiled. So we create a mark on top of it
	err_code = clEnqueueMarkerWithWaitList(
	    C->command_queue(), 1, &_user_event, &event);
	CHECK_OCL_OR_THROW(
	    err_code,
	    std::string("Failure setting the the output event for tool \"") +
	        name() + "\".");
	profiler->end(event);

	return _user_event;
}

void
Reduction::variables()
{
	InputOutput::Variables* vars = CalcServer::singleton()->variables();
	if (!vars->get(_input_name)) {
		std::stringstream msg;
		msg << "The tool \"" << name()
		    << "\" is asking the undeclared input variable \"" << _input_name
		    << "\"." << std::endl;
		LOG(L_ERROR, msg.str());
		throw std::runtime_error("Invalid variable");
	}
	if (!vars->get(_input_name)->isArray()) {
		std::stringstream msg;
		msg << "The tool \"" << name() << "\" is asking the input variable \""
		    << _input_name << "\", which is a scalar." << std::endl;
		LOG(L_ERROR, msg.str());
		throw std::runtime_error("Invalid variable type");
	}
	_input_var = (InputOutput::ArrayVariable*)vars->get(_input_name);
	if (!vars->get(_output_name)) {
		std::stringstream msg;
		msg << "The tool \"" << name()
		    << "\" is asking the undeclared output variable \"" << _output_name
		    << "\"." << std::endl;
		LOG(L_ERROR, msg.str());
		throw std::runtime_error("Invalid variable");
	}
	if (vars->get(_output_name)->isArray()) {
		std::stringstream msg;
		msg << "The tool \"" << name() << "\" is asking the output variable \""
		    << _output_name << "\", which is an array." << std::endl;
		LOG(L_ERROR, msg.str());
		throw std::runtime_error("Invalid variable type");
	}
	_output_var = vars->get(_output_name);
	if (!vars->isSameType(_input_var->type(), _output_var->type())) {
		std::stringstream msg;
		msg << "Mismatching input and output types within the tool \""
		    << name() << "\"." << std::endl;
		LOG(L_ERROR, msg.str());
		msg.str("");
		msg << "\tInput variable \"" << _input_var->name()
		    << "\" is of type \"" << _input_var->type() << "\"." << std::endl;
		LOG0(L_DEBUG, msg.str());
		msg << "\tOutput variable \"" << _output_var->name()
		    << "\" is of type \"" << _output_var->type() << "\"." << std::endl;
		LOG0(L_DEBUG, msg.str());
		throw std::runtime_error("Invalid variable type");
	}

	// The scalar variable event is internally handled by the tool
	setDependencies({ _input_var }, { _output_var });
}

void
Reduction::setupOpenCL()
{
	size_t data_size, local_size, max_local_size;
	cl_int err_code;
	cl_kernel kernel;
	CalcServer* C = CalcServer::singleton();
	InputOutput::Variables* vars = C->variables();

	// Get the elements data size to can allocate local memory later
	data_size = vars->typeToBytes(_input_var->type());

	std::ostringstream source;
	source << REDUCTION_INC << " #define IDENTITY " << _null_val << std::endl;
	source << "T reduce(T a, T b) " << std::endl;
	source << "{ " << std::endl;
	source << "    T c; " << std::endl;
	source << _operation << ";" << std::endl;
	source << "    return c; " << std::endl;
	source << "} " << std::endl;
	source << REDUCTION_SRC;

	// Starts a dummy kernel in order to study the local size that can be used
	local_size = __CL_MAX_LOCALSIZE__;
	kernel = compile_kernel(source.str(), "reduction", flags(local_size));
	err_code = clGetKernelWorkGroupInfo(kernel,
	                                    C->device(),
	                                    CL_KERNEL_WORK_GROUP_SIZE,
	                                    sizeof(size_t),
	                                    &max_local_size,
	                                    NULL);
	CHECK_OCL_OR_THROW(
	    err_code,
	    std::string("Failure querying the work group size in tool \"") +
	        name() + "\".");
	if (max_local_size < __CL_MIN_LOCALSIZE__) {
		LOG(L_ERROR, "insufficient local memory.\n");
		std::stringstream msg;
		msg << "\t" << max_local_size
		    << " local work group size with __CL_MIN_LOCALSIZE__="
		    << __CL_MIN_LOCALSIZE__ << std::endl;
		LOG0(L_DEBUG, msg.str());
		throw std::runtime_error("OpenCL error");
	}
	local_size = max_local_size;
	if (!isPowerOf2(local_size)) {
		local_size = nextPowerOf2(local_size) / 2;
	}
	clReleaseKernel(kernel);

	// Now we can start a loop while the amount of reduced data is greater than
	// one
	size_t n = _n.at(0);
	_n.clear();
	size_t i = 0;
	while (n > 1) {
		// Get work sizes
		_n.push_back(n);
		_local_work_sizes.push_back(local_size);
		_global_work_sizes.push_back(roundUp<size_t>(n, local_size));
		_number_groups.push_back(_global_work_sizes.at(i) /
		                         _local_work_sizes.at(i));
		// Build the output memory object
		cl_mem output = NULL;
		output = clCreateBuffer(C->context(),
		                        CL_MEM_READ_WRITE,
		                        _number_groups.at(i) * data_size,
		                        NULL,
		                        &err_code);
		CHECK_OCL_OR_THROW(
		    err_code,
		    std::string("Failure allocating ") +
		        std::to_string(_number_groups.at(i) * data_size) +
		        " bytes on the device for tool \"" + name() + "\".");
		allocatedMemory(_number_groups.at(i) * data_size + allocatedMemory());
		_mems.push_back(output);
		// Build the kernel
		kernel = compile_kernel(source.str(), "reduction", flags(local_size));
		_kernels.push_back(kernel);

		err_code =
		    clSetKernelArg(kernel, 0, sizeof(cl_mem), (void*)&(_mems.at(i)));
		CHECK_OCL_OR_THROW(
		    err_code,
		    std::string("Failure sending input argument in tool \"") + name() +
		        "\".");
		err_code = clSetKernelArg(
		    kernel, 1, sizeof(cl_mem), (void*)&(_mems.at(i + 1)));
		CHECK_OCL_OR_THROW(
		    err_code,
		    std::string("Failure sending output argument in tool \"") + name() +
		        "\".");
		err_code = C->setKernelSizeArg(kernel, 2, n);
		CHECK_OCL_OR_THROW(
		    err_code,
		    std::string(
		        "Failure sending number of threads argument in tool \"") +
		        name() + "\".");
		err_code = clSetKernelArg(kernel, 3, local_size * data_size, NULL);
		CHECK_OCL_OR_THROW(
		    err_code,
		    std::string("Failure setting local memory in tool \"") + name() +
		        "\".");
		// Setup next step
		std::stringstream msg;
		msg << "\tStep " << i << ", " << n << " elements reduced to "
		    << _number_groups.at(i) << std::endl;
		LOG(L_DEBUG, msg.str());
		n = _number_groups.at(i);
		i++;
	}
}

void
Reduction::setVariables()
{
	cl_int err_code;

	if (_input == *(cl_mem*)_input_var->get()) {
		return;
	}

	err_code = clSetKernelArg(
	    _kernels.at(0), 0, _input_var->typesize(), _input_var->get());
	CHECK_OCL_OR_THROW(err_code,
	                   std::string("Failure setting the input variable \"") +
	                       _input_var->name() + "\" in tool \"" + name() +
	                       "\".");

	_input = *(cl_mem*)_input_var->get();
	_mems.at(0) = _input;
}

const std::string
Reduction::flags(const size_t local_size)
{
	std::ostringstream f;
	if (!_output_var->type().compare("unsigned int")) {
		// Spaces are not a good business into definitions passed as args
		f << " -DT=uint";
	} else if (!_output_var->type().compare("unsigned long")) {
		// Spaces are not a good business into definitions passed as args
		f << " -DT=ulong";
	} else {
		f << " -DT=" << _output_var->type();
	}
	f << " -DLOCAL_WORK_SIZE=" << local_size << "u";

	return f.str();
}

}
} // namespaces
