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
 * @brief Prefix sum
 * (See Aqua::CalcServer::Scan for details)
 * @note Hardcoded versions of the files CalcServer/Scan.cl.in and
 * CalcServer/Scan.hcl.in are internally included as a text array.
 */

#include <sstream>
#include "aquagpusph/AuxiliarMethods.hpp"
#include "aquagpusph/InputOutput/Logger.hpp"
#include "Scan.hpp"
#include "Kernel.hpp"
#include "CalcServer.hpp"

namespace Aqua {
namespace CalcServer {

#ifndef DOXYGEN_SHOULD_SKIP_THIS
#include "aquagpusph/CalcServer/Scan.hcl"
#include "aquagpusph/CalcServer/Scan.cl"
#endif
std::string SCAN_INC = xxd2string(Scan_hcl_in, Scan_hcl_in_len);
std::string SCAN_SRC = xxd2string(Scan_cl_in, Scan_cl_in_len);

Scan::Scan(const std::string name,
                     const std::string var_name,
                     bool once)
  : Tool(name, once)
  , _var_name(var_name)
  , _var(NULL)
  , _init(NULL)
  , _scan(NULL)
  , _sweepup(NULL)
  , _sweepdown(NULL)
{
}

Scan::~Scan()
{
	for (unsigned int i = 0; i < _input.size(); i++) {
		if (i != 0) {
			// The first input data does not belongs to us
			clReleaseMemObject(_input.at(i));
		}
		clReleaseMemObject(_partial.at(i));
		clReleaseMemObject(_flags.at(i));
	}
	if (_init)
		clReleaseKernel(_init);
	if (_scan)
		clReleaseKernel(_scan);
	if (_sweepup)
		clReleaseKernel(_sweepup);
	if (_sweepdown)
		clReleaseKernel(_sweepdown);
}

void
Scan::setup()
{
	std::ostringstream msg;
	msg << "Loading the tool \"" << name() << "\"..." << std::endl;
	LOG(L_INFO, msg.str());

	Tool::setup();
	variables();

	_n = _var->size() / InputOutput::Variables::typeToBytes(_var->type());
	setupOpenCL();

	std::vector<Profile*> profilers;
	profilers.push_back(new EventProfile("init", this));
	for (size_t i = 0; i < _k.size(); i++) {
		std::stringstream reduction_stage;
		reduction_stage << "step " << i + 1 << "/" << _k.size();
		profilers.push_back(new EventProfile(reduction_stage.str(), this));
	}
	Profiler::substages(profilers);
}

cl_event
Scan::_execute(const std::vector<cl_event> events)
{
	cl_event event, out_event;
	cl_int err_code;
	auto C = CalcServer::singleton();

	// We first execute the initialization so have our dummy flags
	const size_t work_group_size = _WGS;
	const size_t global_work_size = roundUp(_n, work_group_size);
	err_code = clEnqueueNDRangeKernel(C->command_queue(),
		                              _init,
		                              1,
		                              NULL,
		                              &global_work_size,
		                              &work_group_size,
		                              0,
		                              NULL,
		                              &event);
	CHECK_OCL_OR_THROW(err_code,
		std::string("Failure executing init in tool \"") + name() + "\".");
	err_code = clFlush(C->command_queue());
	CHECK_OCL_OR_THROW(err_code,
		std::string("Failure flushing the command queue in tool \"") +
			name() + "\" at init");

	auto profiler =
		dynamic_cast<EventProfile*>(Profiler::substages().at(0));
	profiler->start(event);
	profiler->end(event);

	// Now we can execute the scanning recursively
	std::vector<cl_event> wait_events = events;
	wait_events.push_back(event);
	out_event = scan(0, wait_events);

	// Replace the event by the new one
	err_code = clReleaseEvent(event);
	CHECK_OCL_OR_THROW(
		err_code,
		std::string("Failure releasing the init event in tool \"") +
			name() + "\".");
	event = out_event;

	return event;
}

void
Scan::variables()
{
	InputOutput::Variables* vars = CalcServer::singleton()->variables();
	if (!vars->get(_var_name)) {
		std::stringstream msg;
		msg << "The tool \"" << name()
		    << "\" is asking the undeclared input variable \"" << _var_name
		    << "\"." << std::endl;
		LOG(L_ERROR, msg.str());
		throw std::runtime_error("Invalid variable");
	}
	if (!vars->get(_var_name)->isArray()) {
		std::stringstream msg;
		msg << "The tool \"" << name() << "\" is asking the input variable \""
		    << _var_name << "\", which is a scalar." << std::endl;
		LOG(L_ERROR, msg.str());
		throw std::runtime_error("Invalid variable type");
	}
	_var = (InputOutput::ArrayVariable*)vars->get(_var_name);

	setOutputDependencies({ _var });
}

void
Scan::setupOpenCL()
{
	cl_int err_code;
	CalcServer* C = CalcServer::singleton();
	InputOutput::Variables* vars = C->variables();

	// Get the elements data size to can allocate memory later
	cl_ulong data_size = vars->typeToBytes(_var->type());

	std::ostringstream source;
	source << SCAN_INC << std::endl << SCAN_SRC;

	// Check that we have local memory enough to operate
	cl_ulong local_mem;
	err_code = clGetDeviceInfo(C->device(),
	                           CL_DEVICE_LOCAL_MEM_SIZE,
	                           sizeof(cl_ulong),
	                           &local_mem,
	                           NULL);
	CHECK_OCL_OR_THROW(
	    err_code,
	    std::string("Failure querying the available local memory in tool \"") +
	        name() + "\".");
	cl_ulong max_wgs = local_mem / (2 * data_size);
	if (max_wgs < _WGS) {
		LOG(L_ERROR, "insufficient local memory.\n");
		std::stringstream msg;
		msg << "\tA maximum work group size of " << max_wgs
		    << " is tolerated by the device, but "
		    << _WGS << " is required by \"" << name() << "\"" << std::endl;
		LOG0(L_DEBUG, msg.str());
		throw std::runtime_error("OpenCL error");
	}

	// Compile the kernels
	_init = compile_kernel(source.str(), "init", flags());
	_scan = compile_kernel(source.str(), "scan_wrapper", flags());
	_sweepup = compile_kernel(source.str(), "sweep_up_subarrays", flags());
	_sweepdown = compile_kernel(source.str(), "sweep_down_subarrays", flags());

	// Create the inputs of the first stage
	_input.push_back(*(cl_mem*)_var->get());
	cl_mem partial = clCreateBuffer(C->context(),
	                                CL_MEM_READ_WRITE,
	                                _n * sizeof(cl_int),
	                                NULL,
	                                &err_code);
	CHECK_OCL_OR_THROW(
		err_code,
		std::string("Failure allocating ") +
			std::to_string(_n * sizeof(cl_int)) +
			" bytes on the device for partials on  tool \"" + name() + "\".");
	_partial.push_back(partial);
	cl_mem flags = clCreateBuffer(C->context(),
	                              CL_MEM_READ_WRITE,
	                              _n * sizeof(cl_int),
	                              NULL,
	                              &err_code);
	CHECK_OCL_OR_THROW(
		err_code,
		std::string("Failure allocating ") +
			std::to_string(_n * sizeof(cl_int)) +
			" bytes on the device for flags on tool \"" + name() + "\".");
	_flags.push_back(flags);
	allocatedMemory(2 * _n * sizeof(cl_int));
	_k.push_back(_n);

	// Setup the subsequent stages
	const size_t m = 2 * _WGS;
	size_t k = roundUp(_n, m) / m;
	std::stringstream msg;
	msg << "\t" << _n << " elements scanned to " << k << " groups"
	    << std::endl;
	LOG(L_DEBUG, msg.str());
	while (k > 1) {
		cl_mem input_sum, partial_sum, flags_sum;
		input_sum = clCreateBuffer(C->context(),
		                           CL_MEM_READ_WRITE,
		                           k * data_size,
		                           NULL,
		                           &err_code);
		CHECK_OCL_OR_THROW(
			err_code,
			std::string("Failure allocating ") +
				std::to_string(k * data_size) +
				" bytes on the device for input_sum on  tool \"" + name() + "\".");
		partial_sum = clCreateBuffer(C->context(),
		                             CL_MEM_READ_WRITE,
		                             k * sizeof(cl_int),
		                             NULL,
		                             &err_code);
		CHECK_OCL_OR_THROW(
			err_code,
			std::string("Failure allocating ") +
				std::to_string(k * sizeof(cl_int)) +
				" bytes on the device for partials_sum on  tool \"" + name() + "\".");
		flags_sum = clCreateBuffer(C->context(),
		                           CL_MEM_READ_WRITE,
		                           k * sizeof(cl_int),
		                           NULL,
		                           &err_code);
		CHECK_OCL_OR_THROW(
			err_code,
			std::string("Failure allocating ") +
				std::to_string(k * sizeof(cl_int)) +
				" bytes on the device for flags_sum on tool \"" + name() + "\".");
		allocatedMemory(k * (data_size + 2 * sizeof(cl_int)));

		_input.push_back(input_sum);
		_partial.push_back(partial_sum);
		_flags.push_back(flags_sum);
		_k.push_back(k);

		std::stringstream msg;
		msg << "\t" << k << " elements scanned to ";
		k = roundUp(k, m) / m;
		msg << k << " groups" << std::endl;
		LOG(L_DEBUG, msg.str());
	}
	// Add the last output number of groups
	_k.push_back(1);

	// Send the 1st stage input flags to the initializer
	err_code = clSetKernelArg(_init,
	                          0,
	                          sizeof(cl_mem),
	                          (void*)&(_partial.at(0)));
	CHECK_OCL_OR_THROW(
		err_code,
		std::string("Failure sending part argument for init in tool \"") +
			name() + "\".");
	err_code = clSetKernelArg(_init,
	                          1,
	                          sizeof(cl_mem),
	                          (void*)&(_flags.at(0)));
	CHECK_OCL_OR_THROW(
		err_code,
		std::string("Failure sending flag argument for init in tool \"") +
			name() + "\".");
	err_code = C->setKernelSizeArg(_init, 2, _n);
	CHECK_OCL_OR_THROW(
		err_code,
		std::string(
			"Failure sending number n argument for init in tool \"") +
			name() + "\".");

	// Set the local memory to the scanning kernels
	err_code = clSetKernelArg(_scan, 3, 2 * _WGS * data_size, NULL);
	CHECK_OCL_OR_THROW(
		err_code,
		std::string(
			"Failure setting x local memory for scan_wrapper in tool \"") +
			name() + "\".");
	err_code = clSetKernelArg(_scan, 4, 2 * _WGS * sizeof(int), NULL);
	CHECK_OCL_OR_THROW(
		err_code,
		std::string(
			"Failure setting p local memory for scan_wrapper in tool \"") +
			name() + "\".");
	err_code = clSetKernelArg(_scan, 5, 2 * _WGS * sizeof(int), NULL);
	CHECK_OCL_OR_THROW(
		err_code,
		std::string(
			"Failure setting f local memory for scan_wrapper in tool \"") +
			name() + "\".");

	err_code = clSetKernelArg(_sweepup, 6, 2 * _WGS * data_size, NULL);
	CHECK_OCL_OR_THROW(
		err_code,
		std::string(
			"Failure setting x local memory for sweep_up_subarrays in tool \"") +
			name() + "\".");
	err_code = clSetKernelArg(_sweepup, 7, 2 * _WGS * sizeof(int), NULL);
	CHECK_OCL_OR_THROW(
		err_code,
		std::string(
			"Failure setting p local memory for sweep_up_subarrays in tool \"") +
			name() + "\".");
	err_code = clSetKernelArg(_sweepup, 8, 2 * _WGS * sizeof(int), NULL);
	CHECK_OCL_OR_THROW(
		err_code,
		std::string(
			"Failure setting f local memory for sweep_up_subarrays in tool \"") +
			name() + "\".");

	err_code = clSetKernelArg(_sweepdown, 6, 2 * _WGS * data_size, NULL);
	CHECK_OCL_OR_THROW(
		err_code,
		std::string(
			"Failure setting x local memory for sweep_down_subarrays in tool \"") +
			name() + "\".");
	err_code = clSetKernelArg(_sweepdown, 7, 2 * _WGS * sizeof(int), NULL);
	CHECK_OCL_OR_THROW(
		err_code,
		std::string(
			"Failure setting p local memory for sweep_down_subarrays in tool \"") +
			name() + "\".");
	err_code = clSetKernelArg(_sweepdown, 8, 2 * _WGS * sizeof(int), NULL);
	CHECK_OCL_OR_THROW(
		err_code,
		std::string(
			"Failure setting f local memory for sweep_down_subarrays in tool \"") +
			name() + "\".");
}

cl_event
Scan::scan(const unsigned int stage,
           const std::vector<cl_event> events)
{
	// Get the input number of elements and the number of output groups
	const size_t n = _k.at(stage);
	const size_t k = _k.at(stage + 1);
	const size_t work_group_size = _WGS;
	cl_event event;
	cl_int err_code;
	CalcServer* C = CalcServer::singleton();

	if (k == 1) {
		// This is the last step
		err_code = clSetKernelArg(_scan,
		                          0,
		                          sizeof(cl_mem),
		                          (void*)&(_input.at(stage)));
		CHECK_OCL_OR_THROW(
			err_code,
			std::string(
				"Failure sending data argument for scan_wrapper in tool \"") +
				name() + "\".");
		err_code = clSetKernelArg(_scan,
		                          1,
		                          sizeof(cl_mem),
		                          (void*)&(_partial.at(stage)));
		CHECK_OCL_OR_THROW(
			err_code,
			std::string(
				"Failure sending part argument for scan_wrapper in tool \"") +
				name() + "\".");
		err_code = clSetKernelArg(_scan,
		                          2,
		                          sizeof(cl_mem),
		                          (void*)&(_flags.at(stage)));
		CHECK_OCL_OR_THROW(
			err_code,
			std::string(
				"Failure sending flag argument for scan_wrapper in tool \"") +
				name() + "\".");
		err_code = C->setKernelSizeArg(_scan, 6, n);
		CHECK_OCL_OR_THROW(
			err_code,
			std::string(
				"Failure sending number n argument for scan_wrapper in tool \"") +
				name() + "\".");

		err_code = clEnqueueNDRangeKernel(C->command_queue(),
		                                  _scan,
		                                  1,
		                                  NULL,
		                                  &work_group_size,
		                                  &work_group_size,
		                                  events.size(),
		                                  events.data(),
		                                  &event);
		CHECK_OCL_OR_THROW(err_code,
			std::string("Failure executing scan_wrapper in tool \"") +
			            name() + "\".");
		err_code = clFlush(C->command_queue());
		CHECK_OCL_OR_THROW(err_code,
			std::string("Failure flushing the command queue in tool \"") +
				name() + "\" at pass " + std::to_string(stage) + ".");

		auto profiler =
			dynamic_cast<EventProfile*>(Profiler::substages().back());
		profiler->start(event);
		profiler->end(event);
		return event;
	}

	// We are on an intermediate step. Let's start doing the sweep-up
	err_code = clSetKernelArg(_sweepup,
	                          0,
	                          sizeof(cl_mem),
	                          (void*)&(_input.at(stage)));
	CHECK_OCL_OR_THROW(
		err_code,
		std::string(
			"Failure sending data argument for sweep_up_subarrays in tool \"") +
			name() + "\".");
	err_code = clSetKernelArg(_sweepup,
	                          1,
	                          sizeof(cl_mem),
	                          (void*)&(_partial.at(stage)));
	CHECK_OCL_OR_THROW(
		err_code,
		std::string(
			"Failure sending part argument for sweep_up_subarrays in tool \"") +
			name() + "\".");
	err_code = clSetKernelArg(_sweepup,
	                          2,
	                          sizeof(cl_mem),
	                          (void*)&(_flags.at(stage)));
	CHECK_OCL_OR_THROW(
		err_code,
		std::string(
			"Failure sending flag argument for sweep_up_subarrays in tool \"") +
			name() + "\".");
	err_code = clSetKernelArg(_sweepup,
	                          3,
	                          sizeof(cl_mem),
	                          (void*)&(_input.at(stage + 1)));
	CHECK_OCL_OR_THROW(
		err_code,
		std::string(
			"Failure sending data_sum argument for sweep_up_subarrays in tool \"") +
			name() + "\".");
	err_code = clSetKernelArg(_sweepup,
	                          4,
	                          sizeof(cl_mem),
	                          (void*)&(_partial.at(stage + 1)));
	CHECK_OCL_OR_THROW(
		err_code,
		std::string(
			"Failure sending part_sum argument for sweep_up_subarrays in tool \"") +
			name() + "\".");
	err_code = clSetKernelArg(_sweepup,
	                          5,
	                          sizeof(cl_mem),
	                          (void*)&(_flags.at(stage + 1)));
	CHECK_OCL_OR_THROW(
		err_code,
		std::string(
			"Failure sending flag_sum argument for sweep_up_subarrays in tool \"") +
			name() + "\".");
	err_code = C->setKernelSizeArg(_sweepup, 9, n);
	CHECK_OCL_OR_THROW(
		err_code,
		std::string(
			"Failure sending number n argument for sweep_up_subarrays in tool \"") +
			name() + "\".");

	const size_t global_work_size = k * work_group_size;
	err_code = clEnqueueNDRangeKernel(C->command_queue(),
	                                  _sweepup,
	                                  1,
	                                  NULL,
	                                  &global_work_size,
	                                  &work_group_size,
	                                  events.size(),
	                                  events.data(),
	                                  &event);
	CHECK_OCL_OR_THROW(err_code,
		std::string("Failure executing sweep_up_subarrays in tool \"") +
			name() + "\" at pass " + std::to_string(stage) + ".");
	err_code = clFlush(C->command_queue());
	CHECK_OCL_OR_THROW(err_code,
		std::string("Failure flushing the command queue in tool \"") +
			name() + "\" at pass " + std::to_string(stage) + ".");

	auto profiler =
		dynamic_cast<EventProfile*>(Profiler::substages().at(stage + 1));
	profiler->start(event);

	// Now we recursively call to make a new scan stage
	cl_event event_out;
	event_out = scan(stage + 1, {event});

	// Replace the event by the new one
	err_code = clReleaseEvent(event);
	CHECK_OCL_OR_THROW(
		err_code,
		std::string("Failure releasing the subscan event for the step ") +
			std::to_string(stage) + " in tool \"" + name() + "\".");
	event = event_out;

	// And now we can run the sweep-down
	err_code = clSetKernelArg(_sweepdown,
	                          0,
	                          sizeof(cl_mem),
	                          (void*)&(_input.at(stage)));
	CHECK_OCL_OR_THROW(
		err_code,
		std::string(
			"Failure sending data argument for sweep_down_subarrays in tool \"") +
			name() + "\".");
	err_code = clSetKernelArg(_sweepdown,
	                          1,
	                          sizeof(cl_mem),
	                          (void*)&(_partial.at(stage)));
	CHECK_OCL_OR_THROW(
		err_code,
		std::string(
			"Failure sending part argument for sweep_down_subarrays in tool \"") +
			name() + "\".");
	err_code = clSetKernelArg(_sweepdown,
	                          2,
	                          sizeof(cl_mem),
	                          (void*)&(_flags.at(stage)));
	CHECK_OCL_OR_THROW(
		err_code,
		std::string(
			"Failure sending flag argument for sweep_down_subarrays in tool \"") +
			name() + "\".");
	err_code = clSetKernelArg(_sweepdown,
	                          3,
	                          sizeof(cl_mem),
	                          (void*)&(_input.at(stage + 1)));
	CHECK_OCL_OR_THROW(
		err_code,
		std::string(
			"Failure sending data_sum argument for sweep_down_subarrays in tool \"") +
			name() + "\".");
	err_code = clSetKernelArg(_sweepdown,
	                          4,
	                          sizeof(cl_mem),
	                          (void*)&(_partial.at(stage + 1)));
	CHECK_OCL_OR_THROW(
		err_code,
		std::string(
			"Failure sending part_sum argument for sweep_down_subarrays in tool \"") +
			name() + "\".");
	err_code = clSetKernelArg(_sweepdown,
	                          5,
	                          sizeof(cl_mem),
	                          (void*)&(_flags.at(stage + 1)));
	CHECK_OCL_OR_THROW(
		err_code,
		std::string(
			"Failure sending flag_sum argument for sweep_down_subarrays in tool \"") +
			name() + "\".");
	err_code = C->setKernelSizeArg(_sweepdown, 9, n);
	CHECK_OCL_OR_THROW(
		err_code,
		std::string(
			"Failure sending number n argument for sweep_down_subarrays in tool \"") +
			name() + "\".");

	err_code = clEnqueueNDRangeKernel(C->command_queue(),
	                                  _sweepdown,
	                                  1,
	                                  NULL,
	                                  &global_work_size,
	                                  &work_group_size,
	                                  1,
	                                  &event,
	                                  &event_out);
	CHECK_OCL_OR_THROW(err_code,
		std::string("Failure executing sweep_down_subarrays in tool \"") +
			name() + "\" at pass " + std::to_string(stage) + ".");
	err_code = clFlush(C->command_queue());
	CHECK_OCL_OR_THROW(err_code,
		std::string("Failure flushing the command queue in tool \"") +
			name() + "\" at pass " + std::to_string(stage) + ".");

	// Replace the event by the new one
	err_code = clReleaseEvent(event);
	CHECK_OCL_OR_THROW(
		err_code,
		std::string("Failure releasing the subscan event for the step ") +
			std::to_string(stage) + " in tool \"" + name() + "\".");
	event = event_out;
	profiler->end(event);

	return event;
}

const std::string
Scan::flags()
{
	auto var_type = replaceAllCopy(_var->type(), "*", "");
	if (!var_type.compare("unsigned int")) {
		// Spaces are not a good business into definitions passed as args
		return " -DT=uint";
	} else if (!var_type.compare("unsigned long")) {
		// Spaces are not a good business into definitions passed as args
		return " -DT=ulong";
	}
	
	return std::string(" -DT=") + var_type;
}

}
} // namespaces
