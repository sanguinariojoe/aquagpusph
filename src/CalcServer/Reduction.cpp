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

#include <AuxiliarMethods.h>
#include <InputOutput/Logger.h>
#include <CalcServer/Reduction.h>
#include <CalcServer.h>

namespace Aqua {
namespace CalcServer {

#ifndef DOXYGEN_SHOULD_SKIP_THIS
#include "CalcServer/Reduction.hcl"
#include "CalcServer/Reduction.cl"
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
		if (mem && (mem != _mems.front())) // The first element can't be removed
			clReleaseMemObject(mem);
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
}

cl_event
Reduction::_execute(const std::vector<cl_event> events_src)
{
	unsigned int i;
	cl_event event;
	cl_int err_code;
	CalcServer* C = CalcServer::singleton();
	InputOutput::Variables* vars = C->variables();

	setVariables();

	// We must execute several kernel in a sequential way, so we are just adding
	// more events to the wait list.
	std::vector<cl_event> events;
	std::copy(events_src.begin(), events_src.end(), std::back_inserter(events));
	for (i = 0; i < _kernels.size(); i++) {
		cl_uint num_events_in_wait_list = events.size();
		const cl_event* event_wait_list = events.size() ? events.data() : NULL;
		size_t _global_work_size = _global_work_sizes.at(i);
		size_t _local_work_size = _local_work_sizes.at(i);
		err_code = clEnqueueNDRangeKernel(C->command_queue(),
		                                  _kernels.at(i),
		                                  1,
		                                  NULL,
		                                  &_global_work_size,
		                                  &_local_work_size,
		                                  num_events_in_wait_list,
		                                  event_wait_list,
		                                  &event);
		if (err_code != CL_SUCCESS) {
			std::ostringstream msg;
			msg << "Failure executing the step " << i << " within the tool \""
			    << name() << "\"." << std::endl;
			LOG(L_ERROR, msg.str());
			InputOutput::Logger::singleton()->printOpenCLError(err_code);
			throw std::runtime_error("OpenCL execution error");
		}
		events.push_back(event);
	}

	// Get back the result
	cl_uint num_events_in_wait_list = events.size();
	const cl_event* event_wait_list = events.size() ? events.data() : NULL;
	err_code = clEnqueueReadBuffer(C->command_queue(),
	                               _mems.at(_mems.size() - 1),
	                               CL_TRUE,
	                               0,
	                               _output_var->typesize(),
	                               _output_var->get(),
	                               num_events_in_wait_list,
	                               event_wait_list,
	                               &event);
	if (err_code != CL_SUCCESS) {
		std::ostringstream msg;
		msg << "Failure reading back the result within the tool \"" << name()
		    << "\"." << std::endl;
		LOG(L_ERROR, msg.str());
		InputOutput::Logger::singleton()->printOpenCLError(err_code);
		throw std::runtime_error("OpenCL error");
	}

	// Release useless transactional events
	for (auto it = events.begin() + events_src.size(); it < events.end();
	     it++) {
		err_code = clReleaseEvent(*it);
		if (err_code != CL_SUCCESS) {
			std::ostringstream msg;
			msg << "Failure releasing transactional event in the tool \""
			    << name() << "\"." << std::endl;
			LOG(L_ERROR, msg.str());
			InputOutput::Logger::singleton()->printOpenCLError(err_code);
			throw std::runtime_error("OpenCL error");
		}
	}

	// Ensure that the variable is populated, this is in fact a blocking
	// operation
	_output_var->setEvent(event);
	vars->populate(_output_var);

	return event;
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
	if (vars->get(_input_name)->isScalar()) {
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
		msg << "Mismatching input and output types within the tool \"" << name()
		    << "\"." << std::endl;
		LOG(L_ERROR, msg.str());
		msg.str("");
		msg << "\tInput variable \"" << _input_var->name() << "\" is of type \""
		    << _input_var->type() << "\"." << std::endl;
		LOG0(L_DEBUG, msg.str());
		msg << "\tOutput variable \"" << _output_var->name()
		    << "\" is of type \"" << _output_var->type() << "\"." << std::endl;
		LOG0(L_DEBUG, msg.str());
		throw std::runtime_error("Invalid variable type");
	}

	// The scalar variable event is internally handled by the tool
	std::vector<InputOutput::Variable*> deps = { _input_var };
	setDependencies(deps);
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
	if (err_code != CL_SUCCESS) {
		LOG(L_ERROR, "Failure querying the work group size.\n");
		InputOutput::Logger::singleton()->printOpenCLError(err_code);
		clReleaseKernel(kernel);
		throw std::runtime_error("OpenCL error");
	}
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
	unsigned int n = _n.at(0);
	_n.clear();
	unsigned int i = 0;
	while (n > 1) {
		// Get work sizes
		_n.push_back(n);
		_local_work_sizes.push_back(local_size);
		_global_work_sizes.push_back(roundUp(n, local_size));
		_number_groups.push_back(_global_work_sizes.at(i) /
		                         _local_work_sizes.at(i));
		// Build the output memory object
		cl_mem output = NULL;
		output = clCreateBuffer(C->context(),
		                        CL_MEM_READ_WRITE,
		                        _number_groups.at(i) * data_size,
		                        NULL,
		                        &err_code);
		if (err_code != CL_SUCCESS) {
			std::stringstream msg;
			msg << "Failure allocating device memory in the tool \"" << name()
			    << "\"." << std::endl;
			LOG(L_ERROR, msg.str());
			InputOutput::Logger::singleton()->printOpenCLError(err_code);
			throw std::runtime_error("OpenCL allocation error");
		}
		allocatedMemory(_number_groups.at(i) * data_size + allocatedMemory());
		_mems.push_back(output);
		// Build the kernel
		kernel = compile_kernel(source.str(), "reduction", flags(local_size));
		_kernels.push_back(kernel);

		err_code =
		    clSetKernelArg(kernel, 0, sizeof(cl_mem), (void*)&(_mems.at(i)));
		if (err_code != CL_SUCCESS) {
			LOG(L_ERROR, "Failure sending input argument\n");
			InputOutput::Logger::singleton()->printOpenCLError(err_code);
			throw std::runtime_error("OpenCL error");
		}
		err_code = clSetKernelArg(
		    kernel, 1, sizeof(cl_mem), (void*)&(_mems.at(i + 1)));
		if (err_code != CL_SUCCESS) {
			LOG(L_ERROR, "Failure sending output argument\n");
			InputOutput::Logger::singleton()->printOpenCLError(err_code);
			throw std::runtime_error("OpenCL error");
		}
		err_code = clSetKernelArg(kernel, 2, sizeof(cl_uint), (void*)&(n));
		if (err_code != CL_SUCCESS) {
			LOG(L_ERROR, "Failure sending number of threads argument\n");
			InputOutput::Logger::singleton()->printOpenCLError(err_code);
			throw std::runtime_error("OpenCL error");
		}
		err_code = clSetKernelArg(kernel, 3, local_size * data_size, NULL);
		if (err_code != CL_SUCCESS) {
			LOG(L_ERROR, "Failure setting local memory\n");
			InputOutput::Logger::singleton()->printOpenCLError(err_code);
			throw std::runtime_error("OpenCL error");
		}
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
	if (err_code != CL_SUCCESS) {
		std::stringstream msg;
		msg << "Failure setting the input variable \"" << _input_var->name()
		    << "\" to the tool \"" << name() << "\"." << std::endl;
		LOG(L_ERROR, msg.str());
		InputOutput::Logger::singleton()->printOpenCLError(err_code);
		throw std::runtime_error("OpenCL error");
	}

	_input = *(cl_mem*)_input_var->get();
	_mems.at(0) = _input;
}

const std::string
Reduction::flags(const size_t local_size)
{
	std::ostringstream f;
	if (!_output_var->type().compare("unsigned int")) {
		// Spaces are not a good business into definitions passed as args
		f << "-DT=uint";
	} else {
		f << "-DT=" << _output_var->type();
	}
	f << " -DLOCAL_WORK_SIZE=" << local_size << "u";

	return f.str();
}

}
} // namespaces
