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
 * @brief Methods to perform a radix sort using the GPU (or any device
 * supported by OpenCL).
 * (See Aqua::CalcServer::Sort for details)
 * @note Hardcoded versions of the files CalcServer/Sort.cl.in and
 * CalcServer/Sort.hcl.in are internally included as a text array.
 */

#include "aquagpusph/AuxiliarMethods.hpp"
#include "aquagpusph/InputOutput/Logger.hpp"
#include "Sort.hpp"

namespace Aqua {
namespace CalcServer {

#ifndef DOXYGEN_SHOULD_SKIP_THIS
#include "aquagpusph/CalcServer/Sort.hcl"
#include "aquagpusph/CalcServer/Sort.cl"
#endif
std::string SORT_INC = xxd2string(Sort_hcl_in, Sort_hcl_in_len);
std::string SORT_SRC = xxd2string(Sort_cl_in, Sort_cl_in_len);

Sort::Sort(const std::string tool_name,
           const std::string variable,
           const std::string permutations,
           const std::string inv_permutations,
           bool once)
  : Tool(tool_name, once)
  , _var_name(variable)
  , _perms_name(permutations)
  , _inv_perms_name(inv_permutations)
  , _var(NULL)
  , _perms(NULL)
  , _inv_perms(NULL)
  , _n(0)
  , _n_padded(0)
  , _init_kernel(NULL)
  , _start_kernel(NULL)
  , _local_kernel(NULL)
  , _global_kernel(NULL)
  , _inv_perms_kernel(NULL)
  , _vals(NULL)
  , _permut(NULL)
  , _local_work_size(0)
  , _global_work_size(0)
{
	Profiler::substages({ new EventProfile("init", this),
	                      new EventProfile("bitonic", this),
	                      new EventProfile("keys", this),
	                      new EventProfile("values", this),
	                      new EventProfile("inverse keys", this) });
}

Sort::~Sort()
{
	if (_init_kernel)
		clReleaseKernel(_init_kernel);
	if (_start_kernel)
		clReleaseKernel(_start_kernel);
	if (_local_kernel)
		clReleaseKernel(_local_kernel);
	if (_global_kernel)
		clReleaseKernel(_global_kernel);
	if (_inv_perms_kernel)
		clReleaseKernel(_inv_perms_kernel);
	if (_vals)
		clReleaseMemObject(_vals);
	if (_permut)
		clReleaseMemObject(_permut);
}

void
Sort::setup()
{
	std::ostringstream msg;
	msg << "Loading the tool \"" << name() << "\"..." << std::endl;
	LOG(L_INFO, msg.str());

	Tool::setup();

	// Get the variables
	variables();

	// Setup the working tools
	setupOpenCL();
}

cl_event
Sort::_execute(const std::vector<cl_event> events)
{
	cl_int err_code;
	cl_event event, wait_event;

	auto C = CalcServer::singleton();
	auto vars = C->variables();

	// Fisrt we copy everything on our transactional memory objects, which are
	// conveniently padded
	// We also initialize the permutations as the particles id, i.e. each
	// particle is converted on itself.
	wait_event = _var->getWritingEvent();
	err_code = clEnqueueNDRangeKernel(C->command_queue(),
	                                  _init_kernel,
	                                  1,
	                                  NULL,
	                                  &_global_work_size,
	                                  NULL,
	                                  1,
	                                  &wait_event,
	                                  &event);
	CHECK_OCL_OR_THROW(err_code,
	                   std::string("Failure executing \"init\" in tool \"") +
	                       name() + "\".");
	auto init_profiler =
	    dynamic_cast<EventProfile*>(Profiler::substages().at(0));
	init_profiler->start(event);
	init_profiler->end(event);
	wait_event = event;

	// Now we carry ut the bitonic sort on our transactional memory objects
	auto bitonic_profiler =
	    dynamic_cast<EventProfile*>(Profiler::substages().at(1));
	size_t bitonic_gsize = _global_work_size / 2;
	err_code = clEnqueueNDRangeKernel(C->command_queue(),
	                                  _start_kernel,
	                                  1,
	                                  NULL,
	                                  &bitonic_gsize,
	                                  &_local_work_size,
	                                  1,
	                                  &wait_event,
	                                  &event);
	CHECK_OCL_OR_THROW(
	    err_code,
	    std::string("Failure executing \"bitonic_start\" in tool \"") + name() +
	        "\".");
	err_code = clReleaseEvent(wait_event);
	CHECK_OCL_OR_THROW(
	    err_code,
	    std::string("Failure releasing \"init\" event in tool \"") + name() +
	        "\".");
	bitonic_profiler->start(event);
	wait_event = event;

	const unsigned int limit = 2u * _local_work_size;
	for (unsigned int bsize = limit; bsize <= _n_padded; bsize <<= 1) {
		for (unsigned int stride = bsize / 2; stride > 0; stride >>= 1) {
			// Execute either local or global, depending on whether the stride
			// is too big for the local memory or not
			auto kernel = (stride >= limit) ? _global_kernel : _local_kernel;
			std::string kname =
			    (stride >= limit) ? "bitonic_global" : "bitonic_local";
			err_code =
			    clSetKernelArg(kernel, 3, sizeof(cl_uint), (void*)&bsize);
			CHECK_OCL_OR_THROW(
			    err_code,
			    std::string("Failure sending argument 3 to kernel \"") + kname +
			        "\" in tool \"" + name() + "\".");
			err_code =
			    clSetKernelArg(kernel, 4, sizeof(cl_uint), (void*)&stride);
			CHECK_OCL_OR_THROW(
			    err_code,
			    std::string("Failure sending argument 4 to kernel \"") + kname +
			        "\" in tool \"" + name() + "\".");
			err_code = clEnqueueNDRangeKernel(C->command_queue(),
			                                  kernel,
			                                  1,
			                                  NULL,
			                                  &bitonic_gsize,
			                                  &_local_work_size,
			                                  1,
			                                  &wait_event,
			                                  &event);
			CHECK_OCL_OR_THROW(err_code,
			                   std::string("Failure executing kernel \"") +
			                       kname + "\" in tool \"" + name() + "\".");
			err_code = clReleaseEvent(wait_event);
			CHECK_OCL_OR_THROW(err_code,
			                   std::string("Failure releasing \"") + kname +
			                       "\" event in tool \"" + name() + "\".");
			wait_event = event;
		}
	}
	bitonic_profiler->end(wait_event);

	// Copy back the results
	auto inv_perms_events = _inv_perms->getReadingEvents();
	inv_perms_events.push_back(_inv_perms->getWritingEvent());

	auto perms_events = _perms->getReadingEvents();
	perms_events.push_back(_perms->getWritingEvent());
	perms_events.push_back(wait_event);
	err_code = clEnqueueCopyBuffer(C->command_queue(),
	                               _permut,
	                               *(cl_mem*)_perms->get(),
	                               0,
	                               0,
	                               _n * sizeof(cl_uint),
	                               perms_events.size(),
	                               perms_events.data(),
	                               &event);
	CHECK_OCL_OR_THROW(
	    err_code,
	    std::string("Failure copying the permutations in tool \"") + name() +
	        "\".");
	auto keys_profiler =
	    dynamic_cast<EventProfile*>(Profiler::substages().at(2));
	keys_profiler->start(event);
	keys_profiler->end(event);
	inv_perms_events.push_back(event);

	auto vals_events = _var->getReadingEvents();
	vals_events.push_back(_var->getWritingEvent());
	vals_events.push_back(wait_event);
	err_code = clEnqueueCopyBuffer(C->command_queue(),
	                               _vals,
	                               *(cl_mem*)_var->get(),
	                               0,
	                               0,
	                               _n * sizeof(cl_uint),
	                               vals_events.size(),
	                               vals_events.data(),
	                               &event);
	CHECK_OCL_OR_THROW(err_code,
	                   std::string("Failure copying the keys in tool \"") +
	                       name() + "\".");
	auto vals_profiler =
	    dynamic_cast<EventProfile*>(Profiler::substages().at(3));
	vals_profiler->start(event);
	vals_profiler->end(event);
	inv_perms_events.push_back(event);

	err_code = clReleaseEvent(wait_event);
	CHECK_OCL_OR_THROW(
	    err_code,
	    std::string("Failure releasing bitonic event in tool \"") + name() +
	        "\".");

	// And finally we can compute the inverse permutations
	err_code = clEnqueueNDRangeKernel(C->command_queue(),
	                                  _inv_perms_kernel,
	                                  1,
	                                  NULL,
	                                  &_global_work_size,
	                                  NULL,
	                                  inv_perms_events.size(),
	                                  inv_perms_events.data(),
	                                  &event);
	CHECK_OCL_OR_THROW(
	    err_code,
	    std::string("Failure executing \"inverse_keys\" in tool \"") + name() +
	        "\".");
	for (unsigned int i = 0; i < 2; i++) {
		err_code = clReleaseEvent(inv_perms_events.back());
		CHECK_OCL_OR_THROW(
		    err_code,
		    std::string("Failure releasing mem copy event in tool \"") +
		        name() + "\".");
		inv_perms_events.pop_back();
	}
	auto inv_profiler =
	    dynamic_cast<EventProfile*>(Profiler::substages().at(4));
	inv_profiler->start(event);
	inv_profiler->end(event);

	return event;
}

void
Sort::variables()
{
	size_t n;
	auto C = CalcServer::singleton();
	auto vars = C->variables();

	// Check and get the variables
	if (!vars->get(_var_name)) {
		std::ostringstream msg;
		msg << "Tool \"" << name()
		    << "\" is asking for the undeclared variable \"" << _var_name
		    << "\"." << std::endl;
		LOG(L_ERROR, msg.str());
		throw std::runtime_error("Invalid variable");
	}
	if (!vars->get(_var_name)->isArray()) {
		std::ostringstream msg;
		msg << "Tool \"" << name() << "\" is asking for the variable \""
		    << _var_name << "\", which is a scalar." << std::endl;
		LOG(L_ERROR, msg.str());
		throw std::runtime_error("Invalid variable");
	}
	_var_type = vars->get(_var_name)->type();
	if (_var_type.find("vec") != std::string::npos) {
		std::ostringstream msg;
		msg << "Tool \"" << name() << "\" cannot process variable \""
		    << _var_name << "\"." << std::endl;
		LOG(L_ERROR, msg.str());
		msg.str("");
		msg << "\t\"" << _var_type << "\" type is not supported." << std::endl;
		LOG(L_DEBUG, msg.str());
		throw std::runtime_error("Invalid variable type");
	}
	_var = (InputOutput::ArrayVariable*)vars->get(_var_name);
	if (_var_type.find("float") != std::string::npos)
		_var_max = "FLT_MAX";
	else if (_var_type.find("unsigned int") != std::string::npos)
		_var_max = "UINT_MAX";
	else if (_var_type.find("int") != std::string::npos)
		_var_max = "INT_MAX";

	if (!vars->get(_perms_name)) {
		std::ostringstream msg;
		msg << "Tool \"" << name()
		    << "\" is asking for the undeclared permutations variable \""
		    << _perms_name << "\"." << std::endl;
		LOG(L_ERROR, msg.str());
		throw std::runtime_error("Invalid variable");
	}
	if (vars->get(_perms_name)->type().compare("unsigned int*")) {
		std::ostringstream msg;
		msg << "Tool \"" << name()
		    << "\" cannot process permutations variable \"" << _perms_name
		    << "\"." << std::endl;
		LOG(L_ERROR, msg.str());
		msg.str("");
		msg << "\t\"unsigned int*\" type was expected, but \""
		    << vars->get(_perms_name)->type() << "\" has been received."
		    << std::endl;
		LOG(L_DEBUG, msg.str());
		throw std::runtime_error("Invalid variable type");
	}
	_perms = (InputOutput::ArrayVariable*)vars->get(_perms_name);

	if (!vars->get(_inv_perms_name)) {
		std::ostringstream msg;
		msg << "Tool \"" << name()
		    << "\" is asking for the undeclared inverse permutations variable "
		       "\""
		    << _inv_perms_name << "\"." << std::endl;
		LOG(L_ERROR, msg.str());
		throw std::runtime_error("Invalid variable");
	}
	if (vars->get(_inv_perms_name)->type().compare("unsigned int*")) {
		std::ostringstream msg;
		msg << "Tool \"" << name()
		    << "\" cannot process inverse permutations variable \""
		    << _inv_perms_name << "\"." << std::endl;
		LOG(L_ERROR, msg.str());
		msg.str("");
		msg << "\t\"unsigned int*\" type was expected, but \""
		    << vars->get(_inv_perms_name)->type() << "\" has been received."
		    << std::endl;
		LOG(L_DEBUG, msg.str());
		throw std::runtime_error("Invalid variable type");
	}
	_inv_perms = (InputOutput::ArrayVariable*)vars->get(_inv_perms_name);

	// Check the lengths
	_n = _var->size() / vars->typeToBytes(_var->type());
	_n_padded = nextPowerOf2(_n);
	n = _perms->size() / vars->typeToBytes(_perms->type());
	if (n != _n) {
		std::ostringstream msg;
		msg << "Lengths mismatch in tool \"" << name() << "\"." << std::endl;
		LOG(L_ERROR, msg.str());
		msg.str("");
		msg << "\tVariable \"" << _var->name() << "\" has length, n=" << _n
		    << std::endl;
		LOG(L_DEBUG, msg.str());
		msg.str("");
		msg << "\tVariable \"" << _perms->name() << "\" has length, n=" << n
		    << std::endl;
		LOG(L_DEBUG, msg.str());
		throw std::runtime_error("Invalid variable length");
	}
	n = _inv_perms->size() / vars->typeToBytes(_inv_perms->type());
	if (n != _n) {
		std::ostringstream msg;
		msg << "Lengths mismatch in tool \"" << name() << "\"." << std::endl;
		LOG(L_ERROR, msg.str());
		msg.str("");
		msg << "\tVariable \"" << _var->name() << "\" has length, n=" << _n
		    << std::endl;
		LOG(L_DEBUG, msg.str());
		msg.str("");
		msg << "\tVariable \"" << _inv_perms->name() << "\" has length, n=" << n
		    << std::endl;
		LOG(L_DEBUG, msg.str());
		throw std::runtime_error("Invalid variable length");
	}

	setOutputDependencies({ _var, _perms, _inv_perms });
}

void
Sort::setupOpenCL()
{
	cl_int err_code;
	auto C = CalcServer::singleton();
	auto vars = C->variables();

	// Get a first estimation of the local work group size from the available
	// memory
	cl_ulong local_mem_size = 0;
	err_code = clGetDeviceInfo(C->device(),
	                           CL_DEVICE_LOCAL_MEM_SIZE,
	                           sizeof(cl_ulong),
	                           &local_mem_size,
	                           NULL);
	CHECK_OCL_OR_THROW(
	    err_code,
	    std::string("Failure querying CL_DEVICE_LOCAL_MEM_SIZE in tool \"") +
	        name() + "\".");

	cl_ulong thread_size = 2 * (vars->typeToBytes(_var->type()) +
	                            vars->typeToBytes(_perms->type()));
	_local_work_size = (std::min)(local_mem_size / thread_size, 256ul);
	if (!isPowerOf2(_local_work_size)) {
		_local_work_size = nextPowerOf2(_local_work_size) >> 1;
	}

	auto kernels = compileOpenCL();
	bool valid_kernels = true;
	for (auto kernel : kernels) {
		size_t local_work_size;
		err_code = clGetKernelWorkGroupInfo(kernel,
		                                    C->device(),
		                                    CL_KERNEL_WORK_GROUP_SIZE,
		                                    sizeof(size_t),
		                                    &local_work_size,
		                                    NULL);
		CHECK_OCL_OR_THROW(
		    err_code,
		    std::string(
		        "Failure querying CL_KERNEL_WORK_GROUP_SIZE in tool \"") +
		        name() + "\".");
		if (local_work_size < _local_work_size) {
			valid_kernels = false;
			_local_work_size = local_work_size;
		}
	}
	if (!valid_kernels)
		kernels = compileOpenCL();
	_global_work_size = getGlobalWorkSize(_n_padded, _local_work_size);

	// Now
	_init_kernel = kernels.at(0);
	_start_kernel = kernels.at(1);
	_local_kernel = kernels.at(2);
	_global_kernel = kernels.at(3);
	_inv_perms_kernel = kernels.at(4);

	// Setup the memory objects
	setupMems();
	setupArgs();

	std::ostringstream msg;
	msg << "\tn: " << _n << std::endl;
	LOG0(L_DEBUG, msg.str());
	msg.str("");
	msg << "\tn (padded): " << _n_padded << std::endl;
	LOG0(L_DEBUG, msg.str());
	msg.str("");
	msg << "\tlocal work size: " << _local_work_size << std::endl;
	LOG0(L_DEBUG, msg.str());
	msg.str("");
	msg << "\tglobal work size: " << _global_work_size << std::endl;
	LOG0(L_DEBUG, msg.str());
}

std::vector<cl_kernel>
Sort::compileOpenCL()
{
	if (_init_kernel)
		clReleaseKernel(_init_kernel);
	if (_start_kernel)
		clReleaseKernel(_start_kernel);
	if (_local_kernel)
		clReleaseKernel(_local_kernel);
	if (_global_kernel)
		clReleaseKernel(_global_kernel);
	if (_inv_perms_kernel)
		clReleaseKernel(_inv_perms_kernel);

	std::ostringstream source;
	source << SORT_INC << SORT_SRC;

	std::ostringstream flags;
	flags << "-DMAX_LOCAL_SIZE=" << _local_work_size;
	return compile(source.str(),
	               { "init",
	                 "bitonic_start",
	                 "bitonic_local",
	                 "bitonic_global",
	                 "inverse_keys" },
	               flags.str());
}

void
Sort::setupMems()
{
	cl_int err_code;
	auto C = CalcServer::singleton();
	auto vars = C->variables();

	if (_vals)
		clReleaseMemObject(_vals);
	_vals = NULL;
	if (_permut)
		clReleaseMemObject(_permut);
	_permut = NULL;

	_vals = clCreateBuffer(C->context(),
	                       CL_MEM_READ_WRITE,
	                       _n_padded * vars->typeToBytes(_var->type()),
	                       NULL,
	                       &err_code);
	CHECK_OCL_OR_THROW(
	    err_code,
	    std::string("Failure allocating") +
	        std::to_string(_n_padded * vars->typeToBytes(_var->type())) +
	        " bytes for tool \"" + name() + "\".");
	_permut = clCreateBuffer(C->context(),
	                         CL_MEM_READ_WRITE,
	                         _n_padded * vars->typeToBytes(_perms->type()),
	                         NULL,
	                         &err_code);
	CHECK_OCL_OR_THROW(
	    err_code,
	    std::string("Failure allocating") +
	        std::to_string(_n_padded * vars->typeToBytes(_var->type())) +
	        " bytes for tool \"" + name() + "\".");

	allocatedMemory(_n_padded * vars->typeToBytes(_var->type()) +
	                _n_padded * vars->typeToBytes(_perms->type()));
}

void
Sort::setupArgs()
{
	cl_int err_code;
	CalcServer* C = CalcServer::singleton();

	err_code =
	    clSetKernelArg(_init_kernel, 0, sizeof(cl_mem), _var->get_async());
	CHECK_OCL_OR_THROW(
	    err_code,
	    std::string("Failure sending argument 0 to \"init\" in tool \"") +
	        name() + "\".");
	err_code = clSetKernelArg(_init_kernel, 1, sizeof(cl_mem), (void*)&_vals);
	CHECK_OCL_OR_THROW(
	    err_code,
	    std::string("Failure sending argument 1 to \"init\" in tool \"") +
	        name() + "\".");
	err_code = clSetKernelArg(_init_kernel, 2, sizeof(cl_mem), (void*)&_permut);
	CHECK_OCL_OR_THROW(
	    err_code,
	    std::string("Failure sending argument 2 to \"init\" in tool \"") +
	        name() + "\".");
	err_code = clSetKernelArg(_init_kernel, 3, sizeof(cl_uint), (void*)&_n);
	CHECK_OCL_OR_THROW(
	    err_code,
	    std::string("Failure sending argument 3 to \"init\" in tool \"") +
	        name() + "\".");
	err_code =
	    clSetKernelArg(_init_kernel, 4, sizeof(cl_uint), (void*)&_n_padded);
	CHECK_OCL_OR_THROW(
	    err_code,
	    std::string("Failure sending argument 4 to \"init\" in tool \"") +
	        name() + "\".");

	err_code = clSetKernelArg(_start_kernel, 0, sizeof(cl_mem), (void*)&_vals);
	CHECK_OCL_OR_THROW(
	    err_code,
	    std::string(
	        "Failure sending argument 0 to \"bitonic_start\" in tool \"") +
	        name() + "\".");
	err_code =
	    clSetKernelArg(_start_kernel, 1, sizeof(cl_mem), (void*)&_permut);
	CHECK_OCL_OR_THROW(
	    err_code,
	    std::string(
	        "Failure sending argument 1 to \"bitonic_start\" in tool \"") +
	        name() + "\".");

	err_code = clSetKernelArg(_local_kernel, 0, sizeof(cl_mem), (void*)&_vals);
	CHECK_OCL_OR_THROW(
	    err_code,
	    std::string(
	        "Failure sending argument 0 to \"bitonic_local\" in tool \"") +
	        name() + "\".");
	err_code =
	    clSetKernelArg(_local_kernel, 1, sizeof(cl_mem), (void*)&_permut);
	CHECK_OCL_OR_THROW(
	    err_code,
	    std::string(
	        "Failure sending argument 1 to \"bitonic_local\" in tool \"") +
	        name() + "\".");
	err_code =
	    clSetKernelArg(_local_kernel, 2, sizeof(cl_uint), (void*)&_n_padded);
	CHECK_OCL_OR_THROW(
	    err_code,
	    std::string(
	        "Failure sending argument 2 to \"bitonic_local\" in tool \"") +
	        name() + "\".");

	err_code = clSetKernelArg(_global_kernel, 0, sizeof(cl_mem), (void*)&_vals);
	CHECK_OCL_OR_THROW(
	    err_code,
	    std::string(
	        "Failure sending argument 0 to \"bitonic_global\" in tool \"") +
	        name() + "\".");
	err_code =
	    clSetKernelArg(_global_kernel, 1, sizeof(cl_mem), (void*)&_permut);
	CHECK_OCL_OR_THROW(
	    err_code,
	    std::string(
	        "Failure sending argument 1 to \"bitonic_global\" in tool \"") +
	        name() + "\".");
	err_code =
	    clSetKernelArg(_global_kernel, 2, sizeof(cl_uint), (void*)&_n_padded);
	CHECK_OCL_OR_THROW(
	    err_code,
	    std::string(
	        "Failure sending argument 2 to \"bitonic_global\" in tool \"") +
	        name() + "\".");

	err_code = clSetKernelArg(
	    _inv_perms_kernel, 0, sizeof(cl_mem), _perms->get_async());
	CHECK_OCL_OR_THROW(
	    err_code,
	    std::string(
	        "Failure sending argument 0 to \"inverse_keys\" in tool \"") +
	        name() + "\".");
	err_code = clSetKernelArg(
	    _inv_perms_kernel, 1, sizeof(cl_mem), _inv_perms->get_async());
	CHECK_OCL_OR_THROW(
	    err_code,
	    std::string(
	        "Failure sending argument 1 to \"inverse_keys\" in tool \"") +
	        name() + "\".");
	err_code =
	    clSetKernelArg(_inv_perms_kernel, 2, sizeof(cl_uint), (void*)&_n);
	CHECK_OCL_OR_THROW(
	    err_code,
	    std::string(
	        "Failure sending argument 2 to \"inverse_keys\" in tool \"") +
	        name() + "\".");
}

}
} // namespace
