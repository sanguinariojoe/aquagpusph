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
 * @brief OpenCL kernel kernel based tool.
 * (see Aqua::CalcServer::Kernel for details)
 */

#include "aquagpusph/AuxiliarMethods.hpp"
#include "aquagpusph/InputOutput/Logger.hpp"
#include "CalcServer.hpp"
#include "Kernel.hpp"

namespace Aqua {
namespace CalcServer {

/** @brief Callback called when the profiler ending event has finished
 *
 * @param event The Ending event
 * @param cmd_exec_status Event status
 * @param user_data Recasted EventProfile pointer
 */
void CL_CALLBACK
event_profile_cb(cl_event n_event, cl_int cmd_exec_status, void* user_data)
{
	cl_int err_code;
	auto profiler = (EventProfile*)user_data;
	profiler->sample();
}

void
EventProfile::start(cl_event event)
{
#ifdef HAVE_GPUPROFILE
	cl_int err_code;
	_start = event;
	err_code = clRetainEvent(event);
	CHECK_OCL_OR_THROW(
	    err_code,
	    std::string("Failure retaining the event on profiler \"") + name() +
	        "\".");
#endif
}

void
EventProfile::end(cl_event event)
{
#ifdef HAVE_GPUPROFILE
	cl_int err_code;
	_end = event;
	err_code = clRetainEvent(event);
	CHECK_OCL_OR_THROW(
	    err_code,
	    std::string("Failure retaining the event on profiler \"") + name() +
	        "\".");
	err_code = clSetEventCallback(event, CL_COMPLETE, &event_profile_cb, this);
	CHECK_OCL_OR_THROW(
	    err_code,
	    std::string("Failure setting the profiling callback on \"") + name() +
	        "\".");
#endif
}

void
EventProfile::sample()
{
	cl_int err_code;
	cl_ulong start = 0, end = 0;
	err_code = clGetEventProfilingInfo(
	    _start, CL_PROFILING_COMMAND_START, sizeof(cl_ulong), &start, NULL);
	CHECK_OCL_OR_THROW(
	    err_code,
	    std::string("Failure retrieving the profiling start on \"") + name() +
	        "\".");
	err_code = clReleaseEvent(_start);
	CHECK_OCL_OR_THROW(
	    err_code,
	    std::string("Failure releasing the start event on profiler \"") +
	        name() + "\".");
	err_code = clGetEventProfilingInfo(
	    _end, CL_PROFILING_COMMAND_END, sizeof(cl_ulong), &end, NULL);
	CHECK_OCL_OR_THROW(
	    err_code,
	    std::string("Failure retrieving the profiling end on \"") + name() +
	        "\".");
	err_code = clReleaseEvent(_start);
	CHECK_OCL_OR_THROW(
	    err_code,
	    std::string("Failure releasing the end event on profiler \"") +
	        name() + "\".");
	auto C = CalcServer::singleton();
	start += C->device_timer_offset();
	end += C->device_timer_offset();
	Profile::sample(start, end);
}

ArgSetter::ArgSetter(const std::string name,
                     cl_kernel kernel,
                     std::vector<InputOutput::Variable*> vars)
  : Named(name)
  , _kernel(kernel)
  , _vars(vars)
{
	for (unsigned int i = 0; i < vars.size(); i++) {
		Arg* arg = new Arg();
		_args.push_back(Arg());
	}
}

void
ArgSetter::execute()
{
	cl_int err_code;
	for (unsigned int i = 0; i < _vars.size(); i++) {
		auto var = _vars[i];
		if (!var)
			continue;
		auto arg = _args[i];
		if (arg == var)
			continue;
		arg = var;
		// Update the variable
		err_code = clSetKernelArg(_kernel, i, arg.size(), arg.value());
		CHECK_OCL_OR_THROW(err_code,
		                   std::string("Failure setting the variable \"") +
		                       var->name() + "\" (id=" + std::to_string(i) +
		                       ") in tool \"" + name() + "\".");
		_args[i] = arg;
	}
}

bool
ArgSetter::Arg::operator==(InputOutput::Variable* var)
{
	if (!var)
		return false;
	if (var->getWritingEvent() == _event)
		return true;
	if (var->typesize() != _size)
		return false;
	// We have 2 cases here:
	//  - Non-reallocatable array variables, that are never changing so we can
	//    get them asynchronously
	//  - Reallocatable array variables and scalar variables, that we must wait
	//    until they are wrote
	if (!var->isArray() ||
	    ((InputOutput::ArrayVariable*)var)->reallocatable()) {
		auto event = var->getWritingEvent();
		auto err_code = clWaitForEvents(1, &event);
		CHECK_OCL_OR_THROW(err_code,
		                   std::string("Failure waiting for the variable \"") +
		                       var->name() + "\" to be written");
	}
	return !memcmp(var->get_async(), _value, _size);
}

/** @brief Syncing utility to set an user event status to the same status of
 * another OpenCL event
 *
 * This utility shall be used by means of sync_user_event()
 *
 * @param n_event Associated event
 * @param cmd_exec_status Triggering status
 * @param user_data Recasted cl_event pointer
 * @note If cmd_exec_status=CL_COMPLETE clReleaseEvent will be called on top of
 * the user event. Thus, call clRetainEvent if you want to further keep the
 * event
 */
void CL_CALLBACK
cbUserEventSync(cl_event n_event, cl_int cmd_exec_status, void* user_data)
{
	cl_int err_code;
	cl_event event = *(cl_event*)user_data;

	err_code = clSetUserEventStatus(event, cmd_exec_status);
	CHECK_OCL_OR_THROW(err_code, "Failure setting user event status");
	err_code = clReleaseEvent(event);
	CHECK_OCL_OR_THROW(err_code, "Failure releasing user event");
	err_code = clReleaseEvent(n_event);
	CHECK_OCL_OR_THROW(err_code, "Failure releasing the triggering event");
	free(user_data);
}

void
sync_user_event(cl_event user_event, cl_event event)
{
	cl_int err_code;
	void* user_data = malloc(sizeof(cl_event));
	if (!user_data) {
		std::stringstream msg;
		msg << "Failure allocating " << sizeof(cl_event)
		    << " bytes for the user event storage" << std::endl;
		LOG(L_ERROR, msg.str());
		throw std::bad_alloc();
	}

	memcpy(user_data, &user_event, sizeof(cl_event));

	err_code =
	    clSetEventCallback(event, CL_COMPLETE, &cbUserEventSync, user_data);
	CHECK_OCL_OR_THROW(err_code, "Failure setting the events syncing callback");

	// Let's try to flush the command queue
	cl_command_queue queue;
	err_code = clGetEventInfo(event, CL_EVENT_COMMAND_QUEUE,
	                          sizeof(cl_command_queue), &queue, NULL);
	CHECK_OCL_OR_THROW(err_code, "Failure getting the event command queue");
	if (queue) {
		err_code = clFlush(queue);
		CHECK_OCL_OR_THROW(err_code, "Failure flushing the command queue");
	}
}

Kernel::Kernel(const std::string tool_name,
               const std::string kernel_path,
               const std::string entry_point,
               const std::string n,
               bool once)
  : Tool(tool_name, once)
  , _path(kernel_path)
  , _entry_point(entry_point)
  , _n(n)
  , _kernel(NULL)
  , _work_group_size(0)
  , _global_work_size(0)
  , _args_setter(NULL)
{
	Profiler::substages({ new EventProfile("Kernel", this) });
}

Kernel::~Kernel()
{
	if (_kernel)
		clReleaseKernel(_kernel);
	if (_args_setter)
		delete _args_setter;
}

void
Kernel::setup()
{
	std::ostringstream msg;
	msg << "Loading the tool \"" << name() << "\" from the file \"" << path()
	    << "\"..." << std::endl;
	LOG(L_INFO, msg.str());

	Tool::setup();
	make(_entry_point);
	variables(_entry_point);
	computeGlobalWorkSize();
	_args_setter = new ArgSetter(name(), _kernel, _vars);
	if (!_args_setter) {
		std::stringstream msg;
		msg << "Failure creating the arg setter for tool \"" << name()
		    << std::endl;
		LOG(L_ERROR, msg.str());
		throw std::bad_alloc();
	}
	_args_setter->execute();
}

cl_event
Kernel::_execute(const std::vector<cl_event> events)
{
	cl_int err_code;
	auto C = CalcServer::singleton();

	// Asynchronously set the kernel arguments
	_args_setter->execute();

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
Kernel::make(const std::string entry_point,
             const std::string add_flags,
             const std::string header)
{
	unsigned int i;
	cl_program program;
	cl_kernel kernel;
	std::ostringstream source;
	std::ostringstream flags;
	size_t source_length = 0;
	cl_int err_code = CL_SUCCESS;
	size_t work_group_size = 0;
	CalcServer* C = CalcServer::singleton();

	// Read the script file
	std::ifstream script(path());
	if (!script) {
		std::stringstream msg;
		msg << "Failure reading the file \"" << path() << "\"." << std::endl;
		LOG(L_ERROR, msg.str());
		throw std::ifstream::failure(msg.str());
	}
	source << header << script.rdbuf();

	// Setup the default flags
	flags << "-I" << getFolderFromFilePath(path()) << " ";
	if (C->base_path().compare("")) {
		flags << "-I" << C->base_path() << " ";
	}
	// Setup the user registered flags
	for (auto def : C->definitions()) {
		flags << def << " ";
	}
	// Add the additionally specified flags
	flags << add_flags;

	// Try to compile without using local memory
	LOG(L_INFO, "Compiling without local memory... ");
	kernel = compile_kernel(source.str(), entry_point, flags.str());

	// Get the work group size
	err_code = clGetKernelWorkGroupInfo(kernel,
	                                    C->device(),
	                                    CL_KERNEL_WORK_GROUP_SIZE,
	                                    sizeof(size_t),
	                                    &work_group_size,
	                                    NULL);
	CHECK_OCL_OR_THROW(
	    err_code,
	    std::string("Failure querying CL_KERNEL_WORK_GROUP_SIZE in tool \"") +
	        name() + "\".");
	LOG0(L_DEBUG, "OK\n");

	_kernel = kernel;
	_work_group_size = work_group_size;

	// Try to compile with local memory
	LOG(L_INFO, "Compiling with local memory... ");
	flags << " -DLOCAL_MEM_SIZE=" << work_group_size;
	try {
		kernel = compile_kernel(source.str(), entry_point, flags.str());
	} catch (std::runtime_error& e) {
		LOG(L_INFO, "Falling back to no local memory usage.\n");
		return;
	}
	cl_ulong used_local_mem;
	err_code = clGetKernelWorkGroupInfo(kernel,
	                                    C->device(),
	                                    CL_KERNEL_LOCAL_MEM_SIZE,
	                                    sizeof(cl_ulong),
	                                    &used_local_mem,
	                                    NULL);
	if (err_code != CL_SUCCESS) {
		LOG(L_ERROR, "Failure querying the used local memory.\n");
		InputOutput::Logger::singleton()->printOpenCLError(err_code);
		clReleaseKernel(kernel);
		LOG(L_INFO, "Falling back to no local memory usage.\n");
		return;
	}
	cl_ulong available_local_mem;
	err_code = clGetDeviceInfo(C->device(),
	                           CL_DEVICE_LOCAL_MEM_SIZE,
	                           sizeof(cl_ulong),
	                           &available_local_mem,
	                           NULL);
	if (err_code != CL_SUCCESS) {
		LOG(L_ERROR, "Failure querying the device available local memory.\n");
		InputOutput::Logger::singleton()->printOpenCLError(err_code);
		clReleaseKernel(kernel);
		LOG(L_INFO, "Falling back to no local memory usage.\n");
		return;
	}

	if (available_local_mem < used_local_mem) {
		LOG(L_ERROR, "Not enough available local memory.\n");
		InputOutput::Logger::singleton()->printOpenCLError(err_code);
		clReleaseKernel(kernel);
		LOG(L_INFO, "Falling back to no local memory usage.\n");
		return;
	}
	LOG0(L_DEBUG, "OK\n");

	// Swap kernels
	err_code = clReleaseKernel(_kernel);
	CHECK_OCL_OR_THROW(
	    err_code,
	    std::string("Failure releasing non-local memory based kernel ") +
	        "in tool \"" + name() + "\".");
	_kernel = kernel;
}

bool
isKernelArgReadOnly(cl_kernel kernel, cl_uint arg_index)
{
	cl_int err_code;
	cl_kernel_arg_address_qualifier address;
	err_code = clGetKernelArgInfo(kernel,
	                              arg_index,
	                              CL_KERNEL_ARG_ADDRESS_QUALIFIER,
	                              sizeof(cl_kernel_arg_address_qualifier),
	                              &address,
	                              NULL);
	CHECK_OCL_OR_THROW(
	    err_code,
	    std::string("Failure querying CL_KERNEL_ARG_ADDRESS_QUALIFIER ") +
	        " on arg " + std::to_string(arg_index) + ".");
	if (address != CL_KERNEL_ARG_ADDRESS_GLOBAL)
		return true;
	cl_kernel_arg_type_qualifier type;
	err_code = clGetKernelArgInfo(kernel,
	                              arg_index,
	                              CL_KERNEL_ARG_TYPE_QUALIFIER,
	                              sizeof(cl_kernel_arg_type_qualifier),
	                              &type,
	                              NULL);
	CHECK_OCL_OR_THROW(
	    err_code,
	    std::string("Failure querying CL_KERNEL_ARG_TYPE_QUALIFIER ") +
	        " on arg " + std::to_string(arg_index) + ".");
	if (type == CL_KERNEL_ARG_TYPE_CONST)
		return true;
	return false;
}

void
Kernel::variables(const std::string entry_point)
{
	cl_int err_code;

	cl_uint num_args;
	InputOutput::Variables* vars = CalcServer::singleton()->variables();
	err_code = clGetKernelInfo(
	    _kernel, CL_KERNEL_NUM_ARGS, sizeof(cl_uint), &num_args, NULL);
	if (err_code != CL_SUCCESS) {
		LOG(L_WARNING, "Failure asking for CL_KERNEL_NUM_ARGS.\n");
		InputOutput::Logger::singleton()->printOpenCLError(err_code);
	}

	std::vector<InputOutput::Variable*> in_vars, out_vars;
	for (cl_uint i = 0; i < num_args; i++) {
		char* arg_name;
		size_t arg_name_len;
		err_code = clGetKernelArgInfo(
		    _kernel, i, CL_KERNEL_ARG_NAME, 0, NULL, &arg_name_len);
		if (err_code != CL_SUCCESS) {
			LOG(L_WARNING, "Failure asking for CL_KERNEL_ARG_NAME len.\n");
			InputOutput::Logger::singleton()->printOpenCLError(err_code);
		}
		arg_name = new char[arg_name_len + 1];
		if (!arg_name) {
			std::ostringstream msg;
			msg << "Failure allocating " << arg_name_len + 1
			    << " for CL_KERNEL_ARG_NAME" << std::endl;
			LOG(L_ERROR, msg.str());
			throw std::bad_alloc();
		}
		arg_name[arg_name_len] = '\0';
		err_code = clGetKernelArgInfo(
		    _kernel, i, CL_KERNEL_ARG_NAME, arg_name_len, arg_name, NULL);
		if (err_code != CL_SUCCESS) {
			LOG(L_WARNING, "Failure asking for CL_KERNEL_ARG_NAME.\n");
			InputOutput::Logger::singleton()->printOpenCLError(err_code);
		}
		auto var = vars->get(arg_name);
		if (!var) {
			std::stringstream msg;
			msg << "The tool \"" << name()
			    << "\" is asking the undeclared variable \"" << arg_name
			    << "\"." << std::endl;
			LOG(L_ERROR, msg.str());
			throw std::runtime_error("Invalid variable");
		}
		if (isKernelArgReadOnly(_kernel, i))
			in_vars.push_back(var);
		else
			out_vars.push_back(var);
		_vars.push_back(var);
	}

	setDependencies(in_vars, out_vars);
}

void
Kernel::computeGlobalWorkSize()
{
	unsigned int N = 0;
	if (!_work_group_size) {
		LOG(L_ERROR, "Work group size must be greater than 0.\n");
		throw std::runtime_error("Null work group size");
	}
	InputOutput::Variables* vars = CalcServer::singleton()->variables();
	if (_n == "") {
		for (auto var : _vars) {
			if (!var->isArray())
				continue;
			unsigned int var_N = var->size() / vars->typeToBytes(var->type());
			if (var_N > N)
				N = var_N;
		}
	} else {
		try {
			vars->solve("unsigned int", _n, &N);
		} catch (...) {
			LOG(L_ERROR, "Failure evaluating the number of threads.\n");
			throw std::runtime_error("Invalid number of threads");
		}
	}

	_global_work_size = (size_t)roundUp(N, (unsigned int)_work_group_size);
	{
		std::stringstream msg;
		msg << _global_work_size << " threads ("
		    << _global_work_size / _work_group_size << " groups of "
		    << _work_group_size << " threads)" << std::endl;
		LOG(L_INFO, msg.str());
	}
}

}
} // namespace
