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

#include <AuxiliarMethods.h>
#include <InputOutput/Logger.h>
#include <CalcServer.h>
#include <CalcServer/Kernel.h>

namespace Aqua {
namespace CalcServer {

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
{
}

Kernel::~Kernel()
{
	if (_kernel)
		clReleaseKernel(_kernel);
	_kernel = NULL;
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
	setVariables(true);
	computeGlobalWorkSize();
}

/** @brief Callback called when all the dependencies of the
 * Aqua::CalcServer::Kernel tool are fulfilled.
 *
 * This function is just redirecting the work to
 * Aqua::CalcServer::Kernel::setVariables()
 * @param event The triggering event
 * @param event_command_status CL_COMPLETE upon all dependencies successfully
 * fulfilled. A negative integer if one or mor dependencies failed.
 * @param user_data A casted pointer to the Aqua::CalcServer::Kernel
 * tool (or the inherited one)
 */
void CL_CALLBACK
args_setter(cl_event event, cl_int event_command_status, void* user_data)
{
	clReleaseEvent(event);
	auto tool = (Kernel*)user_data;
	if (event_command_status != CL_COMPLETE) {
		std::stringstream msg;
		msg << "Skipping \"" << tool->name() << "\" due to dependency errors."
		    << std::endl;
		LOG(L_WARNING, msg.str());
		clSetUserEventStatus(tool->getUserEvent(), event_command_status);
		clReleaseEvent(tool->getUserEvent());
		return;
	}

	try {
		tool->setVariables();
	} catch (...) {
		std::stringstream msg;
		msg << "Failure setting the \"" << tool->name() << "\" args"
		    << std::endl;
		LOG(L_WARNING, msg.str());
		clSetUserEventStatus(tool->getUserEvent(), -1);
		clReleaseEvent(tool->getUserEvent());
		return;
	}
	clSetUserEventStatus(tool->getUserEvent(), CL_COMPLETE);
	clReleaseEvent(tool->getUserEvent());
}

cl_event
Kernel::_execute(const std::vector<cl_event> events)
{
	cl_int err_code;
	CalcServer* C = CalcServer::singleton();

	// Ask to set the kernel arguments as fast as possible
	cl_event trigger;
	std::vector<cl_event> scalar_deps;
	for (auto var : _scalar_vars) {
		if (var->getWritingEvent())
			scalar_deps.push_back(var->getWritingEvent());
	}
	cl_uint num_scalar_events = scalar_deps.size();
	const cl_event* scalar_events = scalar_deps.size() ? scalar_deps.data() : NULL;
	err_code = clEnqueueMarkerWithWaitList(
	    C->command_queue(), num_scalar_events, scalar_events, &trigger);
	if (err_code != CL_SUCCESS) {
		std::stringstream msg;
		msg << "Failure setting the trigger for tool \"" << name() << "\"."
		    << std::endl;
		LOG(L_ERROR, msg.str());
		InputOutput::Logger::singleton()->printOpenCLError(err_code);
		throw std::runtime_error("OpenCL execution error");
	}

	// Now we create a user event that we will set as completed when we already
	// setted the kernel args
	auto args_event = clCreateUserEvent(C->context(), &err_code);
	if (err_code != CL_SUCCESS) {
		std::stringstream msg;
		msg << "Failure creating the args setter event for tool \"" << name()
		    << "\"." << std::endl;
		LOG(L_ERROR, msg.str());
		InputOutput::Logger::singleton()->printOpenCLError(err_code);
		throw std::runtime_error("OpenCL execution error");
	}
	err_code = clRetainEvent(args_event);
	if (err_code != CL_SUCCESS) {
		std::stringstream msg;
		msg << "Failure retaining the args setter event for tool \"" << name()
		    << "\"." << std::endl;
		LOG(L_ERROR, msg.str());
		InputOutput::Logger::singleton()->printOpenCLError(err_code);
		throw std::runtime_error("OpenCL execution error");
	}
	_user_event = args_event;

	// So it is time to register our callback on our trigger
	err_code = clSetEventCallback(trigger, CL_COMPLETE, args_setter, this);
	if (err_code != CL_SUCCESS) {
		std::stringstream msg;
		msg << "Failure registering the solver callback in tool \"" << name()
		    << "\"." << std::endl;
		LOG(L_ERROR, msg.str());
		InputOutput::Logger::singleton()->printOpenCLError(err_code);
		throw std::runtime_error("OpenCL execution error");
	}

	// And now we can enqueue the kernel execution
	cl_event event;
	std::vector<cl_event> kernel_events = events;
	kernel_events.push_back(args_event);
	err_code = clEnqueueNDRangeKernel(C->command_queue(),
	                                  _kernel,
	                                  1,
	                                  NULL,
	                                  &_global_work_size,
	                                  &_work_group_size,
	                                  kernel_events.size(),
	                                  kernel_events.data(),
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
		msg << "Failure releasing the args setter event for tool \"" << name()
		    << "\"." << std::endl;
		LOG(L_ERROR, msg.str());
		InputOutput::Logger::singleton()->printOpenCLError(err_code);
		throw std::runtime_error("OpenCL execution error");
	}

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
	if (err_code != CL_SUCCESS) {
		LOG(L_ERROR, "Failure querying the work group size.\n");
		InputOutput::Logger::singleton()->printOpenCLError(err_code);
		clReleaseKernel(kernel);
		throw std::runtime_error("OpenCL error");
	}
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
	if (err_code != CL_SUCCESS) {
		LOG(L_WARNING, "Failure releasing non-local memory based kernel.\n");
		InputOutput::Logger::singleton()->printOpenCLError(err_code);
	}
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
	if (err_code != CL_SUCCESS) {
		std::ostringstream msg;
		msg << "Failure asking for CL_KERNEL_ARG_ADDRESS_QUALIFIER on arg "
		    << arg_index << std::endl;
		LOG(L_WARNING, msg.str());
		return false;
	}
	if (address != CL_KERNEL_ARG_ADDRESS_GLOBAL)
		return true;
	cl_kernel_arg_type_qualifier type;
	err_code = clGetKernelArgInfo(kernel,
	                              arg_index,
	                              CL_KERNEL_ARG_TYPE_QUALIFIER,
	                              sizeof(cl_kernel_arg_type_qualifier),
	                              &type,
	                              NULL);
	if (err_code != CL_SUCCESS) {
		std::ostringstream msg;
		msg << "Failure asking for CL_KERNEL_ARG_TYPE_QUALIFIER on arg "
		    << arg_index << std::endl;
		LOG(L_WARNING, msg.str());
		return false;
	}
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
		if (!var->isArray())
			_scalar_vars.push_back(var);
	}

	setDependencies(in_vars, out_vars);
}

void
Kernel::setVariables(bool with_arrays)
{
	unsigned int i;
	cl_int err_code;

	for (unsigned int i = 0; i < _vars.size(); i++) {
		auto var = _vars[i];
		if (var->isArray() && !with_arrays)
			continue;
		// Update the variable
		err_code = clSetKernelArg(_kernel,
		                          i,
		                          var->typesize(),
		                          var->get_async());
		if (err_code != CL_SUCCESS) {
			std::stringstream msg;
			msg << "Failure setting the variable \"" << var->name()
			    << "\" (id=" << i << ") to the tool \"" << name() << "\"."
			    << std::endl;
			LOG(L_ERROR, msg.str());
			InputOutput::Logger::singleton()->printOpenCLError(err_code);
			throw std::runtime_error("OpenCL error");
		}
	}
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
