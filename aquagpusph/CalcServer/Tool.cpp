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
 * @brief Tools virtual environment to allow the user to define/manipulate the
 * tools used to carry out the simulation.
 * (see Aqua::CalcServer::Tool)
 */

#include "Tool.hpp"
#include "CalcServer.hpp"
#include "aquagpusph/InputOutput/Logger.hpp"
#include <queue>
#include <chrono>

namespace Aqua {
namespace CalcServer {

void
Profile::step()
{
	auto C = CalcServer::singleton();
	_step = C->step();
}

void
Profile::sample(cl_ulong start, cl_ulong end)
{
	auto C = CalcServer::singleton();
	C->sample(_step, _tool, this, start, end);
}

Tool::Tool(const std::string tool_name, bool once)
  : Named(tool_name)
  , _once(once)
  , _next_tool(NULL)
  , _prev_tool(NULL)
  , _allocated_memory(0)
  , _n_iters(0)
  , _elapsed_time(0.f)
  , _average_elapsed_time(0.f)
  , _squared_elapsed_time(0.f)
  , _event(NULL)
  , _parent(NULL)
{
}

Tool::~Tool() {}

void
Tool::setup()
{
	std::vector<Tool*> tools = CalcServer::singleton()->tools();

	// Get the next tool in the pipeline
	int i = id_in_pipeline();
	// There are some tools which can be out of the pipeline (e.g. UnSort or
	// reports)
	if ((i >= 0) && (i + 1 < tools.size())) {
		next_tool(tools.at(i + 1));
	}
	return;
}

void CL_CALLBACK
exec_status_check(cl_event event, cl_int event_command_status, void* user_data)
{
	auto tool = (Tool*)user_data;

#ifdef HAVE_GPUPROFILE
	cl_command_type event_type;
	clGetEventInfo(event,
	               CL_EVENT_COMMAND_TYPE,
	               sizeof(cl_command_type),
	               &event_type,
	               NULL);
	if (event_type != CL_COMMAND_USER) {
		cl_ulong start = 0, end = 0;
		clGetEventProfilingInfo(
		    event, CL_PROFILING_COMMAND_START, sizeof(cl_ulong), &start, NULL);
		clGetEventProfilingInfo(
		    event, CL_PROFILING_COMMAND_END, sizeof(cl_ulong), &end, NULL);
		tool->addElapsedTime((cl_double)(end - start) * (cl_double)(1e-06));
	}
#endif

	clReleaseEvent(event);
	if (event_command_status == CL_COMPLETE)
		return;
	std::stringstream msg;
	msg << "Tool \"" << tool->name() << "\" exited with the error code "
	    << event_command_status << "." << std::endl;
	LOG(L_ERROR, msg.str());
	CalcServer::singleton()->raiseSIGINT();
}

bool
need_new_cmd(const Tool* prev_tool, std::vector<cl_event> events)
{
	if (!prev_tool) {
		// First tool, no new cmd required
		return false;
	} else if (!prev_tool->getEvent()) {
		// The tool is not providing an event to check. Better to use the same
		// queue, just in case
		return false;
	}

	// OK, let's check if the event provided by the previous tool is blocking
	// our execution, which would mean that we can get enqueued on the same
	// queue.
	for (auto event : events)
		if (event == prev_tool->getEvent())
			return false;
	// It seems we can get executed before the last tool
	return true;
}

void
Tool::execute()
{
	if (_once && (_n_iters > 0))
		return;

	cl_int err_code;

	// Collect the dependencies
	auto C = CalcServer::singleton();
	if (C->debug_mode()) {
		std::stringstream msg;
		msg << "Executing \"" << name() << "\"..." << std::endl;
		LOG(L_DEBUG, msg.str());
	}
	std::vector<cl_event> events = getEvents();
	if (need_new_cmd(_prev_tool, events))
		C->command_queue(CalcServer::cmd_queue::cmd_queue_new);

	// To avoid a infinite queues we block the execution until the previous
	// execution event is finished.
	if (_event) {
		err_code = clWaitForEvents(1, &_event);
		CHECK_OCL_OR_THROW(err_code,
		                   std::string("Failure waiting for \"") + name() +
		                       "\" tool previous event.");
		err_code = clReleaseEvent(_event);
		CHECK_OCL_OR_THROW(err_code,
		                   std::string("Failure releasing the \"") + name() +
		                       "\" tool previous event.");
	}

	for (auto substage : Profiler::substages())
		substage->step();
	cl_event event = _execute(events);
	_event = event;
	err_code = clFlush(C->command_queue());
	CHECK_OCL_OR_THROW(err_code,
		std::string("Failure flushing the command queue at \"") + name() +
		"\".");

	if (event != NULL) {
		// Add the event to the reading dependencies
		auto in_vars = getInputDependencies();
		for (auto var : in_vars) {
			var->addReadingEvent(event);
		}
		// Replace the writing dependencies event by the new one
		auto out_vars = getOutputDependencies();
		for (auto var : out_vars) {
			var->setEvent(event);
		}

		// Check for errors when the event is marked as completed
		// NOTE: clReleaseEvent() will be called on the callback, so no need to
		// do it here
		err_code = clRetainEvent(event);
		CHECK_OCL_OR_THROW(
		    err_code,
		    std::string("Failure retaining the new event in \"") + name() +
		        "\" tool.");
		err_code =
		    clSetEventCallback(event, CL_COMPLETE, exec_status_check, this);
		CHECK_OCL_OR_THROW(
		    err_code,
		    std::string("Failure registering the callback in tool \"") +
		        name() + "\".");
	}

	// Release the events in the wait list, which were retained by getEvents()
	for (auto e : events) {
		err_code = clReleaseEvent(e);
		CHECK_OCL_OR_THROW(
		    err_code,
		    std::string("Failure releasing a predecessor event in \"") +
		        name() + "\" tool.");
	}

	if (C->debug_mode()) {
		if (event != NULL) {
			err_code = clWaitForEvents(1, &event);
			CHECK_OCL_OR_THROW(err_code,
			                   std::string("Failure waiting for \"") + name() +
			                       "\" tool event.");
		}
		err_code = C->finish();
		CHECK_OCL_OR_THROW(
		    err_code,
		    std::string("Failure calling clFinish() after launching \"") +
		        name() + "\".");
		std::stringstream msg;
		msg << "\"" << name() << "\" finished" << std::endl;
		LOG(L_DEBUG, msg.str());
	}
}

int
Tool::id_in_pipeline()
{
	std::vector<Tool*> tools = CalcServer::singleton()->tools();
	auto it = std::find(tools.begin(), tools.end(), this);
	if (it == tools.end())
		return -1;
	return std::distance(tools.begin(), it);
}

void
Tool::addElapsedTime(float elapsed_time)
{
	_elapsed_time = elapsed_time;
	// Invert the average computation
	_average_elapsed_time *= _n_iters;
	_squared_elapsed_time *= _n_iters;
	// Add the new data
	_average_elapsed_time += elapsed_time;
	_squared_elapsed_time += elapsed_time * elapsed_time;
	// And average it again
	_n_iters++;
	_average_elapsed_time /= _n_iters;
	_squared_elapsed_time /= _n_iters;
}

std::vector<cl_kernel>
Tool::compile(const std::string source,
              const std::vector<std::string> names,
              const std::string additional_flags)
{
	cl_int err_code;
	cl_program program;
	std::vector<cl_kernel> kernels;
	Aqua::CalcServer::CalcServer* C = Aqua::CalcServer::CalcServer::singleton();

	if (!names.size()) {
		LOG(L_WARNING, "compile() called without kernel names\n");
		return kernels;
	}

	std::ostringstream flags;
#ifdef AQUA_DEBUG
	flags << " -DDEBUG ";
#else
	flags << " -DNDEBUG ";
#endif
	flags << " -cl-mad-enable -cl-fast-relaxed-math";
#ifdef HAVE_3D
	flags << " -DHAVE_3D";
#else
	flags << " -DHAVE_2D";
#endif
	flags << " " << additional_flags;

	size_t source_length = source.size();
	const char* source_cstr = source.c_str();
	program = clCreateProgramWithSource(
	    C->context(), 1, &source_cstr, &source_length, &err_code);
	CHECK_OCL_OR_THROW(err_code, "Failure creating the OpenCL program");
	err_code =
	    clBuildProgram(program, 0, NULL, flags.str().c_str(), NULL, NULL);
	if (err_code != CL_SUCCESS) {
		LOG(L_ERROR, "Error compiling the OpenCL script\n");
		InputOutput::Logger::singleton()->printOpenCLError(err_code);
		LOG0(L_ERROR, "--- Build log ---------------------------------\n");
		size_t log_size = 0;
		clGetProgramBuildInfo(
		    program, C->device(), CL_PROGRAM_BUILD_LOG, 0, NULL, &log_size);
		char* log = (char*)malloc(log_size + sizeof(char));
		if (!log) {
			std::stringstream msg;
			msg << "Failure allocating " << log_size
			    << " bytes for the building log" << std::endl;
			LOG0(L_ERROR, msg.str());
			LOG0(L_ERROR, "--------------------------------- Build log ---\n");
			throw std::bad_alloc();
		}
		strcpy(log, "");
		clGetProgramBuildInfo(
		    program, C->device(), CL_PROGRAM_BUILD_LOG, log_size, log, NULL);
		strcat(log, "\n");
		LOG0(L_DEBUG, log);
		LOG0(L_ERROR, "--------------------------------- Build log ---\n");
		free(log);
		log = NULL;
		clReleaseProgram(program);
		throw std::runtime_error("OpenCL compilation error");
	}
	for (auto name : names) {
		kernels.push_back(clCreateKernel(program, name.c_str(), &err_code));
		CHECK_OCL_OR_THROW(err_code,
		                   std::string("Failure creating the \"") + name +
		                       "\" OpenCL kernel.");
	}

	return kernels;
}

cl_kernel
Tool::compile_kernel(const std::string source,
                     const std::string kernel_name,
                     const std::string additional_flags)
{
	return compile(source, { kernel_name }, additional_flags).at(0);
}

const std::vector<cl_event>
Tool::getEvents() const
{
	cl_int err_code;
	std::vector<cl_event> res;
	// First traverse the input variables. We must wait for their writing
	// events
	for (auto var : _in_vars) {
		cl_event event = var->getWritingEvent();
		if (std::find(res.begin(), res.end(), event) != res.end())
			continue;
		res.push_back(event);
	}
	// Now we must parse the output variables. On output variables we must wait
	// for both their writing and their reading events
	for (auto var : _out_vars) {
		auto events = var->getReadingEvents();
		if (var->getWritingEvent())
			events.push_back(var->getWritingEvent());
		for (auto e : events) {
			if (std::find(res.begin(), res.end(), e) != res.end())
				continue;
			res.push_back(e);
		}
	}

	for (auto event : res) {
		// Retain the event until we work with it
		err_code = clRetainEvent(event);
		CHECK_OCL_OR_THROW(err_code,
		                   std::string("Failure retaining an event in \"") +
		                       name() + "\" tool.");
	}

	return res;
}

std::vector<InputOutput::Variable*>
Tool::namesToVars(const std::vector<std::string>& names) const
{
	InputOutput::Variables* vars = CalcServer::singleton()->variables();
	std::vector<InputOutput::Variable*> res;
	for (auto var_name : names) {
		InputOutput::Variable* var = vars->get(var_name);
		if (!var) {
			std::stringstream msg;
			msg << "The tool \"" << name()
			    << "\" is asking the undeclared variable \"" << var_name
			    << "\"." << std::endl;
			LOG(L_ERROR, msg.str());
			throw std::runtime_error("Invalid variable");
		}
		res.push_back(var);
	}
	return res;
}

}
} // namespace
