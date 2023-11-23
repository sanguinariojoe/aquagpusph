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

#ifndef KERNEL_H_INCLUDED
#define KERNEL_H_INCLUDED

#include "aquagpusph/sphPrerequisites.hpp"

#include <vector>
#if __APPLE__
#include <OpenCL/cl.h>
#else
#include <CL/cl.h>
#endif
#include "Tool.hpp"

namespace Aqua {
namespace CalcServer {

/** @brief Synchronize an user event with another OpenCL event
 *
 * clReleaseEvent() will be called on both events. Thus, call
 * clRetainEvent() if you want to further keep them
 * @param user_event User event to be synced
 * @param event Associated OpenCL event
 */
void
sync_user_event(cl_event user_event, cl_event event);

/** @class EventProfile kernel.h CalcServer/Kernel.h
 * @brief Profiler for tools based on OpenCL enqueued commands
 */
class EventProfile : public Aqua::CalcServer::Profile
{
  public:
	/** Constructor
	 * @param name The name
	 * @param tool The owning tool
	 */
	EventProfile(const std::string name, Tool *tool)
	  : Profile(name, tool)
	{
	}

	/** Destructor
	 */
	~EventProfile() {}

	/** @brief Starting event
	 * @note It might be the same event passed to ::end()
	 */
	void start(cl_event event);

	/** @brief Ending event
	 * @note It might be the same event passed to ::start()
	 */
	void end(cl_event event);

	/** @brief Compute the time stamps
	 */
	void sample();

  private:
	/// The starting event
	cl_event _start;

	/// The ending event
	cl_event _end;
};

/** @class ArgsSetter Kernel.h CalcServer/Kernel.h
 * @brief A helper class to asynchronously set the arguments for an OpenCL
 * kernel.
 *
 * This class will provide the caller with an event that shall be waited before
 * the kernel can be enqueued.
 */
class ArgSetter : public Aqua::CalcServer::Named
{
	class Arg;

  public:
	/** Constructor.
	 * @param name A name for logging purposes. Usually the same than the owner
	 * tool
	 * @param kernel The OpenCL kernel.
	 * @param vars The list of variables to be set, on the same order they
	 * appear on @p kernel. NULL can be passed for the arguments that are not
	 * handled by this class.
	 */
	ArgSetter(const std::string name,
	          cl_kernel kernel,
	          std::vector<InputOutput::Variable*> vars);

	/** Destructor
	 */
	~ArgSetter() {}

	/** @brief Get the OpenCL kernel
	 * @return The OpenCL kernel
	 */
	inline cl_kernel getKernel() const { return _kernel; }

	/** @brief Get the list of variables to be set
	 * @return The list of variables
	 */
	inline std::vector<InputOutput::Variable*> getVars() const { return _vars; }

	/** @brief Get the list of variables to be set
	 * @return The list of variables
	 */
	inline std::vector<Arg> getArgs() const { return _args; }

	/** @brief Set the kernel arguments
	 * @return The event that mark when the kernel arguments has been already
	 * set. Remember to call clReleaseEvent() on this.
	 * @warning The arguments are set right after all the variables writing
	 * events are finished. Thus, waiting just for this might not be enough,
	 * but the reading events of the kernel output variables shall be waited
	 * for also
	 */
	void execute();

  private:
	/// OpenCL kernel
	cl_kernel _kernel;

	/// List of variables in the same order they appear as kernel arguments
	std::vector<InputOutput::Variable*> _vars;

	/// List of already set kernel arguments
	std::vector<Arg> _args;
};

/** @class Arg Kernel.h CalcServer/Kernel.h
 * @brief A helper class to check whether a kernel argument has changed or not.
 *
 * This class is reading the variables asynchronously, so if synchronous
 * operations shall be carried out then InputOutput::Variable::sync() must
 * be called before
 */
class ArgSetter::Arg
{
  public:
	/// Constructor.
	Arg()
	  : _size(0)
	  , _value(NULL)
	  , _event(NULL)
	{
	}

	/** @brief Copy Constructor
	 * @param arg The argument to copy
	 */
	Arg(const Arg& arg)
	  : _size(0)
	  , _value(NULL)
	  , _event(NULL)
	{
		set(arg.size(), arg.value(), arg.event());
	}

	/// Destructor
	~Arg()
	{
		if (_value)
			free(_value);
	}

	/** @brief Get the size of the argument
	 * @return The size
	 */
	const size_t size() const { return _size; }

	/** @brief Get the value of the argument
	 * @return The value memory
	 */
	const void* value() const { return _value; }

	/** @brief Get the last recorded writing event
	 * @return The event
	 */
	const cl_event event() const { return _event; }

	/** @brief Check whether the content of a variable is the same than this
	 * argument
	 * @param var Variable to compare with
	 * @return true if the variable still matches the argument, false otherwise
	 */
	inline bool operator==(InputOutput::Variable* var)
	{
		if (!var)
			return false;
		if (var->getWritingEvent() == _event)
			return true;
		if (var->typesize() != _size)
			return false;
		return !memcmp(var->get_async(), _value, _size);
	}

	/** @brief Check whether the content of a variable is the same than this
	 * argument
	 * @param var Variable to compare with
	 * @return false if the variable still matches the argument, true otherwise
	 */
	inline bool operator!=(InputOutput::Variable* var)
	{
		return !(*this == var);
	}

	/** @brief Copy other argument
	 * @param arg Argument to copy
	 * @throw std::bad_alloc If the memory required cannot be allocated
	 */
	inline void operator=(const Arg& arg)
	{
		set(arg.size(), arg.value(), arg.event());
	}

	/** @brief Copy the variable content on the argument
	 * @param var Variable to assign
	 * @throw std::bad_alloc If the memory required cannot be allocated
	 */
	inline void operator=(InputOutput::Variable* var)
	{
		try {
			set(var->typesize(), var->get(), var->getWritingEvent());
		} catch (...) {
			std::stringstream msg;
			msg << "From variable " << var->name() << std::endl;
			LOG0(L_DEBUG, msg.str());
			throw;
		}
	}

  protected:
	/** @brief Set the argument
	 * @param s Size
	 * @param v Value
	 * @param e Event
	 * @throw std::bad_alloc If the memory required cannot be allocated
	 */
	inline void set(const size_t s, const void* v, const cl_event e)
	{
		if (_size != s) {
			if (_value)
				free(_value);
			_size = s;
			_value = malloc(_size);
			if (!_value) {
				std::stringstream msg;
				msg << "Failure allocating " << _size << " bytes" << std::endl;
				LOG(L_ERROR, msg.str());
				throw std::bad_alloc();
			}
		}
		memcpy(_value, v, s);
		_event = e;
	}

  private:
	/// Typesize
	size_t _size;

	/// Value
	void* _value;

	/// Last writing event
	cl_event _event;
};

/** @class Kernel Kernel.h CalcServer/Kernel.h
 * @brief A tool consisting in an OpenCL kernel execution. The variables used
 * in the OpenCL kernel are automatically detected.
 */
class Kernel : public Aqua::CalcServer::Tool
{
  public:
	/** Constructor.
	 * @param tool_name Tool name.
	 * @param kernel_path Kernel path.
	 * @param n Number of threads to launch. An empty string to autocompute it
	 * from the length of the output array variables
	 * @param once Run this tool just once. Useful to make initializations.
	 */
	Kernel(const std::string tool_name,
	       const std::string kernel_path,
	       const std::string entry_point = "entry",
	       const std::string n = "",
	       bool once = false);

	/** Destructor
	 */
	virtual ~Kernel();

	/** Initialize the tool.
	 * @return false if all gone right, true otherwise.
	 */
	void setup();

	/** Get the kernel file path.
	 * @return Tool kernel file path.
	 */
	const std::string path() { return (const std::string)_path; }

	/** Get the work group size
	 * @return Work group size
	 */
	size_t workGroupSize() const { return _work_group_size; }

	/** Get the work group size
	 * @return Work group size
	 */
	size_t globalWorkSize() const { return _global_work_size; }

  protected:
	/** Execute the tool
	 * @param events List of events that shall be waited before safe execution
	 * @return OpenCL event to be waited before accessing the dependencies
	 */
	cl_event _execute(const std::vector<cl_event> events);

	/** Compile the OpenCL program
	 * @param entry_point Program entry point method.
	 * @param flags Compiling additional flags.
	 * @param header Header to be append at the start of the source code.
	 * @return false if all gone right, true otherwise.
	 */
	void make(const std::string entry_point = "entry",
	          const std::string flags = "",
	          const std::string header = "");

	/** Compute the variables required by the program
	 * @param entry_point Program entry point method.
	 * @return false if all gone right, true otherwise.
	 */
	void variables(const std::string entry_point = "main");

	/** Compute the global work size
	 */
	void computeGlobalWorkSize();

  private:
	/// Kernel path
	std::string _path;

	/// Kernel entry point
	std::string _entry_point;

	/// Number of threads expression
	std::string _n;

	/// OpenCL kernel
	cl_kernel _kernel;

	/// work group size
	size_t _work_group_size;

	/// global work size
	size_t _global_work_size;

	/// List of dependencies, in the same order they have as kernel arguments
	std::vector<InputOutput::Variable*> _vars;

	/// Arguments setter for the kernel
	ArgSetter* _args_setter;
};

}
} // namespace

#endif // KERNEL_H_INCLUDED
