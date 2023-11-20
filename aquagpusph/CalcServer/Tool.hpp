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

#ifndef TOOL_H_INCLUDED
#define TOOL_H_INCLUDED

#include "aquagpusph/sphPrerequisites.hpp"
#include "aquagpusph/Variable.hpp"
#include "aquagpusph/AuxiliarMethods.hpp"
#include <math.h>
#include <vector>
#include <tuple>
#include <deque>
#include <algorithm>
#include <numeric>

#ifndef N_PROFILE_SAMPLES
#define N_PROFILE_SAMPLES 100
#endif

namespace Aqua {
namespace CalcServer {

/** @class Named Tool.h CalcServer/Tool.h
 * @brief A helper class for named objects
 */
class Named
{
  public:
	/** Constructor.
	 * @param name The name
	 */
	Named(const std::string name)
	  : _name(name)
	{
	}

	/** Destructor
	 */
	~Named() {}

	/** Set the tool name.
	 * @param tool_name Tool name.
	 */
	inline void name(const std::string tool_name) { _name = tool_name; };

	/** @brief Get the name
	 * @return The name
	 */
	inline const std::string name() const { return _name; }

  private:
	/// Name
	std::string _name;
};

/** @brief Get the average and standard deviation values
 * @param v Vector of values
 * @return the average and standard deviation values
 */
template<typename T>
inline std::tuple<T, T>
stats(const std::deque<T>& v)
{
	if (!v.size())
		return { 0, 0 };

	const T sum = std::accumulate(v.begin(), v.end(), 0.0);
	const T mean = sum / v.size();

	const T sum2 = std::inner_product(v.begin(), v.end(), v.begin(), 0.0);
	const T std = std::sqrt(sum2 / v.size() - mean * mean);

	return { mean, std };
}

/** @class Profile Tool.h CalcServer/Tool.h
 * @brief Profiler subinstance base class
 */
class Profile : public Aqua::CalcServer::Named
{
  public:
	/** Constructor
	 * @param name The name
	 * @param n Number of samples to keep in memory
	 */
	Profile(const std::string name, const unsigned int n = N_PROFILE_SAMPLES)
	  : Named(name)
	  , _n(n)
	{
	}

	/** Destructor
	 */
	virtual ~Profile() {}

	/** @brief Get the elapsed time and its standard deviation
	 * @return The average elapsed time and its standard deviation
	 */
	virtual inline std::tuple<cl_ulong, cl_ulong> elapsed() const
	{

		return stats(elapsed_times());
	}

	/** @brief Get the list of beginnings at each sample
	 * @return The list of beginnings
	 */
	virtual inline std::deque<cl_ulong> begin() const { return _start; }

	/** @brief Get the list of ends at each sample
	 * @return The list of ends
	 */
	virtual inline std::deque<cl_ulong> end() const { return _end; }

  protected:
	/** @brief Add profiling info
	 * @param start Starting timestamp
	 * @param end Ending timestamp
	 */
	inline void sample(cl_ulong start, cl_ulong end)
	{
		_start.push_back(start);
		_end.push_back(end);
		trim();
	}

	/** @brief Get the list of elapsed times
	 * @return The elapsed time for all the took samples
	 */
	inline std::deque<cl_ulong> elapsed_times() const
	{
		std::deque<cl_ulong> v;
		for (unsigned int i = 0; i < _start.size(); i++)
			v.push_back(_end[i] - _start[i]);
		return v;
	}

  private:
	/** @brief Drop the older samples until the required length is obtained
	 */
	inline void trim()
	{
		while (_start.size() > _n) {
			_start.pop_front();
			_end.pop_front();
		}
	}

	/// Number of time stamps to store
	unsigned int _n;

	/// The beginning time stamps in nanoseconds (FIFO list)
	std::deque<cl_ulong> _start;

	/// The ending time stamps in nanoseconds (FIFO list)
	std::deque<cl_ulong> _end;
};

/** @class Profiler Tool.h CalcServer/Tool.h
 * @brief Profiling base class
 *
 * Each tool is a Profiler, which might has several subinstances. Each
 * subinstance counts when each part of the tool computation began and end.
 * This class just collect all those subinstances together
 */
class Profiler
{
  public:
	/** Constructor.
	 * @param instances The list of subinstances that conforms this profiler
	 */
	Profiler() {}

	/** Destructor
	 */
	~Profiler()
	{
		for (auto i : _instances)
			delete i;
	}

	/** @brief Get the subinstances
	 * @return The list of subinstances
	 */
	inline std::vector<Profile*> subinstances() const { return _instances; }

	/** @brief Get the elapsed time as well as its standard deviation
	 *
	 * The elapsed time is the total time took by the tool to compute. Gaps
	 * between instances (the so called overhead) are not included
	 * @return Elapsed time with its standard deviation
	 */
	std::tuple<cl_ulong, cl_ulong> elapsed() const;

	/** @brief Get the list of beginnings at each sample
	 * @return The list of beginnings
	 * @warning It is assumed that the first subinstance is marking the
	 * beginning
	 */
	inline std::deque<cl_ulong> begin() const
	{
		if (!_instances.size())
			return {};
		return _instances.front()->begin();
	}

	/** @brief Get the list of ends at each sample
	 * @return The list of ends
	 * @warning It is assumed that the last subinstance is marking the
	 * end
	 */
	inline std::deque<cl_ulong> end() const
	{
		if (!_instances.size())
			return {};
		return _instances.back()->end();
	}

	/** @brief Get the total time took as well as its standard deviation
	 *
	 * The total time took might differ from the
	 * Aqua::CalcServer::Profiler::elapsed() on the overhead
	 * @return Elapsed time with its standard deviation
	 */
	std::tuple<cl_ulong, cl_ulong> total() const;

	/** @brief Get the delta time
	 * @param t Timer
	 * @param t0 Timer reference to subtract
	 * @return Time delta
	 */
	static inline cl_long delta(const cl_ulong& t, const cl_ulong& t0)
	{
		return (t > t0) ? t - t0 : -(cl_long)(t0 - t);
	}

	/** @brief Get the delta time
	 * @param t List of timers
	 * @param t0 Timer reference to subtract
	 * @return List of time deltas
	 */
	static inline std::deque<cl_long> delta(const std::deque<cl_ulong> t,
	                                        const cl_ulong& t0)
	{
		std::deque<cl_long> dt;
		for (auto val : t)
			dt.push_back(delta(val, t0));
		return dt;
	}

	/** @brief Get the delta time
	 * @param t List of timers
	 * @param t0 List of timer references to subtract
	 * @return List of time deltas
	 */
	static inline std::deque<cl_long> delta(const std::deque<cl_ulong> t,
	                                        const std::deque<cl_ulong> t0)
	{
		std::deque<cl_long> dt;
		for (unsigned int i = 0; i < t.size(); i++)
			dt.push_back(delta(t[i], t0[i]));
		return dt;
	}

  protected:
	/** @brief Set the tool subinstances
	 * @param instances The subinstances
	 */
	inline void subinstances(std::vector<Profile*> instances)
	{
		_instances = instances;
	}

  private:
	/// List of subinstances
	std::vector<Profile*> _instances;
};

/** @class Tool Tool.h CalcServer/Tool.h
 * @brief Tools base class.
 *
 * The way that AQUAgpusph compute each problem is set through a set of tools
 * that are computed sequentially. Several tools can be considered, for
 * instance:
 *   -# A single OpenCL kernel
 *   -# A more complex OpenCL tool like Reductions or LinkList
 *   -# Python scripts
 *   -# Variables set
 */
class Tool
  : public Aqua::CalcServer::Named
  , public Aqua::CalcServer::Profiler
{
  public:
	/** Constructor.
	 * @param tool_name Name of the tool. Useful to identify errors.
	 * @param once Run this tool just once. Useful to make initializations.
	 */
	Tool(const std::string tool_name, bool once = false);

	/** Destructor
	 */
	virtual ~Tool();

	/** Initialize the tool.
	 */
	virtual void setup();

	/** @brief Execute the tool measuring the elapsed time.
	 *
	 * Actually this method is just ensuring that the tool can be executed,
	 * e.g. the tool has been already executed, but it is asked to be ran just
	 * once.
	 * If the tool can be executed, then _execute() method is called, measuring
	 * the time required to carry out the task.
	 * @return false if all gone right, true otherwise.
	 * @note Usually you don't want to overload this method, but the _execute()
	 * protected one.
	 */
	virtual void execute();

	/** Get the next tool to be executed in the pipeline.
	 *
	 * Such tool is usually just the next one in the linearized tools chain.
	 * However, conditional tools may alter the flow
	 * @return Next tool to be executed. NULL if this is the last tool of the
	 * pipeline
	 */
	virtual inline Tool* next_tool() { return _next_tool; }

	/** Get the allocated memory for this tool.
	 * @return allocated memory by this tool.
	 */
	inline size_t allocatedMemory() const { return _allocated_memory; }

	/** Get the number of times that this tool has been called.
	 * @return Number of times this tool has been called.
	 */
	inline unsigned int used_times() const { return _n_iters; }

	/** @brief Add new data to the average and squared elapsed times
	 * @param elapsed_time Elapsed time
	 */
	void addElapsedTime(float elapsed_time);

	/** Get the time consumed by the tool.
	 * @param averaged true if the avergaed time step is required, false
	 * otherwise.
	 * @return time consumed.
	 */
	inline float elapsedTime(bool averaged = true) const
	{
		if (!averaged)
			return _elapsed_time;
		return _average_elapsed_time;
	}

	/** Get the time consumed variance.
	 * @return Time consumed variance.
	 */
	inline float elapsedTimeVariance() const
	{
		return _squared_elapsed_time - pow(_average_elapsed_time, 2);
	}

	/** Get the time consumed standard deviation.
	 * @return Time consumed standard deviation.
	 */
	inline float elapsedTimeDeviation() const
	{
		return sqrt(elapsedTimeVariance());
	}

	/** Get the scope modifier
	 *
	 * Scopes can be used to create groups of tools that can be eventually
	 * enabled/disabled in runtime. This is sueful to create conditions.
	 *
	 * @return 0 if this tool is not modifying the scope, 1 if this tool is
	 * creating a new subscope, and -1 if this tool is closing an already
	 * created subscope, returning to the previous one
	 * @note scopes shall be always balanced
	 */
	virtual const int scope_modifier() { return 0; }

	/** @brief Get the input depedencies of the tool
	 * @return The input and output dependencies
	 */
	inline std::vector<InputOutput::Variable*> getInputDependencies() const
	{
		return _in_vars;
	}

	/** @brief Get the output depedencies of the tool
	 * @return The input and output dependencies
	 */
	inline std::vector<InputOutput::Variable*> getOutputDependencies() const
	{
		return _out_vars;
	}

	/** @brief Get the depedencies of the tool
	 * @return The input and output dependencies
	 */
	inline std::tuple<std::vector<InputOutput::Variable*>,
	                  std::vector<InputOutput::Variable*>>
	getDependencies() const
	{
		return { _in_vars, _out_vars };
	}

	/** @brief Get the tool event
	 *
	 * That event will mark when the tool has already finished
	 * @return The event
	 */
	inline cl_event getEvent() const { return _event; }

  protected:
	/** Get the tool index in the pipeline
	 * @return Index of the tool in the pipeline. -1 if the tool cannot be find
	 */
	int id_in_pipeline();

	/** Set the next tool to be executed in the pipeline.
	 * @param tool Next tool to be executed. NULL if this is the last tool of
	 * the pipeline
	 */
	void next_tool(Tool* tool) { _next_tool = tool; }

	/** Set the allocated memory for this tool.
	 * @param mem_size allocated memory by this tool.
	 */
	void allocatedMemory(size_t mem_size) { _allocated_memory = mem_size; }

	/** Execute the tool
	 * @param events List of events that shall be waited before safe execution
	 * @return OpenCL event to be waited before accessing the dependencies
	 */
	virtual cl_event _execute(const std::vector<cl_event> events)
	{
		return NULL;
	}

	/** @brief Set the depedencies of the tool
	 *
	 * The dependencies are the variables that this tool is either reading or
	 * writing.
	 *
	 * @param inputs Names of the input variables
	 * @param outputs Names of the output variables
	 */
	inline void setDependencies(std::vector<std::string> inputs,
	                            std::vector<std::string> outputs)
	{
		_in_vars = namesToVars(inputs);
		_out_vars = namesToVars(outputs);
	}

	/** @brief Set the depedencies of the tool
	 *
	 * The dependencies are the variables that this tool is either reading or
	 * writing.
	 *
	 * @param inputs Input variables
	 * @param outputs Output variables
	 */
	inline void setDependencies(std::vector<InputOutput::Variable*> inputs,
	                            std::vector<InputOutput::Variable*> outputs)
	{
		_in_vars = inputs;
		_out_vars = outputs;
	}

	/** @brief Set the depedencies of the tool
	 *
	 * The dependencies are the variables that this tool is either reading or
	 * writing.
	 * @warning This function will treat all the dependencies as outputs
	 * @param vars Output variable names
	 */
	inline void setDependencies(std::vector<std::string> vars)
	{
		std::vector<std::string> inputs;
		setDependencies(inputs, vars);
	}

	/** @brief Set the reading depedencies of the tool
	 *
	 * The dependencies are the variables that this tool is either reading or
	 * writing.
	 *
	 * @param inputs Names of the input variables
	 */
	inline void setInputDependencies(std::vector<std::string> vars)
	{
		_in_vars = namesToVars(vars);
	}

	/** @brief Set the reading depedencies of the tool
	 *
	 * The dependencies are the variables that this tool is either reading or
	 * writing.
	 *
	 * @param inputs Input variables
	 */
	inline void setInputDependencies(std::vector<InputOutput::Variable*> vars)
	{
		_in_vars = vars;
	}

	/** @brief Set the writing depedencies of the tool
	 *
	 * The dependencies are the variables that this tool is either reading or
	 * writing.
	 *
	 * @param inputs Names of the output variables
	 */
	inline void setOutputDependencies(std::vector<std::string> vars)
	{
		_out_vars = namesToVars(vars);
	}

	/** @brief Set the writing depedencies of the tool
	 *
	 * The dependencies are the variables that this tool is either reading or
	 * writing.
	 *
	 * @param inputs Output variables
	 */
	inline void setOutputDependencies(std::vector<InputOutput::Variable*> vars)
	{
		_out_vars = vars;
	}

	/** @brief Set the depedencies of the tool
	 *
	 * The dependencies are the variables that this tool is either reading or
	 * writing.
	 * @warning This function will treat all the dependencies as outputs
	 * @param vars Output variables
	 */
	inline void setDependencies(std::vector<InputOutput::Variable*> vars)
	{
		std::vector<InputOutput::Variable*> inputs;
		setDependencies(inputs, vars);
	}

	/** @brief Compile an OpenCL source code and generate the corresponding
	 * kernel
	 *
	 * With this method several operations are carried out at the same time.
	 * First the program is compiled and linked. Afterwards, the required
	 * kernels are extracted, and the program object is released
	 *
	 * @param source Source code to be compiled
	 * @param names Function names to be extracted in the kernel
	 * @param flags Additional compilation flags. Some flags are used by
	 * default:
	 *   - -DDEBUG/-DNDEBUG depending on whether DEBUG mode is enabled or not
	 *   - -cl-mad-enable -cl-fast-relaxed-math
	 *   - -DHAVE_2D/-DHAVE_3D depending on whether 2D or 3D is considered
	 * @return Kernel instances
	 */
	static std::vector<cl_kernel> compile(const std::string source,
	                                      const std::vector<std::string> names,
	                                      const std::string flags = "");

	/** @brief Compile an OpenCL source code and generate the corresponding
	 * kernel
	 *
	 * With this method several operations are carried out at the same time.
	 * First the program is compiled and linked. Afterwards, the required kernel
	 * is extracted, and the program object is released
	 *
	 * @param source Source code to be compiled
	 * @param kernel_name Function name to be extracted in the kernel
	 * @param flags Additional compilation flags. Some flags are used by
	 * default:
	 *   - -DDEBUG/-DNDEBUG depending on whether DEBUG mode is enabled or not
	 *   - -cl-mad-enable -cl-fast-relaxed-math
	 *   - -DHAVE_2D/-DHAVE_3D depending on whether 2D or 3D is considered
	 * @return Kernel instance
	 */
	static cl_kernel compile_kernel(const std::string source,
	                                const std::string kernel_name,
	                                const std::string flags = "");

  private:
	/** @brief Get the list of events that this tool shall wait for
	 * @return List of events
	 * @note The events returned have been retained, so call clReleaseEvent()
	 * after using them
	 */
	const std::vector<cl_event> getEvents() const;

	/** @brief Get a list of variables from their names
	 * @param names List of names
	 * @return List of varaibles
	 */
	std::vector<InputOutput::Variable*> namesToVars(
	    const std::vector<std::string>& names) const;

	/// true if the tool shall be run just once, false otherwise
	bool _once;

	/// Next tool in the execution pipeline
	Tool* _next_tool;

	/// Total auxiliar memory allocated in the device
	size_t _allocated_memory;

	/// Times that this tool has been called
	unsigned int _n_iters;

	/// Average elapsed time
	float _elapsed_time;

	/// Average elapsed time
	float _average_elapsed_time;

	/// Average squared elapsed time
	float _squared_elapsed_time;

	/// List of input dependencies
	std::vector<InputOutput::Variable*> _in_vars;

	/// List of output dependencies
	std::vector<InputOutput::Variable*> _out_vars;

	/// Output event of the tool
	cl_event _event;
};

}
} // namespace

#endif // TOOL_H_INCLUDED
