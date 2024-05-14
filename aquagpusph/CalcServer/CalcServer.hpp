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
 * @brief The calculation main entry point.
 * (See Aqua::CalcServer::CalcServer for details)
 */

#ifndef CALCSERVER_H_INCLUDED
#define CALCSERVER_H_INCLUDED

#include <vector>
#include <deque>
#include <map>
#include <string>
#include <iterator>
#include <chrono>

#include "aquagpusph/sphPrerequisites.hpp"
#include "aquagpusph/ProblemSetup.hpp"
#include "aquagpusph/TimeManager.hpp"
#include "aquagpusph/Variable.hpp"
#include "aquagpusph/Singleton.hpp"
#include "Tool.hpp"

#ifndef N_PROFILING_SNAPSHOTS
#define N_PROFILING_SNAPSHOTS 2
#endif

#define CHECK_OCL_OR_THROW(err_code, msg)                                      \
	if (err_code != CL_SUCCESS) {                                              \
		LOG(L_ERROR, std::string(msg) + "\n");                                 \
		Aqua::InputOutput::Logger::singleton()->printOpenCLError(err_code);    \
		throw std::runtime_error("OpenCL error");                              \
	}

namespace Aqua {
/// @namespace Aqua::CalcServer Calculation server name space.
namespace CalcServer {

class UnSort;
typedef InputOutput::ProblemSetup::sphSettings::debug_opts debug_opts;

/** @class CalcServer CalcServer.h CalcServer.h
 * @brief Exception raised when the user manually interrupts the simulation.
 *
 * The target of this exception is handling the users interrumptions without
 * returning back a core dump in the command line, which might be wrongly
 * interpreted as an error.
 */
class user_interruption : public std::runtime_error
{
  public:
	/// Constructor
	user_interruption(const std::string msg)
	  : std::runtime_error(msg){};
};

/// Profiling sample
typedef struct _ProfilingSample
{
	/// The triggering tool
	Tool* tool;
	/// The name of the sample
	std::string name;
	/// The starting timer
	cl_ulong start;
	/// The ending timer
	cl_ulong end;
} ProfilingSample;

/// Profiling snapshot
typedef struct _ProfilingSnapshot
{
	/// The profiling step
	cl_ulong step;
	/// The list of samples
	std::vector<ProfilingSample> samples;
} ProfilingSnapshot;

/** @class ProfilingInfo CalcServer.h CalcServer.h
 * @brief A FIFO list of profiling snapshots
 */
class ProfilingInfo
{
  public:
	/** @brief Constructor.
	 * @param n Number of snapshots to keep
	 */
	ProfilingInfo(const cl_uint n = N_PROFILING_SNAPSHOTS)
	  : _step(0)
	  , _n(n)
	{
	}

	/// Destructor
	~ProfilingInfo() {}

	/** @brief Get the number of times the tools pack has been executed
	 *
	 * This might not match InputOutput::TimeManager::step(). This is actually
	 * meant for profiling purposes
	 */
	inline cl_ulong step() const { return _step; }

	/** @brief Add a new sample to an specific step
	 * @param step The profiling step
	 * @param tool Triggering tool
	 * @param name Sample name
	 * @param start Starting timer
	 * @param end Ending timer
	 */
	void sample(cl_ulong step,
	            Tool* tool,
	            std::string name,
	            cl_ulong start,
	            cl_ulong end);

	/** @brief Add a new sample to an specific step
	 * @param step The profiling step
	 * @param tool Triggering tool
	 * @param substage Tool computation substage
	 * @param start Starting timer
	 * @param end Ending timer
	 */
	inline void sample(cl_ulong step,
	                   Tool* tool,
	                   Profile* substage,
	                   cl_ulong start,
	                   cl_ulong end)
	{
		std::string name = tool->name() + "::" + substage->name();
		auto parent = tool;
		while (parent = parent->parent()) {
			name = parent->name() + "::" + name;
		}
		sample(step, tool, name, start, end);
	}

	/** @brief Get the stored snapshots
	 * @return The list of stored snapshots
	 */
	std::deque<ProfilingSnapshot> get() const;

	/** @brief Get the delta time
	 *
	 * Whereas the timers are big unsigned integers, the delta is conversely
	 * a relatively small signed integer, although we keep the original 64 bits
	 * size
	 * @param t Timer
	 * @param t0 Timer reference to subtract
	 * @return Time delta
	 */
	static inline cl_long delta(const cl_ulong& t, const cl_ulong& t0)
	{
		return (t > t0) ? t - t0 : -(cl_long)(t0 - t);
	}

  protected:
	/** @brief Let the profiler know that a new step of tools execution is
	 * about to start
	 */
	void newStep();

  private:
	/// A counter on the times the tools pack have been called
	cl_ulong _step;
	/// Maximum number of snapshots
	cl_uint _n;
	/// List of snaphsots
	std::deque<ProfilingSnapshot> _snapshots;
};

/** @class CalcServer CalcServer.h CalcServer.h
 * @brief Entity that perform the main work of the simulation.
 *
 * In the Aqua::CalcServer::CalcServer a time subloop is performed where the
 * SPH simulation is performed while no output files should be updated.
 * @note Updating output files require to download data from the server, which
 * due to the low bandwidth asigned is usually a bottleneck, hence letting the
 * Aqua::CalcServer::CalcServer works without interrumptions is the best
 * optimization technique.
 * @remarks Some output files are managed internally by this class, like the
 * log file, the energy file, pressure sensors file, or bounds file.
 */
class CalcServer
  : public Aqua::Singleton<Aqua::CalcServer::CalcServer>
  , public Aqua::CalcServer::ProfilingInfo
{
  public:
	/** @brief Constructor.
	 * @param sim_data Simulation data read from XML files
	 */
	CalcServer(const Aqua::InputOutput::ProblemSetup& sim_data);

	/// Destructor
	~CalcServer();

	/** @brief Raise a SIGINT/SIGTERM signal
	 *
	 * This can be useful for instance to asynchronously interrupt the
	 * execution after detecting errors
	 * @note This function will never force the execution termination. If the
	 * signal has been already rised somewhere else, this function will just
	 * gently wait for the execution to be finished
	 */
	void raiseSIGINT();

	/** @brief Internal time loop.
	 *
	 * Calculation server will be iterating while no output files should be
	 * updated (or even the simulation is finished).
	 * @param t_manager Time manager to let the calculation server when shall
	 * stop the internal loop.
	 */
	void update(InputOutput::TimeManager& t_manager);

	/// Setup some additional simulation data.
	/** Even thought this work is associated with the constructor CalcServer(),
	 * when something may fail it is preferable to let it to a separated method
	 * that could report errors, allowing the program to deal with them.
	 */
	void setup();

	/** Get the variables manager
	 * @return Variables manager
	 */
	inline InputOutput::Variables* variables() { return &_vars; }

	/** Get the definitions registered.
	 * @return List of definitions.
	 */
	inline std::vector<std::string> definitions() const { return _definitions; }

	/** Get the tools registered.
	 * @return List of tools.
	 */
	inline std::vector<Tool*> tools() const { return _tools; }

	/** Get the active context
	 * @return OpenCL context
	 */
	inline cl_context context() const { return _context; }

	/** Get the platform
	 * @return OpenCL platform
	 */
	inline cl_platform_id platform() const { return _platform; }

	/** Get the device
	 * @return OpenCL device
	 */
	inline cl_device_id device() const { return _device; }

	enum cmd_queue
	{
		cmd_queue_current,
		cmd_queue_new,
		cmd_queue_mpi
	};

	/** Get a command queue
	 *
	 * To allow the asynchronous "As-Fast-As-Possible" execution of commands
	 * they can be enqueued on different queues. The queues are always reseted
	 * when a new time step is executed.
	 * @param which cmd_queue::cmd_queue_current to enqueue the command on the
	 * current command queue, cmd_queue::cmd_queue_new to raise a new one where
	 * commands can be executed in parallel, cmd_queue::cmd_queue_mpi to get
	 * the command queue designed to exchange info with other processes
	 * @return The command queue
	 * @raise std::runtime_error if the command queue cannot be created
	 */
	cl_command_queue command_queue(
	    cmd_queue which = cmd_queue::cmd_queue_current);

	/** Call clFinish() on top of all known command queues
	 * @return CL_SUCCESS if clFinish() is correctly called on all queues, an
	 * error code otherwise
	 */
	inline cl_int finish() const
	{
		cl_int err_code = CL_SUCCESS;
		for (auto queue : _command_queues)
			err_code |= clFinish(queue);
		err_code |= clFinish(_command_queue_parallel);
		return err_code;
	}

	/** Create a marker to be set as completed when all the enqued commnands on
	 * all the command queues are finished.
	 * @return the master event
	 */
	cl_event marker() const;

	/** @brief Get the offset between the device timer and the host one
	 * @return The time offset, in nanoseconds
	 */
	inline cl_ulong device_timer_offset() const { return _device_timer_offset; }

	/** @brief Get the host timer
	 * @return The host timer in nanoseconds
	 */
	static inline cl_ulong host_timer()
	{
		auto now = std::chrono::system_clock::now();
		return (cl_ulong)std::chrono::duration_cast<std::chrono::nanoseconds>(
		           now.time_since_epoch())
		    .count();
	}

	/** @brief Download a unsorted variable from the device.
	 * @param var_name Variable to unsort and download.
	 * @param offset The offset in bytes in the memory object to read from.
	 * @param cb The size in bytes of data being downloaded.
	 * @param ptr The host memory where the data should be copied
	 * @return The data download event, NULL if errors are detected.
	 * @note The caller must wait for the events (clWaitForEvents) before
	 * accessing the downloaded data.
	 * @remarks The caller must call clReleaseEvent to destroy the event.
	 * Otherwise a memory leak can be expected.
	 */
	cl_event getUnsortedMem(const std::string var_name,
	                        size_t offset,
	                        size_t cb,
	                        void* ptr);

	/** @brief Get the AQUAgpusph root path.
	 * @return AQUAgpusph root path
	 */
	const std::string base_path() const { return _base_path.c_str(); }

	/** @brief Report if any tools debug option is enabled.
	 * @see Aqua::InputOutput::ProblemSetup::sphSettings::debug_tools
	 * @return true if any tools debug mode option is enabled, false otherwise
	 */
	inline bool debug_mode() const
	{
		return _sim_data.settings.debug_tools;
	}

	/** @brief Report if the tools debug mode option is enabled.
	 * @param opt The options to check
	 * @see Aqua::InputOutput::ProblemSetup::sphSettings::debug_tools
	 * @return true if the tools debug mode option is enabled, false otherwise
	 */
	inline bool debug_mode(debug_opts opt) const
	{
		return _sim_data.settings.debug_tools & opt;
	}

	/** @brief Get the variables associated to an event, if any
	 * @return The list of variables that depends on the event, both for
	 * writing and reading
	 */
	std::tuple<std::vector<InputOutput::Variable*>,
	           std::vector<InputOutput::Variable*>>
	eventVars(cl_event event) const;

	/** @brief Get the tool that emitted an event, if any
	 * @return The tool that created the event. Null if no tool is found
	 */
	Tool* eventTool(cl_event event) const;

  private:
	/** Setup the OpenCL stuff.
	 */
	void setupOpenCL();
	/** Prints all the available platforms and devices returned by OpenCL.
	 */
	void queryOpenCL();
	/** Get a platform from the available ones.
	 */
	void setupPlatform();
	/** Get the available devices in the selected platform.
	 */
	void setupDevices();

	/// Number of available OpenCL platforms
	cl_uint _num_platforms;
	/// List of OpenCL platforms
	cl_platform_id* _platforms;
	/// Number of devices
	cl_uint _num_devices;
	/// List of devices
	cl_device_id* _devices;
	/// OpenCL context
	cl_context _context;
	/// Selected platform
	cl_platform_id _platform;
	/// Selected device
	cl_device_id _device;
	/// Current main command queue
	unsigned int _command_queue_current;
	/// Main command queues
	std::vector<cl_command_queue> _command_queues;
	/** Command queue for task carried out in parallel tasks (e.g. events
	 * callbacks).
	 */
	cl_command_queue _command_queue_parallel;

	/// The difference between the device timer and the system clock
	cl_ulong _device_timer_offset;

	/// User registered variables
	InputOutput::Variables _vars;

	/// User registered definitions
	std::vector<std::string> _definitions;

	/// User registered tools
	std::vector<Tool*> _tools;

	/** @brief AQUAgpusph root path.
	 *
	 * This path is added to the OpenCL include paths.
	 */
	std::string _base_path;

	/** @brief Currently executed tool/report.
	 *
	 * Useful to can report runtime OpenCL implementation errors (see
	 * clCreateContext).
	 *
	 * @note Since this data is automatically readed by OpenCL at the time of
	 * reporting errors, old-style C string is more convenient
	 */
	char* _current_tool_name;

	/** Map with the unsorter for each variable. Storing the unsorters should
	 * dramatically reduce the saving files overhead in some platforms
	 */
	std::map<std::string, UnSort*> unsorters;

	/// Simulation data read from XML files
	Aqua::InputOutput::ProblemSetup _sim_data;
};

}
} // namespace

#endif // CALCSERVER_H_INCLUDED
