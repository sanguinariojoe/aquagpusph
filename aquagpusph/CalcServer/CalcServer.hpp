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

#ifndef _ITEMS
/** @def _ITEMS
 * @brief Number of items in a group
 * @note Must be power of 2, and in some devices greather than 32.
 * @see Aqua::CalcServer::RadixSort
 */
#define _ITEMS 128
#endif

#ifndef _GROUPS
/** @def _GROUPS
 * @brief Number of groups
 * @note Must be power of 2, and the amount of data should be divisible by
 * _ITEMS*_GROUPS
 * @see Aqua::CalcServer::RadixSort
 */
#define _GROUPS 32
#endif

namespace Aqua {
/// @namespace Aqua::CalcServer Calculation server name space.
namespace CalcServer {

class UnSort;

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
class CalcServer : public Aqua::Singleton<Aqua::CalcServer::CalcServer>
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

	/** Get the command queue
	 *
	 * Two different command queues are available. Indeed, OpenCL specification
	 * allows to create the command queues with the
	 * CL_QUEUE_OUT_OF_ORDER_EXEC_MODE_ENABLE option, which allows that a
	 * subsequent kernel/task execution begins before the preceding one has
	 * already finished. However, the specification does not allows, by any
	 * means, that a preceding kernel/task execution can be triggered after
	 * a subsequent one. Therefore, parallel tasks shall use a different queue
	 * to ensure the correct order within each queue
	 *
	 * @param parallel true if the command queue for task executed in parallel
	 * is queried, false otherwise
	 * @return The command queue
	 */
	inline cl_command_queue command_queue(bool parallel = false) const
	{
		if (parallel)
			return _command_queue_parallel;
		return _command_queue;
	}

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
			now.time_since_epoch()).count();
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

	/** @brief Report if the tools debug mode is enabled.
	 * @see Aqua::InputOutput::ProblemSetup::sphSettings::debug_tools
	 * @return true if the tools debug mode is enabled, false otherwise
	 */
	inline bool debug_mode() const { return _sim_data.settings.debug_tools; }

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
	/// Main command queue
	cl_command_queue _command_queue;
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

  private:
	/// Simulation data read from XML files
	Aqua::InputOutput::ProblemSetup _sim_data;
};

}
} // namespace

#endif // CALCSERVER_H_INCLUDED
