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

#include <stdlib.h>
#include <dlfcn.h>
#include <limits>
#include <string>
#include <stack>
#include <assert.h>
#include <signal.h>
#include <tuple>
#include <mutex>

#include "CalcServer.hpp"
#include "aquagpusph/AuxiliarMethods.hpp"
#include "aquagpusph/InputOutput/Logger.hpp"
#include "Assert.hpp"
#include "Conditional.hpp"
#include "Copy.hpp"
#include "Kernel.hpp"
#include "LinkList.hpp"
#include "Python.hpp"
#include "RadixSort.hpp"
#include "Sort.hpp"
#include "Reduction.hpp"
#include "Set.hpp"
#include "SetScalar.hpp"
#include "UnSort.hpp"
#include "Reports/Performance.hpp"
#include "Reports/Screen.hpp"
#include "Reports/TabFile.hpp"
#include "Reports/SetTabFile.hpp"

#ifdef HAVE_MPI
#include <mpi.h>
#include "MPISync.hpp"
#endif

namespace Aqua {
namespace CalcServer {

/// @brief Have been a SIGINT already registered?
static bool sigint_received = false;

/** @brief Handle SIGINT signals
 *
 * The first time a SIGINT is received, Aqua::CalcServer::sigint_received is set
 * to true, such that, at the end of the current time step the simulation will
 * stop, the last output will be printed, and the resources will be correctly
 * released.
 *
 * If SIGINT is received twice, then this handler will enforce the inmediate
 * program exit.
 *
 * @param s Recevied signal, SIGINT
 */
void
sigint_handler(int s)
{
	// Log the reception, and afterwards the processing. That way, in case of
	// MPI jobs we can know if some uncoordinated processes have failed to
	// correctly finish the job
	LOG(L_WARNING, "SIGINT/SIGTERM received\n");
	if (sigint_received) {
		// The user asked more than once to stop the simulation, force it
		LOG(L_ERROR, "Forced program exit (SIGINT/SIGTERM received twice)\n");
		exit(EXIT_FAILURE);
	}
	sigint_received = true;
}

/** @brief A mutex to avoid several threads writing profiling smaples
 * simultaneously
 */
std::mutex profiling_mutex;

std::deque<ProfilingSnapshot>
ProfilingInfo::get() const
{
	const std::lock_guard<std::mutex> lock(profiling_mutex);
	return _snapshots;
}

void
ProfilingInfo::sample(cl_ulong step,
                      Tool* tool,
                      std::string name,
                      cl_ulong start,
                      cl_ulong end)
{
	const std::lock_guard<std::mutex> lock(profiling_mutex);
	/// Traverse the snapshots looking for the step one
	for (auto it = _snapshots.rbegin();  it != _snapshots.rend(); ++it ) {
		if ((*it).step == step) {
			(*it).samples.push_back({tool, name, start, end});
			return;
		}
	}
}

void
ProfilingInfo::newStep()
{
	const std::lock_guard<std::mutex> lock(profiling_mutex);
	_step++;
	ProfilingSnapshot snapshot;
	snapshot.step = _step;
	_snapshots.push_back(snapshot);
	if (_snapshots.size() > _n)
		_snapshots.pop_front();
}


CalcServer::CalcServer(const Aqua::InputOutput::ProblemSetup& sim_data)
  : _num_platforms(0)
  , _platforms(NULL)
  , _num_devices(0)
  , _devices(NULL)
  , _context(NULL)
  , _platform(NULL)
  , _device(NULL)
  , _command_queue(NULL)
  , _command_queue_parallel(NULL)
  , _device_timer_offset(0)
  , _current_tool_name(NULL)
  , _sim_data(sim_data)
{
	unsigned int i, j;

	setupOpenCL();

	_base_path = _sim_data.settings.base_path;
	_current_tool_name = new char[256];
	memset(_current_tool_name, '\0', 256 * sizeof(char));

	unsigned int num_sets = _sim_data.sets.size();
	unsigned int N = 0;
	for (auto set : _sim_data.sets) {
		N += set->n();
	}

	const unsigned int num_icell = nextPowerOf2(roundUp(N, _ITEMS * _GROUPS));

	// Register default scalars
	std::ostringstream valstr;
	int mpi_rank = 0, mpi_size = 1;
#ifdef HAVE_MPI
	try {
		mpi_rank = MPI::COMM_WORLD.Get_rank();
		mpi_size = MPI::COMM_WORLD.Get_size();
	} catch (MPI::Exception e) {
		std::ostringstream msg;
		msg << "Error getting MPI rank and size. " << std::endl
		    << e.Get_error_code() << ": " << e.Get_error_string() << std::endl;
		LOG(L_ERROR, msg.str());
		throw;
	}
	assert(mpi_rank >= 0);
	assert(mpi_size > 0);
#endif
	valstr.str("");
	valstr << mpi_rank;
	_vars.registerVariable("mpi_rank", "unsigned int", "", valstr.str());
	valstr.str("");
	valstr << mpi_size;
	_vars.registerVariable("mpi_size", "unsigned int", "", valstr.str());
#ifdef HAVE_3D
	unsigned int dims = 3;
#else
	unsigned int dims = 2;
#endif
	valstr.str("");
	valstr << dims;
	_vars.registerVariable("dims", "unsigned int", "", valstr.str());
	_vars.registerVariable("t", "float", "", "0");
	_vars.registerVariable("dt", "float", "", "0");
	_vars.registerVariable("iter", "unsigned int", "", "0");
	_vars.registerVariable("frame", "unsigned int", "", "0");
	valstr.str("");
	valstr << std::numeric_limits<float>::max();
	_vars.registerVariable("end_t", "float", "", valstr.str());
	valstr.str("");
	valstr << std::numeric_limits<unsigned int>::max();
	_vars.registerVariable("end_iter", "unsigned int", "", valstr.str());
	_vars.registerVariable("end_frame", "unsigned int", "", valstr.str());
	valstr.str("");
	valstr << N;
	_vars.registerVariable("N", "unsigned int", "", valstr.str());
	valstr.str("");
	valstr << _sim_data.sets.size();
	_vars.registerVariable("n_sets", "unsigned int", "", valstr.str());
	valstr.str("");
	valstr << num_icell;
	_vars.registerVariable("n_radix", "unsigned int", "", valstr.str());
	// Number of cells in x, y, z directions, and the total (n_x * n_y * n_z)
	_vars.registerVariable("n_cells", "uivec4", "", "1, 1, 1, 1");
	// Kernel support
	_vars.registerVariable("support", "float", "", "2");

	// Register default arrays
	valstr.str("");
	valstr << N;
	_vars.registerVariable("id", "unsigned int*", valstr.str(), "");
	_vars.registerVariable("r", "vec*", valstr.str(), "");
	_vars.registerVariable("iset", "unsigned int*", valstr.str(), "");
	valstr.str("");
	valstr << num_icell;
	_vars.registerVariable("id_sorted", "unsigned int*", valstr.str(), "");
	_vars.registerVariable("id_unsorted", "unsigned int*", valstr.str(), "");
	_vars.registerVariable("icell", "unsigned int*", valstr.str(), "");
	_vars.registerVariable("ihoc", "unsigned int*", "n_cells_w", "");

	// Register the user variables and arrays
	for (i = 0; i < _sim_data.variables.names.size(); i++) {
		_vars.registerVariable(_sim_data.variables.names.at(i),
		                       _sim_data.variables.types.at(i),
		                       _sim_data.variables.lengths.at(i),
		                       _sim_data.variables.values.at(i));
	}

	// Register the user definitions
	for (i = 0; i < _sim_data.definitions.names.size(); i++) {
		valstr.str("");
		valstr << "-D" << _sim_data.definitions.names.at(i);
		if (_sim_data.definitions.values.at(i).compare("")) {
			if (!_sim_data.definitions.evaluations.at(i)) {
				valstr << "=" << _sim_data.definitions.values.at(i);
			} else {
				float defval = 0.f;
				_vars.solve(
				    "float", _sim_data.definitions.values.at(i), &defval);
				// We need to specify the format, to ensure that a decimal
				// number (i.e. a number with decimal point) is retrieved, so
				// C++ streamer can't be applied here
				char defvalstr[128];
				snprintf(defvalstr, sizeof(defvalstr), "%#G", defval);
				valstr << "=" << defvalstr << "f";
			}
		}
		_definitions.push_back(valstr.str());
	}

	// Register the tools
	for (auto t : _sim_data.tools) {
		bool once = false;
		if (!t->get("once").compare("true")) {
			once = true;
		}
		// Computation tools
		if (!t->get("type").compare("kernel")) {
			std::string tool_path = t->get("path");
			if (!isFile(tool_path) && isFile(_base_path + "/" + tool_path)) {
				tool_path = _base_path + "/" + tool_path;
			}
			Kernel* tool = new Kernel(t->get("name"),
			                          tool_path,
			                          t->get("entry_point"),
			                          t->get("n"),
			                          once);
			_tools.push_back(tool);
		} else if (!t->get("type").compare("copy")) {
			Copy* tool =
			    new Copy(t->get("name"), t->get("in"), t->get("out"), once);
			_tools.push_back(tool);
		} else if (!t->get("type").compare("python")) {
			Python* tool = new Python(t->get("name"), t->get("path"), once);
			_tools.push_back(tool);
		} else if (!t->get("type").compare("set")) {
			Set* tool =
			    new Set(t->get("name"), t->get("in"), t->get("value"), once);
			_tools.push_back(tool);
		} else if (!t->get("type").compare("set_scalar")) {
			SetScalar* tool = new SetScalar(
			    t->get("name"), t->get("in"), t->get("value"), once);
			_tools.push_back(tool);
		} else if (!t->get("type").compare("reduction")) {
			Reduction* tool = new Reduction(t->get("name"),
			                                t->get("in"),
			                                t->get("out"),
			                                t->get("operation"),
			                                t->get("null"),
			                                once);
			_tools.push_back(tool);
		} else if (!t->get("type").compare("link-list")) {
			LinkList* tool = new LinkList(t->get("name"), t->get("in"));
			_tools.push_back(tool);
		} else if (!t->get("type").compare("radix-sort")) {
			RadixSort* tool = new RadixSort(t->get("name"),
			                                t->get("in"),
			                                t->get("perm"),
			                                t->get("inv_perm"),
			                                once);
			_tools.push_back(tool);
		} else if (!t->get("type").compare("sort")) {
			Sort* tool = new Sort(t->get("name"),
			                      t->get("in"),
			                      t->get("perm"),
			                      t->get("inv_perm"),
			                      once);
			_tools.push_back(tool);
		} else if (!t->get("type").compare("assert")) {
			Assert* tool =
			    new Assert(t->get("name"), t->get("condition"), once);
			_tools.push_back(tool);
		} else if (!t->get("type").compare("if")) {
			If* tool = new If(t->get("name"), t->get("condition"), once);
			_tools.push_back(tool);
		} else if (!t->get("type").compare("while")) {
			While* tool = new While(t->get("name"), t->get("condition"), once);
			_tools.push_back(tool);
		} else if (!t->get("type").compare("endif") ||
		           !t->get("type").compare("end")) {
			End* tool = new End(t->get("name"), once);
			_tools.push_back(tool);
		} else if (!t->get("type").compare("mpi-sync")) {
#ifdef HAVE_MPI
			std::vector<std::string> fields =
			    split(replaceAllCopy(t->get("fields"), " ", ""), ',');
			std::vector<unsigned int> procs;
			if (t->get("processes").compare("")) {
				unsigned int val;
				std::vector<std::string> exprs =
				    split_formulae(t->get("processes"));
				for (auto expr : exprs) {
					_vars.solve("unsigned int", expr, (void*)(&val));
					procs.push_back(val);
				}
			}
			MPISync* tool = new MPISync(
			    t->get("name"), t->get("mask"), fields, procs, once);
			_tools.push_back(tool);
#else
			LOG(L_ERROR, "AQUAgpusph has been compiled without MPI support\n");
			LOG0(L_DEBUG, "Tools of type 'mpi-sync' cannot be considered");
#endif
		} else if (!t->get("type").compare("installable")) {
			void* handle = dlopen(t->get("path").c_str(), RTLD_LAZY);
			if (!handle) {
				std::ostringstream msg;
				msg << "Installable tool \"" << t->get("name")
				    << "\" failed loading \"" << t->get("path") << "\" library."
				    << std::endl;
				LOG(L_ERROR, msg.str());
				LOG0(L_DEBUG, dlerror());
				throw std::runtime_error("failure loading library");
			}

			Tool* (*maker)(const std::string, bool);
			maker = (Tool * (*)(const std::string, bool))
			    dlsym(handle, "create_object");
			if (!maker) {
				std::ostringstream msg;
				msg << "Installable tool \"" << t->get("name")
				    << "\" failed loading \"create_tool\" symbol." << std::endl;
				LOG(L_ERROR, msg.str());
				LOG0(L_DEBUG, dlerror());
				throw std::runtime_error("failure loading library");
			}
			Tool* tool = (Tool*)maker(t->get("name"), once);
			_tools.push_back(tool);
		} else if (!t->get("type").compare("dummy")) {
			Tool* tool = new Tool(t->get("name"), once);
			_tools.push_back(tool);
		}
		// Reports
		else if (!t->get("type").compare("report_screen")) {
			bool bold = false;
			if (!t->get("bold").compare("true") ||
			    !t->get("bold").compare("True")) {
				bold = true;
			}
			Reports::Screen* tool = new Reports::Screen(
			    t->get("name"), t->get("fields"), t->get("color"), bold);
			_tools.push_back(tool);
		} else if (!t->get("type").compare("report_file")) {
			Reports::TabFile* tool = new Reports::TabFile(
			    t->get("name"), t->get("fields"), t->get("path"));
			_tools.push_back(tool);
		} else if (!t->get("type").compare("report_particles")) {
			// Get the first particle associated to this set
			unsigned int set_id = std::stoi(t->get("set"));
			if (set_id >= _sim_data.sets.size()) {
				std::ostringstream msg;
				msg << "Report \"" << t->get("name")
				    << "\" requested the particles set " << set_id
				    << " but just " << _sim_data.sets.size() << " can be found."
				    << std::endl;
				LOG(L_ERROR, msg.str());
				throw std::runtime_error("particles set out of bounds");
			}
			unsigned int first = 0;
			for (j = 0; j < set_id; j++) {
				first += _sim_data.sets.at(j)->n();
			}

			// And the ipf and fps
			unsigned int ipf = std::stoi(t->get("ipf"));
			float fps = std::stof(t->get("fps"));

			Reports::SetTabFile* tool =
			    new Reports::SetTabFile(t->get("name"),
			                            t->get("fields"),
			                            first,
			                            _sim_data.sets.at(set_id)->n(),
			                            t->get("path"),
			                            ipf,
			                            fps);
			_tools.push_back(tool);
		} else if (!t->get("type").compare("report_performance")) {
			bool bold = false;
			if (!t->get("bold").compare("true") ||
			    !t->get("bold").compare("True")) {
				bold = true;
			}
			Reports::Performance* tool = new Reports::Performance(
			    t->get("name"), t->get("color"), bold, t->get("path"));
			_tools.push_back(tool);
		}
		// Error
		else {
			std::ostringstream msg;
			msg << "Unrecognized tool type \"" << t->get("type")
			    << "\" when parsing the tool \"" << t->get("name") << "\"."
			    << std::endl;
			LOG(L_ERROR, msg.str());
			throw std::runtime_error("Invalid tool type");
		}
	}

	// Register the reporters
	for (auto r : _sim_data.reports) {
		if (!r->get("type").compare("screen")) {
			bool bold = false;
			if (!r->get("bold").compare("true") ||
			    !r->get("bold").compare("True")) {
				bold = true;
			}
			Reports::Screen* tool = new Reports::Screen(
			    r->get("name"), r->get("fields"), r->get("color"), bold);
			_tools.push_back(tool);
		} else if (!r->get("type").compare("file")) {
			Reports::TabFile* tool = new Reports::TabFile(
			    r->get("name"), r->get("fields"), r->get("path"));
			_tools.push_back(tool);
		} else if (!r->get("type").compare("particles")) {
			// Get the first particle associated to this set
			unsigned int set_id = std::stoi(r->get("set"));
			if (set_id >= _sim_data.sets.size()) {
				std::ostringstream msg;
				msg << "Report \"" << r->get("name")
				    << "\" requested the particles set " << set_id
				    << " but just " << _sim_data.sets.size() << " can be found."
				    << std::endl;
				LOG(L_ERROR, msg.str());
				throw std::runtime_error("particles set out of bounds");
			}
			unsigned int first = 0;
			for (j = 0; j < set_id; j++) {
				first += _sim_data.sets.at(j)->n();
			}

			// And the ipf and fps
			unsigned int ipf = std::stoi(r->get("ipf"));
			float fps = std::stof(r->get("fps"));

			Reports::SetTabFile* tool =
			    new Reports::SetTabFile(r->get("name"),
			                            r->get("fields"),
			                            first,
			                            _sim_data.sets.at(set_id)->n(),
			                            r->get("path"),
			                            ipf,
			                            fps);
			_tools.push_back(tool);
		} else if (!r->get("type").compare("performance")) {
			bool bold = false;
			if (!r->get("bold").compare("true") ||
			    !r->get("bold").compare("True")) {
				bold = true;
			}
			Reports::Performance* tool = new Reports::Performance(
			    r->get("name"), r->get("color"), bold, r->get("path"));
			_tools.push_back(tool);
		} else {
			std::ostringstream msg;
			msg << "Unrecognized report type \"" << r->get("type")
			    << "\" when parsing the report \"" << r->get("name") << "\"."
			    << std::endl;
			LOG(L_ERROR, msg.str());
			throw std::runtime_error("Invalid report type");
		}
	}

	std::ostringstream msg;
	msg << "Allocated memory = " << _vars.allocatedMemory() << " bytes"
	    << std::endl;
	LOG(L_INFO, msg.str());

	setup();

	signal(SIGINT, sigint_handler);
	signal(SIGTERM, sigint_handler);
}

CalcServer::~CalcServer()
{
	unsigned int i;
	delete[] _current_tool_name;

	if (_command_queue)
		clReleaseCommandQueue(_command_queue);
	if (_command_queue_parallel)
		clReleaseCommandQueue(_command_queue_parallel);
	if (_context)
		clReleaseContext(_context);
	_context = NULL;

	if (_platforms)
		delete[] _platforms;
	_platforms = NULL;
	if (_devices)
		delete[] _devices;
	_devices = NULL;

	for (auto tool : _tools) {
		delete tool;
	}
	_tools.clear();

	for (auto& unsorter : unsorters) {
		delete unsorter.second;
	}
}

void
CalcServer::raiseSIGINT()
{
	if (sigint_received)
		return;
	sigint_handler(SIGINT);
}

void
CalcServer::update(InputOutput::TimeManager& t_manager)
{
	unsigned int i;
	while (!t_manager.mustPrintOutput() && !t_manager.mustStop()) {
#ifdef HAVE_MPI
		try {
			MPI::COMM_WORLD.Barrier();
		} catch (MPI::Exception e) {
			LOG(L_INFO, "MPI error while syncing at the beggining\n");
			std::ostringstream msg;
			msg << e.Get_error_code() << ": " << e.Get_error_string()
			    << std::endl;
			LOG0(L_DEBUG, msg.str());
			MPI::COMM_WORLD.Abort(-1);
			throw;
		}
#endif
		InputOutput::Logger::singleton()->initFrame();

		// Execute the tools
		ProfilingInfo::newStep();
		Tool* tool = _tools.front();
		while (tool) {
			try {
				strncpy(_current_tool_name, tool->name().c_str(), 255);
				tool->execute();
				tool = tool->next_tool();
			} catch (std::runtime_error& e) {
				sleep(__ERROR_SHOW_TIME__);
				throw;
			}
		}
		strcpy(_current_tool_name, "__post execution__");

		// Key events
		if (sigint_received) {
			LOG(L_WARNING, "Interrumption request (SIGINT/SIGTERM) detected\n");
			sleep(__ERROR_SHOW_TIME__);
			throw user_interruption("Simulation interrupted by the user");
		}

		InputOutput::Logger::singleton()->endFrame();
	}
}

cl_event
CalcServer::getUnsortedMem(const std::string var_name,
                           size_t offset,
                           size_t cb,
                           void* ptr)
{
	cl_int err_code;

	// Generate the unsorted if it does not exist yet
	UnSort* unsorter = NULL;
	if (unsorters.find(var_name) == unsorters.end()) {
		unsorter = new UnSort("__unsort " + var_name, var_name.c_str());
		try {
			unsorter->setup();
		} catch (std::runtime_error& e) {
			delete unsorter;
			return NULL;
		}
		unsorters.insert(std::make_pair(var_name, unsorter));
	}
	// Get the unsorter
	unsorter = unsorters[var_name];
	try {
		unsorter->execute();
	} catch (std::runtime_error& e) {
		return NULL;
	}
	cl_mem mem = unsorter->output();
	cl_event event, event_wait = unsorter->input()->getReadingEvents().back();
	err_code = clEnqueueReadBuffer(command_queue(),
	                               mem,
	                               CL_FALSE,
	                               offset,
	                               cb,
	                               ptr,
	                               1,
	                               &event_wait,
	                               &event);
	if (err_code != CL_SUCCESS) {
		std::ostringstream msg;
		msg << "Failure receiving the variable \"" << var_name
		    << "\" from server." << std::endl;
		LOG(L_ERROR, msg.str());
		InputOutput::Logger::singleton()->printOpenCLError(err_code);
		return NULL;
	}

	return event;
}

void
CalcServer::setupOpenCL()
{
	LOG(L_INFO, "Initializating OpenCL...\n");
	queryOpenCL();
	setupPlatform();
	setupDevices();
	LOG(L_INFO, "OpenCL is ready to work!\n");
}

/** @brief Get the OpenCL major and minor version from the string returned by
 * clGetPlatformInfo() and clGetDeviceInfo()
 * @param version The version string
 * @return The major and minor version numbers. 0.0 if the string is
 * ill-formatted
 */
std::tuple<unsigned int, unsigned int>
opencl_version(std::string version)
{
	auto words = split(version);
	if ((words.size() < 2) || (words.at(0) != "OpenCL")) {
		return { 0u, 0u };
	}
	auto version_nums = split(words.at(1), '.');
	if (version_nums.size() < 2) {
		return { 0u, 0u };
	}
	return { std::stoi(version_nums.at(0)), std::stoi(version_nums.at(1)) };
}

/** @brief Get the paltform's OpenCL major and minor version
 * @param platform The platform identifier
 * @return The major and minor version numbers. 0.0 if the version cannot be
 * recovered
 */
std::tuple<unsigned int, unsigned int>
opencl_version(cl_platform_id platform)
{
	cl_int err_code;
	char* opencl_version_str;
	size_t opencl_version_str_len;
	unsigned int opencl_version_major, opencl_version_minor;
	err_code = clGetPlatformInfo(
	    platform, CL_PLATFORM_VERSION, 0, NULL, &opencl_version_str_len);
	if (err_code != CL_SUCCESS) {
		return { 0u, 0u };
	}
	opencl_version_str = new char[opencl_version_str_len];
	if (!opencl_version_str) {
		std::ostringstream msg;
		msg << "Failure allocating " << opencl_version_str_len
		    << " bytes for the CL_PLATFORM_VERSION string" << std::endl;
		LOG(L_ERROR, msg.str());
		throw std::bad_alloc();
	}
	err_code = clGetPlatformInfo(platform,
	                             CL_PLATFORM_VERSION,
	                             opencl_version_str_len,
	                             opencl_version_str,
	                             NULL);
	if (err_code != CL_SUCCESS) {
		return { 0u, 0u };
	}
	auto [version_major, version_minor] = opencl_version(opencl_version_str);
	delete[] opencl_version_str;
	return { version_major, version_minor };
}

/** @brief Get the device's OpenCL major and minor version
 * @param platform The platform identifier
 * @return The major and minor version numbers. 0.0 if the version cannot be
 * recovered
 */
std::tuple<unsigned int, unsigned int>
opencl_version(cl_device_id device)
{
	cl_int err_code;
	char* opencl_version_str;
	size_t opencl_version_str_len;
	unsigned int opencl_version_major, opencl_version_minor;
	err_code = clGetDeviceInfo(
	    device, CL_DEVICE_VERSION, 0, NULL, &opencl_version_str_len);
	if (err_code != CL_SUCCESS) {
		return { 0u, 0u };
	}
	opencl_version_str = new char[opencl_version_str_len];
	if (!opencl_version_str) {
		std::ostringstream msg;
		msg << "Failure allocating " << opencl_version_str_len
		    << " bytes for the CL_DEVICE_VERSION string" << std::endl;
		LOG(L_ERROR, msg.str());
		throw std::bad_alloc();
	}
	err_code = clGetDeviceInfo(device,
	                           CL_DEVICE_VERSION,
	                           opencl_version_str_len,
	                           opencl_version_str,
	                           NULL);
	if (err_code != CL_SUCCESS) {
		return { 0u, 0u };
	}
	auto [version_major, version_minor] = opencl_version(opencl_version_str);
	delete[] opencl_version_str;
	return { version_major, version_minor };
}

/** @brief Wrapper to clGetDeviceInfo() for strings related operations
 * @param device The device
 * @param param_name The parameter
 * @return The returned string
 */
std::string
getDeviceInfoStr(cl_device_id device, cl_device_info param_name)
{
	cl_int err_code;
	char* param;
	size_t param_len;
	err_code = clGetDeviceInfo(device, param_name, 0, NULL, &param_len);
	if (err_code != CL_SUCCESS) {
		LOG(L_ERROR, "Failure getting the device param len.\n");
		InputOutput::Logger::singleton()->printOpenCLError(err_code);
		throw std::runtime_error("OpenCL error");
	}
	param = new char[param_len + 1];
	if (!param) {
		std::ostringstream msg;
		msg << "Failure allocating " << param_len + 1
		    << " bytes for the device's param" << std::endl;
		LOG(L_ERROR, msg.str());
		throw std::bad_alloc();
	}
	param[param_len] = '\0';
	err_code = clGetDeviceInfo(device, param_name, param_len, param, NULL);
	if (err_code != CL_SUCCESS) {
		LOG(L_ERROR, "Failure getting the device param.\n");
		InputOutput::Logger::singleton()->printOpenCLError(err_code);
		throw std::runtime_error("OpenCL error");
	}
	std::string res(param);
	delete[] param;
	return res;
}

void
CalcServer::queryOpenCL()
{
	cl_int err_code;
	cl_uint i, j, num_devices = 0;
	cl_device_id* devices;
	char* name;
	size_t name_len;
	_platforms = NULL;
	// Gets the total number of platforms
	err_code = clGetPlatformIDs(0, NULL, &_num_platforms);
	if (err_code != CL_SUCCESS) {
		LOG(L_ERROR, "Failure getting the number of platforms.\n");
		InputOutput::Logger::singleton()->printOpenCLError(err_code);
		throw std::runtime_error("OpenCL error");
	}
	// Get the array of platforms
	_platforms = new cl_platform_id[_num_platforms];
	if (!_platforms) {
		std::ostringstream msg;
		msg << "Failure allocating " << _num_platforms * sizeof(cl_platform_id)
		    << " bytes for OpenCL platforms array." << std::endl;
		LOG(L_ERROR, msg.str());
		throw std::bad_alloc();
	}
	err_code = clGetPlatformIDs(_num_platforms, _platforms, NULL);
	if (err_code != CL_SUCCESS) {
		LOG(L_ERROR, "Failure getting the platforms list.\n");
		InputOutput::Logger::singleton()->printOpenCLError(err_code);
		throw std::runtime_error("OpenCL error");
	}
	for (i = 0; i < _num_platforms; i++) {
		auto [version_major, version_minor] = opencl_version(_platforms[i]);
		std::ostringstream info;
		info << "Platform " << i << "..." << std::endl;
		LOG0(L_INFO, info.str());

		err_code = clGetPlatformInfo(
		    _platforms[i], CL_PLATFORM_NAME, 0, NULL, &name_len);
		if (err_code != CL_SUCCESS) {
			LOG(L_ERROR, "Failure getting the platform name length.\n");
			InputOutput::Logger::singleton()->printOpenCLError(err_code);
			throw std::runtime_error("OpenCL error");
		}
		name = new char[name_len + 1];
		if (!name) {
			std::ostringstream msg;
			msg << "Failure allocating " << name_len + 1
			    << " bytes for the platform's name" << std::endl;
			LOG(L_ERROR, msg.str());
			throw std::bad_alloc();
		}
		name[name_len] = '\0';
		err_code = clGetPlatformInfo(
		    _platforms[i], CL_PLATFORM_NAME, name_len, name, NULL);
		if (err_code != CL_SUCCESS) {
			LOG(L_ERROR, "Failure getting the platform name.\n");
			InputOutput::Logger::singleton()->printOpenCLError(err_code);
			throw std::runtime_error("OpenCL error");
		}
		info.str("");
		info << "\tNAME: " << name << std::endl;
		LOG0(L_DEBUG, info.str());
		delete[] name;

		info.str("");
		info << "\tOPENCL: " << version_major << "." << version_minor
		     << std::endl;
		LOG0(L_DEBUG, info.str());

		// Get the number of devices
		err_code = clGetDeviceIDs(
		    _platforms[i], CL_DEVICE_TYPE_ALL, 0, NULL, &num_devices);
		if (err_code != CL_SUCCESS) {
			std::ostringstream msg;
			msg << "Failure getting the number of devices (platform " << i
			    << ")." << std::endl;
			LOG(L_ERROR, msg.str());
			InputOutput::Logger::singleton()->printOpenCLError(err_code);
			throw std::runtime_error("OpenCL error");
		}
		// Gets the devices array
		devices = new cl_device_id[num_devices];
		if (!devices) {
			std::ostringstream msg;
			msg << "Failure allocating " << num_devices * sizeof(cl_device_id)
			    << " bytes for OpenCL devices array (platform " << i << ")."
			    << std::endl;
			LOG(L_ERROR, msg.str());
			throw std::bad_alloc();
		}
		err_code = clGetDeviceIDs(
		    _platforms[i], CL_DEVICE_TYPE_ALL, num_devices, devices, NULL);
		if (err_code != CL_SUCCESS) {
			std::ostringstream msg;
			msg << "Failure getting the devices list (platform " << i << ")."
			    << std::endl;
			LOG(L_ERROR, msg.str());
			InputOutput::Logger::singleton()->printOpenCLError(err_code);
			throw std::runtime_error("OpenCL error");
		}
		// Shows device arrays
		for (j = 0; j < num_devices; j++) {
			info.str("");
			info << "\tDevice " << j << "..." << std::endl;
			LOG0(L_INFO, info.str());

			info.str("");
			info << "\t\tNAME: " << getDeviceInfoStr(devices[j], CL_DEVICE_NAME)
			     << std::endl;
			LOG0(L_DEBUG, info.str());

			info.str("");
			info << "\t\tVENDOR: "
			     << getDeviceInfoStr(devices[j], CL_DEVICE_VENDOR) << std::endl;
			LOG0(L_DEBUG, info.str());

			cl_device_type dType;
			err_code = clGetDeviceInfo(devices[j],
			                           CL_DEVICE_TYPE,
			                           sizeof(cl_device_type),
			                           &dType,
			                           NULL);
			if (err_code != CL_SUCCESS) {
				LOG(L_ERROR, "Failure getting the device type.\n");
				InputOutput::Logger::singleton()->printOpenCLError(err_code);
				throw std::runtime_error("OpenCL error");
			}
			if (dType == CL_DEVICE_TYPE_CPU)
				LOG0(L_DEBUG, "\t\tTYPE: CL_DEVICE_TYPE_CPU\n");
			else if (dType == CL_DEVICE_TYPE_GPU)
				LOG0(L_DEBUG, "\t\tTYPE: CL_DEVICE_TYPE_GPU\n");
			else if (dType == CL_DEVICE_TYPE_ACCELERATOR)
				LOG0(L_DEBUG, "\t\tTYPE: CL_DEVICE_TYPE_ACCELERATOR\n");
			else if (dType == CL_DEVICE_TYPE_DEFAULT)
				LOG0(L_DEBUG, "\t\tTYPE: CL_DEVICE_TYPE_DEFAULT\n");
#ifdef CL_DEVICE_TYPE_CUSTOM
			else if (dType == CL_DEVICE_TYPE_CUSTOM)
				LOG0(L_DEBUG, "\t\tTYPE: CL_DEVICE_TYPE_CUSTOM\n");
#endif
			else {
				info.str("");
				info << "\t\tTYPE: " << dType << std::endl;
				LOG0(L_DEBUG, info.str());
			}
			auto [version_major, version_minor] = opencl_version(devices[j]);
			info.str("");
			info << "\t\tOPENCL: " << version_major << "." << version_minor
			     << std::endl;
			LOG0(L_DEBUG, info.str());
		}
		delete[] devices;
		devices = NULL;
	}
}

void
CalcServer::setupPlatform()
{
	int rank = 0;
#ifdef HAVE_MPI
	try {
		rank = MPI::COMM_WORLD.Get_rank();
	} catch (MPI::Exception e) {
		std::ostringstream msg;
		msg << "Error getting MPI rank. " << std::endl
		    << e.Get_error_code() << ": " << e.Get_error_string() << std::endl;
		LOG(L_ERROR, msg.str());
		throw;
	}
#endif
	if (rank >= _sim_data.settings.devices.size()) {
		std::ostringstream msg;
		msg << "Process " << rank << " has not an OpenCL declared device ("
		    << _sim_data.settings.devices.size() << " devices declared)"
		    << std::endl;
		LOG(L_ERROR, msg.str());
		throw std::runtime_error("Out of bounds");
	}
	const unsigned int platform_id =
	    _sim_data.settings.devices.at(rank).platform_id;
	if (platform_id >= _num_platforms) {
		std::ostringstream msg;
		LOG(L_ERROR, "The requested OpenCL platform can't be used.\n");
		msg << "Platform " << platform_id << " has been selected, but just "
		    << _num_platforms << " are available." << std::endl;
		LOG0(L_DEBUG, msg.str());
		throw std::runtime_error("Out of bounds");
	}
	_platform = _platforms[platform_id];
	auto [version_major, version_minor] = opencl_version(_platform);
	if ((version_major < 1) || ((version_major == 1) && (version_minor < 2))) {
		std::ostringstream msg;
		LOG(L_ERROR, "The requested OpenCL platform can't be used.\n");
		msg << "Platform " << platform_id << " supports just OpenCL "
		    << version_major << "." << version_minor << std::endl;
		LOG0(L_DEBUG, msg.str());
		throw std::runtime_error("Insufficient OpenCL support");
	}
}

/** @brief Create an OpenCL command queue
 *
 * This function is just a wrapper for backguard compatibility. More
 * specifically, clCreateCommandQueue is deprecated since OpenCL 1.2
 *
 * @param context OpenCL context
 * @param device OpenCL device
 * @param errcode_ret Returning error code
 * @see
 * https://www.khronos.org/registry/OpenCL/sdk/2.0/docs/man/xhtml/deprecated.html
 * @see
 * https://www.khronos.org/registry/OpenCL/sdk/2.0/docs/man/xhtml/clCreateCommandQueueWithProperties.html
 */
cl_command_queue
create_command_queue(cl_context context,
                     cl_device_id device,
                     cl_int* errcode_ret)
{
#if (OPENCL_PLATFORM_MAJOR > 1) ||                                             \
    ((OPENCL_PLATFORM_MAJOR == 1) && (OPENCL_PLATFORM_MINOR > 1))
	const cl_queue_properties properties[4] = {
		CL_QUEUE_PROPERTIES,
		CL_QUEUE_OUT_OF_ORDER_EXEC_MODE_ENABLE | CL_QUEUE_PROFILING_ENABLE,
		0
	};
	return clCreateCommandQueueWithProperties(
	    context, device, properties, errcode_ret);
#else
	cl_command_queue_properties properties =
	    CL_QUEUE_OUT_OF_ORDER_EXEC_MODE_ENABLE | CL_QUEUE_PROFILING_ENABLE;
	return clCreateCommandQueue(context, device, properties, errcode_ret);
#endif
}

/** @brief Runtime error reporting tool
 *
 * Errors reported in this way directly depends on the implementation.
 * @param errinfo is a pointer to an error string.
 * @param private_info pointer to binary data that is returned by the OpenCL
 * implementation that can be used to log additional information helpful in
 * debugging the error.
 * @param cb Size of the binary data, #private_info.
 * @param user_data Current tool name.
 */
void CL_CALLBACK
context_error_notify(const char* errinfo,
                     const void* private_info,
                     size_t cb,
                     void* user_data)
{
	const char* current_tool_name = (const char*)user_data;

	std::ostringstream msg;
	msg << "OpenCL implementation reported a runtime error";
	if (strlen(current_tool_name)) {
		msg << "(" << current_tool_name << ")";
	}
	msg << ":" << std::endl;

	LOG(L_ERROR, msg.str());
	msg.str("");
	msg << errinfo << std::endl;
	LOG0(L_DEBUG, msg.str());
}

/** @brief Sample both the system clock and the devie timer to compute the
 * offset.
 * @param event The triggering event
 * @param event_command_status CL_COMPLETE
 * @param user_data A casted pointer to the
 * Aqua::CalcServer::CalcServer::_device_timer_offset value
 */
void CL_CALLBACK
device_timer_sampler(cl_event event,
                     cl_int event_command_status,
                     void* user_data)
{
	cl_int err_code;
	cl_ulong device_timer;

	auto host_timer = CalcServer::host_timer();

	err_code = clGetEventProfilingInfo(
	    event, CL_PROFILING_COMMAND_END, sizeof(cl_ulong), &device_timer, NULL);
	if (err_code != CL_SUCCESS) {
		LOG(L_ERROR, "Failure sampling the device timer.\n");
		InputOutput::Logger::singleton()->printOpenCLError(err_code);
		return;
	}

	cl_ulong* device_timer_offset = (cl_ulong*)user_data;
	*device_timer_offset = host_timer - device_timer;
}

void
CalcServer::setupDevices()
{
	cl_int err_code;
	cl_uint i;
	_devices = NULL;

	// Get the selected device index and type
	int rank = 0;
#ifdef HAVE_MPI
	try {
		rank = MPI::COMM_WORLD.Get_rank();
	} catch (MPI::Exception e) {
		std::ostringstream msg;
		msg << "Error getting MPI rank. " << std::endl
		    << e.Get_error_code() << ": " << e.Get_error_string() << std::endl;
		LOG(L_ERROR, msg.str());
		throw;
	}
#endif
	if (rank >= _sim_data.settings.devices.size()) {
		std::ostringstream msg;
		msg << "Process " << rank << " has not an OpenCL declared device ("
		    << _sim_data.settings.devices.size() << " devices declared)"
		    << std::endl;
		LOG(L_ERROR, msg.str());
		throw std::runtime_error("Out of bounds");
	}
	const unsigned int device_id =
	    _sim_data.settings.devices.at(rank).device_id;
	const cl_device_type device_type =
	    _sim_data.settings.devices.at(rank).device_type;

	// Gets the number of valid devices
	err_code = clGetDeviceIDs(_platform, device_type, 0, NULL, &_num_devices);
	if (err_code != CL_SUCCESS) {
		LOG(L_ERROR, "Failure getting the number of devices.\n");
		InputOutput::Logger::singleton()->printOpenCLError(err_code);
		throw std::runtime_error("OpenCL error");
	}
	if (device_id >= _num_devices) {
		LOG(L_ERROR, "The selected device can't be used.\n");
		std::ostringstream msg;
		msg << "Device " << device_id << " has been selected, but just "
		    << _num_devices << " devices are available." << std::endl;
		LOG0(L_DEBUG, msg.str());
		if (device_type == CL_DEVICE_TYPE_ALL)
			LOG0(L_DEBUG, "\t\tCL_DEVICE_TYPE_ALL filter activated\n");
		else if (device_type == CL_DEVICE_TYPE_CPU)
			LOG0(L_DEBUG, "\t\tCL_DEVICE_TYPE_CPU filter activated\n");
		else if (device_type == CL_DEVICE_TYPE_GPU)
			LOG0(L_DEBUG, "\t\tCL_DEVICE_TYPE_GPU filter activated\n");
		else if (device_type == CL_DEVICE_TYPE_ACCELERATOR)
			LOG0(L_DEBUG, "\t\tCL_DEVICE_TYPE_ACCELERATOR filter activated\n");
		else if (device_type == CL_DEVICE_TYPE_DEFAULT)
			LOG0(L_DEBUG, "\t\tCL_DEVICE_TYPE_DEFAULT filter activated\n");

		throw std::runtime_error("Out of bounds");
	}
	// Gets the devices array
	_devices = new cl_device_id[_num_devices];
	if (!_devices) {
		std::ostringstream msg;
		msg << "Failure allocating " << _num_devices * sizeof(cl_device_id)
		    << " bytes for the selected devices array." << std::endl;
		LOG(L_ERROR, msg.str());
		throw std::bad_alloc();
	}
	err_code = clGetDeviceIDs(
	    _platform, device_type, _num_devices, _devices, &_num_devices);
	if (err_code != CL_SUCCESS) {
		LOG(L_ERROR, "Failure getting the devices list.\n");
		InputOutput::Logger::singleton()->printOpenCLError(err_code);
		throw std::runtime_error("OpenCL error");
	}
	// Create a devices context
	_context = clCreateContext(0,
	                           _num_devices,
	                           _devices,
	                           context_error_notify,
	                           _current_tool_name,
	                           &err_code);
	if (err_code != CL_SUCCESS) {
		LOG(L_ERROR, "Failure creating an OpenCL context.\n");
		InputOutput::Logger::singleton()->printOpenCLError(err_code);
		throw std::runtime_error("OpenCL error");
	}

	// Select the appropriate device
	_device = _devices[device_id];

	auto [version_major, version_minor] = opencl_version(_device);
	if ((version_major < 1) || ((version_major == 1) && (version_minor < 2))) {
		std::ostringstream msg;
		LOG(L_ERROR, "The requested OpenCL device can't be used.\n");
		msg << "Device " << device_id << " supports just OpenCL "
		    << version_major << "." << version_minor << std::endl;
		LOG0(L_DEBUG, msg.str());
		throw std::runtime_error("Insufficient OpenCL support");
	}

	// Create the command queues
	_command_queue = create_command_queue(_context, _device, &err_code);
	if (err_code != CL_SUCCESS) {
		LOG(L_ERROR, "Failure generating the main command queue\n");
		InputOutput::Logger::singleton()->printOpenCLError(err_code);
		throw std::runtime_error("OpenCL error");
	}
	_command_queue_parallel =
	    create_command_queue(_context, _device, &err_code);
	if (err_code != CL_SUCCESS) {
		LOG(L_ERROR, "Failure generating the parallel command queue\n");
		InputOutput::Logger::singleton()->printOpenCLError(err_code);
		throw std::runtime_error("OpenCL error");
	}

	// Let's sample the device timer
	cl_event sampler;
	auto trigger = clCreateUserEvent(_context, &err_code);
	if (err_code != CL_SUCCESS) {
		LOG(L_ERROR, "Failure generating the trigger event\n");
		InputOutput::Logger::singleton()->printOpenCLError(err_code);
		throw std::runtime_error("OpenCL error");
	}
	err_code =
	    clEnqueueMarkerWithWaitList(_command_queue, 1, &trigger, &sampler);
	if (err_code != CL_SUCCESS) {
		LOG(L_ERROR, "Failure generating the sampler event\n");
		InputOutput::Logger::singleton()->printOpenCLError(err_code);
		throw std::runtime_error("OpenCL error");
	}
	err_code = clSetEventCallback(
	    sampler, CL_COMPLETE, device_timer_sampler, &_device_timer_offset);
	if (err_code != CL_SUCCESS) {
		LOG(L_ERROR, "Failure registering the timer sampling callback\n");
		InputOutput::Logger::singleton()->printOpenCLError(err_code);
		throw std::runtime_error("OpenCL error");
	}
	err_code = clSetUserEventStatus(trigger, CL_COMPLETE);
	if (err_code != CL_SUCCESS) {
		LOG(L_ERROR, "Failure triggering the timer sampler\n");
		InputOutput::Logger::singleton()->printOpenCLError(err_code);
		throw std::runtime_error("OpenCL error");
	}
	err_code = clWaitForEvents(1, &sampler);
	if (err_code != CL_SUCCESS) {
		LOG(L_ERROR, "Failure waiting for the timer sampler\n");
		InputOutput::Logger::singleton()->printOpenCLError(err_code);
		throw std::runtime_error("OpenCL error");
	}
	for (auto e : { trigger, sampler }) {
		err_code = clReleaseEvent(e);
		if (err_code != CL_SUCCESS) {
			LOG(L_ERROR, "Failure releasing the timer sampler events\n");
			InputOutput::Logger::singleton()->printOpenCLError(err_code);
			throw std::runtime_error("OpenCL error");
		}
	}

	{
		std::ostringstream msg;
		msg << "Device timer offset = " << _device_timer_offset << " ns"
		    << std::endl;
		LOG(L_INFO, msg.str());
	}
}

void
CalcServer::setup()
{
	unsigned int i, j;
	cl_uint err_code = 0;

	// Check for the required variables that must be defined by the user
	if (!_vars.get("h")) {
		LOG(L_ERROR, "Missing kernel length \"h\".\n");
		throw std::runtime_error("Undeclared kernel length variable");
	}
	if (_vars.get("h")->type().compare("float")) {
		std::ostringstream msg;
		msg << "Kernel length \"h\" must be of type \"float\", not \""
		    << _vars.get("h")->type() << "\"" << std::endl;
		LOG(L_ERROR, msg.str());
		throw std::runtime_error("Invalid kernel length variable type");
	}
	float h = *(float*)_vars.get("h")->get();
	if (h <= 0.f) {
		std::ostringstream msg;
		msg << "Kernel length \"h\" must be a positive value (h=" << h << ")"
		    << std::endl;
		LOG(L_ERROR, msg.str());
		throw std::runtime_error("Invalid kernel length value");
	}

	std::ostringstream msg;
	msg << "Found kernel length, h = " << h << " [m]" << std::endl;
	LOG(L_INFO, msg.str());

	// Setup the scalar data per particles sets
	for (i = 0; i < _sim_data.sets.size(); i++) {
		InputOutput::ProblemSetup::sphParticlesSet* set = _sim_data.sets.at(i);
		for (j = 0; j < set->scalarNames().size(); j++) {
			std::string name = set->scalarNames().at(j);
			std::string val = set->scalarValues().at(j);
			if (!_vars.get(name)) {
				std::ostringstream msg;
				msg << "Particles set " << i
				    << " asks for the undeclared variable \"" << name << "\"."
				    << std::endl;
				LOG(L_ERROR, msg.str());
				throw std::runtime_error("Invalid variable");
			}
			if (_vars.get(name)->type().find('*') == std::string::npos) {
				std::ostringstream msg;
				msg << "Particles set " << i
				    << " can't set the value of a scalar variable, \"" << name
				    << "\"." << std::endl;
				LOG(L_ERROR, msg.str());
				throw std::runtime_error("Invalid variable type");
			}
			InputOutput::ArrayVariable* var =
			    (InputOutput::ArrayVariable*)_vars.get(name);
			size_t typesize = _vars.typeToBytes(_vars.get(name)->type());
			size_t len = _vars.get(name)->size() / typesize;
			if (len != _sim_data.sets.size()) {
				std::ostringstream msg;
				msg << "Particles set " << i
				    << " can't set the value of the array variable, \"" << name
				    << "\", because its length, " << len
				    << "differs from the number of sets, "
				    << _sim_data.sets.size() << std::endl;
				LOG(L_ERROR, msg.str());
				throw std::runtime_error("Invalid variable length");
			}
			void* data = malloc(typesize);
			try {
				_vars.solve(_vars.get(name)->type(), val, data);
			} catch (...) {
				std::ostringstream msg;
				msg << "Particles set " << i
				    << " failed evaluating variable, \"" << name << "\"."
				    << std::endl;
				LOG(L_ERROR, msg.str());
				throw;
			}
			cl_mem mem = *(cl_mem*)_vars.get(name.c_str())->get();
			err_code = clEnqueueWriteBuffer(_command_queue,
			                                mem,
			                                CL_TRUE,
			                                i * typesize,
			                                typesize,
			                                data,
			                                0,
			                                NULL,
			                                NULL);
			free(data);
			data = NULL;
			if (err_code != CL_SUCCESS) {
				std::ostringstream msg;
				msg << "Particles set " << i << " failed sending variable, \""
				    << name << "\" to the computational device." << std::endl;
				LOG(L_ERROR, msg.str());
				InputOutput::Logger::singleton()->printOpenCLError(err_code);
				throw std::runtime_error("OpenCL error");
			}
		}
	}

	// Setup the tools
	for (auto tool : _tools) {
		tool->setup();
	}
}

}
} // namespace
