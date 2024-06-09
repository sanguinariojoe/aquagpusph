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
 * @brief Simulation configuration data structures.
 * (See Aqua::InputOutput::ProblemSetup for details)
 */

#ifndef PROBLEMSETUP_H_INCLUDED
#define PROBLEMSETUP_H_INCLUDED

/** \def __NO_OUTPUT_MODE__
 * No print output is selected.
 */
#define __NO_OUTPUT_MODE__ 0
/** \def __FPS_MODE__
 * FPS (frames per second) mode is selected
 * @remarks could not be combinable with IPF mode
 */
#define __FPS_MODE__ 1 << 0
/** \def __IPF_MODE__
 * IPF (iterations per frame) mode is selected
 * @remarks could not be combinable with FPS mode
 */
#define __IPF_MODE__ 1 << 1
/** \def __TIME_MODE__
 * Time control mode is selected
 * @remarks Is equal to __FPS_MODE__
 */
#define __TIME_MODE__ 1 << 0
/** \def __ITER_MODE__
 * Iterations control mode is selected
 * @remarks Is equal to __IPF_MODE__
 */
#define __ITER_MODE__ 1 << 1
/** \def __FRAME_MODE__
 * Frames control mode is selected
 * @remarks Is equal to __IPF_MODE__
 */
#define __FRAME_MODE__ 1 << 2

#include "sphPrerequisites.hpp"

#include <math.h>
#include <string>
#include <vector>
#include <vector>
#include <map>

#if __APPLE__
#include <OpenCL/cl.h>
#else
#include <CL/cl.h>
#endif

namespace Aqua {
namespace InputOutput {

/** @class ProblemSetup ProblemSetup.h ProblemSetup.h
 * @brief Simulation configuration data.
 *
 * All the XML input configuration files should have the following form:
 * @code{.xml}
    <?xml version="1.0" ?>
    <sphInput>
        ...
    </sphInput>
 * @endcode
 * Where `"..."` is the configuration data.
 *
 * See the user manual chapter 4.
 *
 * @see Aqua::InputOutput::State
 */
class ProblemSetup
{
  public:
	/** @brief Constructor.
	 */
	ProblemSetup();

	/// Copy constructor
	ProblemSetup(const ProblemSetup& p);

	/** @brief Destructor.
	 */
	~ProblemSetup();

	/** @brief Set the number of dimensions
	 * @param n The number of dimensions, either 2 or 3
	 */
	inline void dims(unsigned int n) { _dims = n; }

	/** @brief Get the the number of dimensions
	 * @return The number of dimensions, either 2 or 3
	 */
	inline unsigned int dims() const { return _dims; }

	/** @brief Set the maximum number of main OpenCL queues to be spawned
	 * @param n The number of command queues
	 * @note The parallel MPI command queue used for communications is not
	 * took into account here
	 */
	inline void nCmdQueues(unsigned int n) { _n_cmd_queues = n; }

	/** @brief Get the maximum number of main OpenCL queues to be spawned
	 * @return The number of command queues
	 * @note The parallel MPI command queue used for communications is not
	 * took into account here
	 */
	inline unsigned int nCmdQueues() const { return _n_cmd_queues; }

	/** @brief General program settings.
	 *
	 * These setting are set between the following XML tags:
	 * @code{.xml}
	    <Settings>
	    </Settings>
	 * @endcode
	 *
	 * @see Aqua::InputOutput::ProblemSetup
	 */
	class sphSettings
	{
	  public:
		/// Constructor.
		sphSettings();
		/// Destructor
		~sphSettings(){};

		/** @brief Save the output in case of failure
		 *
		 * If true, the particles sets and simulation state will be saved in
		 * case of failure/fatal error. Otherwise the simulation will just stop.
		 *
		 * You can disable save on fail with the following tag (save on fail is
		 * enabled by default):
		 * `<SaveOnFail value="false" />`
		 */
		bool save_on_fail;

		/// @brief Possible debugging modes
		typedef enum _debug {
			/// Do not debug
			NO_DEBUG = 0,
			/// Print information about the events and tools execution
			DEBUG_EVENTS = 1<<0,
			/// Print information about kernel arguments setting
			DEBUG_ARGS = 1<<1,
			/// Wait for a tool dependencies to finish after executing it
			DEBUG_DEPS_SYNC = 1<<2,
			/// Wait for a tool to finish right after executing it
			DEBUG_SYNC = 1<<3
		} debug_opts;

		/** @brief Debugging mode
		 *
		 * Several options can be activated (see ::debug_opts).
		 * Some options are just printing additional information, while some
		 * others are forcing tools to be executed one by one, while
		 * they are logged.
		 * This obviously has a quite negative impact on the performance, but
		 * can be of use to find the tool which is causing issues.
		 *
		 * You can enable the tools debugging with the following tag:
		 * `<DebugTools events="true" args="true" sync="true" deps="true" />`
		 *
		 * Please, note that when sync attribute is true, deps attribute
		 * becomes useless
		 */
		int debug_tools;

		/** @brief AQUAgpusph root path.
		 *
		 * Usually this option is automatically set by the basic module, using
		 * the tag RootPath.
		 * This path is added to the OpenCL include paths.
		 */
		std::string base_path;

		/** @brief General program settings.
		*
		* These setting are set between the following XML tags:
		* @code{.xml}
		    <Settings>
		    </Settings>
		* @endcode
		*
		* @see Aqua::InputOutput::ProblemSetup
		*/
		class device
		{
		  public:
			/// @brief Possible patch enabling states
			typedef enum _patch {
				/// Let AQUAgpusph decide
				AUTO = 0,
				/// Enabled
				ENABLED = 1,
				/// Disabled
				DISABLED = 2
			} patch_state;

			/** @brief Constructor
			 *
			 * @param platform_index Index of the OpenCL platform to use
			 * @param device_index Index of the OpenCL device to use from the
			 * platform
			 * @param t Type of OpenCL device
			 * @param bits The address bits of the device. 0 to consider
			 * CL_DEVICE_ADDRESS_BITS
			 * @see https://registry.khronos.org/OpenCL/sdk/3.0/docs/man/html/clGetDeviceInfo.html
			 */
			device(const unsigned int platform_index,
			       const unsigned int device_index,
			       const cl_device_type t = CL_DEVICE_TYPE_ALL,
			       const unsigned int bits = 32)
			  : platform_id(platform_index)
			  , device_id(device_index)
			  , device_type(t)
			  , addr_bits(bits)
			  , patches({{"nvidia_#4665567", patch_state::AUTO}})
			{};

			/// Destructor
			~device(){};

			/** @brief Index of the OpenCL platform to use.
			 *
			 * AQUAgpusph is providing the available OpenCL platforms, and the
			 * devices into them.
			 *
			 * This field can be set with the tag `Device`, for instance:
			 * `<Device platform="0" device="0" type="GPU" addr_bits="0" />`
			 */
			unsigned int platform_id;

			/** @brief Index of the OpenCL device to use in the platform
			 * #platform_id.
			 *
			 * AQUAgpusph is providing the available OpenCL platforms, and the
			 * devices into them.
			 *
			 * This field can be set with the tag `Device`, for instance:
			 * `<Device platform="0" device="0" type="GPU" addr_bits="0" />`
			 *
			 * @remarks The index of the device is refered to the available ones
			 * compatibles with the selected type #device_type.
			 */
			unsigned int device_id;

			/** @brief Type of devices that will be considered in the platform
			 * #platform_id.
			 *
			 * This field can be set with the tag `Device`, for instance:
			 * `<Device platform="0" device="0" type="GPU" addr_bits="0" />`
			 *
			 * @see #device_id.
			 */
			cl_device_type device_type;

			/** @brief Address bits of the device.
			 *
			 * If not provided, 32 will be considered, i.e. the
			 * CL_DEVICE_ADDRESS_BITS will be queried
			 *
			 * This field can be set with the tag `Device`, for instance:
			 * `<Device platform="0" device="0" type="GPU" addr_bits="0" />`
			 * @see https://registry.khronos.org/OpenCL/sdk/3.0/docs/man/html/clGetDeviceInfo.html
			 */
			unsigned int addr_bits;

			/** @brief List of patchs to be forcibly enabled/disabled
			 *
			 * The several patches currently exist:
			 *
			 *  - nvidia_#4665567 : NVIDIA CUDA platform bug #4665567.
			 *                      clEnqueueMarkerWithWaitList() does not
			 *                      properly works
			 *
			 * If a patch is not forcibly enabled/disabled, AQUAgpusph will
			 * automatically decide whether it is required or not.
			 *
			 * This field can be set with the tag `Device`, using the
			 * attributes enable_patch and disable_patch. The patchs are
			 * separated by commas. If a patch is found simultaneously on
			 * enable_patch and disable_patch attributes, it will be disabled.
			 * For instance:
			 * `<Device platform="0" device="0" type="GPU" enable_patch="nvidia_#4665567" />`
			 * @see https://registry.khronos.org/OpenCL/sdk/3.0/docs/man/html/clGetDeviceInfo.html
			 */
			std::map<std::string, int> patches;

			/** @brief Is the patch enabled?
			 * @param patch Patch name
			 * @return true if the patch is enabled, false if it is either
			 * disabled or automatically decided by AQUAgpusph
			 */
			inline bool isPatchEnabled(const std::string name) const
			{
				return patches.at(name) == patch_state::ENABLED;
			}

			/** @brief Is the patch disabled?
			 * @param patch Patch name
			 * @return true if the patch is disabled, false if it is either
			 * enabled or automatically decided by AQUAgpusph
			 */
			inline bool isPatchDisabled(const std::string name) const
			{
				return patches.at(name) == patch_state::DISABLED;
			}
		};

		/** @brief The list of devices to be considered
		 *
		 * Devices can be added with the tag `Device`, for instance:
		 * `<Device platform="0" device="0" type="GPU" />`
		 *
		 * In case MPI is enabled, devices are assigned to the processes in the
		 * same order they are declared in the XML file. Thus, the XML file and
		 * the MPI machine file must hold a strict relation
		 */
		std::vector<ProblemSetup::sphSettings::device> devices;

	};

	/// Stored settings
	sphSettings settings;

	/** @struct sphVariables
	 * @brief Simulation variables registered.
	 *
	 * In order to make AQUAgpusph more extensible, variables can be registered
	 * in runtime, specifing them in the input XML files.
	 *
	 * The variables can be either scalars or arrays:
	 *    - int, unsigned int, float.
	 *    - ivec2, ivec3, ivec4.
	 *    - uivec2, uivec3, uivec4.
	 *    - vec2, vec3, vec4.
	 *    - ivec, uivec, vec.
	 *    - int*, unsigned int*, float*.
	 *    - ivec2*, ivec3*, ivec4*.
	 *    - uivec2*, uivec3*, uivec4*.
	 *    - vec2*, vec3*, vec4*.
	 *    - ivec*, uivec*, vec*.
	 *
	 * These setting are set between the following XML tags:
	 * @code{.xml}
	    <Variables>
	    </Variables>
	 * @endcode
	 *
	 * @see Aqua::InputOutput::Variables
	 * @see Aqua::InputOutput::ProblemSetup
	 */
	class sphVariables
	{
	  public:
		/// Constructor.
		sphVariables(){};
		/// Destructor
		~sphVariables(){};

		/// Name of the variables
		std::vector<std::string> names;
		/// Type of variables
		std::vector<std::string> types;
		/// Lengths
		std::vector<std::string> lengths;
		/// Values
		std::vector<std::string> values;

		/** @brief Add a new variable.
		 *
		 * It can be repeated, such that the variable will be overwritten later
		 * by the last instance found.
		 *
		 * @param name Name of the variable.
		 * @param type Type of the varaible.
		 * @param length Array length, 1 for scalars, 0 for arrays that will
		 * not be allocated at the start (for instance the heads of chains,
		 * which requires the number of cells).
		 * @param value Variable value, NULL for arrays. It is optional for
		 * scalar variables.
		 */
		void registerVariable(std::string name,
		                      std::string type,
		                      std::string length,
		                      std::string value);
	};

	/// Variables storage
	sphVariables variables;

	/** @struct sphDefinitions
	 * @brief OpenCL kernels compilation definitions.
	 *
	 * In order to make AQUAgpusph more extensible, preprocessor definitions
	 * can be set, such that they will be defined for all the compiled kernels.
	 *
	 * 3 types of definitions can be practised:
	 *    - Named definitions, like -DHAVE_3D
	 *    - Value definitions, like -DNEIGH_CELLS=9
	 *    - Float evaluated definitions, like -DH=h, where h will be replaced
	 *      by the variable value
	 *
	 * The definitions are set between the following XML tags:
	 * @code{.xml}
	    <Definitions>
	    </Definitions>
	 * @endcode
	 *
	 * @see Aqua::InputOutput::ProblemSetup
	 * @warning The definitions are made right after the variables setup, and
	 * they are set just one time, so you cannot use them to pass changing
	 * variables to the OpenCL kernels.
	 */
	struct sphDefinitions
	{
		/// Constructor.
		sphDefinitions(){};
		/// Destructor
		~sphDefinitions(){};

		/// Name of the definition
		std::vector<std::string> names;
		/// Value of the definition, empty for named definitions.
		std::vector<std::string> values;
		/** True if the value should be evaluated as a math expression, false
		 * otherwise
		 */
		std::vector<bool> evaluations;

		/** @brief Add a new definition.
		 *
		 * It can be repeated, such that the definition will be overwritten
		 * later by the last instance found.
		 *
		 * @param name Name of the definition.
		 * @param value Value of the definition, empty for named definitions.
		 * @param evaluate True if the value should be evaluated as a math
		 * expression, false otherwise.
		 */
		void define(const std::string name,
		            const std::string value,
		            const bool evaluate);

		/** @brief Reports if a the required name has been already defined.
		 * @param name Name of the definition.
		 * @return true if it is already defined, false otherwise.
		 */
		bool isDefined(std::string name);

		/** @brief Undefine a registered definition
		 * @param name Name of the definition.
		 */
		void undefine(std::string name);
	};

	/// Definitions storage
	sphDefinitions definitions;

	/** @class sphTool ProblemSetup.h ProblemSetup.h
	 * @brief Tool to be executed.
	 */
	class sphTool
	{
	  public:
		/** Constructor
		 */
		sphTool(){};

		/** Destructor
		 */
		~sphTool(){};

		/** Add new data value
		 * @param name Name of the data value
		 * @param value Value
		 */
		void set(const std::string name, const std::string value);

		/** Get a data value
		 * @param name Name of the data value
		 * @return Value, NULL if the variable does not exist
		 */
		const std::string get(const std::string name);

		/** Get a data value
		 * @param index Index of the data value
		 * @return Value, NULL if the variable does not exist
		 */
		const std::string get(unsigned int index);

		/** Get a data value name
		 * @param index Index of the data value
		 * @return Name, NULL if the variable does not exist
		 */
		const std::string getName(unsigned int index);

		/** Get the number of data values
		 * @return Number of data values
		 */
		unsigned int n() const { return _data.size(); }

	  private:
		/** Check if exist a data value
		 * @param name Name of the data value
		 * @return true if the data value has been found, false otherwise
		 */
		bool has(const std::string name);

		/// Registered data
		std::map<std::string, std::string> _data;
	};

	/// Array of tools
	std::vector<sphTool*> tools;

	/** @brief Helper function to get the number of already defined instances of
	 * the same tool
	 *
	 * Since the same tool can be added several times (using wildcards), at the
	 * time of removing it the user should be careful before of checking that it
	 * is the last remaining instance.
	 *
	 * @param tool Tool instance to look for
	 * @return The number of instances found
	 */
	unsigned int toolInstances(ProblemSetup::sphTool* tool);

	/** @brief Array of reports.
	 *
	 * Reports are a some kind of special tools dedicated to generate summary
	 * outputs.
	 */
	std::vector<sphTool*> reports;

	/** @struct sphTimingParameters
	 * @brief Simulation time flow options.
	 *
	 * This options are used to control the simulation time step, the output
	 * files frequency or the simulation finish criteria.
	 *
	 * These setting are set between the following XML tags:
	 * @code{.xml}
	    <Timing>
	    </Timing>
	 * @endcode
	 *
	 * @see Aqua::InputOutput::ProblemSetup
	 * @see Aqua::InputOutput::TimeManager
	 */
	struct sphTimingParameters
	{
		/** @brief Simulation finish criteria to apply.
		 *
		 * Must a combination of the following options:
		 *   - #__TIME_MODE__ : Maximum simulated time.
		 *   - #__ITER_MODE__ : Maximum number os time steps.
		 *   - #__FRAME_MODE__ : Maximum number of output frames printed.
		 *
		 * This field can be set with the tag `Option`, for instance:
		 *   - `<Option name="End" type="Time" value="1.0" />`
		 *   - `<Option name="End" type="Steps" value="10000" />`
		 *   - `<Option name="End" type="Frames" value="10" />`
		 *
		 * @warning If no simulation stop criteria is selected, the simulation
		 * will not stop until an error is reached, or it is manually stopped.
		 */
		unsigned int sim_end_mode;

		/** @brief Simulation finish time instant.
		 *
		 * If #__TIME_MODE__ is set in #sim_end_mode, the simulation will be
		 * stopped when this simulation time instant is reached.
		 *
		 * This field can be set with the tag `Option`, for instance:
		 * `<Option name="End" type="Time" value="1.0" />`
		 *
		 * @see #sim_end_mode
		 */
		float sim_end_time;

		/** @brief Total number of time steps to compute.
		 *
		 * If #__ITER_MODE__ is set in #sim_end_mode, the simulation will be
		 * stopped when this simulation time steps are computed.
		 *
		 * This field can be set with the tag `Option`, for instance:
		 * `<Option name="End" type="Steps" value="10000" />`
		 *
		 * @see #sim_end_mode
		 */
		int sim_end_step;

		/** @brief Total number of output frames to write.
		 *
		 * If #__FRAME_MODE__ is set in #sim_end_mode, the simulation will be
		 * stopped when this output frames are written.
		 *
		 * This field can be set with the tag `Option`, for instance:
		 * `<Option name="End" type="Frames" value="10" />`
		 *
		 * @see #sim_end_mode
		 */
		int sim_end_frame;

		/** @brief Particles output updating criteria to apply.
		 *
		 * The particles output may be hard disk heavily demanding, hardly
		 * affecting the general program performance as well, therefore it is
		 * strongly recommended to try to write it as less often as possible
		 *
		 * Must a combination of the following options:
		 *   - #__NO_OUTPUT_MODE__ : The log file will not be never updated
		 *     (default value).
		 *   - #__FPS_MODE__ : Frames per second.
		 *   - #__IPF_MODE__ : Iterations (time steps) per frame.
		 *
		 * This field can be set with the tag `Option`, for instance:
		 *   - `<Option name="Output" type="FPS" value="120" />`
		 *   - `<Option name="Output" type="IPF" value="10" />`
		 *
		 * @see Aqua::InputOutput::Particles
		 */
		unsigned int output_mode;

		/** @brief Particles output updating rate.
		 *
		 * If #__FPS_MODE__ is set in #output_mode, the particles output will be
		 * updated this value times per second of simulations.
		 *
		 * This field can be set with the tag `Option`, for instance:
		 * `<Option name="Output" type="FPS" value="120" />`
		 *
		 * @see #output_mode
		 */
		float output_fps;

		/** @brief Particles output updating rate.
		 *
		 * If #__IPF_MODE__ is set in #output_mode, the particles output will be
		 * updated every time that this value of time steps is computed.
		 *
		 * This field can be set with the tag `Option`, for instance:
		 * `<Option name="Output" type="IPF" value="10" />`
		 *
		 * @see #output_mode
		 */
		int output_ipf;

		/** @brief Time step \f$ \Delta t \f$ computation method.
		 *
		 * The time step \f$ \Delta t \f$ can be computed in 3 ways (they
		 * cannot be combined):
		 *   - #__DT_VARIABLE__ : Time step \f$ \Delta t \f$ is recomputed each
		 *     iteration.
		 *   - #__DT_FIXCALCULATED__ : Time step \f$ \Delta t \f$ is computed at
		 *     the start of the simulation, being constant after that.
		 *   - #__DT_FIX__ : Time step \f$ \Delta t \f$ is manually specified by
		 *     the user.
		 *
		 * This field can be set with the tag `Option`, for instance:
		 *   - `<Option name="TimeStep" value="Variable" />`
		 *   - `<Option name="TimeStep" value="Fixed" />`
		 *   - `<Option name="TimeStep" value="0.1" />`
		 *
		 * @see Aqua::InputOutput::TimeManager
		 * @see Aqua::CalcServer::TimeStep
		 */
		unsigned int dt_mode;
	} time_opts;

	/** @class sphParticlesSet ProblemSetup.h ProblemSetup.h
	 * @brief Particles set data.
	 *
	 * In AQUAgpusph the particles can be deivided in sets.
	 * Each set may have different physic properties, or it can be used just to
	 * separate the particles by type (i.e.
	 * Fluid particles/Boundaries/Sensors).
	 *
	 * Each particles set should be loaded from a different files, as well as
	 * saved in different files.
	 *
	 * These setting are set between the following XML tags:
	 * @code{.xml}
	    <ParticlesSet n="100000">
	    </ParticlesSet>
	 * @endcode
	 * Where @paramname{n} is the number of particles for this set of
	 * particles.
	 *
	 * @see Aqua::InputOutput::Particles
	 */
	class sphParticlesSet
	{
	  public:
		/// Constructor
		sphParticlesSet();

		/// Destructor
		~sphParticlesSet();

		/** @brief Set the number of particles
		 * @param N Number of particles.
		 */
		void n(size_t N) { _n = N; }

		/** @brief Get the number of particles
		 * @return Number of particles.
		 */
		size_t n() const { return _n; }

		/** @brief Add a scalar property for this particles set.
		 *
		 * This field can be set with the tag `Scalar`, for instance:
		 * `<Scalar name="delta" value="1.0" />`
		 */
		void addScalar(std::string name, std::string value);

		/** @brief Get the scalar names list
		 * @return scalar names.
		 */
		std::vector<std::string> scalarNames() const { return _snames; }

		/** @brief Get the scalar values list
		 * @return scalar values.
		 */
		std::vector<std::string> scalarValues() const { return _svalues; }

		/** @brief Set the file path from which the particles should be read.
		 *
		 * This field can be set with the tag `Load`, for instance:
		 * `<Load format="ASCII" file="./Fluid.dat" fields="r,normal" />`
		 *
		 * Several constants can be used in the file name. See
		 * Aqua::setStrConstants()
		 *
		 * @param path File path.
		 * @param format File format.
		 * @param fields Fields to be loaded from the file
		 * @see Aqua::InputOutput::Particles
		 */
		void input(std::string path, std::string format, std::string fields);

		/** @brief Get the file path from the particles would be read.
		 * @return File path.
		 * @see input()
		 */
		const std::string inputPath() const { return _in_path; }

		/** @brief Get the input file format
		 * @return File format.
		 * @see input()
		 */
		const std::string inputFormat() const { return _in_format; }

		/** @brief Get the input file fields
		 * @return Fields to be loaded list.
		 * @see input()
		 */
		std::vector<std::string> inputFields() const { return _in_fields; }

		/** @brief Set the file path where the particles should be written.
		 *
		 * This field can be set with the tag `Save`, for instance:
		 * `<Save format="VTK" file="output.proc{mpi_rank}.{index}.vtk"
		 * fields="r,normal" />`
		 *
		 * Several scaper strings can be used to let AQUAgpusph compute the
		 * file path, at the appropriate Aqua::InputOutput::Particles instance.
		 * See Aqua::newFilePath()
		 *
		 * @param path File path.
		 * @param format File format.
		 * @param fields Fields to be loaded from the file
		 * @see Aqua::InputOutput::Particles
		 * @note In case `"%d"`/`"{index}"` is not provided, it will be appended
		 * at the end of the file name, as well as the appropriate file
		 * extension (backward compatibility)
		 */
		void output(std::string path, std::string format, std::string fields);

		/** @brief Get the output file path
		 * @return File path.
		 * @see output()
		 */
		const std::string outputPath() const { return _out_path; }

		/** @brief Get the output file format
		 * @return File format.
		 * @see output()
		 */
		const std::string outputFormat() const { return _out_format; }

		/** @brief Get the output file fields
		 * @return Fields to be written list.
		 * @see output()
		 */
		std::vector<std::string> outputFields() const { return _out_fields; }

	  private:
		/// Number of particles
		size_t _n;

		/// Scalars names
		std::vector<std::string> _snames;
		/// Scalar values
		std::vector<std::string> _svalues;

		/// Particles data file to load
		std::string _in_path;

		/// Format of the particles data file to load
		std::string _in_format;

		/// Fields to load from the file
		std::vector<std::string> _in_fields;

		/// Fluid particles data file to write
		std::string _out_path;

		/// Fluid particles data file to write format
		std::string _out_format;

		/// Fields to write in the file
		std::vector<std::string> _out_fields;
	};

	/// Array of particles sets
	std::vector<sphParticlesSet*> sets;

  private:
	/// The number of dimensions, either 2 or 3
	unsigned int _dims;

	/** Maximum number of main OpenCL command queues (there is always a
	 * parallel one for MPI communications)
	 */
	unsigned int _n_cmd_queues;

	/// Hack to avoid local copies try to remove the class
	bool _copy;
};

}
} // namespace

#endif // PROBLEMSETUP_H_INCLUDED
