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

/** \def __H5Part__
 * H5Part output is selected.
 */
#define __H5Part__ 1 << 0
/** \def __VTK__
 * VTK output is selected.
 */
#define __VTK__ 1 << 1
/** \def __TECPLOT__
 * TECPLOT output is selected.
 */
#define __TECPLOT__ 1 << 2

/** \def __DT_VARIABLE__
 * Time step will be calculated every iteration. This method guarantee the
 * Courant number conservation.
 */
#define __DT_VARIABLE__ 0
/** \def __DT_FIXCALCULATED__
 * Time step will be calculated at the start of simulation, and conserved along
 * it.
 */
#define __DT_FIXCALCULATED__ 1
/** \def __DT_FIX__
 * Time step will be provided by user.
 */
#define __DT_FIX__ 2

/** \def __BOUNDARY_ELASTIC__
 * Elastic bounce will be selected as the boundary condition.
 */
#define __BOUNDARY_ELASTIC__ 0
/** \def __BOUNDARY_FIXED__
 * Fixed particles will be selected as the boundary condition.
 */
#define __BOUNDARY_FIXED__ 1
/** \def __BOUNDARY_DELEFFE__
 * Boundary integrals will be selected as the boundary condition.
 */
#define __BOUNDARY_DELEFFE__ 2

#include <sphPrerequisites.h>

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <vector>
#include <deque>
#include <map>
#include <string>

#include <CL/cl.h>

#include <Singleton.h>

namespace Aqua{ namespace InputOutput{

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
class ProblemSetup : public Aqua::Singleton<Aqua::InputOutput::ProblemSetup>
{
public:
    /** @brief Constructor.
     *
     * In this method some initial values will be assigned, however it is
     * expected that Aqua::InputOutput::FileManager is overwritting them.
     *
     * perform() method should be called after Aqua::InputOutput::FileManager
     * was setup the simulation configuration data.
     */
    ProblemSetup();

    /** @brief Destructor.
     */
    ~ProblemSetup();

    /** @brief Compute the kernel length \f$ h \f$ and the corrected dynamic
     * viscosity \f$ \mu \f$ (due to the artificial viscosity factor
     * \f$ \alpha \f$).
     *
     * This method must be called after Aqua::InputOutput::FileManager was set
     * the simulation configuration.
     *
     * @return false if all gone right, true otherwise.
     */
    bool perform();

    /** @struct sphSettings
     * @brief General program settings.
     *
     * These setting are set between the following XML tags:
     * @code{.xml}
        <Settings>
        </Settings>
     * @endcode
     *
     * @see Aqua::InputOutput::ProblemSetup
     */
    struct sphSettings
    {
        /// Constructor.
        void init();
        /// Destructor
        void destroy();

        /** @brief Verbose level.
         *
         *   - 2 = Slowest alternative, prints relevant data every timestep.
         *   - 1 = Prints relevant data every frame (i.e. when output data
         *     files are updated).
         *   - 0 = Fastest alternative, prints minimum data.
         *
         * This field can be set with the tag `Verbose`, for instance:
         * `<Verbose level="1" />`
         */
        int verbose_level;

        /** @brief Index of the OpenCL platform to use.
         *
         * AQUAgpusph is providing the available OpenCL platforms, and the
         * devices into them.
         *
         * This field can be set with the tag `Device`, for instance:
         * `<Device platform="0" device="0" type="GPU" />`
         */
        unsigned int platform_id;

        /** @brief Index of the OpenCL device to use in the platform
         * #platform_id.
         *
         * AQUAgpusph is providing the available OpenCL platforms, and the
         * devices into them.
         *
         * This field can be set with the tag `Device`, for instance:
         * `<Device platform="0" device="0" type="GPU" />`
         *
         * @remarks The index of the device is refered to the available ones
         * compatibles with the selected type #device_type.
         */
        unsigned int device_id;

        /** @brief Type of devices that will be considered in the platform
         * #platform_id.
         *
         * This field can be set with the tag `Device`, for instance:
         * `<Device platform="0" device="0" type="GPU" />`
         *
         * @see #device_id.
         */
        cl_device_type device_type;
    }settings;

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
    struct sphVariables
    {
        /// Name of the variables
        std::deque<char*> names;
        /// Type of variables
        std::deque<char*> types;
        /// Lengths
        std::deque<char*> lengths;
        /// Values
        std::deque<char*> values;

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
         * @param save true if the variable should be saved, false otherwise.
         */
        void registerVariable(const char* name,
                              const char* type,
                              const char* length,
                              const char* value);

        /// Remove all the stored variables.
        void destroy();
    }variables;

    /** @struct sphDefinitions
     * @brief OpenCL kernels compilation definitions.
     *
     * In order to make AQUAgpusph more extensible, preprocessor definitions
     * can be set, such that they will be defined for all the compiled kernels.
     *
     * 3 types of definitions can be practised:
     *    - Named definitions, like -DHAVE_3D
     *    - Value definitions, like -DNEIGH_CELLS=9
     *    - Float evaluated definitions, like -DH=h, where h will be replaced by
     *      the variable value
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
        /// Name of the definition
        std::deque<char*> names;
        /// Value of the definition, empty for named definitions.
        std::deque<char*> values;
        /** True if the value should be evaluated as a math expression, false
         * otherwise
         */
        std::deque<bool> evaluations;

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
        void registerDefinition(const char* name,
                                const char* value,
                                const bool evaluate);

        /// Remove all the stored definitions.
        void destroy();
    }definitions;

    /** @class sphTool ProblemSetup.h ProblemSetup.h
     * @brief Tool to be executed.
     */
    class sphTool
    {
    public:
        /** Constructor
         */
        sphTool();

        /** Destructor
         */
        ~sphTool();

        /** Add new data value
         * @param name Name of the data value
         * @param value Value
         */
        void set(const char* name, const char* value);

        /** Get a data value
         * @param name Name of the data value
         * @return Value, NULL if the variable does not exist
         */
        const char* get(const char* name);

        /** Get a data value
         * @param index Index of the data value
         * @return Value, NULL if the variable does not exist
         */
        const char* get(unsigned int index);

        /** Get a data value name
         * @param index Index of the data value
         * @return Name, NULL if the variable does not exist
         */
        const char* getName(unsigned int index);

        /** Get the number of data values
         * @return Number of data values
         */
        unsigned int n() const {return _data.size();}

    private:
        /** Check if exist a data value
         * @param name Name of the data value
         * @return true if the data value has been found, false otherwise
         */
        bool has(const char* name);

        /// Registered data
        std::map<std::string, std::string> _data;
    };

    /// Array of tools
    std::deque<sphTool*> tools;

    /** @brief Array of reports.
     *
     * Reports are a some kind of special tools dedicated to generate summary
     * outputs.
     */
    std::deque<sphTool*> reports;

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

        /** @brief Log file updating criteria to apply.
         *
         * Must a combination of the following options:
         *   - #__NO_OUTPUT_MODE__ : The log file will not be never updated
         *     (default value).
         *   - #__FPS_MODE__ : Frames per second.
         *   - #__IPF_MODE__ : Iterations (time steps) per frame.
         *
         * This field can be set with the tag `Option`, for instance:
         *   - `<Option name="LogFile" type="FPS" value="120" />`
         *   - `<Option name="LogFile" type="IPF" value="10" />`
         *
         * @see Aqua::InputOutput::Log
         */
        unsigned int log_mode;

        /** @brief Log file updating rate.
         *
         * If #__FPS_MODE__ is set in #log_mode, the log file will be updated
         * this value times per second of simulations.
         *
         * This field can be set with the tag `Option`, for instance:
         * `<Option name="LogFile" type="FPS" value="120" />`
         *
         * @see #log_mode
         */
        float log_fps;

        /** @brief Log file updating rate.
         *
         * If #__IPF_MODE__ is set in #log_mode, the log file will be updated
         * every time that this value of time steps is computed.
         *
         * This field can be set with the tag `Option`, for instance:
         * `<Option name="LogFile" type="IPF" value="10" />`
         *
         * @see #log_mode
         */
        int log_ipf;

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
    }time_opts;

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
        void n(unsigned int N){_n = N;}

        /** @brief Get the number of particles
         * @return Number of particles.
         */
        unsigned int n() const {return _n;}

        /** @brief Add a scalar property for this particles set.
         *
         * This field can be set with the tag `Scalar`, for instance:
         * `<Scalar name="gamma" value="1.0" />`
         */
        void addScalar(const char* name, const char* value);

        /** @brief Get the scalar names list
         * @return scalar names.
         */
        std::deque<char*> scalarNames() const {return _snames;}

        /** @brief Get the scalar values list
         * @return scalar values.
         */
        std::deque<char*> scalarValues() const {return _svalues;}

        /** @brief Set the file path from which the particles should be read.
         *
         * This field can be set with the tag `Option`, for instance:
         * `<Load format="ASCII" file="./Fluid.dat" fields="pos,normal" />`
         *
         * @param path File path.
         * @param format File format.
         * @param fields Fields to be loaded from the file
         * @see Aqua::InputOutput::Particles
         */
        void input(const char* path, const char* format, const char* fields);

        /** @brief Get the file path from the particles would be read.
         * @return File path.
         * @see input()
         */
        const char* inputPath() const {return _in_path;}

        /** @brief Get the input file format
         * @return File format.
         * @see input()
         */
        const char* inputFormat() const {return _in_format;}

        /** @brief Get the input file fields
         * @return Fields to be loaded list.
         * @see input()
         */
        std::deque<char*> inputFields() const {return _in_fields;}

        /** @brief Set the file path where the particles should be written.
         *
         * This field can be set with the tag `Option`, for instance:
         * `<Save format="VTK" file="output" fields="pos,normal" />`
         *
         * @param path File path.
         * @param format File format.
         * @param fields Fields to be loaded from the file
         * @see Aqua::InputOutput::Particles
         */
        void output(const char* path, const char* format, const char* fields);

        /** @brief Get the output file path
         * @return File path.
         * @see output()
         */
        const char* outputPath() const {return _out_path;}

        /** @brief Get the output file format
         * @return File format.
         * @see output()
         */
        const char* outputFormat() const {return _out_format;}

        /** @brief Get the output file fields
         * @return Fields to be written list.
         * @see output()
         */
        std::deque<char*> outputFields() const {return _out_fields;}
    private:
        /// Number of particles
        unsigned int _n;

        /// Scalars names
        std::deque<char*> _snames;
        /// Scalar values
        std::deque<char*> _svalues;

        /// Particles data file to load
        char *_in_path;

        /// Format of the particles data file to load
        char *_in_format;

        /// Fields to load from the file
        std::deque<char*> _in_fields;

        /// Fluid particles data file to write
        char *_out_path;

        /// Fluid particles data file to write format
        char *_out_format;

        /// Fields to write in the file
        std::deque<char*> _out_fields;
    };

    /// Array of particles sets
    std::deque<sphParticlesSet*> sets;
};

}}  // namespace

#endif // PROBLEMSETUP_H_INCLUDED
