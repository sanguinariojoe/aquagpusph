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
     * @brief Simulation variable registered.
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
        /// Save flags
        std::deque<bool> saves;

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
                              const char* value,
                              const bool save);

        /// Remove all the stored variables.
        void destroy();
    }variables;

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

    /** @struct sphOpenCLKernels
     * @brief OpenCL kernels source code files that should be loaded for each
     * tool of Aqua::CalcServer::CalcServer.
     *
     * These setting are set between the following XML tags:
     * @code{.xml}
        <OpenCL>
        </OpenCL>
     * @endcode
     *
     * @note AQUAgpusph provide the file OpenCLMain.xml in the folder resources
     * that will select a set of kernels that should be enough for the almost
     * simulations.
     * @see Aqua::InputOutput::ProblemSetup
     */
    struct sphOpenCLKernels
    {
        /// Constructor.
        void init();
        /// Destructor.
        void destroy();

        /** @brief Predictor time integration stage.
         *
         * This field can be set with the tag `Predictor`, for instance:
         * `<Predictor file="OpenCL/Predictor" />`
         *
         * @note The `.cl` extension will be automatically added.
         */
        char *predictor;
        /** @brief Particles data permutation using the #radix_sort data.
         *
         * This field can be set with the tag `Permutate`, for instance:
         * `<Permutate file="OpenCL/Permutate" />`
         *
         * @note The `.cl` extension will be automatically added.
         */
        char *link_list;
        /** @brief LinkList computation.
         *
         * This field can be set with the tag `LinkList`, for instance:
         * `<LinkList file="OpenCL/LinkList" />`
         *
         * @note The `.cl` extension will be automatically added.
         */
        char *permutate;
        /** @brief Fluid particles interactions.
         *
         * This field can be set with the tag `Rates`, for instance:
         * `<Rates file="OpenCL/Rates" />`
         *
         * @note The `.cl` extension will be automatically added.
         */
        char *rates;
        /** @brief Corrector time integration stage.
         *
         * This field can be set with the tag `Corrector`, for instance:
         * `<Corrector file="OpenCL/Corrector" />`
         *
         * @note The `.cl` extension will be automatically added.
         */
        char *corrector;
        /** @brief Time step computation.
         *
         * This field can be set with the tag `TimeStep`, for instance:
         * `<TimeStep file="OpenCL/TimeStep" />`
         *
         * @note The `.cl` extension will be automatically added.
         */
        char *time_step;

        /** @brief Reduction, i.e. Prefix sums, maximum/minimum value, etc...
         *
         * This field can be set with the tag `Reduction`, for instance:
         * `<Reduction file="OpenCL/Reduction" />`
         *
         * @note The `.cl` extension will be automatically added.
         */
        char *reduction;
        /** @brief RadixSort, i.e. Particles sort by the cell index.
         *
         * This field can be set with the tag `RadixSort`, for instance:
         * `<RadixSort file="OpenCL/RadixSort" />`
         *
         * @note The `.cl` extension will be automatically added.
         */
        char *radix_sort;

        /** @brief Density interpolation (formerly density reinitialization).
         *
         * The density reinitialization imply that the density resulting from
         * the mass conservation equation is replaced by \f$ \rho(\mathbf{x}) =
           \int_{\Omega} \rho(\mathbf{y}) W(\mathbf{y} - \mathbf{x})
           \mathrm{d}\mathbf{x} \f$
         *
         * This field can be set with the tag `DensInterpolation`, for instance:
         * `<DensInterpolation file="OpenCL/DensInterpolation" />`
         *
         * @note The `.cl` extension will be automatically added.
         */
        char *dens_int;
        /** @brief 0th order correction (formerly Shepard correction).
         *
         * The resulting variation rates will be renormalized dividing by the
         * Shepard factor \f$ \gamma(\mathbf{x}) = \int_{\Omega}
           W(\mathbf{y} - \mathbf{x}) \mathrm{d}\mathbf{x} \f$, such that 0th
         * order fields will be consistently computed.
         *
         * This field can be set with the tag `Shepard`, for instance:
         * `<Shepard file="OpenCL/Shepard" />`
         *
         * @note The `.cl` extension will be automatically added.
         */
        char *shepard;
        /** @brief Simplest boundary condition consisting in a elastic bounce.
         *
         * This field can be set with the tag `ElasticBounce`, for instance:
         * `<ElasticBounce file="OpenCL/ElasticBounce" />`
         *
         * @note The `.cl` extension will be automatically added.
         */
        char *elastic_bounce;
        /** @brief Boundary integrals based boundary condition.
         *
         * This field can be set with the tag `DeLeffe`, for instance:
         * `<DeLeffe file="OpenCL/DeLeffe" />`
         *
         * @note The `.cl` extension will be automatically added.
         */
        char *de_Leffe;
        /** @brief Ghost particles based boundary condition.
         *
         * This field can be set with the tag `GhostParticles`, for instance:
         * `<GhostParticles file="OpenCL/GhostParticles" />`
         *
         * @note The `.cl` extension will be automatically added.
         */
        char *ghost;
        /** @brief Forces and moments of the fluid.
         *
         * This field can be set with the tag `Torque`, for instance:
         * `<Torque file="OpenCL/Torque" />`
         *
         * @note The `.cl` extension will be automatically added.
         */
        char *torque;
        /** @brief Energy of the fluid.
         *
         * This field can be set with the tag `Energy`, for instance:
         * `<Energy file="OpenCL/Energy" />`
         *
         * @note The `.cl` extension will be automatically added.
         */
        char *energy;
        /** @brief Fluid bound box, and minimum and maximum velocities.
         *
         * This field can be set with the tag `Bounds`, for instance:
         * `<Bounds file="OpenCL/Bounds" />`
         *
         * @note The `.cl` extension will be automatically added.
         */
        char *bounds;
        /** @brief Domain check in order to exclude the particles far away.
         *
         * Particles running far away may break the simulations due to the
         * increasing number of cells, so it is strongly recommended to set
         * a domain to destroy the problematic particles.
         *
         * This field can be set with the tag `Domain`, for instance:
         * `<Domain file="OpenCL/Domain" />`
         *
         * @note The `.cl` extension will be automatically added.
         */
        char *domain;
        /// Deprecated field.
        char *portal;
    }OpenCL_kernels;

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
        /** @brief Starting time instant.
         *
         * If a negative value is provided a stabilization period will be
         * considered.
         *
         * This field can be set with the tag `Option`, for instance:
         * `<Option name="Start" value="0.0" />`
         *
         * @note During the stabilization stage motions will not be computed,
         * and particles velocity will be vanished (infinite viscosity model).
         */
        float t0;

        /** @brief Time step at the start of the simulation.
         *
         * This option is designed for simulations continuation, for other
         * simulations this options should be equal to 0.
         *
         * This field should not be manually set.
         */
        float dt0;

        /** @brief Starting time step index.
         *
         * This option is designed for simulations continuation, for other
         * simulations this options should be equal to 0.
         *
         * This field should not be manually set.
         */
        unsigned int step0;

        /** @brief Starting output frame index.
         *
         * This option is designed for simulations continuation, for other
         * simulations this options should be equal to 0.
         *
         * This field should not be manually set.
         */
        unsigned int frame0;

        /** @brief Time step \f$ \Delta t \f$ manually set by the user.
         *
         * This field can be set with the tag `Option`, for instance to set
         * \f$ \Delta t = 0.1 \mathrm{s} \f$:
         * `<Option name="TimeStep" value="0.1"/>`
         *
         * @see #dt_mode
         */
        float dt;

        /** @brief Minimum time step.
         *
         * This option will be useful just if variable time step is set.
         * When the computed time step value below this options value is
         * detected, it will be clamped.
         *
         * This field can be set with the tag `Option`, for instance:
         * `<Option name="MinTimeStep" value="0.0"/>`
         *
         * @see #dt_mode
         */
        float dt_min;

        /** @brief Courant factor.
         *
         * This option will be useful just if variable time step is set.
         * The computed required time step will be multiplied by the Courant
         * factor.
         * Low Courant factor values implies higher computational costs,
         * however high values may lead instabilities.
         *
         * This field can be set with the tag `Option`, for instance:
         * `<Option name="Courant" value="0.25"/>`
         *
         * @see #dt_mode
         */
        float courant;

        /** @brief Velocity clamping to match with the minimum time step
         * #dt_min.
         *
         * The velocity of the particles will be conveniently clamped to ensure
         * that the computed time step will not be lower than #dt_min.
         *
         * This field can be set with the tag `Option`, for instance:
         * `<Option name="ClampVel" value="False"/>`
         *
         * @see #dt_min
         */
        bool velocity_clamp;

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

        /** @brief Energy report file updating criteria to apply.
         *
         * Must a combination of the following options:
         *   - #__NO_OUTPUT_MODE__ : The log file will not be never updated
         *     (default value).
         *   - #__FPS_MODE__ : Frames per second.
         *   - #__IPF_MODE__ : Iterations (time steps) per frame.
         *
         * This field can be set with the tag `Option`, for instance:
         *   - `<Option name="EnFile" type="FPS" value="120" />`
         *   - `<Option name="EnFile" type="IPF" value="10" />`
         *
         * @see Aqua::InputOutput::Energy
         */
        unsigned int energy_mode;

        /** @brief Energy report file updating rate.
         *
         * If #__FPS_MODE__ is set in #energy_mode, the energy file will be
         * updated this value times per second of simulations.
         *
         * This field can be set with the tag `Option`, for instance:
         * `<Option name="EnFile" type="FPS" value="120" />`
         *
         * @see #energy_mode
         */
        float energy_fps;

        /** @brief Energy report file updating rate.
         *
         * If #__IPF_MODE__ is set in #energy_mode, the energy file will be
         * updated every time that this value of time steps is computed.
         *
         * This field can be set with the tag `Option`, for instance:
         * `<Option name="EnFile" type="IPF" value="10" />`
         *
         * @see #energy_mode
         */
        int energy_ipf;

        /** @brief Fluid bounds file updating criteria to apply.
         *
         * Must a combination of the following options:
         *   - #__NO_OUTPUT_MODE__ : The log file will not be never updated
         *     (default value).
         *   - #__FPS_MODE__ : Frames per second.
         *   - #__IPF_MODE__ : Iterations (time steps) per frame.
         *
         * This field can be set with the tag `Option`, for instance:
         *   - `<Option name="BoundsFile" type="FPS" value="120" />`
         *   - `<Option name="BoundsFile" type="IPF" value="10" />`
         *
         * @see Aqua::InputOutput::Bounds
         */
        unsigned int bounds_mode;

        /** @brief Fluid bounds file updating rate.
         *
         * If #__FPS_MODE__ is set in #bounds_mode, the bounds file will be
         * updated this value times per second of simulations.
         *
         * This field can be set with the tag `Option`, for instance:
         * `<Option name="BoundsFile" type="FPS" value="120" />`
         *
         * @see #bounds_mode
         */
        float bounds_fps;

        /** @brief Fluid bounds file updating rate.
         *
         * If #__IPF_MODE__ is set in #bounds_mode, the bounds file will be
         * updated every time that this value of time steps is computed.
         *
         * This field can be set with the tag `Option`, for instance:
         * `<Option name="BoundsFile" type="IPF" value="10" />`
         *
         * @see #bounds_mode
         */
        int bounds_ipf;

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

    /** @struct sphSPHParameters
     * @brief Simulation SPH related parameters.
     *
     * This options are used to control the physics simulation parameters.
     *
     * These setting are set between the following XML tags:
     * @code{.xml}
        <SPH>
        </SPH>
     * @endcode
     *
     * @see Aqua::InputOutput::ProblemSetup
     */
    struct sphSPHParameters
    {
        /** @brief Gravity acceleration \f$ \mathbf{g} \f$, present in the
         * momentum equation.
         *
         * This field can be set with the tag `Option`, for instance:
         * `<Option name="g" x="0.0" y="0.0" z="-9.81" />`
         *
         * @warning In 2D simulations the field \f$ z \f$ does not exist.
         */
        vec g;
        /** @brief Kernel height factor \f$ \frac{h}{\Delta \mathbf{r}} \f$.
         *
         * This field can be set with the tag `Option`, for instance:
         * `<Option name="hfac" value="4.0" />`
         */
        float hfac;
        /** @brief Typical particle distance \f$ \Delta \mathbf{r} \f$.
         *
         * The particles distance may (and probably will) be non constant.
         * Since this value controls the kernel height, if several particles
         * distances can be selected, the largest one will be the best option.
         *
         * This field can be set with the tag `Option`, for instance:
         * `<Option name="deltar" x="0.001" y="0.001" z="0.001" />`
         *
         * @warning In 2D simulations the field \f$ z \f$ does not exist.
         */
        vec deltar;
        /** @brief Kernel height \f$ h \f$.
         *
         * This value cannot be manually set but it is computed as #hfac
         * \f$ \cdot \f$ #deltar.
         */
        float h;
        /** @brief Characteristic sound speed \f$ c_s \f$.
         *
         * The sound speed is controlling the fluid compressibility such that
         * higher values means lower fluid compressibility.
         * However, increasing the sound speed will result in lower time step,
         * and therefore in higher computational costs.
         * Usually 10 times the maximum expected velocity should be enough.
         *
         * This field can be set with the tag `Option`, for instance:
         * `<Option name="cs" value="15.0" />`
         */
        float cs;

        /** @brief Frequency of Link-List computation (in time steps between
         * computation events).
         *
         * The Link-List process imply a radix-sort computation which can be
         * computationally demanding, so reducing the times that it will be
         * computed may lead some performance improvements.
         * On the other hand, increasing the Link-List validity period results
         * in a cells length grow, and therefore an increment of neighs as well,
         * which can cause a performance lost.
         *
         * This field can be set with the tag `Option`, for instance:
         * `<Option name="LLSteps" value="1" />`
         *
         * @see Aqua::CalcServer::LinkList
         */
        unsigned int link_list_steps;

        /** @brief Frequency of density interpolations (in time steps between
         * computation events).
         *
         * The density interpolation may smooth the density field, however
         * its consistency has not been demonstrated, and it may result in some
         * instabilities, so it is recommended to let it disabled (0 value).
         *
         * This field can be set with the tag `Option`, for instance:
         * `<Option name="DensSteps" value="0" />`
         *
         * @see Aqua::CalcServer::DensityInterpolation
         */
        unsigned int dens_int_steps;
        /** @brief Minimum value for the density.
         *
         * If the density of a particle becomes lower than this one, it will be
         * conveniently clamped.
         * This trick may help to avoid that a punctual instability can
         * critically affect the rest of the simulation.
         *
         * This field can be set with the tag `Option`, for instance:
         * `<Option name="DensBounds" min="800.0" max="1200.0" />`
         *
         * @see #rho_max
         */
        float rho_min;
        /** @brief Maximum value for the density.
         *
         * If the density of a particle becomes greater than this one, it will
         * be conveniently clamped.
         * This trick may help to avoid that a punctual instability can
         * critically affect the rest of the simulation.
         *
         * This field can be set with the tag `Option`, for instance:
         * `<Option name="DensBounds" min="800.0" max="1200.0" />`
         *
         * @see #rho_min
         */
        float rho_max;

        /** @brief Boundary condition technique to be applied.
         *
         * It should be one of the following values:
         *    - #__BOUNDARY_ELASTIC__ : Elastic bounce.
         *    - #__BOUNDARY_FIXED__ : Fixed particles.
         *    - #__BOUNDARY_DELEFFE__ : Boundary integrals.
         *
         * This field can be set with the tag `Option`, for instance:
         *     - `<Option name="Boundary" value="ElasticBounce" />`
         *     - `<Option name="Boundary" value="FixedParticles" />`
         *     - `<Option name="Boundary" value="DeLeffe" />`
         *
         * @note Ghost particles boundary condition are managed in a different
         * way.
         * @see Aqua::CalcServer::Boundary::ElasticBounce
         * @see Aqua::CalcServer::Boundary::DeLeffe
         * @see Aqua::CalcServer::Rates
         */
        unsigned int boundary_type;
        /** @brief Slip boundary condition.
         *
         * Event though no-slip boundary conditions may be desirable, sometimes
         * the Reynolds number of the simulation is not large enough to
         * adequately solve the boundary layer, and in this case free-slip
         * boundary condition should be applied.
         *
         * It should be one of the following values:
         *     - 0 = No Slip condition.
         *     - 1 = Free Slip condition.
         *
         * This field can be set with the tag `Option`, for instance:
         *     - `<Option name="SlipCondition" value="FreeSlip" />`
         *     - `<Option name="SlipCondition" value="NoSlip" />`
         */
        unsigned int slip_condition;
        /** @brief Elastic factor for the elastic bounce boundary condition.
         *
         * If #boundary_type is #__BOUNDARY_ELASTIC__, this value will set the
         * amount of kinetic energy conserved in the interaction. A factor of 1
         * imply that the velocity of the particle will be preserved (except for
         * the direction), while a factor of 0 imply that the particle will loss
         * all its normal to the boundary velocity.
         *
         * The tangential velocity is not affected.
         *
         * This field can be set with the tag `Option`, for instance:
         * `<Option name="BoundElasticFactor" value="0.0" />`
         *
         * @see #boundary_type
         */
        float elastic_factor;
        /** @brief Interaction distance for the elastic bounce boundary
         * condition.
         *
         * If #boundary_type is #__BOUNDARY_ELASTIC__, this value will set the
         * distance where the interaction with the boundary will be computed.
         * This distance is set as a factor of the kernel height \f$ h \f$
         *
         * This field can be set with the tag `Option`, for instance:
         * `<Option name="BoundDist" value="0.1" />`
         *
         * @see #hfac
         */
        float elastic_dist;

        /** @brief 0th order correction.
         *
         * Also known as Shepard correction.
         *
         * The computed rates of change will be renormalized such that constant
         * fields will be consistently interpolated.
         * The Shepard term is computed as
         * \f$ \gamma(\mathbf{x}) = \int_{\Omega}
         *     W(\mathbf{y} - \mathbf{x})
         *     \mathrm{d}\mathbf{x} \f$.
         *
         * Must be one of the following values:
         *   - 0 = No 0th order correction.
         *   - 1 = It will be applied to the computed forces.
         *   - 2 = It will be applied to computed density rate.
         *   - 3 = It will be applied to both of them.
         *
         * This field can be set with the tag `Option`, for instance:
         *   - `<Option name="Shepard" value="None" />`
         *   - `<Option name="Shepard" value="Force" />`
         *   - `<Option name="Shepard" value="Dens" />`
         *   - `<Option name="Shepard" value="ForceDens" />`
         *
         * @remarks 0th order correction (formerly Shepard correction) is
         * mandatory when Boundary Integrals based boundary condition is
         * imposed, but since Shepard term can carry some instabilities, it
         * is strongly recommended to do not apply it if it is not required by
         * the formulation.
         */
        unsigned int has_shepard;

        /** @brief Existence of a computational domain.
         *
         * If true, domain bounds will be established such that the particles
         * out of the domain will be replaced by zero mass fixed particles.
         *
         * Doing this it can be avoided that a far particle, in a far Link-List
         * cell, cause the program fail due to the large number of required
         * cells.
         * The particle will not be deleted.
         *
         * This field can be set with the tag `Option`, for instance:
         * `<Option name="Domain" x="0.0" y="0.0" z="0.0" l="1.0" d="1.0" h="1.0" move="true" />`
         *
         * @warning In the 2D version just the fields \f$ x \f$, \f$ y \f$,
         * \f$ l \f$ and \f$ h \f$ exists.
         * @see #domain_min
         * @see #domain_max
         * @see #domain_motion
         */
        bool has_domain;
        /** @brief Computational domain lower corner.
         *
         * @see #has_domain
         */
        vec domain_min;
        /** @brief Computational domain higher corner.
         *
         * @see #has_domain
         */
        vec domain_max;
        /** @brief Computational domain motion.
         *
         * If true, the domain will follow the motions. It is useful for closed
         * solids, like tanks, to assert that the domain follows the fluid.
         *
         * @see #has_domain
         */
        bool domain_motion;
    }SPH_opts;

    /** @struct sphFluidParameters
     * @brief Fluid specie physic parameters.
     *
     * In the simulations several fluid species may coexist (multiphase
     * simulations for instance).
     * In this structure the physic parameters of each fluid specie are stored.
     *
     * These setting are set between the following XML tags:
     * @code{.xml}
        <Fluid n="100000">
        </Fluid>
     * @endcode
     * Where @paramname{n} is the number of particles for this fluid specie.
     *
     * @see Aqua::InputOutput::FluidManager
     */
    struct sphFluidParameters
    {
        /// Constructor
        void init();

        /// Destructor
        void destroy();

        /// Number of particles
        unsigned int n;

        /** @brief Gamma exponent \f$ \gamma \f$
         *
         * The EOS \f$ p = \frac{c_s^2 \rho_0}{\gamma}
             \left(
                \left( \frac{\rho}{\rho_0} \right)^\gamma - 1
             \right) \f$ is considered.
         *
         * As it can be appreciated, higher \f$ \gamma \f$ values result in
         * lower incompressible flows.
         *
         * This field can be set with the tag `Option`, for instance:
         * `<Option name="gamma" value="1.0" />`
         */
        float gamma;

        /** @brief Density of reference \f$ \rho_0 \f$.
         *
         * This field can be set with the tag `Option`, for instance:
         * `<Option name="refd" value="998.0" />`
         */
        float refd;

        /** @brief Dynamic viscosity \f$ \mu \f$.
         *
         * This field can be set with the tag `Option`, for instance:
         * `<Option name="Viscdyn" value="0.000894" />`
         */
        float visc_dyn;

        /** @brief Kinematic viscosity \f$ \nu \f$.
         * Resulting from the division \f$ \frac{\mu}{\rho_0} \f$
         */
        float visc_kin;

        /** @brief Artificial viscosity \f$\alpha\f$ value.
         *
         * \f$\alpha\f$ define the artificial viscosity of the fluid, such that
         * the dynamic viscosity will be corrected as
         * \f$ \mu \geq \frac{\alpha}{8} \rho c_s h \f$.
         *
         * This field can be set with the tag `Option`, for instance:
         * `<Option name="alpha" value="0.01" />`
         *
         * @note \f$\alpha \geq 0.01\f$ is suggested if \f$\delta\f$ is null.
         * @see #delta
         */
        float alpha;

        /** @brief Continuity equation diffusive term factor \f$\delta\f$.
         *
         * Formulation known as \f$\delta\f$-SPH.
         *
         * This term may recover the simulation stability with lower entropy
         * effect than the artificial viscosity, and therefore you may consider
         * setting \f$\alpha = 0\f$ and \f$\delta = 1\f$.
         *
         * However the positive sign of the generated entropy is not granted by
         * \f$\delta\f$-SPH (Second Law of Thermodynamics is not fulfilled).
         *
         * This field can be set with the tag `Option`, for instance:
         * `<Option name="delta" value="0.0" />`
         *
         * @see #alpha
         */
        float delta;

        /** @brief Corrected dynamic viscosity.
         *
         * @see #alpha
         */
        float visc_dyn_corrected;

        /** @brief File path which should be read to load the particles
         *
         * This field can be set with the tag `Option`, for instance:
         * `<Load format="ASCII" file="./Fluid.dat" />`
         *
         * @see Aqua::InputOutput::Particles
         * @see #in_format
         */
        char *in_path;

        /** @brief Format of the tile to be read to load the particles
         *
         * This field can be set with the tag `Option`, for instance:
         * `<Load format="ASCII" file="./Fluid.dat" />`
         *
         * @see Aqua::InputOutput::Particles
         * @see #in_path
         */
        char *in_format;

        /** @brief File path which should be written to save the particles
         *
         * This field can be set with the tag `Option`, for instance:
         * `<Save format="VTK" file="./output" />`
         *
         * The extension will be automatically set depending on the selected
         * output format.
         *
         * @see Aqua::InputOutput::Particles
         * @see #out_format
         */
        char *out_path;

        /** @brief Format of the tile to be read to load the particles
         *
         * This field can be set with the tag `Option`, for instance:
         * `<Save format="VTK" file="./output" />`
         *
         * @see Aqua::InputOutput::Particles
         * @see #out_path
         */
        char *out_format;
    }*fluids;

    /// Number of fluids
    unsigned int n_fluids;
    /// Add a new fluid to the list #fluids
    void addFluid();

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

    /** @class sphMoveParameters ProblemSetup.h ProblemSetup.h
     * @brief Motion definition parameters.
     *
     * Several motions can be simultaneously applied in AQUAgpusph, which can be
     * used to impose different motions to each Fluid, or to create motions by
     * composing several sub-motions.
     *
     * These setting are set between the following XML tags:
     * @code{.xml}
        <Movements>
        </Movements>
     * @endcode
     *
     * @see Aqua::CalcServer::Movement::Movement
     */
    class sphMoveParameters
    {
    public:
        /// Constructor
        sphMoveParameters();
        /// Destructor
        ~sphMoveParameters();

        /// Type of motion.
        /**
         *   - -1 = No motion.
         *   - 0 = Quaternion. No external data is loaded, so it can be only
         *     controlled by the OpenCL script.
         *   - 1 = LIQuaternion. Linearly interpolated quaternion data.
         *   - 2 = C1Quaternion. C1 continuous interpolated quaternion data.
         *   - 3 = ScriptQuaternion. Python controlled quaternion data.
         */
        int type;

        /// Motions definition XML file
        char* path;
    };

    /// List of motions
    std::deque<sphMoveParameters*> motions;

    /** @class sphSensorsParameters ProblemSetup.h ProblemSetup.h
     * @brief Sensors definition.
     *
     * The sensors are a special type of particles with null mass, where some
     * fields are interpolated.
     *
     * These setting are set between the following XML tags:
     * @code{.xml}
        <Sensors>
        </Sensors>
     * @endcode
     *
     * @see Aqua::CalcServer::Sensors
     */
    class sphSensorsParameters
    {
    public:
        /// Constructor
        sphSensorsParameters();
        /// Destructor
        ~sphSensorsParameters();

        /// Output frequency
        float fps;

        /// OpenCL script
        char *script;

        /// List of positions
        std::deque<vec> pos;

        /** @brief Add a new sensor to the list.
         *
         * @param position Position of sensor.
         * @return false if all gone right, true otherwise.
         */
        bool add(vec position);
    }SensorsParameters;

    /** \struct sphPortal
     * @brief Outdated data, just ignore it.
     */
    struct sphPortal
    {
        /** \struct Portal
         * @brief Outdated data, just ignore it.
         */
        struct Portal
        {
            /// Start corner
            vec corner;
            /// Up vector
            vec up;
            /// Side vector
            vec side;
            /// Normal
            vec normal;
        }in,out;
    };
    /// List of portal pairs
    std::deque<sphPortal*> portals;

    /** @class sphGhostParticles ProblemSetup.h ProblemSetup.h
     * @brief Ghost particles mirroring boundaries definition.
     *
     * The Ghost particles based boundary condition is similar to the fixed
     * particles one, except because the fluid extension particles are generated
     * by a mirroring process with respect to the boundary.
     *
     * In AQUAgpusph a GPU oriented version called "Virtual Ghost Particles" is
     * used, where the mirrored particles are not stored but are computed in the
     * interactions computation stage.
     *
     * Therefore the boundaries should be defined here to can perform the
     * mirroring later.
     *
     * These setting are set between the following XML tags:
     * @code{.xml}
        <GhostParticles>
        </GhostParticles>
     * @endcode
     *
     * @see Aqua::CalcServer::Boundary::GhostParticles
     */
    class sphGhostParticles
    {
    public:
        #ifdef HAVE_3D
            /** \struct Wall
             * @brief Wall defined by 4 points.
             */
            struct Wall
            {
                /// 1st Corner
                vec p1;
                /// 2nd Corner
                vec p2;
                /// 3rd Corner
                vec p3;
                /// 4th Corner
                vec p4;
                /// Normal to face
                vec n;
                /// 1st Corner velocity
                vec v1;
                /// 2nd Corner velocity
                vec v2;
                /// 3rd Corner velocity
                vec v3;
                /// 4th Corner velocity
                vec v4;
            };
            /// Stored walls
            std::deque<Wall*> walls;

            /** @brief Add a quadrangular wall.
             * @param p1 1st corner of wall.
             * @param p2 2nd corner of wall.
             * @param p3 3rd corner of wall.
             * @param p4 4th corner of wall.
             * @param v1 1st corner of wall velocity.
             * @param v2 2nd corner of wall velocity.
             * @param v3 3rd corner of wall velocity.
             * @param v4 4th corner of wall velocity.
             * @return false if all gone right, true otherwise.
             * @remarks Normal will be computed using corners data,
             * in order to get best results is strongly recommended
             * try to use planar faces.
             * @warning Program will assume corners connected as
             * \f$p_1 \rightarrow p_2 \rightarrow p_3 \rightarrow p_4
             * \rightarrow p_1\f$.
             */
            bool add(vec p1,
                     vec p2,
                     vec p3,
                     vec p4,
                     vec v1,
                     vec v2,
                     vec v3,
                     vec v4);

            /** @brief Add a triangular wall.
             * @param p1 1st corner of wall.
             * @param p2 2nd corner of wall.
             * @param p3 3rd corner of wall.
             * @param v1 1st corner of wall velocity.
             * @param v2 2nd corner of wall velocity.
             * @param v3 3rd corner of wall velocity.
             * @return false if all gone right, true otherwise.
             */
            bool add(vec p1,
                     vec p2,
                     vec p3,
                     vec v1,
                     vec v2,
                     vec v3);
        #else
            /** \struct Wall
             * @brief Wall defined by 2 points.
             */
            struct Wall
            {
                /// 1st Corner
                vec p1;
                /// 2nd Corner
                vec p2;
                /// Normal to face
                vec n;
                /// 1st Corner velocity
                vec v1;
                /// 2nd Corner velocity
                vec v2;
            };
            /// Stored walls
            std::deque<Wall*> walls;

            /** @brief Add wall.
             * @param p1 1st corner of wall.
             * @param p2 2nd corner of wall.
             * @param v1 1st corner of wall velocity.
             * @param v2 2nd corner of wall velocity.
             * @return false if all gone right, true otherwise.
             */
            bool add(vec p1,
                     vec p2,
                     vec v1,
                     vec v2);
        #endif

        /** @brief Pressure extension model.
         *     - 0 = ASM (Antisymmetric)
         *     - 1 = SSM (Symmetric)
         *     - 2 = Takeda (Variation of ASM).
         * @note Hydrostatic pressure will be added.
         */
        uint p_extension;

        /** @brief Normal velocity extension model.
         *     - 0 = ASM (Antisymmetric)
         *     - 1 = SSM (Symmetric)
         *     - 2 = Takeda (Variation of ASM).
         *     - 3 = U0M (No velocity will considered).
         * @note Velocity extension is composed with the wall velocity.
         */
        uint vn_extension;

        /** Tangent velocity extension model.
         *     - 0 = ASM (Antisymmetric)
         *     - 1 = SSM (Symmetric)
         *     - 2 = Takeda (Variation of ASM).
         *     - 3 = U0M (No velocity will considered).
         * @note Velocity extension is composed with the wall velocity.
         */
        uint vt_extension;

    }ghost_particles;

};

}}  // namespace

#endif // PROBLEMSETUP_H_INCLUDED
