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

#include <CL/cl.h>

#include <vector>
#include <map>
#include <string>
#include <iterator>

#include <sphPrerequisites.h>
#include <ProblemSetup.h>
#include <TimeManager.h>
#include <Variable.h>
#include <Singleton.h>
#include <CalcServer/Tool.h>

#ifndef _ITEMS
    /** @def _ITEMS
     * @brief Number of items in a group
     * @note Must be power of 2, and in some devices greather than 32.
     * @see Aqua::CalcServer::RadixSort
     */
    #define _ITEMS  128
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


namespace Aqua{
/// @namespace Aqua::CalcServer Calculation server name space.
namespace CalcServer{

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
    user_interruption(const std::string msg) : std::runtime_error(msg) {};
};

/** @class CalcServer CalcServer.h CalcServer.h
 * @brief Entity devoted to perform the main work of the simulation.
 *
 * In the Aqua::CalcServer::CalcServer a time subloop is performed where the
 * SPH simulation is performed while no output files should be updated.
 * @note Updating output files require to download data from the server, which
 * due to the low bandwidth asigned is usually a bottleneck, hence letting the
 * Aqua::CalcServer::CalcServer works without interrumptions is a great
 * optimization technique.
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

    /** @brief Internal time loop.
     * 
     * Calculation server will be iterating while no output files should be
     * updated (or even the simulation is finished).
     * @param t_manager Time manager to let the calculation server when shall
     * stop the internal loop.
     */
    void update(InputOutput::TimeManager& t_manager);

    /** @brief Setup some additional simulation data.
     * Even thought this work is associated with the constructor CalcServer(),
     * when something may fail it is preferable to let it to a separated method
     * that could report errors, allowing the program to deal with them.
     */
    void setup();

    /** @brief Download a unsorted variable from the device
     *
     * AQUAgpusph is specifically designed to allow particles be sorted by their
     * cell index, in such a way the list of particles in a specific cell can
     * be computed by means of an increasing index. However, at the time of
     * saving information it is usually required to get the information sorted
     * again.
     *
     * @param var_name Variable to unsort and download
     * @param offset The offset in bytes in the memory object to read from
     * @param cb The size in bytes of data being downloaded
     * @param ptr The host memory where the data should be copied
     * @return The data download event, NULL if errors are detected
     * @note The caller must wait for the event be complete (clWaitForEvents)
     * before accessing the downloaded data
     * @note The caller must call clReleaseEvent() on the returned event
     * when finished
     */
    cl_event getUnsortedMem(const std::string& var_name,
                            const size_t& offset,
                            const size_t& cb,
                            void *ptr);

    /** @brief Get the variables manager
     * @return Variables manager
     */
    inline InputOutput::Variables* variables() {return &_vars;}

    /** @brief Get the registered definitions.
     * @return List of definitions.
     */
    inline const std::vector<std::string>& definitions() const {
        return _definitions;
    }

    /** @brief Get the registered tools
     * @return List of tools.
     */
    inline const std::vector<Tool*>& tools() const {return _tools;}

    /** @brief Get the current OpenCL context
     * @return OpenCL context
     */
    inline const cl_context& context() const {return _context;}

    /** @brief Get the platform
     * @return OpenCL platform
     */
    inline const cl_platform_id& platform() const {return _platform;}

    /** @brief Get the device
     * @return OpenCL device
     */
    inline const cl_device_id& device() const {return _device;}

    /** @brief Get the command queue
     * @return OpenCL command queue
     */
    inline const cl_command_queue& command_queue() const {
        return _command_queue;
    }

    /** @brief Get the AQUAgpusph root path.
     * @return AQUAgpusph root path
     */
    inline const std::string& base_path() const {return _base_path;}

private:
    /** @brief Setup the OpenCL stuff
     */
    void setupOpenCL();
    /** @brief Prints all the available OpenCL platforms and devices
     */
    void queryOpenCL();
    /** @brief Get the selected platform from the available ones
     */
    void setupPlatform();
    /** @brief Compute the available devices in the selected platform.
     */
    void setupDevices();

    /// Simulation data read from XML files
    Aqua::InputOutput::ProblemSetup _sim_data;

    /// Number of available OpenCL platforms
    cl_uint _num_platforms;
    /// List of OpenCL platforms
    cl_platform_id *_platforms;
    /// Number of devices
    cl_uint _num_devices;
    /// List of devices
    cl_device_id *_devices;
    /// OpenCL context
    cl_context _context;
    /// OpenCL command queue
    cl_command_queue *_command_queues;
    /// Selected platform
    cl_platform_id _platform;
    /// Selected device
    cl_device_id _device;
    /// Selected command queue
    cl_command_queue _command_queue;

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
};

}}  // namespace

#endif // CALCSERVER_H_INCLUDED
