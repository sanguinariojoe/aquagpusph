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

#include <deque>

#include <sphPrerequisites.h>
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

/** @class CalcServer CalcServer.h CalcServer.h
 * @brief Entity that perform the main work of the simulation.
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
    /// Constructor.
    /** The constructor will just initialize default values, however setup()
     * should be called after that.
     */
    CalcServer();

    /// Destructor
    ~CalcServer();

    /// Internal time loop.
    /** Calculation server will be iterating while no output files should be
     * updated (or the simulation is finished).
     *
     * @return false if all gone right, true otherwise.
     */
    bool update();

    /// Setup some additional simulation data.
    /** Even thought this work is associated with the constructor CalcServer(),
     * when something may fail it is preferable to let it to a separated method
     * that could report errors, allowing the program to deal with them.
     * @return false if the calculation server has been successfully setup,
     * true otherwise
     */
    bool setup();

    /** Get the variables manager
     * @return Variables manager
     */
    InputOutput::Variables* variables() const {return _vars;}

    /** Get the tools registered.
     * @return List of tools.
     */
    std::deque<Tool*> tools() const {return _tools;}

    /** Get the active context
     * @return OpenCL context
     */
    cl_context context() const{return _context;}

    /** Get the platform
     * @return OpenCL platform
     */
    cl_platform_id platform() const{return _platform;}

    /** Get the device
     * @return OpenCL device
     */
    cl_device_id device() const{return _device;}

    /** Get the command queue
     * @return OpenCL command queue
     */
    cl_command_queue command_queue() const{return _command_queue;}

    /** Download a unsorted variable from the device.
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
    cl_event getUnsortedMem(const char* var_name,
                            size_t offset,
                            size_t cb,
                            void *ptr);
private:
    /// Setup the OpenCL stuff.
    /**
     * @return false if the OpenCL environment has been succesfully built,
     * true otherwise
     */
    bool setupOpenCL();
    /// Prints all the available platforms and devices returned by OpenCL.
    /**
     * @return false if the OpenCL environment can be succesfully built,
     * true otherwise
     */
    bool queryOpenCL();
    /// Get a platform from the available ones.
    /**
     * @return false if a platform could be obtained, true otherwise
     */
    bool setupPlatform();
    /// Get the available devices in the selected platform.
    /**
     * @return false if the devices have been succesfully obtained, true
     * otherwise
     */
    bool setupDevices();

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
    InputOutput::Variables *_vars;

    /// User registered tools
    std::deque<Tool*> _tools;
};

}}  // namespace

#endif // CALCSERVER_H_INCLUDED
