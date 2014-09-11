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

#include <sphPrerequisites.h>
#include <Variable.h>
#include <Singleton.h>

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

    /// Transfer the data from the computational device to the host
    /** Transferring data from/to the computational device (managed by
     * Aqua::CalcServer::CalcServer) is a slow operation, so try to avoid
     * calling this method too often.
     *
     * @param dest Array where the data should be copied.
     * @param orig Computational device allocated data to copy.
     * @param size Size of the data to copy.
     * @param offset Offset into the array to start reading.
     * @return false if the data has been successfully copied, true otherwise.
     */
    bool getData(void *dest, cl_mem orig, size_t size, size_t offset=0);

    /// Transfer the data from the host device to the computational device.
    /** Transferring data from/to the computational device (managed by
     * Aqua::CalcServer::CalcServer) is a slow operation, so try to avoid
     * calling this method too often.
     *
     * @param dest Identifier of the destination in the server.
     * @param orig Array of data to read.
     * @param size Size of the data to copy.
     * @return false if the data has been successfully copied, true otherwise.
     */
    bool sendData(cl_mem dest, void* orig, size_t size);

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
private:
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
};

}}  // namespace

#endif // CALCSERVER_H_INCLUDED
