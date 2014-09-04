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
 * @brief Aqua::CalcServer tools base class.
 * (See Aqua::CalcServer::Kernel for details)
 */

#ifndef KERNEL_H_INCLUDED
#define KERNEL_H_INCLUDED

#include <sphPrerequisites.h>

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <CL/cl.h>

#include <AuxiliarMethods.h>

namespace Aqua{ namespace CalcServer{

/** @class Kernel Kernel.h CalcServer/Kernel.h
 * @brief Base class for all the tools of the calculation server
 * Aqua::CalcServer.
 */
class Kernel
{
public:
    /** @brief Constructor.
     * @param kernel_name Kernel name.
     * The kernel name will be used later to refer to the results of the tool.
     */
    Kernel(const char* kernel_name);

    /// Destructor
    virtual ~Kernel();

    /** @brief Maximum allowed local work size for the selected device.
     * @param n Amount of data to solve.
     * If 0 provided, Aqua::CalcServer::N will be used.
     * @param queue Command queue of the computational device. If it is NULL,
     * the first command queue present in Aqua::CalcServer will be selected.
     * @return The maximum local work size.
     */
    virtual size_t localWorkSize(unsigned int n=0,cl_command_queue queue=NULL);

    /** @brief Get the global work size.
     *
     * The global work size is associated with the applied local work size.
     * @param size Local work size.
     * @param n Amount of data to solve. If 0 is provided, CalcServer number
     * of particles will be selected.
     * @return Global work size.
     * @see Aqua::getGlobalWorkSize
     */
    virtual size_t globalWorkSize(size_t size, unsigned int n=0);

    /** @brief Set the kernel name.
     * @param kernel_name Kernel name.
     */
    void name(const char* kernel_name);
    /** @brief Get the kernel name.
     * @return Kernel name.
     */
    const char* name(){return (const char*)_name;}

    #ifdef HAVE_GPUPROFILE
        /** @brief Set the kernel time consumed.
         * @param t Kernel time consumed.
         */
        void profileTime(float t){_time = t;}
        /** @brief Get the kernel time consumed.
         * @return Kernel time consumed.
         */
        float profileTime(){return _time;}
    #endif

private:
    /// Kernel name
    char* _name;
    #ifdef HAVE_GPUPROFILE
        /// Kernel real time consumed
        float _time;
    #endif
};

}}  // namespace

#endif // KERNEL_H_INCLUDED
