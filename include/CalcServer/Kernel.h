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

#ifndef KERNEL_H_INCLUDED
#define KERNEL_H_INCLUDED

// ----------------------------------------------------------------------------
// Include Prerequisites
// ----------------------------------------------------------------------------
#include <sphPrerequisites.h>

// ----------------------------------------------------------------------------
// Include standar libraries
// ----------------------------------------------------------------------------
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

// ----------------------------------------------------------------------------
// Include OpenCL libraries
// ----------------------------------------------------------------------------
#include <CL/cl.h>

// ----------------------------------------------------------------------------
// Include auxiliar methods
// ----------------------------------------------------------------------------
#include <AuxiliarMethods.h>

namespace Aqua{ namespace CalcServer{

/** @class Kernel Kernel.h CalcServer/Kernel.h
 * @brief All calculation subclasses of CalcServer are called Kernel,
 * as OpenCL nomenclature natural extension. See Inherited diagram.
 */
class Kernel
{
public:
	/** Constructor.
	 * @param kernelName Kernel name.
	 */
	Kernel(const char* kernelName);

	/** Destructor
	 */
	~Kernel();

	/** Returns maximum local work size for selected device.
	 * @param n Amount of data to solve. If 0 provided,
	 * CalcServer number of particles will selected.
	 * @param queue Command queue where device is allocated.
	 * If NULL, first command queue on CalcServer will
	 * selected.
	 * @return Global work size.
	 * @note Probably you don't need overload this method.
	 */
	virtual size_t localWorkSize(unsigned int n=0,cl_command_queue queue=NULL);

	/** Returns global work size from designed local size.
	 * @param size Local work size.
	 * @param n Amount of data to solve. If 0 provided,
	 * CalcServer number of particles will selected.
	 * @return Global work size.
	 * @note Probably you don't need overload this method.
	 */
	virtual size_t globalWorkSize(size_t size, unsigned int n=0);

	/** Set the kernel name.
	 * @param kernelName Kernel name.
	 */
	void name(const char* kernelName);
	/** Get the kernel name.
	 * @return Kernel name.
	 */
	const char* name(){return (const char*)mName;}

	#ifdef HAVE_GPUPROFILE
	    /** Set the kernel time consumed.
	     * @param t Kernel time consumed.
	     */
	    void profileTime(float t){_time = t;}
	    /** Get the kernel time consumed.
	     * @return Kernel time consumed.
	     */
	    float profileTime(){return _time;}
	#endif

private:
	/// Kernel name
	char* mName;
	#ifdef HAVE_GPUPROFILE
	    /// Kernel real time consumed
	    float _time;
	#endif
};

}}  // namespace

#endif // KERNEL_H_INCLUDED
