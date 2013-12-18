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

#ifndef BOUNDS_H_INCLUDED
#define BOUNDS_H_INCLUDED

// ----------------------------------------------------------------------------
// Include Generic kernel
// ----------------------------------------------------------------------------
#include <CalcServer/Kernel.h>

// ----------------------------------------------------------------------------
// Include Reduction tool
// ----------------------------------------------------------------------------
#include <CalcServer/Reduction.h>

#ifndef __BOUNDS_COORDS_MAX_OP__
#define __BOUNDS_COORDS_MAX_OP__ 1
#endif
#ifndef __BOUNDS_COORDS_MIN_OP__
#define __BOUNDS_COORDS_MIN_OP__ 2
#endif
#ifndef __BOUNDS_VEL_MAX_OP__
#define __BOUNDS_VEL_MAX_OP__ 4
#endif
#ifndef __BOUNDS_VEL_MIN_OP__
#define __BOUNDS_VEL_MIN_OP__ 8
#endif

namespace Aqua{ namespace CalcServer{

/** @class Bounds Bounds.h CalcServer/Bounds.h
 * @brief Computes fluid particles bounds, including coordinates
 * bounds and maximum and minimum velocities.
 */
class Bounds : public Aqua::CalcServer::Kernel
{
public:
	/** Constructor.
	 */
	Bounds();

	/** Destructor.
	 */
	~Bounds();

	/** Get maximum fluid particles coordinates.
	 * @return Maximum coordinates [m].
	 * @warning Remember call execute before use this method.
	 */
	vec maxCoords(){return mMaxCoords;}

	/** Get minimum fluid particles coordinates.
	 * @return Minimum coordinates [m].
	 * @warning Remember call execute before use this method.
	 */
	vec minCoords(){return mMinCoords;}

	/** Get maximum fluid particles velocity.
	 * @return Maximum velocity [m/s].
	 * @warning Remember call execute before use this method.
	 */
	vec maxVel(){return mMaxVel;}

	/** Get minimum fluid particles velocity.
	 * @return Minimum velocity [m/s].
	 * @warning Remember call execute before use this method.
	 */
	vec minVel(){return mMinVel;}

	/** Compute the bounds.
	 * @return false if all gone right. \n true otherwise.
	 */
	bool execute();

private:
	/** Compute the maximum or the minimum desired value.
	 * @param output Output computed value.
	 * @param data Input data array.
	 * @param op Operation to compute, __BOUNDS_COORDS_MAX_OP__,
	 * __BOUNDS_COORDS_MIN_OP__, __BOUNDS_VEL_MAX_OP__ or __BOUNDS_VEL_MIN_OP__.
	 * @return false if all gone right. \n true otherwise.
	 */
	bool execute(vec *output, int op);

	/** Setup Bounds OpenCL stuff.
	 * @return false if all gone right. \n true otherwise.
	 */
	bool setupBounds();

	/** Setup Reductions stuff.
	 * @return false if all gone right. \n true otherwise.
	 */
	bool setupReduction();

	/// Server allocated auxiliar memory.
	cl_mem mDevMem;
	/// Host allocated maximum coordinates.
	vec mMaxCoords;
	/// Host allocated minimum coordinates.
	vec mMinCoords;
	/// Host allocated maximum velocity.
	vec mMaxVel;
	/// Host allocated minimum velocity.
	vec mMinVel;
	/// Kernel path
	char *mPath;
	/// OpenCL program
	cl_program program;
	/// OpenCL maximum coordinates setup kernel
	cl_kernel clMaxCoordsKernel;
	/// OpenCL minimum coordinates setup kernel
	cl_kernel clMinCoordsKernel;
	/// OpenCL maximum velocity setup kernel
	cl_kernel clMaxVelKernel;
	/// OpenCL minimum velocity setup kernel
	cl_kernel clMinVelKernel;
	/// Global work size
	size_t global_work_size;
	/// Local work size
	size_t local_work_size;
    /// Maximum coordiantes reduction tool
    Reduction *maxCoordsReduction;
    /// Minimum coordinates reduction tool
    Reduction *minCoordsReduction;
    /// Maximum velocity reduction tool
    Reduction *maxVelReduction;
    /// Minimum velocity reduction tool
    Reduction *minVelReduction;
};

}}  // namespace

#endif // BOUNDS_H_INCLUDED
