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
 * @brief Tool to compute the fluid bounds box, as well as the minimum and
 * maximum velocities.
 * (See Aqua::CalcServer::Bounds for details)
 */

#ifndef BOUNDS_H_INCLUDED
#define BOUNDS_H_INCLUDED

#include <CalcServer/Kernel.h>
#include <CalcServer/Reduction.h>

#ifndef __BOUNDS_COORDS_MAX_OP__
/** @def __BOUNDS_COORDS_MAX_OP__
 * @brief Maximum bounds box coordinates computation.
 */
#define __BOUNDS_COORDS_MAX_OP__ 1
#endif
#ifndef __BOUNDS_COORDS_MIN_OP__
/** @def __BOUNDS_COORDS_MIN_OP__
 * @brief Minimum bounds box coordinates computation.
 */
#define __BOUNDS_COORDS_MIN_OP__ 2
#endif
#ifndef __BOUNDS_VEL_MAX_OP__
/** @def __BOUNDS_COORDS_MIN_OP__
 * @brief Maximum velocity computation.
 */
#define __BOUNDS_VEL_MAX_OP__ 4
#endif
#ifndef __BOUNDS_VEL_MIN_OP__
/** @def __BOUNDS_COORDS_MIN_OP__
 * @brief Minimum velocity computation.
 */
#define __BOUNDS_VEL_MIN_OP__ 8
#endif

namespace Aqua{ namespace CalcServer{

/** @class Bounds Bounds.h CalcServer/Bounds.h
 * @brief Computes the fluid particles bounds box, and the maximum and minimum
 * velocities.
 *
 * The bounds box is defined as the smallest box where all the fluid particles
 * are included inside.
 *
 * To do it this tool is working as follows:
 *   -# The velocity and position fields in the shorted space are copied in a
 *      helper memory buffer.
 *   -# The values of the non fluid particles (boundary or sensor ones) are
 *      filtered, i.e. For the minimum value computation INFINITY values are set
 *      while for the maximum components -INFINITY is used.
 *   -# The corresponding reduction is processed to get the minimum/maximum
 *      value.
 *
 * @see Bounds.cl
 * @see Aqua::InputOutput::Bounds
 * @see Aqua::CalcServer::Reduction
 */
class Bounds : public Aqua::CalcServer::Kernel
{
public:
    /// Constructor.
    Bounds();

    /// Destructor.
    ~Bounds();

    /** @brief Get maximum coordinates of the fluid bounds box.
     * @return Maximum coordinates [m].
     */
    vec maxCoords(){return _pos_max;}

    /** @brief Get minimum coordinates of the fluid bounds box.
     * @return Minimum coordinates [m].
     */
    vec minCoords(){return _pos_min;}

    /** @brief Get maximum computed fluid particles velocity.
     * @return Maximum velocity [m/s].
     */
    vec maxVel(){return _vel_max;}

    /** @brief Get minimum computed fluid particles velocity.
     * @return Minimum velocity [m/s].
     */
    vec minVel(){return _vel_min;}

    /** @brief Perform the work.
     * @return false if all gone right, true otherwise.
     */
    bool execute();

private:
    /** @brief Compute the maximum or the minimum desired value.
     * @param output Output computed value.
     * @param op Operation to compute:
     *   - #__BOUNDS_COORDS_MAX_OP__ for maximum coordinate
     *   - #__BOUNDS_COORDS_MIN_OP__ for minimum coordinate
     *   - #__BOUNDS_VEL_MAX_OP__ for maximum velocity
     *   - #__BOUNDS_VEL_MIN_OP__ for minimum velocity
     * @return false if all gone right, true otherwise.
     */
    bool execute(vec *output, int op);

    /** @brief Setup the OpenCL stuff.
     * @return false if all gone right, true otherwise.
     */
    bool setupBounds();

    /** @brief Setup the reduction tools.
     * @return false if all gone right, true otherwise.
     */
    bool setupReduction();

    /// Server allocated helper memory to perform the reduction.
    cl_mem _device_mem;
    /// Host allocated maximum coordinates.
    vec _pos_max;
    /// Host allocated minimum coordinates.
    vec _pos_min;
    /// Host allocated maximum velocity.
    vec _vel_max;
    /// Host allocated minimum velocity.
    vec _vel_min;
    /// Kernel path
    char *_path;
    /// OpenCL program
    cl_program _program;
    /// Maximum coordinates filter.
    cl_kernel _pos_max_kernel;
    /// Minimum coordinates filter.
    cl_kernel _pos_min_kernel;
    /// Maximum velocity filter.
    cl_kernel _vel_max_kernel;
    /// Minimum velocity filter.
    cl_kernel _vel_min_kernel;
    /// Global work size
    size_t _global_work_size;
    /// Local work size
    size_t _local_work_size;
    /// Maximum coordinates reduction tool
    Reduction *_pos_max_reduction;
    /// Minimum coordinates reduction tool
    Reduction *_pos_min_reduction;
    /// Maximum velocity reduction tool
    Reduction *_vel_max_reduction;
    /// Minimum velocity reduction tool
    Reduction *_vel_min_reduction;
};

}}  // namespace

#endif // BOUNDS_H_INCLUDED
