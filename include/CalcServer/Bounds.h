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

#include <CalcServer/Kernel.h>
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
 * @brief Computes the fluid particles bounds, which includes the coordinates
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

    /** Get maximum computed fluid particles coordinates.
     * @return Maximum coordinates [m].
     */
    vec maxCoords(){return _pos_max;}

    /** Get minimum computed fluid particles coordinates.
     * @return Minimum coordinates [m].
     */
    vec minCoords(){return _pos_min;}

    /** Get maximum computed fluid particles velocity.
     * @return Maximum velocity [m/s].
     */
    vec maxVel(){return _vel_max;}

    /** Get minimum computed fluid particles velocity.
     * @return Minimum velocity [m/s].
     */
    vec minVel(){return _vel_min;}

    /** Compute the bounds.
     * @return false if all gone right, true otherwise.
     */
    bool execute();

private:
    /** Compute the maximum or the minimum desired value.
     * @param output Output computed value.
     * @param data Input data array.
     * @param op Operation to compute:
     *   - __BOUNDS_COORDS_MAX_OP__ for maximum coordinate
     *   - __BOUNDS_COORDS_MIN_OP__ for minimum coordinate
     *   - __BOUNDS_VEL_MAX_OP__ for maximum velocity
     *   - __BOUNDS_VEL_MIN_OP__ for minimum velocity
     * @return false if all gone right, true otherwise.
     */
    bool execute(vec *output, int op);

    /** Setup the OpenCL stuff.
     * @return false if all gone right, true otherwise.
     */
    bool setupBounds();

    /** Setup the reduction tool.
     * @return false if all gone right, true otherwise.
     */
    bool setupReduction();

    /// Server allocated auxiliar memory to perform the reduction.
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
    /** Maximum coordinates per particle kernel. The maximum coordinates for
     * each particle will be the position of the particle for fluid particles,
     * and the vector (vec)(-INFINITY,-INFINITY,-INFINITY,0.f) for the other
     * types of particles.
     */
    cl_kernel _pos_max_kernel;
    /** Minimum coordinates per particle kernel. The minimum coordinates for
     * each particle will be the position of the particle for fluid particles,
     * and the vector (vec)(INFINITY,INFINITY,INFINITY,0.f) for the other
     * types of particles.
     */
    cl_kernel _pos_min_kernel;
    /** Maximum velocity per particle kernel. The maximum velocity for
     * each particle will be the velocity of the particle for fluid particles,
     * and the vector (vec)(0.f,0.f,0.f,0.f) for the other
     * types of particles.
     */
    cl_kernel _vel_max_kernel;
    /** Minimum velocity per particle kernel. The minimum velocity for
     * each particle will be the velocity of the particle for fluid particles,
     * and the vector (vec)(INFINITY,INFINITY,INFINITY,0.f) for the other
     * types of particles.
     */
    cl_kernel _vel_min_kernel;
    /// Global work size
    size_t _global_work_size;
    /// Local work size
    size_t _local_work_size;
    /// Maximum coordiantes reduction tool
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
