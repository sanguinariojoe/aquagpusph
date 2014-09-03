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

#ifndef ELASTICBOUNCE_H_INCLUDED
#define ELASTICBOUNCE_H_INCLUDED

// ----------------------------------------------------------------------------
// Include Generic kernel
// ----------------------------------------------------------------------------
#include <CalcServer/Kernel.h>

namespace Aqua{ namespace CalcServer{ namespace Boundary{

/** @class ElasticBounce ElasticBounce.h CalcServer/Boundary/ElasticBounce.h
 * @brief The simplest boundary condition, which look for the particles which
 * will trespass a wall, performing an elastic bounce (with a configurable
 * elastic factor, shuch that the particle trayectory can be mirrored with
 * respect to the wall or the particle can be blocked just changing the
 * elastic factor).
 * @note This boundary condition is ussually used in a combination with one of
 * the other boundary conditions just in order to grant that the particles
 * cannot trespass the solid walls.
 */
class ElasticBounce : public Aqua::CalcServer::Kernel
{
public:
    /** Constructor.
     */
    ElasticBounce();

    /** Destructor.
     */
    ~ElasticBounce();

    /** Compute the boundary condition.
     * @return false if all gone right, true otherwise.
     */
    bool execute();

private:
    /** Setup the OpenCL stuff
     * @return false if all gone right, true otherwise.
     */
    bool setupOpenCL();

    /// OpenCL script path
    char* _path;

    /// OpenCL program
    cl_program _program;
    /// OpenCL vertices set kernel
    cl_kernel _kernel;
    /// Wall element radius (it is guessed as a circular element)
    cl_float _r;
    /// Global work size (calculated with local_work_size).
    size_t _global_work_size;
    /// Local work size (default value = 256)
    size_t _local_work_size;
};

}}}  // namespace

#endif // ELASTICBOUNCE_H_INCLUDED
