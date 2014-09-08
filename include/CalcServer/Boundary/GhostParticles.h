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
 * @brief Ghost particles fluid extension based boundary condition.
 * (See Aqua::CalcServer::Boundary::GhostParticles for details)
 */

#ifndef GHOSTPARTICLES_H_INCLUDED
#define GHOSTPARTICLES_H_INCLUDED

#include <CalcServer/Kernel.h>

#include <vector>
#include <deque>

namespace Aqua{ namespace CalcServer{ namespace Boundary{

/** @class GhostParticles GhostParticles.h CalcServer/Boundary/GhostParticles.h
 * @brief Boundary condition based on extending the fluid with a mirroring
 * process of the particles with respect to the wall.
 *
 * The fields on the extended fluid are conveniently set.
 *
 * In order to improve the consistency the model used for the pressure, normal
 * velocity, and tangent velocity can be modified, but the following model is
 * set by default (Free slip):
 *   - pressModel = "SSM"
 *   - normalUModel = "ASM"
 *   - tangentUModel = "SSM"
 *
 * @see GhostParticles.cl
 */
class GhostParticles : public Aqua::CalcServer::Kernel
{
public:
    /// Constructor.
    GhostParticles();

    /// Destructor.
    ~GhostParticles();

    /** @brief Perform the work.
     * @return false if all gone right, true otherwise.
     */
    bool execute();

private:
    /** @brief Setup the OpenCL stuff
     * @return false if all gone right, true otherwise.
     */
    bool setupOpenCL();

    /// OpenCL script path
    char* _path;

    /// OpenCL program
    cl_program _program;
    /// OpenCL vertices set kernel
    cl_kernel _kernel;
    /// Global work size (calculated with local_work_size).
    size_t _global_work_size;
    /// Local work size (default value = 256)
    size_t _local_work_size;
    /// true if \f$delta\f$-SPH (cont. eq. diffusive term) must be applied.
    bool _is_delta;
};

}}}  // namespace

#endif // GHOSTPARTICLES_H_INCLUDED
