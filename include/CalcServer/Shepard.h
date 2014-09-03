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

#ifndef SHEPARD_H_INCLUDED
#define SHEPARD_H_INCLUDED

// ----------------------------------------------------------------------------
// Include Generic kernel
// ----------------------------------------------------------------------------
#include <CalcServer/Kernel.h>

namespace Aqua{ namespace CalcServer{

/** @class Shepard Shepard.h CalcServer/Shepard.h
 * @brief Renormalization factor applied to forces and density variation
 * rates.
 * @remarks 0th order correction (formerly Shepard correction) is mandatory
 * when DeLeffe boundary condition (formerly boundary integrals) is imposed,
 * but since Shepard term can be a little bit unstable, it is strongly
 * recommended to don't apply it if it is not required.
 */
class Shepard : public Aqua::CalcServer::Kernel
{
public:
    /** Constructor.
     */
    Shepard();

    /** Destructor.
     */
    ~Shepard();

    /** Shepard boundary computation.
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
    /// OpenCL shepard term kernel
    cl_kernel _kernel;
    /// Global work size
    size_t _global_work_size;
    /// Local work size
    size_t _local_work_size;
};

}}  // namespace

#endif // SHEPARD_H_INCLUDED
