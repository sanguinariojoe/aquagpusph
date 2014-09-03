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

#ifndef TIMESTEP_H_INCLUDED
#define TIMESTEP_H_INCLUDED

#include <CalcServer/Kernel.h>
#include <CalcServer/Reduction.h>

namespace Aqua{ namespace CalcServer{

/** @class TimeStep TimeStep.h CalcServer/TimeStep.h
 * @brief Time step computation. The time step could be set in 3 different
 * ways:
 *   -# Provided fixed time step: The used set a fixed time step.
 *   -# Computed fixed time step: The time step will be computed at \f$t=0\f$
 *   -# Variable time step: The time step is recomputed each iteration.
 * To compute a time step a courant condition is imposed such that a particle
 * cannot move more than \f$ 0.1 h \f$ per time step, or a pressure wave
 * cannot be transported more than \f$ h \f$ in a time step
 * (\f$ t \leq \frac{h}{\mathrm{max}(c_s, 10 \vert\mathbf{u}\vert_{max})} \f$)
 */
class TimeStep : public Aqua::CalcServer::Kernel
{
public:
    /** Constructor.
     */
    TimeStep();

    /** Destructor.
     */
    ~TimeStep();

    /** Time step computation.
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
    /// OpenCL kernel
    cl_kernel _kernel;
    /// Global work size (calculated with local_work_size).
    size_t _global_work_size;
    /// Local work size (default value = 256)
    size_t _local_work_size;

    /// Convective time step reduction tool
    Reduction *_reduction;

    /// Main time step
    float _dt;

    /// Time step clamped flag (in order to avoid the annoying clamping error)
    bool _is_dt_clamp;
};

}}  // namespace

#endif // TIMESTEP_H_INCLUDED
