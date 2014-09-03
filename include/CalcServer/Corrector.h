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

#ifndef CORRECTOR_H_INCLUDED
#define CORRECTOR_H_INCLUDED

#include <CalcServer/Kernel.h>

namespace Aqua{ namespace CalcServer{

/** @class Corrector Corrector.h CalcServer/Corrector.h
 * @brief Corrector stage. Time integration uses a quasi-second order
 * Predictor-Corrector integration scheme:
 *   - \f$ \mathbf{u}_{n+1} = \mathbf{u}_{n} + \Delta t \left(
        \mathbf{g} +
        \left. \frac{\mathrm{d}\mathbf{u}}{\mathrm{d}t} \right\vert_{n+1/2}
     \right)
     + \frac{\Delta t}{2} \left(
        \left. \frac{\mathrm{d}\mathbf{u}}{\mathrm{d}t} \right\vert_{n + 1/2} -
        \left. \frac{\mathrm{d}\mathbf{u}}{\mathrm{d}t} \right\vert_{n - 1/2}
     \right)
     \f$
 *   - \f$ \mathbf{r}_{n+1} = \mathbf{r}_{n} + \Delta t \, \mathbf{u}_{n}
     + \frac{\Delta t^2}{2} \left(
        \mathbf{g} +
        \left. \frac{\mathrm{d}\mathbf{u}}{\mathrm{d}t} \right\vert_{n+1/2}
     \right)
     \f$
 *   - \f$ \rho_{n+1} = \rho_{n} + \Delta t
        \left. \frac{\mathrm{d}\rho}{\mathrm{d}t} \right\vert_{n+1/2}
     + \frac{\Delta t}{2} \left(
        \left. \frac{\mathrm{d}\rho}{\mathrm{d}t} \right\vert_{n + 1/2} -
        \left. \frac{\mathrm{d}\rho}{\mathrm{d}t} \right\vert_{n - 1/2}
     \right)
     \f$
 */
struct Corrector : public Aqua::CalcServer::Kernel
{
    /** Constructor.
     */
    Corrector();

    /** Destructor
     */
    ~Corrector();

    /** Executes the time integration corrector stage.
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
    /// OpenCL main kernel
    cl_kernel _kernel;
    /// Velocity clamping kernel
    cl_kernel _vel_clamp_kernel;
    /// Global work size
    size_t _global_work_size;
    /// Local work size
    size_t _local_work_size;
};

}}  // namespace

#endif // CORRECTOR_H_INCLUDED
