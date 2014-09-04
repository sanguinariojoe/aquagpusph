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
 * @brief Density field geometrical interpolation.
 * (See Aqua::CalcServer::DensityInterpolation for details)
 */

#ifndef DENSITYINTERPOLATION_H_INCLUDED
#define DENSITYINTERPOLATION_H_INCLUDED

#include <CalcServer/Kernel.h>

namespace Aqua{ namespace CalcServer{

/** @class DensityInterpolation DensityInterpolation.h CalcServer/DensityInterpolation.h
 * @brief Density interpolation is a smoothing density field technique.
 *
 * When the \f$\alpha\f$ and \f$\delta\f$ parameters are not large enough a
 * characteristic noise may be appreciated in the pressure field.
 *
 * In order to partially fix it some authors proposed to interpolate the density
 * field, replacing the value resulting from the continuity equation.
 *
 * The density field is interpolated with the usual SPH expression:
 * \f$ \rho_i = \frac{1}{\gamma_i} \sum_j W \left(
        \mathbf{r}_j - \mathbf{r}_i
   \right) m_j \f$,
 * where \f$ \gamma_i \f$ is the Shepard factor.
 * @see DensInt.cl
 * @see Aqua::InputOutput::ProblemSetup::sphSPHParameters::dens_int_steps
 * @warning This model may cause some instabilities, and its consistency is not
 * formalized, hence it is strongly recommended to don't apply it if you don't
 * know what are you doing.
 */
struct DensityInterpolation : public Aqua::CalcServer::Kernel
{
    /// Constructor.
    DensityInterpolation();

    /// Destructor.
    ~DensityInterpolation();

    /** @brief Perform the work
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
    /// OpenCL kernel
    cl_kernel _kernel;
    /// Global work size
    size_t _global_work_size;
    /// Local work size
    size_t _local_work_size;
};

}}  // namespace

#endif // DENSITYINTERPOLATION_H_INCLUDED
