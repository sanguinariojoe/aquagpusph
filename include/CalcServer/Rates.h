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
 * @brief Fluid interactions computation.
 * (See Aqua::CalcServer::Rates for details)
 */

#ifndef RATES_H_INCLUDED
#define RATES_H_INCLUDED

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

/** @class Rates Rates.h CalcServer/Rates.h
 * @brief Fluid interactions computation.
 *
 * The following data is computed:
 *    -# \f$ \frac{\mathrm{d} \mathbf{u}}{\mathrm{d} t} \f$
 *    -# \f$ \frac{\mathrm{d} \rho}{\mathrm{d} t} \f$
 *    -# \f$ \gamma = \sum_i \frac{W \left(
            \mathbf{r}_j - \mathbf{r}_i
       \right)}{\rho_j} m_j \f$.
 *    -# \f$ \nabla \gamma = \sum_i \frac{\nabla W \left(
            \mathbf{r}_j - \mathbf{r}_i
       \right)}{\rho_j} m_j \f$
 *
 * In this process the fluid particles interactions are computed, as well as
 * the interactions with the fixed particles.
 *
 * For optimization purposes, during the same process, some data is interpolated
 * at the boundary elements and sensors.
 *
 * @see Rates.cl
 */
class Rates : public Aqua::CalcServer::Kernel
{
public:
    /// Constructor.
    Rates();

    /// Destructor.
    ~Rates();

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
    char *_path;

    /// OpenCL program
    cl_program _program;
    /// Particles interaction kernel
    cl_kernel _kernel;
    /// Global work size.
    size_t _global_work_size;
    /// Local work size
    size_t _local_work_size;
    /// true if \f$delta\f$-SPH (cont. eq. diffusive term) must be applied.
    bool _is_delta;
};

}}  // namespace

#endif // RATES_H_INCLUDED
