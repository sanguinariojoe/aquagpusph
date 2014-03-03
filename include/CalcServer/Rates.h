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

#ifndef RATES_H_INCLUDED
#define RATES_H_INCLUDED

#include <CalcServer/Kernel.h>

namespace Aqua{ namespace CalcServer{

/** @class Rates Rates.h CalcServer/Rates.h
 * @brief Compute the particles interactions to get the following data:
 *    -# \f$ \frac{\mathrm{d} \mathbf{u}}{\mathrm{d} t} \f$
 *    -# \f$ \frac{\mathrm{d} \rho}{\mathrm{d} t} \f$
 *    -# \f$ \gamma = \sum_i \frac{W \left(
            \mathbf{r}_j - \mathbf{r}_i
       \right)}{\rho_j} m_j \f$.
 *    -# \f$ \nabla \gamma = \sum_i \frac{\nabla W \left(
            \mathbf{r}_j - \mathbf{r}_i
       \right)}{\rho_j} m_j \f$
 * In this process the interaction between fluid particles, and the effect of
 * the fixed particles over the fluid particles is included.
 * Some data for the boundary elements and for the sensors is computed as
 * well, like the interpolated pressure and density.
 */
class Rates : public Aqua::CalcServer::Kernel
{
public:
	/** Constructor.
	 */
	Rates();

	/** Destructor.
	 */
	~Rates();

	/** Compute the particles interaction.
	 * @return false if all gone right, true otherwise.
	 */
	bool execute();

private:
	/** Setup the OpenCL stuff
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
