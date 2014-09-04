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
 * @brief Particles out of domain filter.
 * (See Aqua::CalcServer::Domain for details)
 */

#ifndef DOMAIN_H_INCLUDED
#define DOMAIN_H_INCLUDED

#include <CalcServer/Kernel.h>

namespace Aqua{ namespace CalcServer{

/** @class Domain Domain.h CalcServer/Domain.h
 * @brief Check for particles out of the domain.
 *
 * When a particle out of bounds is detected it is converted in a fixed zero
 * mass particle such that it is not taken into account in the Link-List process
 * anymore.
 *
 * If this filter is not applied exist a risk about singular particles that may
 * scape from the domain, forcing the generation of a huge list of cells, and
 * eventually causing the simulation fail.
 *
 * The total mass lost due to the particles out of the domain will be computed
 * along the simulation.
 * @see Domain.cl
 * @see Aqua::InputOutput::ProblemSetup::sphSPHParameters::has_domain
 */
class Domain : public Aqua::CalcServer::Kernel
{
public:
    /// Constructor.
    Domain();

    /// Destructor.
    ~Domain();

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
    /// OpenCL kernel
    cl_kernel _kernel;
    /// Global work size
    size_t _global_work_size;
    /// Local work size
    size_t _local_work_size;
};

}}  // namespace

#endif // DOMAIN_H_INCLUDED
