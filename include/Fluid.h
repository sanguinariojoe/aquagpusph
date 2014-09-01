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
 * @brief Fluid particles host allocated data.
 * (See Aqua::InputOutput::Fluid for details)
 */

#ifndef FLUID_H_INCLUDED
#define FLUID_H_INCLUDED

#include <sphPrerequisites.h>

#include <Singleton.h>

namespace Aqua{ namespace InputOutput{

/** @class Fluid Fluid.h Fluid.h
 * @brief Fluid particles data allocated in the host.
 * This class is a intermediate layer to download/uploaded the data between the
 * host (see ::main) and the calculation server (see
 * Aqua::CalcServer::CalcServer).
 *
 * More specifically, the particles data loaded from the input files will be
 * stored here, and then transfered to the server.
 *
 * Later, during the simulation, the particles data will be downloaded from the
 * server and stored here to be printed in the output files.
 *
 * @see Aqua::InputOutput::FileManager
 * @see Aqua::InputOutput::Particles
 */
class Fluid : public Aqua::Singleton<Aqua::InputOutput::Fluid>
{
public:
    /** Constructor.
     */
    Fluid();

    /** Destructor.
     */
    ~Fluid();

    /** Retrieves all data from the calculation server.
     *
     * @return false if all gone right, true otherwise.
     */
    bool retrieveData();

    /** Get the number of particles in this fluid.
     *
     * @return Number of particles.
     */
    unsigned int n(){return num_particles;}

    /** Get the number of particles in this fluid.
     *
     * @return Number of particles.
     * @see n()
     */
    unsigned int nParticles(){return n();}

    /** Get the number of fluids.
     *
     * @return Number of fuid species.
     */
    unsigned int nFluids(){return num_fluids;}

    /** Fixed/sensor/regular part flag:
     *   - imove > 0 for regular fluid particles.
     *   - imove = 0 for sensors.
     *   - imove < 0 for boundary elements/particles.
     */
    int *imove;
    /// Fluid identifier (Several Fluid instances may exist).
    int *ifluid;
    /// Position \f$ \mathbf{r} \f$.
    vec *pos;
    /// Normal \f$ \mathbf{n} \f$.
    vec *normal;
    /// Velocity \f$ \mathbf{u} \f$.
    vec *v;
    /// Density \f$ \rho \f$.
    float *dens;
    /// Density variation rate \f$ \frac{\mathrm{d}\rho}{\mathrm{d}t} \f$.
    float *drdt;
    /// Acceleration \f$ \frac{\mathrm{d}\mathbf{u}}{\mathrm{d}t} \f$.
    vec *f;
    /// Mass \f$ m \f$.
    float *mass;
    /// Pressure \f$ p \f$.
    float *press;
    /** Shepard term \f$ \gamma(\mathbf{x}) = \int_{\Omega}
           W(\mathbf{y} - \mathbf{x})
           \mathrm{d}\mathbf{x} \f$.
     */
    float *shepard;
private:
    /// Number of particles
    unsigned int num_particles;
    /// Number of fluids
    unsigned int num_fluids;
};

}}  // namespace

#endif // FLUID_H_INCLUDED
