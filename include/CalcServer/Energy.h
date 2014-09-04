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
 * @brief Tool to compute the fluid energy components.
 * (See Aqua::CalcServer::Energy for details)
 */

#ifndef ENERGY_H_INCLUDED
#define ENERGY_H_INCLUDED

#include <CalcServer/Kernel.h>
#include <CalcServer/Reduction.h>

namespace Aqua{ namespace CalcServer{

/** @class Energy Energy.h CalcServer/Energy.h
 * @brief Computes the fluid energy components.
 *
 * The following energy components are considered:
 *   -# Potential energy: \f$ E_{pot} = - \sum_i m_i
     \mathbf{g} \cdot \mathbf{r}_i \f$.
 *   -# Kinetic energy: \f$ E_{kin} = \sum_i \frac{1}{2} m_i
     \vert \mathbf{u}_i \vert^2 \f$.
 *   -# Internal energy: \f$ U = \int_0^t \sum_i \frac{p_i}{\rho_i^2}
     \left(
        \frac{\mathrm{d} \rho_i}{\mathrm{d} t}
        - \left. \frac{\mathrm{d} \rho_i}{\mathrm{d} t} \right\vert_F
     \right) m_i \mathrm{d}t \f$.
 *   -# Enthalpy: \f$ H = \int_0^t \sum_i \frac{p_i}{\rho_i^2}
     \frac{\mathrm{d} \rho_i}{\mathrm{d} t} m_i \mathrm{d}t \f$.
 *   -# Entropy: \f$ TS = U - H \f$.
 *   -# Total energy: \f$ E = U + E_{kin} \f$.
 *
 * This tool is computing each of the aforementioned energy components for each
 * particle, reducing it later with a prefix sum.
 *
 * @remarks Since some energy components must be integrated in time, a low
 * energy file output/update frequency may imply poor results.
 * @see Energy.cl
 * @see Aqua::InputOutput::Energy
 * @see Aqua::CalcServer::Reduction
 */
class Energy : public Aqua::CalcServer::Kernel
{
public:
    /// Constructor.
    Energy();

    /// Destructor.
    ~Energy();

    /** @brief Get the total energy.
     * @return Total energy: \f$ E = U + E_{kin} \f$.
     * @warning The viscous dissipation is not implemented yet.
     */
    float energy(){return _E.x + _E.w;}

    /** @brief Get the internal energy.
     * @return Internal energy: \f$ U = \int_0^t \sum_i \frac{p_i}{\rho_i^2}
         \left(
            \frac{\mathrm{d} \rho_i}{\mathrm{d} t}
            - \left. \frac{\mathrm{d} \rho_i}{\mathrm{d} t} \right\vert_F
         \right) m_i \mathrm{d}t \f$.
     * @warning The viscous dissipation is not implemented yet.
     */
    float internalEnergy(){return _E.x;}

    /** @brief Get the enthalpy.
     * @return Enthalpy: \f$ H = \int_0^t \sum_i \frac{p_i}{\rho_i^2}
         \frac{\mathrm{d} \rho_i}{\mathrm{d} t} m_i \mathrm{d}t \f$.
     * @warning The viscous dissipation is not implemented yet.
     */
    float enthalpy(){return _E.y;}

    /** @brief Get the entropy.
     * @return Entropy: \f$ TS = U - H \f$.
     * @warning The viscous dissipation is not implemented yet.
     */
    float entropy(){return _E.x - _E.y;}

    /** @brief Get the potential energy.
     * @return Potential energy: \f$ E_{pot} = - \sum_i m_i
         \mathbf{g} \cdot \mathbf{r}_i \f$.
     */
    float potentialEnergy(){return _E.z;}

    /** @brief Get the total kinetic energy.
     * @return Kinetic energy: \f$ E_{kin} = \sum_i \frac{1}{2} m_i
         \vert \mathbf{u}_i \vert^2 \f$.
     */
    float kineticEnergy(){return _E.w;}

    /** @brief Perform the work.
     * @return false if all gone right, true otherwise.
     */
    bool execute();

private:
    /** @brief Setup the OpenCL stuff.
     * @return false if all gone right, true otherwise.
     */
    bool setupEnergy();

    /** @brief Setup the reduction tools.
     * @return false if all gone right, true otherwise.
     */
    bool setupReduction();

    /// Server allocated energy array.
    cl_mem _device_energy;
    /// Host allocated energy rate of change.
    vec4 _dEdt;
    /// Host allocated energy.
    vec4 _E;
    /// Last time when the energy was computed
    float _time;
    /// Kernel path
    char *_path;
    /// OpenCL program
    cl_program _program;
    /// OpenCL kernel
    cl_kernel _kernel;
    /// Global work size
    size_t _global_work_size;
    /// Local work size
    size_t _local_work_size;
    /// Energy values reduction tool
    Reduction *_reduction;
};

}}  // namespace

#endif // ENERGY_H_INCLUDED
