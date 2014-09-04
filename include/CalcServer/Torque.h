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
 * @brief Tool to compute the fluid global force and moment.
 * (See Aqua::CalcServer::Torque for details)
 */

#ifndef TORQUE_H_INCLUDED
#define TORQUE_H_INCLUDED

#include <CalcServer/Kernel.h>
#include <CalcServer/Reduction.h>

namespace Aqua{ namespace CalcServer{

/** @class Torque Torque.h CalcServer/Torque.h
 * @brief Computes fluid global force and moment.
 *
 * The moment is computed respect to a provided point called Center Of Rotation
 * @paramname{COR}.
 *
 * This tool is computing the force and the moment for each particle,
 * conveniently reducing it later with a prefix sum.
 *
 * @see Torque.cl
 * @see Aqua::CalcServer::Reduction
 */
class Torque : public Aqua::CalcServer::Kernel
{
public:
    /// Constructor.
    Torque();

    /// Destructor.
    ~Torque();

    /** @brief Set the Center of rotation.
     * @param COR Center of rotation.
     */
    void cor(vec COR){_cor = COR;}

    /** @brief Get the Center of rotation.
     * @return Center of rotation.
     */
    vec cor(){return _cor;}

    /** @brief Get the moment respect to the center of rotation.
     * @return Fluid torque.
     * @see cor().
     */
    vec torque(){return _torque;}

    /** @brief Get the force.
     * @return Fluid force.
     */
    vec force(){return _force;}

    /** @brief Perform the work.
     * @return false if all gone right, true otherwise.
     */
    bool execute();

private:
    /** @brief Setup the OpenCL stuff.
     * @return false if all gone right, true otherwise.
     */
    bool setupTorque();

    /** @brief Setup the reduction tool
     * @return false if all gone right, true otherwise.
     */
    bool setupReduction();

    /// Center of rotation
    vec _cor;

    /// Server allocated torque.
    cl_mem _device_torque;
    /// Host allocated torque
    vec _torque;
    /// Server allocated torque.
    cl_mem _device_force;
    /// Host allocated torque
    vec _force;
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
    /// Torque value reduction tool
    Reduction *_torque_reduction;
    /// Force value reduction tool
    Reduction *_force_reduction;
};

}}  // namespace

#endif // TORQUE_H_INCLUDED
