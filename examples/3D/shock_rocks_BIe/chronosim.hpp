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
 * @brief The simulation of a st of rocks defined as a rigid body, using
 * https://projectchrono.org/
 */

#ifndef INSTALLABLEDEMO_H_INCLUDED
#define INSTALLABLEDEMO_H_INCLUDED

#include <aquagpusph/CalcServer/Tool.hpp>
#include <chrono/physics/ChSystemNSC.h>
#include <chrono/physics/ChBody.h>
#include <chrono/physics/ChForce.h>

namespace Aqua{ namespace CalcServer{

/** @class RocksSim RocksSim.h
 * @brief Rocks on a flow simulation
 */
class RocksSim : public Aqua::CalcServer::Tool
{
public:
    /** Constructor, with same arguments than Tool.
     * @param name Tool name.
     * @param once Run this tool just once. Useful to make initializations.
     */
    RocksSim(const std::string tool_name, bool once);

    /** Destructor.
     */
    ~RocksSim();

    /** Initialize the tool.
     */
    void setup();

protected:
    /** Execute the tool.
     */
    cl_event _execute(const std::vector<cl_event> events);

private:
    /// Chrono simulation
    std::shared_ptr<chrono::ChSystemNSC> _sys;

    /// The ground
    std::shared_ptr<chrono::ChBody> _ground;

    /// The rocks
    std::vector<std::shared_ptr<chrono::ChBody>> _rocks;

    /// The hydrodynamic forces
    std::vector<std::shared_ptr<chrono::ChForce>> _forces;

    /// The hydrodynamic torques
    std::vector<std::shared_ptr<chrono::ChForce>> _torques;

    /// Velocity at the begginning of the time step
    std::vector<std::shared_ptr<chrono::ChVector3d>> _drdt0;

    /// Acceleration at the begginning of the time step
    std::vector<std::shared_ptr<chrono::ChVector3d>> _ddrddt0;

    /// Angular velocity at the begginning of the time step
    std::vector<std::shared_ptr<chrono::ChVector3d>> _dadt0;

    /// Angular acceleration at the begginning of the time step
    std::vector<std::shared_ptr<chrono::ChVector3d>> _ddaddt0;
};

}}  // namespace

#endif // INSTALLABLEDEMO_H_INCLUDED 
