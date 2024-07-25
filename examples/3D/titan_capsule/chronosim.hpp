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
 * @brief The simulation of the Titan capsule rigid body, using
 * https://projectchrono.org/
 */

#ifndef INSTALLABLEDEMO_H_INCLUDED
#define INSTALLABLEDEMO_H_INCLUDED

#include <aquagpusph/CalcServer/Tool.hpp>
#include <chrono/physics/ChSystemNSC.h>
#include <chrono/physics/ChBody.h>
#include <chrono/physics/ChForce.h>

namespace Aqua{ namespace CalcServer{

/** @class TitanSim TitanSim.h
 * @brief Titan capsule landing simulation
 */
class TitanSim : public Aqua::CalcServer::Tool
{
public:
    /** Constructor, with same arguments than Tool.
     * @param name Tool name.
     * @param once Run this tool just once. Useful to make initializations.
     */
    TitanSim(const std::string tool_name, bool once);

    /** Destructor.
     */
    ~TitanSim();

    /** Initialize the tool.
     */
    void setup();

protected:
    /** Execute the tool.
     */
    cl_event _execute(const std::vector<cl_event> events);

private:
    /// The pitch angle of the simulation
    float _pitch;

    /// The starting vertical velocity
    float _vel;

    /// The initial position of the center of gravity
    float _cogz;

    /// Chrono simulation
    std::shared_ptr<chrono::ChSystemNSC> _sys;

    /// The titan capsule chrono entity
    std::shared_ptr<chrono::ChBody> _titan;

    /// The force from AQUAgpusph to Chrono
    std::shared_ptr<chrono::ChForce> _force;

    /// The moment from AQUAgpusph to Chrono
    std::shared_ptr<chrono::ChForce> _moment;

};

}}  // namespace

#endif // INSTALLABLEDEMO_H_INCLUDED 
