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
 * @brief The simulation of the Apollo capsule rigid body, using
 * https://projectchrono.org/
 */

#include "chronosim.hpp"
#include <aquagpusph/CalcServer/CalcServer.hpp>
#include <aquagpusph/InputOutput/Logger.hpp>
#include <cmath>

extern "C" Aqua::CalcServer::ApolloSim* create_object(
    const std::string name, bool once)
{
    return new Aqua::CalcServer::ApolloSim(name, once);
}

namespace Aqua{ namespace CalcServer{

ApolloSim::ApolloSim(const std::string name, bool once)
    : Tool(name, once)
{
}

ApolloSim::~ApolloSim()
{
}

#define LB2KG 0.4535924
#define IN2M 0.0254

void
ApolloSim::setup()
{
    Tool::setup();

    // Get the configuration variables
    auto vars = CalcServer::singleton()->variables();
    _pitch = *((float*)vars->get("pitch")->get(true));
    _vel = *((float*)vars->get("u0")->get(true));
    _cogz = *((float*)vars->get("cogz")->get(true));
    _pitch *= M_PI / 180.0;

    _sys = chrono_types::make_shared<chrono::ChSystemNSC>();
    _apollo = chrono_types::make_shared<chrono::ChBody>();
    _sys->AddBody(_apollo);
    _force = chrono_types::make_shared<chrono::ChForce>();
    _moment = chrono_types::make_shared<chrono::ChForce>();
    _apollo->AddForce(_force);
    _apollo->AddForce(_moment);

    _apollo->SetName("Apollo");
    _sys->SetGravitationalAcceleration(chrono::ChVector3d(0, 0, -9.81));
    // Simulating Space Capsule Water Landing with Explicit Finite Element
    // Method
    // _apollo->SetMass(16200 * LB2KG);
    // _apollo->SetInertiaXX(chrono::ChVector3d(
    //     66169823 * LB2KG * IN2M * IN2M,
    //     71179525 * LB2KG * IN2M * IN2M,
    //     80721115 * LB2KG * IN2M * IN2M
    // ));
    // Pitching Angle on Space Capsule Water Landing Using Smooth Particle
    // Hydrodynamic Method
    _apollo->SetMass(3900);
    _apollo->SetInertiaXX(chrono::ChVector3d(5560, 5270, 4180));

    _force->SetMode(chrono::ChForce::FORCE);
    _force->SetFrame(chrono::ChForce::BODY);
    _force->SetAlign(chrono::ChForce::WORLD_DIR);
    _force->SetVrelpoint(chrono::ChVector3d(0, 0, 0));
    _moment->SetMode(chrono::ChForce::TORQUE);
    _moment->SetFrame(chrono::ChForce::BODY);
    _moment->SetAlign(chrono::ChForce::WORLD_DIR);
    _moment->SetVrelpoint(chrono::ChVector3d(0, 0, 0));

    _apollo->SetPos(chrono::ChVector3d(0, 0, _cogz));
    _apollo->SetLinVel(chrono::ChVector3d(0, 0, -_vel));
    // On NWU coordinates:
    //    x : Positive moment = positive roll = portside goes up
    //    y : Positive moment = positive pitch = bow goes down
    //    z : Positive moment = positive yaw = bow goes to the portside
    chrono::ChQuaternion<double> R;
    R.SetFromCardanAnglesXYZ(chrono::ChVector3d(0, _pitch, 0));
    _apollo->SetRot(R);

    setInputDependencies({"dt", "Force_p", "Moment_p"});
    setOutputDependencies({"motion_r", "motion_drdt", "motion_ddrddt",
                           "motion_a", "motion_dadt", "motion_ddaddt",
                           "forces_r"});
};

void
setForce(std::shared_ptr<chrono::ChForce> var, vec4 value)
{
    chrono::ChVector3d v(value.x, value.y, value.z);
    if (v.IsNull()) {
        var->SetMforce(0.0);
    } else {
        const auto norm = v.Length();
        v.Scale(1.0 / norm);
        var->SetVrelpoint(chrono::ChVector3d(0, 0, 0));
        var->SetDir(v);
        var->SetMforce(norm);
    }
}

void
setVec(Aqua::InputOutput::Variable* var, chrono::ChVector3d value)
{
    vec4 v;
    v.x = value.x();
    v.y = value.y();
    v.z = value.z();
    v.w = 0.f;
    var->set(&v, true);
}

cl_event
ApolloSim::_execute(const std::vector<cl_event> UNUSED_PARAM events)
{
    auto vars = CalcServer::singleton()->variables();
    // Get the forces from AQUAgpusph
    float dt = *((float*)vars->get("dt")->get(true));
    const vec4 F = *((vec4*)vars->get("Force_p_iset")->get(true));
    const vec4 M = *((vec4*)vars->get("Moment_p_iset")->get(true));
    // Simulate the body motion
    setForce(_force, F);
    setForce(_moment, M);
    _sys->DoStepDynamics(dt);
    // Get the new position and angle
    const chrono::ChVector3d r = _apollo->GetPos();
    const chrono::ChVector3d drdt = _apollo->GetLinVel();
    const chrono::ChVector3d ddrddt = _apollo->GetLinAcc();
    const chrono::ChVector3d a = _apollo->GetRot().GetCardanAnglesXYZ();
    const chrono::ChVector3d dadt = _apollo->GetAngVelLocal();
    const chrono::ChVector3d ddaddt = _apollo->GetAngAccLocal();
    setVec(vars->get("motion_r"), r);
    setVec(vars->get("motion_drdt"), drdt);
    setVec(vars->get("motion_ddrddt"), ddrddt);
    setVec(vars->get("motion_a"), a);
    setVec(vars->get("motion_dadt"), dadt);
    setVec(vars->get("motion_ddaddt"), ddaddt);
    setVec(vars->get("forces_r"), r);
    return NULL;
}

}}  // namespaces
 
