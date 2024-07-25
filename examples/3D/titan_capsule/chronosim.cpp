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

#include "chronosim.hpp"
#include <aquagpusph/CalcServer/CalcServer.hpp>
#include <aquagpusph/InputOutput/Logger.hpp>
#include <cmath>

extern "C" Aqua::CalcServer::TitanSim* create_object(
    const std::string name, bool once)
{
    return new Aqua::CalcServer::TitanSim(name, once);
}

namespace Aqua{ namespace CalcServer{

TitanSim::TitanSim(const std::string name, bool once)
    : Tool(name, once)
{
}

TitanSim::~TitanSim()
{
}

#define LB2KG 0.4535924
#define IN2M 0.0254
#define VEL 9.66
#define COGz 1.2366120307925486

void
TitanSim::setup()
{
    Tool::setup();

    // Get the configuration variables
    auto vars = CalcServer::singleton()->variables();
    _pitch = *((float*)vars->get("pitch")->get(true));
    _vel = *((float*)vars->get("u0")->get(true));
    _cogz = *((float*)vars->get("cogz")->get(true));
    _pitch *= M_PI / 180.0;

    _sys = chrono_types::make_shared<chrono::ChSystemNSC>();
    _titan = chrono_types::make_shared<chrono::ChBody>();
    _sys->AddBody(_titan);
    _force = chrono_types::make_shared<chrono::ChForce>();
    _moment = chrono_types::make_shared<chrono::ChForce>();
    _titan->AddForce(_force);
    _titan->AddForce(_moment);

    _titan->SetName("Titan");
    _sys->SetGravitationalAcceleration(chrono::ChVector3d(0, 0, -9.81));
    // Simulating Space Capsule Water Landing with Explicit Finite Element
    // Method
    // _titan->SetMass(16200 * LB2KG);
    // _titan->SetInertiaXX(chrono::ChVector3d(
    //     66169823 * LB2KG * IN2M * IN2M,
    //     71179525 * LB2KG * IN2M * IN2M,
    //     80721115 * LB2KG * IN2M * IN2M
    // ));
    // Pitching Angle on Space Capsule Water Landing Using Smooth Particle
    // Hydrodynamic Method
    _titan->SetMass(3900);
    _titan->SetInertiaXX(chrono::ChVector3d(4180, 5270, 5560));

    _force->SetMode(chrono::ChForce::FORCE);
    _force->SetFrame(chrono::ChForce::BODY);
    _force->SetAlign(chrono::ChForce::WORLD_DIR);
    _moment->SetMode(chrono::ChForce::TORQUE);
    _moment->SetFrame(chrono::ChForce::BODY);
    _moment->SetAlign(chrono::ChForce::WORLD_DIR);

    _titan->SetPos(chrono::ChVector3d(0, 0, _cogz));
    _titan->SetLinVel(chrono::ChVector3d(0, 0, -_vel));
    // On NWU coordinates:
    //    x : Positive moment = negative roll = portside goes down
    //    y : Positive moment = negative pitch = bow goes up
    //    z : Positive moment = negative yaw = bow goes to the boardside
    // Thus we are inverting the angles to match the AQUAgpusph criteria
    chrono::ChMatrix33<double> R;
    R.SetFromCardanAnglesXYZ(chrono::ChVector3d(0, -_pitch, 0));
    _titan->SetRot(R);

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
TitanSim::_execute(const std::vector<cl_event> UNUSED_PARAM events)
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
    const chrono::ChVector3d r = _titan->GetPos();
    const chrono::ChVector3d drdt = _titan->GetLinVel();
    const chrono::ChVector3d ddrddt = _titan->GetLinAcc();
    const chrono::ChVector3d a = -_titan->GetRotMat().GetCardanAnglesXYZ();
    vec4 a_prev = *((vec4*)vars->get("motion_a")->get(true));
    std::cout << "DT " << dt << std::endl;
    const chrono::ChVector3d dadt = dt > 0.f ?
        (a - chrono::ChVector3d(a_prev.x, a_prev.y, a_prev.z)) / dt :
        chrono::ChVector3d(0, 0, 0);
    vec4 dadt_prev = *((vec4*)vars->get("motion_dadt")->get(true));
    const chrono::ChVector3d ddaddt = dt > 0.f ?
        (dadt - chrono::ChVector3d(dadt_prev.x, dadt_prev.y, dadt_prev.z)) / dt :
        chrono::ChVector3d(0, 0, 0);
    std::cout << "dadt " << dadt.y() << std::endl;
    std::cout << "dadt_prev " << dadt_prev.y << std::endl;
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
 
