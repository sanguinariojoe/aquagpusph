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
 * @brief The simulation of a Ball as rigid body, using
 * https://projectchrono.org/
 */

#define _USE_MATH_DEFINES
#include "chronosim.hpp"
#include <chrono/physics/ChBodyEasy.h>
#include <aquagpusph/CalcServer/CalcServer.hpp>
#include <aquagpusph/InputOutput/Logger.hpp>
#include <cmath>
#include <filesystem>
#include <stdexcept>

/** @brief Function called by AQUAgpusph to retrieve a Aqua::CalcServer::Tool
 * derived object that will be treated as any other AQUAgpusoh tool
 * @return BritanSim object
*/
extern "C" Aqua::CalcServer::BritanSim* create_object(
    const std::string name, bool once)
{
    return new Aqua::CalcServer::BritanSim(name, once);
}

namespace Aqua{ namespace CalcServer{

BritanSim::BritanSim(const std::string name, bool once)
    : Tool(name, once)
{
}

// Destructor
BritanSim::~BritanSim()
{
}

void
BritanSim::setup()
{
    Tool::setup();

    // Get the configuration variables from AQUAgpusph (set on SPH.xml)
    auto vars = CalcServer::singleton()->variables();
    const float rho = *((float*)vars->get("RHOP")->get(true));
    const float R = *((float*)vars->get("R")->get(true));

    _sys = chrono_types::make_shared<chrono::ChSystemNSC>();
    _ball = chrono_types::make_shared<chrono::ChBody>();
    _sys->AddBody(_ball);
    _force = chrono_types::make_shared<chrono::ChForce>();
    _moment = chrono_types::make_shared<chrono::ChForce>();
    _ball->AddForce(_force);
    _ball->AddForce(_moment);

    _ball->SetName("Ball");    
    _sys->SetGravitationalAcceleration(chrono::ChVector3d(0, 0, 0));

    const double V = 4.0 / 3.0 * M_PI * pow(R, 3.0);
    _ball->SetMass(V * rho);
    const double Ir = 2.0 / 5.0 * pow(R, 2.0);
    _ball->SetInertiaXX(chrono::ChVector3d(rho * V * Ir));

    _force->SetMode(chrono::ChForce::FORCE);
    _force->SetFrame(chrono::ChForce::BODY);
    _force->SetAlign(chrono::ChForce::WORLD_DIR);
    _force->SetVrelpoint(chrono::ChVector3d(0, 0, 0));
    _moment->SetMode(chrono::ChForce::TORQUE);
    _moment->SetFrame(chrono::ChForce::BODY);
    _moment->SetAlign(chrono::ChForce::WORLD_DIR);
    _moment->SetVrelpoint(chrono::ChVector3d(0, 0, 0));

    _ball->SetPos(chrono::ChVector3d(0, 0, 0));
    _ball->SetLinVel(chrono::ChVector3d(0, 0, 0));

    _sys->SetTimestepperType(chrono::ChTimestepper::Type::EULER_EXPLICIT);
    _sys->Setup();

    setInputDependencies({"dt", "iter_midpoint", "iter_midpoint_max",
                          "Force_p_iset", "Moment_p_iset"});
    setOutputDependencies({"motion_r", "motion_drdt", "motion_ddrddt",
                           "motion_a", "motion_dadt", "motion_ddaddt",
                           "forces_r"});
};

/**
 * @brief Set the Chrono Force value from an AQUAgpusph 3D vector
 * @param var The output Chrono force
 * @param value The input value
 */
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

/** @brief Convert a Chrono vector into an AQUAgpusph one, and set it into an
 * AQUAgpusph variable
 * @param var The output AQUAgpusph variable
 * @param value The input Chrono vector
 */
void
setVec(Aqua::InputOutput::Variable* var, chrono::ChVector3d value)
{
    //create a vec4 vector, populate it
    vec4 v;
    v.x = value.x();
    v.y = value.y();
    v.z = value.z();
    v.w = 0.f;

    //set the aqua variable
    var->set(&v, true);
}

cl_event
BritanSim::_execute(const std::vector<cl_event> UNUSED_PARAM events)
{
    auto vars = CalcServer::singleton()->variables();
    // Check whether we are on the midpoint, or at the final iteration
    const unsigned int iter =
        *((unsigned int*)vars->get("iter_midpoint")->get(true));
    const unsigned int iter_max =
        *((unsigned int*)vars->get("iter_midpoint_max")->get(true));
    const bool is_midpoint = iter < iter_max;
    // Get the forces from AQUAgpusph
    float dt = *((float*)vars->get("dt")->get(true));
    const vec4 F = *((vec4*)vars->get("Force_p_iset")->get(true));
    const vec4 M = *((vec4*)vars->get("Moment_p_iset")->get(true));

    if (iter == 0) {
        // At the beggining of the time step we must copy the results from
        // the other instance of this solver
        vec4 data;
        chrono::ChQuaternion<double> R;
        data = *((vec4*)vars->get("motion_r")->get(true));
        _ball->SetPos(chrono::ChVector3d(data.x, data.y, data.z));
        data = *((vec4*)vars->get("motion_drdt")->get(true));
        _ball->SetLinVel(chrono::ChVector3d(data.x, data.y, data.z));
        data = *((vec4*)vars->get("motion_ddrddt")->get(true));
        _ball->SetLinAcc(chrono::ChVector3d(data.x, data.y, data.z));
        data = *((vec4*)vars->get("motion_a")->get(true));
        R.SetFromCardanAnglesXYZ(chrono::ChVector3d(data.x, data.y, data.z));
        _ball->SetRot(R);
        data = *((vec4*)vars->get("motion_dadt")->get(true));
        _ball->SetAngVelLocal(chrono::ChVector3d(data.x, data.y, data.z));
        data = *((vec4*)vars->get("motion_ddaddt")->get(true));
        _ball->SetAngAccLocal(chrono::ChVector3d(data.x, data.y, data.z));
    }

    // Book-keeping, so we can restore the state within the midpoint iterations
    // This is actually needed also at the final iterator because of the Euler
    // explicit (see below)
    double T = _sys->GetChTime();
    chrono::ChState X(_sys->GetNumCoordsPosLevel(), _sys.get());
    chrono::ChStateDelta V(_sys->GetNumCoordsVelLevel(), _sys.get());
    chrono::ChStateDelta A(_sys->GetNumCoordsVelLevel(), _sys.get());
    chrono::ChVectorDynamic<> L(_sys->GetNumConstraints());
    const chrono::ChVector3d drdt0 = _ball->GetLinVel();
    const chrono::ChVector3d ddrddt0 = _ball->GetLinAcc();
    const chrono::ChVector3d dadt0 = _ball->GetAngVelLocal();
    const chrono::ChVector3d ddaddt0 = _ball->GetAngAccLocal();
    _sys->StateGather(X, V, T);
    _sys->StateGatherAcceleration(A);
    _sys->StateGatherReactions(L);

    // Compute the dynamics
    setForce(_force, F);
    setForce(_moment, M);
    _sys->DoStepDynamics(dt);

    // Explicit Euler is not considering the force we setted, so we must rewind
    // and repeat
    if (!is_midpoint) {
        _sys->SetChTime(T);
        _sys->StateScatterReactions(L);
        _sys->StateScatterAcceleration(A);
        _sys->StateScatter(X, V, T, true);
        setForce(_force, F);
        setForce(_moment, M);
        _sys->DoStepDynamics(dt);
    }

    // Get the new position and angle
    chrono::ChVector3d r = _ball->GetPos();
    const chrono::ChVector3d drdt = _ball->GetLinVel();
    const chrono::ChVector3d ddrddt = _ball->GetLinAcc();
    const chrono::ChVector3d a = _ball->GetRot().GetCardanAnglesXYZ();
    const chrono::ChVector3d dadt = _ball->GetAngVelLocal();
    const chrono::ChVector3d ddaddt = _ball->GetAngAccLocal();

    // On the explicit Euler scheme the position is integrated directly from
    // the velocity at the beggining. We want to use the midpoint velocity
    // instead
    r = r + ddrddt * (0.5 * dt * dt);
    _ball->SetPos(r);

    // Update AQUAgpusph
    if (is_midpoint) {
        _sys->SetChTime(T);
        _sys->StateScatterReactions(L);
        _sys->StateScatterAcceleration(A);
        _sys->StateScatter(X, V, T, true);
        setVec(vars->get("motion_drdt"), 0.5 * (drdt + drdt0));
        setVec(vars->get("motion_ddrddt"), 0.5 * (ddrddt + ddrddt0));
        setVec(vars->get("motion_dadt"), 0.5 * (dadt + dadt0));
        setVec(vars->get("motion_ddaddt"), 0.5 * (ddaddt + ddaddt0));
        vars->populate("motion_drdt");
        vars->populate("motion_ddrddt");
        vars->populate("motion_dadt");
        vars->populate("motion_ddaddt");
    } else {
        setVec(vars->get("motion_r"), r);
        setVec(vars->get("motion_a"), a);
        setVec(vars->get("forces_r"), r);
        setVec(vars->get("motion_drdt"), drdt);
        setVec(vars->get("motion_ddrddt"), ddrddt);
        setVec(vars->get("motion_dadt"), dadt);
        setVec(vars->get("motion_ddaddt"), ddaddt);
        vars->populate("motion_r");
        vars->populate("motion_a");
        vars->populate("forces_r");
        vars->populate("motion_drdt");
        vars->populate("motion_ddrddt");
        vars->populate("motion_dadt");
        vars->populate("motion_ddaddt");
    }
    return NULL;
}

}}  // namespaces
 
