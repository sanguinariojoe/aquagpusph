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
#include <chrono/physics/ChBodyEasy.h>
#include <aquagpusph/CalcServer/CalcServer.hpp>
#include <aquagpusph/InputOutput/Logger.hpp>
#include <cmath>

// The density of the rock material
#define ROCK_DENSITY 89.4
// The envelope size, that should match the blender setup
#define ENVELOPE_SIZE 0.001

extern "C" Aqua::CalcServer::RocksSim* create_object(
    const std::string name, bool once)
{
    return new Aqua::CalcServer::RocksSim(name, once);
}

namespace Aqua{ namespace CalcServer{

/** @brief Get the number of digits of a number
 * @param n the Number
 * @return The number of digits
 */
unsigned int num_digits(unsigned int n)
{
    unsigned int digits = 1;
    while (n /= 10)
        digits++;
    return digits;
}

/** @brief Convert an unsigned integer to a string, appending leading zeroes
 * @param n the Number
 * @param digits The minimum number of digits of the resulting string
 * @return The string
 */
std::string int2string(unsigned int n, unsigned int digits)
{
    std::string str = std::to_string(n);
    if (str.length() < digits)
        str.insert(0, digits - str.length(), '0');

    return str;
}

RocksSim::RocksSim(const std::string name, bool once)
    : Tool(name, once)
{
}

RocksSim::~RocksSim()
{
}

void
RocksSim::setup()
{
    Tool::setup();

    // Get the configuration variables
    auto vars = CalcServer::singleton()->variables();
    const unsigned int n_solids =
        *((unsigned int*)vars->get("n_solids")->get(true));
    const float L = *((float*)vars->get("L")->get(true));
    //const float rho = *((float*)vars->get("REFD")->get(true));
    
    printf("Chronosim: Number of solids %u", n_solids);
    printf("Chronosim: Number of solids %f", L);

    // Setup the chrono system
    _sys = chrono_types::make_shared<chrono::ChSystemNSC>();
    _sys->SetCollisionSystemType(chrono::ChCollisionSystem::Type::BULLET);
    _sys->SetGravitationalAcceleration(chrono::ChVector3d(0, 0, 0));

    // Setup the floor
    auto ground_mat =
        chrono_types::make_shared<chrono::ChContactMaterialNSC>();
    ground_mat->SetFriction(1.0);

    _ground = chrono_types::make_shared<chrono::ChBodyEasyBox>(
        2 * L, 2 * L, 2 * L,            // Box size (we choose the position later)
        1e6,                            // Density (fixed, it does not matter)
        false,                          // No visual needed
        true,                           // Collisions enabled
        ground_mat);                    // Material
    _ground->SetName("g");
    _ground->GetCollisionModel()->SetEnvelope(ENVELOPE_SIZE);
    _sys->AddBody(_ground);
    _ground->SetPos(chrono::ChVector3d(0.0, 0.0, -(L + ENVELOPE_SIZE)));
    _ground->SetFixed(true);

    // Setup the rocks
    const unsigned int digits = num_digits(n_solids);


    for (unsigned int i=0; i < n_solids; i++) {
        auto trimesh = chrono::ChTriangleMeshConnected::CreateFromSTLFile(
            std::string("ball.") + int2string(i, digits) + ".subdivided.stl");
        double vol;
        chrono::ChVector3d cog;
        chrono::ChMatrix33<> inertia;
        trimesh->ComputeMassProperties(true, vol, cog, inertia);
        trimesh->Transform(-cog, chrono::ChMatrix33<>(1));

        // Setup the rock body
        auto rock = chrono_types::make_shared<chrono::ChBody>();
        _rocks.push_back(rock);
        _sys->Add(rock);
        rock->SetName(std::string("ball.") + std::to_string(i));
        rock->SetMass(vol * ROCK_DENSITY);
        rock->SetInertia(inertia * ROCK_DENSITY);

        // Setup the collision model
        // The rock is somewhere, so we are displacing the mesh to the origin,
        // and then we are moving the body to the COG
        auto coll_model =
            chrono_types::make_shared<chrono::ChCollisionModel>();
        coll_model->SetSafeMargin(0.01);
        coll_model->SetEnvelope(ENVELOPE_SIZE);
        auto rock_mat =
            chrono_types::make_shared<chrono::ChContactMaterialNSC>();
        rock_mat->SetFriction(0.005);
        rock_mat->SetDampingF(0.01);
        auto coll_shape =
            chrono_types::make_shared<chrono::ChCollisionShapeTriangleMesh>(
                rock_mat, trimesh, false, false, 0.010f);
        coll_model->AddShape(coll_shape, chrono::ChFrame<>(
            chrono::ChVector3d(0, 0, 0), chrono::QUNIT));
        rock->AddCollisionModel(coll_model);
        rock->EnableCollision(true);
        rock->SetPos(cog);

        // Add the forces
        /*auto bouyancy = chrono_types::make_shared<chrono::ChForce>();
        rock->AddForce(bouyancy);
        bouyancy->SetMode(chrono::ChForce::FORCE);
        bouyancy->SetFrame(chrono::ChForce::BODY);
        bouyancy->SetAlign(chrono::ChForce::WORLD_DIR);
        bouyancy->SetVrelpoint(chrono::ChVector3d(0, 0, 0));
        bouyancy->SetDir(chrono::ChVector3d(0, 0, 1));
        bouyancy->SetMforce(vol * rho * 9.81);*/

        auto force = chrono_types::make_shared<chrono::ChForce>();
        rock->AddForce(force);
        _forces.push_back(force);
        force->SetMode(chrono::ChForce::FORCE);
        force->SetFrame(chrono::ChForce::BODY);
        force->SetAlign(chrono::ChForce::WORLD_DIR);
        force->SetVrelpoint(chrono::ChVector3d(0, 0, 0));
        auto torque = chrono_types::make_shared<chrono::ChForce>();
        rock->AddForce(torque);
        _torques.push_back(torque);
        torque->SetMode(chrono::ChForce::TORQUE);
        torque->SetFrame(chrono::ChForce::BODY);
        torque->SetAlign(chrono::ChForce::WORLD_DIR);
        torque->SetVrelpoint(chrono::ChVector3d(0, 0, 0));
    }

    
    // _sys->SetTimestepperType(chrono::ChTimestepper::Type::EULER_IMPLICIT);
    _sys->Setup();

    std::vector<std::string> indeps({"dt"}), outdeps;
    for (unsigned int i=0; i < n_solids; i++) {
        indeps.push_back(
            std::string("ball_") + int2string(i, digits) + "_Force_p");
        indeps.push_back(
            std::string("ball_") + int2string(i, digits) + "_Moment_p");
        outdeps.push_back(
            std::string("ball_") + int2string(i, digits) + "_forces_r");
        outdeps.push_back(
            std::string("ball_") + int2string(i, digits) + "_motion_r");
        outdeps.push_back(
            std::string("ball_") + int2string(i, digits) + "_motion_drdt");
        outdeps.push_back(
            std::string("ball_") + int2string(i, digits) + "_motion_ddrddt");
        outdeps.push_back(
            std::string("ball_") + int2string(i, digits) + "_motion_a");
        outdeps.push_back(
            std::string("ball_") + int2string(i, digits) + "_motion_dadt");
        outdeps.push_back(
            std::string("ball_") + int2string(i, digits) + "_motion_ddaddt");
    }
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
RocksSim::_execute(const std::vector<cl_event> UNUSED_PARAM events)
{
    auto vars = CalcServer::singleton()->variables();
    float dt = *((float*)vars->get("dt")->get(true));

    // Apply the forces to the rocks
    unsigned int n_solids = _rocks.size();
    const unsigned int digits = num_digits(n_solids);
    for (unsigned int i = 0; i < n_solids; i++) {
        std::string prefix = std::string("ball_") + int2string(i, digits);
        const vec4 F = *((vec4*)vars->get(prefix + "_Force_p")->get(true));
        const vec4 M = *((vec4*)vars->get(prefix + "_Moment_p")->get(true));

        setForce(_forces[i], F);
        setForce(_torques[i], M);
    }

    // Compute the dynamics
    _sys->DoStepDynamics(dt);


    // Get the new positions and angles
    for (unsigned int i = 0; i < n_solids; i++) {
        std::string prefix = std::string("ball_") + int2string(i, digits);
        auto rock = _rocks[i];
        const chrono::ChVector3d r = rock->GetPos();
        const chrono::ChVector3d drdt = rock->GetLinVel();
        const chrono::ChVector3d ddrddt = rock->GetLinAcc();
        const chrono::ChVector3d a = rock->GetRot().GetCardanAnglesXYZ();
        const chrono::ChVector3d dadt = rock->GetAngVelLocal();
        const chrono::ChVector3d ddaddt = rock->GetAngAccLocal();
        setVec(vars->get(prefix + "_forces_r"), r);
        setVec(vars->get(prefix + "_motion_r"), r);
        setVec(vars->get(prefix + "_motion_drdt"), drdt);
        setVec(vars->get(prefix + "_motion_ddrddt"), ddrddt);
        setVec(vars->get(prefix + "_motion_a"), a);
        setVec(vars->get(prefix + "_motion_dadt"), dadt);
        setVec(vars->get(prefix + "_motion_ddaddt"), ddaddt);
        vars->populate(prefix + "_forces_r");
        vars->populate(prefix + "_motion_r");
        vars->populate(prefix + "_motion_drdt");
        vars->populate(prefix + "_motion_ddrddt");
        vars->populate(prefix + "_motion_a");
        vars->populate(prefix + "_motion_dadt");
        vars->populate(prefix + "_motion_ddaddt");
    }
    return NULL;
}

}}  // namespaces
 
