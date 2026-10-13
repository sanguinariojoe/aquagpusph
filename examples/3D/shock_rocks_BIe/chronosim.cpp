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
#define ROCK_DENSITY 2.0
// The envelope size, that should match the blender setup
#define ENVELOPE_SIZE 0.001

/** @brief Creates a function "create_object" that returns a 
           RocksSim object. RocksSim is the derived class 
           from Aqua::CalcServer::Tool
           Extern "C" makes this C namestyle
    @return RocksSim object
*/
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

// initiates a tool element, Aqua type
RocksSim::RocksSim(const std::string name, bool once)
    : Tool(name, once)
{
}

// Destructor
RocksSim::~RocksSim()
{
}

void
RocksSim::setup()
{
    //Tool is from Aqua
    Tool::setup();

    // Get the configuration variables
    auto vars = CalcServer::singleton()->variables();

    // get number of solids
    const unsigned int n_solids =
        *((unsigned int*)vars->get("n_solids")->get(true));
    //get the Length L
    const float L = *((float*)vars->get("L")->get(true));
    //const float rho = *((float*)vars->get("REFD")->get(true));
    
    // Setup the chrono system
    // ChSystemNSC Non-Smooth Contact 
    // Rigid collision and friction
    _sys = chrono_types::make_shared<chrono::ChSystemNSC>();

    // Set Bullet collision system 
    // the alternative is default Chrono system
    // bullet was previous default
    _sys->SetCollisionSystemType(chrono::ChCollisionSystem::Type::BULLET);

    // Set gravity for the simmulation
    _sys->SetGravitationalAcceleration(chrono::ChVector3d(0, 0, -9.81));

    // Setup the floor
    // Uses impulses & Coulomb friction
    auto ground_mat =
        chrono_types::make_shared<chrono::ChContactMaterialNSC>();
    
    // very high friction
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
    // Setup the format of the objects
    // the rock files 
    // are called rock.0.stl rock.00.stl rock.000.stl
    // depending on the number of rocks
    // how may digits to have
    const unsigned int digits = num_digits(n_solids);
    
    // setup every rock
    for (unsigned int i=0; i < n_solids; i++) {

        auto trimesh = chrono::ChTriangleMeshConnected::CreateFromSTLFile(
            std::string("rock.") + int2string(i, digits) + ".subdivided.stl");
        
        double vol;
        chrono::ChVector3d cog;// cog center of gravity?
        chrono::ChMatrix33<> inertia;

        // ComputeMassProperties gets density volume center of gravity and innertia
        // here True in first parameter means full body
        // and not only surface.
        // oputput is per unit of density
        // that is everything must be multiplied by density afterwards
        trimesh->ComputeMassProperties(true, vol, cog, inertia);

        // move the center of the ball to the center of gravity 
        // and apply there the inertia
        // chrono::ChMatrix33<>(1) means do not rotate the ball
        trimesh->Transform(-cog, chrono::ChMatrix33<>(1));

        // Setup the rock body
        auto rock = chrono_types::make_shared<chrono::ChBody>();
        _rocks.push_back(rock);//method of std::vector
        _sys->Add(rock);// Add to simulation
        
        rock->SetName(std::string("rock.") + std::to_string(i));
        
        // note multiplication of previously per unit of density magnitudes
        rock->SetMass(vol * ROCK_DENSITY);
        rock->SetInertia(inertia * ROCK_DENSITY);

        // Setup the collision model
        // The rock is somewhere, so we are displacing the mesh to the origin,
        // and then we are moving the body to the COG
        auto coll_model =
            chrono_types::make_shared<chrono::ChCollisionModel>();
        coll_model->SetSafeMargin(0.01);
        coll_model->SetEnvelope(ENVELOPE_SIZE);

        // material of the objects
        auto rock_mat =
            chrono_types::make_shared<chrono::ChContactMaterialNSC>();
        
        // Material uses impulses & Coulomb friction
        rock_mat->SetFriction(0.005);
        rock_mat->SetDampingF(0.01);
        
        // Triangles for the collision model of chrono
        // paramters are
        // Defined material, STL mess, is_static, is_convex, 
        // radius thinkness of the skin of the triangles
        auto coll_shape =
            chrono_types::make_shared<chrono::ChCollisionShapeTriangleMesh>(
                rock_mat, trimesh, false, false, 0.010f);
        
        // ChFrame<> The ChFrame defines the relative position and orientation
        // of the shape with respect to the body's Center of Gravity (COG).
        // chrono::ChVector3d(0, 0, 0): This means the center of the STL mesh coincides with the body's COG.
        // This only works correctly because of the previously performed trimesh->Transform(-cog, ...)
        // chrono::QUNIT: This is the "Identity Quaternion." It represents zero rotation. The mesh will be oriented exactly as it was designed in the STL file.
        // AddShape(shape, frame)
        coll_model->AddShape(coll_shape, chrono::ChFrame<>(
            chrono::ChVector3d(0, 0, 0), chrono::QUNIT));

        rock->AddCollisionModel(coll_model);
        rock->EnableCollision(true);// enable collison model
        rock->SetPos(cog);// moving the ball back to its cog

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
        force->SetMode(chrono::ChForce::FORCE);// force and not torque
        force->SetFrame(chrono::ChForce::BODY);// The force is attached to the ball's "nose." 
        //If the ball flips, the force flips with it (like a jet engine). 
        // position of the forces set via SetVrelpoint() is interpreted in Local Coordinates of the ball.
        // It is set at the point (0,0,0), the force is "attached" to the center of the ball. 
        // Even as the ball flies across the map, the force stays perfectly centered on the ball's mass.
        force->SetAlign(chrono::ChForce::WORLD_DIR);
        // This defines the Direction where the vector points. The direction of the force is fixed relative to the Inertial World Frame (the X, Y, Z axes of the universe).
        // Even if the ball starts tumbling or spinning at high speeds after the impact, the force will always point in the same direction (e.g., always pushing "East").
        force->SetVrelpoint(chrono::ChVector3d(0, 0, 0));//see previous comment


        //set the torque
        auto torque = chrono_types::make_shared<chrono::ChForce>();
        rock->AddForce(torque);
        _torques.push_back(torque);
        torque->SetMode(chrono::ChForce::TORQUE);
        torque->SetFrame(chrono::ChForce::BODY);
        torque->SetAlign(chrono::ChForce::WORLD_DIR);
        torque->SetVrelpoint(chrono::ChVector3d(0, 0, 0));

    }//loop on particles finish here

    
    // _sys->SetTimestepperType(chrono::ChTimestepper::Type::EULER_IMPLICIT);
    _sys->Setup();//final step before one enters the simulation loop

    std::vector<std::string> indeps({"dt"}), outdeps;// creates two
    // vectors of strings, indeps and outdeps, indeps initialized with 
    // list {"dt"} and outdeps non initalized

    for (unsigned int i=0; i < n_solids; i++) {
        // adds members to the vector, 
        // custom co-simulation wrappers
        // setting which variables must be updated or satisfied
        // before the next calculation step (inchrono) can proceed
        indeps.push_back(
            std::string("rock_") + int2string(i, digits) + "_Force_p");
        indeps.push_back(
            std::string("rock_") + int2string(i, digits) + "_Moment_p");

        // add members to the vector outdeps
        // outdeps 
        //defines what the system writes or provides to other 
        //modules after a calculation is finished.
        outdeps.push_back(
            std::string("rock_") + int2string(i, digits) + "_forces_r");
        outdeps.push_back(
            std::string("rock_") + int2string(i, digits) + "_motion_r");
        outdeps.push_back(
            std::string("rock_") + int2string(i, digits) + "_motion_drdt");
        outdeps.push_back(
            std::string("rock_") + int2string(i, digits) + "_motion_ddrddt");
        outdeps.push_back(
            std::string("rock_") + int2string(i, digits) + "_motion_a");
        outdeps.push_back(
            std::string("rock_") + int2string(i, digits) + "_motion_dadt");
        outdeps.push_back(
            std::string("rock_") + int2string(i, digits) + "_motion_ddaddt");
    }
};

/**
 * @brief Set the Force object * 
 * @param var The Chrono project force pointer one
 * whises to get
 * @param value the Opencl vector that is already known
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

/** @brief get a vector from Chrono and set an Aqua variable
 * @param var The aqua variable that will be set
 * @param value The known Chrono vector
 */
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


/**
 * @brief Calculate time step integration
 * for every moving object
 * get force and momentum
 * calculate movement in the time step dt
 * actualize position, velocity, etc...
 * of each object
 * @param
 * @return This munction returns a cl_event that 
 * is fixed to NULL
 */
cl_event
RocksSim::_execute(const std::vector<cl_event> UNUSED_PARAM events)
{
    //  type of aqua
    auto vars = CalcServer::singleton()->variables();

    //get the timestep, convert to float* deindiriction
    float dt = *((float*)vars->get("dt")->get(true));

    // Apply the forces to the rocks
    unsigned int n_solids = _rocks.size();//amount of moving objects
    const unsigned int digits = num_digits(n_solids);// format
    for (unsigned int i = 0; i < n_solids; i++) {
        //get force and momentum from aqua
        // forces and moments are know at this stage
        // that is the singleton object know them
        std::string prefix = std::string("rock_") + int2string(i, digits);
        const vec4 F = *((vec4*)vars->get(prefix + "_Force_p")->get(true));
        const vec4 M = *((vec4*)vars->get(prefix + "_Moment_p")->get(true));

        //set the forces for Chrono
        setForce(_forces[i], F);
        setForce(_torques[i], M);
    }

    // Compute the dynamics
    _sys->DoStepDynamics(dt);


    // Get the new positions and angles
    for (unsigned int i = 0; i < n_solids; i++) {
        std::string prefix = std::string("rock_") + int2string(i, digits);
        
        //_rocks defined in hpp
        // here the new possition, vel, etc is got from Chrono
        auto rock = _rocks[i];
        const chrono::ChVector3d r = rock->GetPos(); // Possition
        const chrono::ChVector3d drdt = rock->GetLinVel(); //modern Chrono GetPos_dt
        const chrono::ChVector3d ddrddt = rock->GetLinAcc(); //modern Chrono GetPos_dtdt
        const chrono::ChVector3d a = rock->GetRot().GetCardanAnglesXYZ(); //GetRot quaternion GetCardanAnglesXYZ angles
        const chrono::ChVector3d dadt = rock->GetAngVelLocal(); // Angular velocity
        const chrono::ChVector3d ddaddt = rock->GetAngAccLocal(); // Angular Acceleration

        // new things are transferred
        setVec(vars->get(prefix + "_forces_r"), r);
        setVec(vars->get(prefix + "_motion_r"), r);
        setVec(vars->get(prefix + "_motion_drdt"), drdt);
        setVec(vars->get(prefix + "_motion_ddrddt"), ddrddt);
        setVec(vars->get(prefix + "_motion_a"), a);
        setVec(vars->get(prefix + "_motion_dadt"), dadt);
        setVec(vars->get(prefix + "_motion_ddaddt"), ddaddt);

        // do not know that is this probably just putting thigs inside of aqua
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
 
