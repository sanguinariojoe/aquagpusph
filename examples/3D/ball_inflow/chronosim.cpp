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

#include "chronosim.hpp"
#include <chrono/physics/ChBodyEasy.h>
#include <aquagpusph/CalcServer/CalcServer.hpp>
#include <aquagpusph/InputOutput/Logger.hpp>
#include <cmath>
#include <filesystem>
#include <stdexcept>

// The density of the rock material
#define ROCK_DENSITY 89.4
//#define ROCK_DENSITY 0.894

//#define TO_KG_MM3 1.0e-9
//#define TO_M2 1.0e-6

// The envelope size, that should match the blender setup
// is this value correct?
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
        // for example 7/=10 is 0
        // problems source in my opision
        // python script does not follow
        // same logic
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
        // insert at the beguining 0 in first possitional argument
        str.insert(0, digits - str.length(), '0');

    return str;
}

// Arcane meaning for non-initiated
// initiates a tool element. What is this?
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
    printf("Chronosim: Starting the setup\n");
    //Tool is Aqua I believe
    Tool::setup();

    // Get the configuration variables
    auto vars = CalcServer::singleton()->variables();

    // unsafe cast to (unsigned int*) of vars sold true
    // dereference with *()
    // very arcane operation that clearly get number of solids
    const unsigned int n_solids =
        *((unsigned int*)vars->get("n_solids")->get(true));
    //simmilar arcane operation to get the Length L
    printf("Chronosim: getting L\n");
    const float L = *((float*)vars->get("L")->get(true));
    //const float rho = *((float*)vars->get("REFD")->get(true));
    printf("Chronosim: getting B\n");
    const float B = *((float*)vars->get("B")->get(true));
    printf("Chronosim: getting H\n");
    const float H = *((float*)vars->get("H")->get(true));
    


    printf("Chronosim: Number of solids %u\n", n_solids);

    printf("Chronosim: Value of L %f\n", L);
    printf("Chronosim: Value of B %f\n", B);
    printf("Chronosim: Value of H %f\n", H);

    // Setup the chrono system
    // ChSystemNSC Non-Smooth Contact 
    // Rigid collision and friction
    _sys = chrono_types::make_shared<chrono::ChSystemNSC>(); 

    // Set Bullet collision system 
    // the alternative is default Chrono system
    // bullet was previous default
    _sys->SetCollisionSystemType(chrono::ChCollisionSystem::Type::BULLET);

    // Set zero gravity for the ball simmulation
    _sys->SetGravitationalAcceleration(chrono::ChVector3d(0, 0, 0));

    // Setup the floor
    // Uses impulses & Coulomb friction
    auto ground_mat =
        chrono_types::make_shared<chrono::ChContactMaterialNSC>();

    // very high friction
    ground_mat->SetFriction(1.0);

    _ground = chrono_types::make_shared<chrono::ChBodyEasyBox>(
        L, B, H,                        // Box size (we choose the position later)
        1e6,                            // Density (fixed, it does not matter)
        false,                          // No visual needed
        true,                           // Collisions enabled
        ground_mat);                    // Material
    _ground->SetName("g");
    _ground->GetCollisionModel()->SetEnvelope(ENVELOPE_SIZE);
    _sys->AddBody(_ground);
    _ground->SetPos(chrono::ChVector3d(0.0, 0.0, -(H + ENVELOPE_SIZE)));
    //_ground->SetPos(chrono::ChVector3d(0.0, 0.0, -(L + ENVELOPE_SIZE)));
    // Grund is static
    _ground->SetFixed(true);
    printf("Ground defined\n");

    // Setup the format of the objects
    // the ball files 
    // are called bal.0.stl bal.00.stl bal.000.stl
    // how may digits to have
    const unsigned int digits = num_digits(n_solids);


    // setup every rock
    for (unsigned int i=0; i < n_solids; i++) {

        printf("Loading object %d\n", i);
        auto nameballfile = std::string("ball.") + int2string(i, digits) + ".subdivided.stl";
        std::cout << "from file " << nameballfile << std::endl;

        if (std::filesystem::exists(nameballfile)) {
            std::cout << "File exists! Reading file" << std::endl;
        } else {
            std::cout << "File not found." << std::endl;
            throw std::runtime_error("Required file not found: " + nameballfile);
        }
                
        auto trimesh = chrono::ChTriangleMeshConnected::CreateFromSTLFile(nameballfile);
        double scale_factor = 0.001;

        // 3. Apply the transformation "The Clean Way"
        // Parameters: (Translation vector, Rotation/Scaling matrix)
        trimesh->Transform(chrono::ChVector3d(0, 0, 0), chrono::ChMatrix33<>(scale_factor));
        //chrono::ChVector3d min_v, max_v;
        //trimesh->GetBoundingBox(min_v, max_v);

        auto aabb = trimesh->GetBoundingBox(); 
        // std::cout << std::format("Mesh Size: X={:.3e}, Y={:.3e}, Z={:.3e} meters\n", 
        //      max_v.x() - min_v.x(), 
        //      max_v.y() - min_v.y(), 
        //      max_v.z() - min_v.z());
        std::cout << "Ball size Min: " << aabb.min.x() << " Max: " << aabb.max.x() << std::endl;
        //trimesh is a ChTriangleMeshConnected object

        // auto trimesh = chrono::ChTriangleMeshConnected::CreateFromSTLFile(
        //     std::string("ball.") + int2string(i, digits) + ".subdivided.stl");

        double vol;
        chrono::ChVector3d cog; // cog center of gravity?
        chrono::ChMatrix33<> inertia;

        // ComputeMassProperties gets density volume center of gravity and innertia
        // here True in first parameter means full body
        // and not only surface.
        // oputput is per unit of density
        // that is everything must be multiplied by density afterwards
        trimesh->ComputeMassProperties(true, vol, cog, inertia);
        printf("Volume of object %d is %e", i, vol);
        // move the center of the ball to the center of gravity 
        // and apply there the inertia
        // chrono::ChMatrix33<>(1) means do not rotate the ball
        trimesh->Transform(-cog, chrono::ChMatrix33<>(1));

        printf("Ball %d readed. Inertia calculated.\n", i);

        // Setup the rock body
        auto rock = chrono_types::make_shared<chrono::ChBody>();
        _rocks.push_back(rock); //method of std::vector
        _sys->Add(rock); // Add to simulation

        rock->SetName(std::string("ball.") + std::to_string(i)); // why not int2string local method?

        // note multiplication of previously per unit of density magnitudes 
        rock->SetMass(vol * ROCK_DENSITY);        
        rock->SetInertia(inertia * ROCK_DENSITY);

        // Setup the collision model
        // The rock is somewhere, so we are displacing the mesh to the origin,
        // and then we are moving the body to the COG
        auto coll_model =
            chrono_types::make_shared<chrono::ChCollisionModel>();
        coll_model->SetSafeMargin(0.0001); // margin to avoid fusion of objects
        coll_model->SetEnvelope(ENVELOPE_SIZE);

        // material of the objects
        auto rock_mat =
            chrono_types::make_shared<chrono::ChContactMaterialNSC>(); // MAterial uses impulses & Coulomb friction
        rock_mat->SetFriction(0.5);
        rock_mat->SetDampingF(0.1);

        // Triangles for the collision model of chrono
        // paramters are
        // Defined material, STL mess, is_static, is_convex, 
        // radius thinkness of the skin of the triangles
        auto coll_shape =
            chrono_types::make_shared<chrono::ChCollisionShapeTriangleMesh>(
                rock_mat, trimesh, false, true, 0.0001f);

        // ChFrame<> The ChFrame defines the relative position and orientation
        // of the shape with respect to the body's Center of Gravity (COG).
        // chrono::ChVector3d(0, 0, 0): This means the center of the STL mesh coincides with the body's COG.
        // This only works correctly because of the previously performed trimesh->Transform(-cog, ...)
        // chrono::QUNIT: This is the "Identity Quaternion." It represents zero rotation. The mesh will be oriented exactly as it was designed in the STL file.
        auto coll_frame = chrono::ChFrame<>(chrono::ChVector3d(0, 0, 0), chrono::QUNIT);

        // AddShape(shape, frame)
        coll_model->AddShape(coll_shape, coll_frame);
        
        rock->AddCollisionModel(coll_model); // outdated
        rock->EnableCollision(true); // enable collison model
        rock->SetPos(cog); // moving the ball back to its cog

        // Add the forces
        /*auto bouyancy = chrono_types::make_shared<chrono::ChForce>();
        rock->AddForce(bouyancy);
        bouyancy->SetMode(chrono::ChForce::FORCE);
        bouyancy->SetFrame(chrono::ChForce::BODY);
        bouyancy->SetAlign(chrono::ChForce::WORLD_DIR);
        bouyancy->SetVrelpoint(chrono::ChVector3d(0, 0, 0));
        bouyancy->SetDir(chrono::ChVector3d(0, 0, 1));
        bouyancy->SetMforce(vol * rho * 9.81);*/
        
        // Add the force
        auto force = chrono_types::make_shared<chrono::ChForce>();
        rock->AddForce(force);
        _forces.push_back(force);
        force->SetMode(chrono::ChForce::FORCE); // force and not torque
        force->SetFrame(chrono::ChForce::BODY); // The force is attached to the ball's "nose." 
        //If the ball flips, the force flips with it (like a jet engine). not understand what this means
        // To me ChForce::AlignmentFrame::WORLD (force always points "North" in the map, regardless of how the ball spins)
        //looks more natural
        // position of the forces set via SetVrelpoint() is interpreted in Local Coordinates of the ball.
        // IT is set at the point (0,0,0), the force is "attached" to the center of the ball. 
        // Even as the ball flies across the map, the force stays perfectly centered on the ball's mass.
        // World, force is applied to the same position in the space        
        force->SetAlign(chrono::ChForce::WORLD_DIR);
        // This defines the Direction where the vector points. The direction of the force is fixed relative to the Inertial World Frame (the X, Y, Z axes of the universe).
        // Even if the ball starts tumbling or spinning at high speeds after the impact, the force will always point in the same direction (e.g., always pushing "East").
        force->SetVrelpoint(chrono::ChVector3d(0, 0, 0)); //see previous comment


        //set the torque
        auto torque = chrono_types::make_shared<chrono::ChForce>();
        rock->AddForce(torque);
        _torques.push_back(torque);
        torque->SetMode(chrono::ChForce::TORQUE);
        torque->SetFrame(chrono::ChForce::BODY);
        torque->SetAlign(chrono::ChForce::WORLD_DIR);
        torque->SetVrelpoint(chrono::ChVector3d(0, 0, 0));

    } //loop on particles finish here

    
    // _sys->SetTimestepperType(chrono::ChTimestepper::Type::EULER_IMPLICIT);
    _sys->Setup();//final step before one enters the simulation loop

    std::vector<std::string> indeps({"dt"}), outdeps; // creates two
    // vectors of strings, indeps and outdeps, indeps initialized with 
    // list {"dt"} and outdeps non initalized

    for (unsigned int i=0; i < n_solids; i++) {
        
        // adds members to the vector, 
        // indeps independent variables ?
        indeps.push_back(
            std::string("ball_") + int2string(i, digits) + "_Force_p");
        indeps.push_back(
            std::string("ball_") + int2string(i, digits) + "_Moment_p");

        // add members to the vector outdeps
        // outdeps 
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
    //create a vec4 vector, populate it
    vec4 v;
    v.x = value.x();
    v.y = value.y();
    v.z = value.z();
    v.w = 0.f;

    //set the aqua variable
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
 * I do not know why this is ok
 * I do not know why not a nullptr is used
 */

cl_event
RocksSim::_execute(const std::vector<cl_event> UNUSED_PARAM events)
{
    // arcane type of aqua
    auto vars = CalcServer::singleton()->variables();
    
    //get the timestep, convert to float* deindiriction
    float dt = *((float*)vars->get("dt")->get(true));

    // Apply the forces to the rocks
    unsigned int n_solids = _rocks.size(); //amount of moving objects
    const unsigned int digits = num_digits(n_solids); // format
    for (unsigned int i = 0; i < n_solids; i++) {

        //get force and momentum from aqua
        // forces and moments are know at this stage
        // that is the singleton object know them
        std::string prefix = std::string("ball_") + int2string(i, digits);
        const vec4 F = *((vec4*)vars->get(prefix + "_Force_p")->get(true));
        const vec4 M = *((vec4*)vars->get(prefix + "_Moment_p")->get(true));
        // printf("x force of ball %e\n", F.x);
        // printf("y force of ball %e\n", F.y);
        // printf("z force of ball %e\n", F.z);
        // printf("w force of ball %e\n", F.w);
        //set the forces for Chrono
        setForce(_forces[i], F);
        setForce(_torques[i], M);
    }

    // Compute the dynamics
    _sys->DoStepDynamics(dt);

//    auto debval = _rocks[0]->GetLinVel().x();
//    auto debval = _rocks[0]->GetLinVel().Length();
//    printf("x velocity of ball %e\n", debval);

    // Get the new positions and angles
    for (unsigned int i = 0; i < n_solids; i++) {
        std::string prefix = std::string("ball_") + int2string(i, digits);
        
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
    return nullptr;
}

}}  // namespaces
 
