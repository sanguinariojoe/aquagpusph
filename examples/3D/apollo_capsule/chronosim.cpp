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

/**
 * @brief C API entry for AQUAgpusph.
 * AQUAgpusph will call this function to receive an 
 * Aqua::CalcServer::Tool derived object.
 * The later, in turn, 
 * will be afterwards considered by any other tool.
 * 
 */

extern "C" Aqua::CalcServer::ApolloSim* create_object(
    const std::string name, bool once)
{
    return new Aqua::CalcServer::ApolloSim(name, once);
}

namespace Aqua{ namespace CalcServer{

// initiates a tool element. Aqua inner 
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
    printf("Chronosim: Starting the setup\n");
    
    //Tool is Aqua 
    Tool::setup();

    // Get the configuration variables
    auto vars = CalcServer::singleton()->variables();

    _pitch = *((float*)vars->get("pitch")->get(true));
    _vel = *((float*)vars->get("u0")->get(true));
    _cogz = *((float*)vars->get("cogz")->get(true)); 
    _pitch *= std::numbers::pi / 180.0; 

    // Setup the chrono system
    // ChSystemNSC Non-Smooth Contact 
    // Rigid collision and friction
    _sys = chrono_types::make_shared<chrono::ChSystemNSC>();

    //Create a generic new body, without properties
    _apollo = chrono_types::make_shared<chrono::ChBody>();

    //Add or Addbody is almost the same
    _sys->AddBody(_apollo);

    // Add a generic force and momentum
    _force = chrono_types::make_shared<chrono::ChForce>();
    _moment = chrono_types::make_shared<chrono::ChForce>();
    // Add them to the body
    _apollo->AddForce(_force);
    _apollo->AddForce(_moment);

    // add an string identifiyer to the apollo
    _apollo->SetName("Apollo");

    //Add gravity in z direction
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

    //mass direcly introduced as known given value
    _apollo->SetMass(3900);

    // Set the diagonal moments of innertia $I_{xx}$, $I_{yy}$, $I_{zz}$
    _apollo->SetInertiaXX(chrono::ChVector3d(5560, 5270, 4180));

    //define the parameters of the force

    // force and not torque
    _force->SetMode(chrono::ChForce::FORCE);
    // Even as the body flies across the map, the force stays perfectly centered on the body's mass.
    _force->SetFrame(chrono::ChForce::BODY);
    // This defines the Direction where the vector points. 
    // The direction of the force is fixed relative to the 
    // Inertial World Frame (the X, Y, Z axes of the universe).
    // Even if the ball starts tumbling or spinning at high speeds 
    // after the impact, the force will always point in the same 
    // direction (e.g., always pushing "East").
    _force->SetAlign(chrono::ChForce::WORLD_DIR);
    
    // It is set at the point (0,0,0), the force is "attached" to the center of the body.
    _force->SetVrelpoint(chrono::ChVector3d(0, 0, 0));
    
    
    //define the paramter of the momentum
    // torque and not force
    _moment->SetMode(chrono::ChForce::TORQUE);
    // see above
    _moment->SetFrame(chrono::ChForce::BODY);
    // see above
    _moment->SetAlign(chrono::ChForce::WORLD_DIR);
    // see above
    _moment->SetVrelpoint(chrono::ChVector3d(0, 0, 0));

    // moving the body back to its cog? maybe in z direction?
    _apollo->SetPos(chrono::ChVector3d(0, 0, _cogz));
    // initial velocity
    _apollo->SetLinVel(chrono::ChVector3d(0, 0, -_vel));

    // On NWU coordinates:
    //    x : Positive moment = positive roll = portside goes up
    //    y : Positive moment = positive pitch = bow goes down
    //    z : Positive moment = positive yaw = bow goes to the portside

    //create a quaternion
    chrono::ChQuaternion<double> R;
    // cardan angles known, set them inside
    R.SetFromCardanAnglesXYZ(chrono::ChVector3d(0, _pitch, 0));
    // give initial turn angles
    _apollo->SetRot(R);

    _sys->SetTimestepperType(chrono::ChTimestepper::Type::EULER_EXPLICIT);//define time integrator

    _sys->Setup();//final step before one enters the simulation loop

    // custom co-simulation wrappers
    // setting which variables must be updated or satisfied
    // before the next calculation step can proceed

    //maybe defined by the use of midpoint integrator
    //in aqua
    setInputDependencies({"dt", "iter_midpoint", "iter_midpoint_max",
                          "Force_p_iset", "Moment_p_iset"});
    
    // defines what the system writes or provides to other 
    // modules after a calculation is finished.
    setOutputDependencies({"motion_r", "motion_drdt", "motion_ddrddt",
                           "motion_a", "motion_dadt", "motion_ddaddt",
                           "forces_r"});

};//setup finalizes here


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
 * get force and momentum
 * calculate movement in the time step dt
 * actualize position, velocity, etc...
 * of each object
 * @param
 * @return This munction returns a cl_event that 
 * is fixed to NULL
 */
cl_event
ApolloSim::_execute(const std::vector<cl_event> UNUSED_PARAM events)
{   
    // type of aqua
    auto vars = CalcServer::singleton()->variables();

    // Check whether we are on the midpoint, or at the final iteration
    const unsigned int iter =
        *((unsigned int*)vars->get("iter_midpoint")->get(true));
    const unsigned int iter_max =
        *((unsigned int*)vars->get("iter_midpoint_max")->get(true));
    const bool is_midpoint = iter < iter_max;

    // Get the forces from AQUAgpusph
    // get the timestep, convert to float* deindiriction
    float dt = *((float*)vars->get("dt")->get(true));

    // get force and momentum from aqua
    // forces and moments are know at this stage
    // that is the singleton object know them
    const vec4 F = *((vec4*)vars->get("Force_p_iset")->get(true));
    const vec4 M = *((vec4*)vars->get("Moment_p_iset")->get(true));

    if (iter == 0) {
        // At the beggining of the time step we must copy the results from
        // the other instance of this solver 
        vec4 data;
        chrono::ChQuaternion<double> R;

        data = *((vec4*)vars->get("motion_r")->get(true));
        _apollo->SetPos(chrono::ChVector3d(data.x, data.y, data.z));
        
        data = *((vec4*)vars->get("motion_drdt")->get(true));
        _apollo->SetLinVel(chrono::ChVector3d(data.x, data.y, data.z));//modern Chrono SetPos_dt
        
        data = *((vec4*)vars->get("motion_ddrddt")->get(true));
        _apollo->SetLinAcc(chrono::ChVector3d(data.x, data.y, data.z));//modern Chrono SetPos_dt
        
        data = *((vec4*)vars->get("motion_a")->get(true));
        R.SetFromCardanAnglesXYZ(chrono::ChVector3d(data.x, data.y, data.z));
        _apollo->SetRot(R);

        data = *((vec4*)vars->get("motion_dadt")->get(true));
        _apollo->SetAngVelLocal(chrono::ChVector3d(data.x, data.y, data.z));
        
        data = *((vec4*)vars->get("motion_ddaddt")->get(true));
        _apollo->SetAngAccLocal(chrono::ChVector3d(data.x, data.y, data.z));
    }

    // Book-keeping, so we can restore the state within the midpoint iterations
    // This is actually needed also at the final iterator because of the Euler
    // explicit (see below)
    double T = _sys->GetChTime();//current simulation time

    // GetNumCoordsPosLevel: integer representing the total number of 
    // coordinates (degrees of freedom) at the position level currently in your simulation
    // _sys.get() returns the raw pointer to the memory address where the system lives
    // chrono::ChState initializing a chrono::ChState object. 
    // This is a specialized container used by Chrono to store the entire "Position Level" 
    // state of the system—the x,y,z coordinates and quaternions for every moving body
    chrono::ChState X(_sys->GetNumCoordsPosLevel(), _sys.get());
    
    //store the velocities and accelerations of your system.
    chrono::ChStateDelta V(_sys->GetNumCoordsVelLevel(), _sys.get());
    chrono::ChStateDelta A(_sys->GetNumCoordsVelLevel(), _sys.get());

    // _sys->GetNumConstraints() returns an integer representing the 
    // total number of scalar constraint equations currently active in your system 
    // L(sys->GetNumConstraints()), creates a container specifically
    // designed to hold the Lagrange Multipliers (λ) for simulation.
    // applied to constrains
    chrono::ChVectorDynamic<> L(_sys->GetNumConstraints());

    //get vels and accelerations
    const chrono::ChVector3d drdt0 = _apollo->GetLinVel();
    const chrono::ChVector3d ddrddt0 = _apollo->GetLinAcc();
    const chrono::ChVector3d dadt0 = _apollo->GetAngVelLocal();
    const chrono::ChVector3d ddaddt0 = _apollo->GetAngAccLocal();

    // "synchronize" mathematical containers with the current physical state of the simulation.
    // It pulls the raw data out of the ChBody objects and packs them into the 
    // ChState and ChStateDelta vectors initialized earlier (X and V).
    _sys->StateGather(X, V, T);
    _sys->StateGatherAcceleration(A);//acceleration 6 dimensions
    _sys->StateGatherReactions(L);//constrains 

    // Compute the dynamics
    // introduce in chrono the soliciations
    setForce(_force, F);
    setForce(_moment, M);

    //calculate
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
    chrono::ChVector3d r = _apollo->GetPos();
    const chrono::ChVector3d drdt = _apollo->GetLinVel();
    const chrono::ChVector3d ddrddt = _apollo->GetLinAcc();
    const chrono::ChVector3d a = _apollo->GetRot().GetCardanAnglesXYZ();
    const chrono::ChVector3d dadt = _apollo->GetAngVelLocal();
    const chrono::ChVector3d ddaddt = _apollo->GetAngAccLocal();

    // On the explicit Euler scheme the position is integrated directly from
    // the velocity at the beggining. We want to use the midpoint velocity
    // instead
    r = r + ddrddt * (0.5 * dt * dt);
    _apollo->SetPos(r);

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
        // here the new possition, vel, etc is got from Chrono
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
    return nullptr;
}

}}  // namespaces
 
