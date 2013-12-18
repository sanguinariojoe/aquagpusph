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

#ifndef CALCSERVER_H_INCLUDED
#define CALCSERVER_H_INCLUDED

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include <CL/cl.h>

#include <sphPrerequisites.h>
#include <CalcServer/Predictor.h>
#include <CalcServer/Grid.h>
#include <CalcServer/LinkList.h>
#include <CalcServer/RadixSort.h>
#include <CalcServer/Rates.h>
#include <CalcServer/Boundary/ElasticBounce.h>
#include <CalcServer/Boundary/DeLeffe.h>
#include <CalcServer/Boundary/GhostParticles.h>
#include <CalcServer/Shepard.h>
#include <CalcServer/Corrector.h>
#include <CalcServer/TimeStep.h>
#include <CalcServer/DensityInterpolation.h>
#include <CalcServer/Torque.h>
#include <CalcServer/Domain.h>
#include <CalcServer/Movements/Movement.h>
#include <CalcServer/Movements/LIQuaternion.h>
#include <CalcServer/Movements/C1Quaternion.h>
#include <CalcServer/Movements/ScriptQuaternion.h>
#include <CalcServer/Sensors.h>
#include <CalcServer/Energy.h>
#include <CalcServer/Bounds.h>
#include <CalcServer/Portal/Portal.h>
#include <Singleton.h>

namespace Aqua{
/** @namespace CalcServer Calculation server name space.
 */
namespace CalcServer{

/** @class CalcServer CalcServer.h CalcServer.h
 * @brief Entity that perform the main work of the simulation. Therefore this
 * class have an internal loop which iterates while no data must be sent to
 * the host for an output file writing process.
 * @remarks Some output files are managed internally by this class, like the
 * log file, the energy file, pressure sensors file, or bounds file.
 */
class CalcServer : public Aqua::Singleton<Aqua::CalcServer::CalcServer>
{
public:
	/** Constructor
	 */
	CalcServer();

	/** Destructor
	 */
	~CalcServer();

	/** Internal server loop. Calculation server will be iterating while the
	 * next output file event is reached (or the simulation is finished).
	 * @return false if all gone right, true otherwise.
	 */
	bool update();

	/** Transfer the data from the computational device to the host.
	 * @param dest Array where the data should be copied.
	 * @param orig Computational device allocated data to copy.
	 * @param size Size of the data to copy.
	 * @param offset Offset into the array to start reading.
	 * @return false if the data has been succesfully copied, true otherwise.
	 */
	bool getData(void *dest, cl_mem orig, size_t size, size_t offset=0);

	/** Send data to a computational device array.
	 * @param dest Identifier of the destination in the server.
	 * @param orig Array of data to read.
	 * @param size Size of the data to copy.
	 * @return false if the data has been succesfully copied, true otherwise.
	 */
	bool sendData(cl_mem dest, void* orig, size_t size);

	/** Allocate memory in the server.
	 * @param size Amount of memory in bytes to allocate.
	 * @return The new memory object, NULL if the memory could not been
	 * allocated.
	 */
	cl_mem allocMemory(size_t size);

	/** Fills the fluid using an OpenCL script
	 * @return false if the fluid has been succesfully set, true otherwise
	 * @warning The OpenCL script initialization way is deprecated.
	 */
	bool fillFluid();

	/** Setup the calculation server with the data recopilated by the host during
	 * the simulation setup process.
	 * @return false if the caluclation server has been succesfully setup, true otherwise
	 */
	bool setup();

	/** Prints the log file.
	 */
	void printLog();

	/** Prints the energy report. The energy components are:
     *     -# Internal energy: \f$ E = E_pot + E_kin \f$
     *     -# Potential energy: \f$ E_pot = m \frac{p}{\rho} \f$
     *     -# Kinetic energy: \f$ E_kin = \frac{1}{2} m \vert v \vert^2 \f$
	 */
	void printEnergy();

	/** Prints the bounds report, including the bounds of the fluid particles
	 * the maximum and minimum velocities.
	 */
	void printBounds();

	/** Perform the energy calculation.
	 */
	void energy();

	/** Perform the bounds calculation.
	 */
	void bounds();

	/// Verbose level
	int verbose_level;

	/// Device type
	cl_device_type clDeviceType;
	/// Number of available platforms
	cl_uint clNPlatforms;
	/// Array of platforms
	cl_platform_id *clPlatforms;
	/// Number of devices
	cl_uint clNDevices;
	/// Array of devices
	cl_device_id *clDevices;
	/// OpenCL context
	cl_context clContext;
	/// OpenCL command queue
	cl_command_queue *clComQueues;
	/// Selected platform
	cl_platform_id clPlatform;
	/// Selected device
	cl_device_id clDevice;
	/// Selected command queue
	cl_command_queue clComQueue;

	/// Memory allocated in bytes
	size_t AllocatedMem;
	/// Total number of particles (including sensors)
	unsigned int N;
	/// Number of fluid particles (including boundary areas and fixed particles)
	unsigned int n;
	/// Number of sensors
	unsigned int nSensors;
	/// Number of fluids
	unsigned int nfluid;
	/// Total number of cells allocated
	unsigned int lxydim;
	/// Number of cells in each direction
	uivec lvec;
	/// Total number of cells
	unsigned int lxy;
	/// Length of the cells
	float rdist;
	/// Minimum particle position (fluid or fixed, it is used for the Link-List definition)
	vec posmin;
	/// Maximum particle position (fluid or fixed, it is used for the Link-List definition)
	vec posmax;

	/// Maximum kernel height value, to control the linklist grid
	float h;
	/// Gravity force
	vec g;
	/// kernel height and particles distance ratio.
	float hfac;
	/// Time step divisor (The inverse of the Courant number).
	float dt_divisor;
	/** LinkList time steps validity. If this value is grater than 1, the link-list process
	 * will be avoided sometimes, but the cells size should be increased to ensure that the
	 * neighbours list still being valid, increasing the number of neighbours per particles.
	 * Sometimes some performance can be gained increasing this value.
	 */
	unsigned int link_list_steps;
	/// LinkList main step, to control if a new link-list is required.
	unsigned int mLLStep;
	/// Cell size increasing factor (resulting from link_list_steps, deltar, hfac and dt_divisor).
	float CellFac;
	/** Density field interpolation steps. The density field is usually computed as an
	 * evolution process fulfilling the continuity equation, but in order to reduce the
	 * pressure noise sometimes it can be interpolated through the SPH model, being
	 * therefore a geometric result.
	 * It is not recommended to use this trick.
	 */
	unsigned int dens_int_steps;
	/// Density interpolation main step, to control if the density field should be interpolated.
	unsigned int mDensStep;

	/// Time step
	float dt;
	/// Sound speed
	float cs;

	/** Fixed/fluid particles flag.
	 *   - > 0 for every fluid. \n
	 *   - = 0 for sensors. \n
	 *   - < 0 for fixed particles or boundary elements.
	 */
	cl_mem imove;
	/// Fixed/fluid particles flag (predictor-corrector backup variable).
	cl_mem imovein;
	/// Fluid identifier
	cl_mem ifluid;
	/// Fluid identifier (predictor-corrector backup variable)
	cl_mem ifluidin;
	/// Position
	cl_mem pos;
	/// Normal, for boundary particles/elements
	cl_mem normal;
	/// Velocity
	cl_mem v;
	/// Acceleration
	cl_mem f;
	/// Density
	cl_mem dens;
	/// Density change rate
	cl_mem drdt;
	/** Density change rate (restricted to the numerical diffusive term of
     * \f$\delta\f$-SPH technique)
     */
	cl_mem drdt_F;
	/// Mass
	cl_mem mass;
	/// Kernel height
	cl_mem hp;
	/// Pressure
	cl_mem press;
	/// Position (predictor-corrector backup variable)
	cl_mem posin;
	/// Normal (sorted-unsorted backup variable)
	cl_mem normalin;
	/// Velocity (predictor-corrector backup variable)
	cl_mem vin;
	/// Acceleration (predictor-corrector backup variable)
	cl_mem fin;
	/// Density (predictor-corrector backup variable)
	cl_mem densin;
	/// Density change rate (predictor-corrector backup variable)
	cl_mem drdtin;
	/// Mass (predictor-corrector backup variable)
	cl_mem massin;
	/// Kernel height (predictor-corrector backup variable)
	cl_mem hpin;
	/// Pressure (predictor-corrector backup variable)
	cl_mem pressin;
	/// Viscous timestep term
	cl_mem sigma;
	/// Convective timestep term
	cl_mem dtconv;
	/// Shepard term (0th correction)
	cl_mem shepard;
	/// Shepard term gradient (0th correction)
	cl_mem gradShepard;
	/** Permutations vector, that gives for each particle index in
	 * the sorted space, their index in the unsorted space.
	 * @note See also lcell variable.
	 */
	cl_mem permutation;
	/** Inversed permutations vector, that gives for each particle
	 * index in the unsorted space, their index in the sorted space.
	 * @note See also lcell variable.
	 */
	cl_mem reversePermutation;
	/** Head of chain for each cell in the sorted space, meaning the
	 * first particle in each cell.
	 * @note See also lcell variable.
	 */
	cl_mem ihoc;
	/** Mark about cells which contains at least one fluid particle.
	 * @note See also lcell variable.
	 */
	cl_mem isValidCell;
	/** Cell where each particle is placed, in the sorted space.
	 * @note To improve the use of multicore platform the particles
	 * are sorted by their cell position at the linklist phase. Therefore
	 * two identifier spaces must be considered, the original one, and
	 * the sorted one.
	 */
	cl_mem lcell;

	/** lcell vector dimension, that must be power of two in order to
	 * use the radix sort process.
	 * This dimension is applied to the permutation vectors as well.
	 */
	unsigned int nLcell;
	/// EOS gamma value for each fluid
	cl_mem gamma;
	/// Reference density for each fluid
	cl_mem refd;
	/// Dynamic viscosity for each fluid
	cl_mem visc_dyn;
	/// Kinetic viscosity for each fluid
	cl_mem visc_kin;
	/// Alpha corrected dynamic viscosity for each fluid.
	cl_mem visc_dyn_corrected;
	/// Continuity equation diffusive term multiplier.
	cl_mem delta;
	/// Minimum value of the time step (it is a reduction process result)
	cl_mem DT;

	/// Sensors type for each sensor
	cl_mem sensorMode;

	/** Platform identifiers:
	 *   - -1 = Any valid platform (Generic one)
	 *   - 0 = nVIDIA
	 *   - 1 = AMD
	 *   - 2 = Intel
	 */
	int PlatformID;

	/** Internal energy: \f$ U = \int_0^t sum_i \frac{p_i}{\rho_i^2}
	 *   \left(
     *      \frac{\mathrm{d} \rho_i}{\mathrm{d} t}
     *      - \left. \frac{\mathrm{d} \rho_i}{\mathrm{d} t} \right\vert_F
	 *   \right) m_i \mathrm{d}t \f$.
	 * @warning The viscous dissipation is not implemented yet.
	 */
	float eint;
    /** Kinetic energy: \f$ E_{kin} = \sum_i \frac{1}{2} m_i
	 *   \vert \mathbf{u}_i \vert^2 \f$.
	 * @warning The viscous dissipation is not implemented yet.
	 */
	float ekin;
	/** Total energy: \f$ E = U + E_{kin} \f$
	 * @warning The viscous dissipation is not implemented yet.
	 */
	float etot;
	/// Minimum X,Y,Z coordinates for the fluid
	vec minCoords;
	/// Maximum X,Y,Z coordinates for the fluid
	vec maxCoords;
	/// Minimum velocity for the fluid
	float minV;
	/// Maximum velocity for the fluid
	float maxV;
	/// true if the energy data has been calculated, false otherwise
	bool energyPerformed;
	/// true if the bounds data has been calculated, false otherwise
	bool boundsPerformed;
	/// Total fluid mass
	float fluidMass;

	// --------------------------------------------
	// Kernels
	// --------------------------------------------
	/// Predictor stage.
	Predictor *mPredictor;
	/// Grid, that determine the number of cells at each direction.
	Grid *mGrid;
	/// Link-list that allocate each particle in a cell to know the neighbours list.
	LinkList *mLinkList;
	/// Rates stage, where the SPH interactions are performed.
	Rates *mRates;
	/// ElasticBounce boundary condition, used only if it is selected.
	Boundary::ElasticBounce *mElasticBounce;
	/// DeLeffe boundary condition, used only if it is selected.
	Boundary::DeLeffe *mDeLeffe;
	/// Ghost particles, used only if it is selected.
	Boundary::GhostParticles *mGhost;
	/// 0th order correction stage, used only if it is selected.
	Shepard *mShepard;
	/// Corrector stage.
	Corrector *mCorrector;
	/// Domain bounds test.
	Domain *mDomain;
	/// Time-step calculation for the next iteration.
	TimeStep *mTimeStep;
	/// Density reinitialization process.
	DensityInterpolation *mDensInt;
	/// Motions array.
	std::deque<Movement::Movement*> mMoves;
	/// Sensors measurement.
	Sensors *mSensors;
	/// Energy computation.
	Energy *mEnergy;
	/// Bounds computation.
	Bounds *mBounds;
	/// Portals array.
	std::deque<Portal::Portal*> mPortals;
private:
	/** Setup the OpenCL stuff.
	 * @return false if the OpenCL environment has been succesfully built, true otherwise
	 */
	bool setupOpenCL();
	/** Prints all the available platforms and devices returned by OpenCL.
	 * @return false if the OpenCL environment can be succesfully built, true otherwise
	 */
	bool queryOpenCL();
	/** Get a platform from the available ones.
	 * @return false if a platform could be obtained, true otherwise
	 */
	bool getPlatform();
	/** Get the available devices in the selected platform.
	 * @return false if the devices have been succesfully obtained, true otherwise
	 */
	bool getDevices();

};

}}  // namespace

#endif // CALCSERVER_H_INCLUDED
