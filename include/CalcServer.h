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
	 * @warning This initialization way is deprecated and it will be removed in
	 * future releases.
	 */
	bool fillFluid();

	/** Setup the calculation server with the data recopilated by the host
	 * during the initialization process.
	 * @return false if the caluclation server has been succesfully setup,
     * true otherwise
	 */
	bool setup();

	/** Updates the log file.
	 */
	void printLog();

	/** Updates the energy report. The energy components printed are:
     *   -# Potential energy: \f$ E_{pot} = - \sum_i m_i
	     \mathbf{g} \cdot \mathbf{r}_i \f$.
     *   -# Kinetic energy: \f$ E_{kin} = \sum_i \frac{1}{2} m_i
	     \vert \mathbf{u}_i \vert^2 \f$.
	 *   -# Internal energy: \f$ U = \int_0^t sum_i \frac{p_i}{\rho_i^2}
	     \left(
            \frac{\mathrm{d} \rho_i}{\mathrm{d} t}
            - \left. \frac{\mathrm{d} \rho_i}{\mathrm{d} t} \right\vert_F
	     \right) m_i \mathrm{d}t \f$.
	 *   -# Enthalpy: \f$ H = \int_0^t sum_i \frac{p_i}{\rho_i^2}
         \frac{\mathrm{d} \rho_i}{\mathrm{d} t} m_i \mathrm{d}t \f$.
     *   -# Entropy: \f$ TS = U - H \f$.
     *   -# Total energy: \f$ E = U + E_{kin} \f$.
	 */
	void printEnergy();

	/** Updates the bounds report, where the bounding box of the fluid
	 * particles and the maximum and minimum velocities is printed.
	 */
	void printBounds();

	/** Performs the energy calculation. When this method is called a flag
	 * is set in order to don't recompute it until the next time step.
	 */
	void energy();

	/** Perform the bounds calculation. When this method is called a flag
	 * is set in order to don't recompute it until the next time step.
	 */
	void bounds();

	/// Verbose level
	int verbose_level;

	/// Number of available platforms
	cl_uint num_platforms;
	/// Array of platforms
	cl_platform_id *platforms;
	/// Number of devices
	cl_uint num_devices;
	/// Array of devices
	cl_device_id *devices;
	/// OpenCL context
	cl_context context;
	/// OpenCL command queue
	cl_command_queue *command_queues;
	/// Selected platform
	cl_platform_id platform;
	/// Selected device
	cl_device_id device;
	/// Selected command queue
	cl_command_queue command_queue;

	/// Memory allocated in bytes
	size_t allocated_mem;
	/// Total number of particles (including sensors)
	unsigned int N;
	/// Number of fluid particles (including boundary areas and fixed particles)
	unsigned int n;
	/// Number of sensors
	unsigned int num_sensors;

	/// Number of fluids
	unsigned int num_fluids;
	/// Total number of cells allocated
	unsigned int num_cells_allocated;

	/** Minimum particle position (fluid and fixed particles are taken into
     * account)
     */
	vec pos_min;
	/** Maximum particle position (fluid and fixed particles are taken into
     * account)
     */
	vec pos_max;

	/// Number of cells in each direction
	uivec num_cells_vec;
	/// Total number of cells
	unsigned int num_cells;
	/// Length of the cells
	float cell_length;
	/// Maximum kernel height value, to control the linklist grid
	float h;

	/// Gravity force \f$ \mathbf{g} \f$
	vec g;
	/// Kernel height factor \f$ \frac{h}{\Delta r} \f$.
	float hfac;
	/// Inverse of the Courant factor \f$ \frac{1}{Cf} \f$.
	float dt_divisor;
	/** LinkList time steps validity. If this value is grater than 1, the
	 * link-list computation process will be avoided some time steps, but the
	 * cells size should be increased to ensure that the neighbours list still
	 * being valid, increasing the number of neighbours per particles.
	 * Sometimes some performance can be gained increasing this value.
	 */
	unsigned int link_list_steps;
	/// LinkList main step, to control if a new link-list is required.
	unsigned int link_list_step;
	/** Cell size increasing factor (resulting from link_list_steps, deltar,
     * hfac and dt_divisor).
     */
	float cell_length_factor;
	/** Density field interpolation steps. The density field is usually
	 * computed as an evolution process fulfilling the continuity equation, but
	 * in order to reduce the pressure noise sometimes it can be interpolated
	 * through the SPH model, being therefore a geometrical result.
	 * It is strongly not recommended to apply this technique.
	 */
	unsigned int dens_int_steps;
	/** Density interpolation main step, to control if the density field should
	 * be interpolated.
	 */
	unsigned int dens_int_step;

	/// Time step \f$ \Delta t \f$.
	float dt;
	/// Sound speed \f$ c_s \f$.
	float cs;

	/** Fixed/fluid particles flag.
	 *   - imove > 0 for every fluid.
	 *   - imove = 0 for sensors.
	 *   - imove < 0 for fixed particles or boundary elements.
	 */
	cl_mem imove;
	/// Fixed/fluid particles flag (predictor-corrector backup variable).
	cl_mem imovein;
	/// Fluid identifier.
	cl_mem ifluid;
	/// Fluid identifier (predictor-corrector backup variable).
	cl_mem ifluidin;
	/// Position \f$ \mathbf{r} \f$.
	cl_mem pos;
	/// Normal, for boundary particles/elements \f$ \mathbf{n} \f$.
	cl_mem normal;
	/// Velocity \f$ \mathbf{u} \f$.
	cl_mem v;
	/// Acceleration \f$ \frac{\mathrm{d}\mathbf{u}}{\mathrm{d}t} \f$.
	cl_mem f;
	/// Density \f$ \rho \f$.
	cl_mem dens;
	/// Density change rate \f$ \frac{\mathrm{d}\rho}{\mathrm{d}t} \f$.
	cl_mem drdt;
	/** Density change rate (restricted to the numerical diffusive term of
     * \f$\delta\f$-SPH technique)
     */
	cl_mem drdt_F;
	/// Mass \f$ m \f$.
	cl_mem mass;
	/// Kernel height \f$ h \f$.
	cl_mem hp;
	/// Pressure \f$ p \f$.
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
	/// Viscous timestep term \f$ \Delta t_{visc} \f$.
	cl_mem sigma;
	/// Convective timestep term \f$ \Delta t_{conv} \f$.
	cl_mem dtconv;
	/** Shepard term (0th correction) \f$ \gamma(\mathbf{x}) = \int_{Omega}
           W(\mathbf{y} - \mathbf{x})
           \mathrm{d}\mathbf{x} \f$.
     */
	cl_mem shepard;
	/** Shepard term gradient (0th correction) \f$ \gamma(\mathbf{x}) = \int_{Omega}
           \nabla W(\mathbf{y} - \mathbf{x})
           \mathrm{d}\mathbf{x} \f$.
     */
	cl_mem shepard_gradient;
	/** Permutations vector, that gives for each particle index in the sorted
	 * space, their index in the unsorted space.
	 */
	cl_mem permutation;
	/** Inversed permutations vector, that gives for each particle index in the
	 * unsorted space, their index in the sorted space.
	 */
	cl_mem permutation_inverse;
	/** Head of chain for each cell in the sorted space, meaning the
	 * first particle in each cell.
	 */
	cl_mem ihoc;
	/** Mark about cells which contains at least one fluid particle.
	 */
	cl_mem cell_has_particles;
	/** Cell where each particle is placed, in the sorted space.
	 * @note To improve the use of multicore platforms the particles are sorted
	 * by their cell position at the linklist phase. Therefore two identifier
	 * spaces must be considered, the original one (formerly unsorted space),
	 * and the sorted one.
	 */
	cl_mem icell;

	/** icell array dimension, that must be power of two in order to use the
	 * radix sort process. This dimension is applied to the permutation vectors
	 * as well.
	 */
	unsigned int num_icell;
	/** EOS gamma exponent \f$ \gamma \f$.
     * \f$ p = \frac{c_s^2 \rho_0}{\gamma}
		 \left(
            \left( \frac{\rho}{\rho_0} \right)^\gamma - 1
         \right) \f$.
     */
	cl_mem gamma;
	/// Density of reference \f$ \rho_0 \f$.
	cl_mem refd;
	/// Dynamic viscosity \f$ \mu \f$.
	cl_mem visc_dyn;
	/// Kinetic viscosity \f$ \nu \f$.
	cl_mem visc_kin;
	/// Alpha corrected dynamic viscosity.
	cl_mem visc_dyn_corrected;
	/// Continuity equation diffusive term factor \f$\delta\f$.
	cl_mem delta;
	/// Minimum value of the time step (it is a reduction process result)
	cl_mem DT;

	/// Sensors types
	cl_mem sensor_mode;

	/** Internal energy: \f$ U = \int_0^t sum_i \frac{p_i}{\rho_i^2}
	     \left(
            \frac{\mathrm{d} \rho_i}{\mathrm{d} t}
            - \left. \frac{\mathrm{d} \rho_i}{\mathrm{d} t} \right\vert_F
	     \right) m_i \mathrm{d}t \f$.
	 * @warning The viscous dissipation is not implemented yet.
	 */
	float eint;
    /** Kinetic energy: \f$ E_{kin} = \sum_i \frac{1}{2} m_i
	    \vert \mathbf{u}_i \vert^2 \f$.
	 * @warning The viscous dissipation is not implemented yet.
	 */
	float ekin;
	/** Total energy: \f$ E = U + E_{kin} \f$
	 * @warning The viscous dissipation is not implemented yet.
	 */
	float etot;
	/// Minimum X,Y,Z coordinates for the fluid
	vec min_fluid_bound;
	/// Maximum X,Y,Z coordinates for the fluid
	vec max_fluid_bound;
	/// Minimum velocity for the fluid
	float min_v;
	/// Maximum velocity for the fluid
	float max_v;
	/// true if the energy data has been calculated, false otherwise
	bool energy_computed;
	/// true if the bounds data has been calculated, false otherwise
	bool bounds_computed;
	/// Total fluid mass
	float fluid_mass;

	// --------------------------------------------
	// Kernels
	// --------------------------------------------
	/// Predictor stage.
	Predictor *predictor;
	/// Grid, that determine the number of cells at each direction.
	Grid *grid;
	/// Link-list that allocate each particle in a cell to know the neighbours list.
	LinkList *link_list;
	/// Rates stage, where the SPH interactions are performed.
	Rates *rates;
	/// ElasticBounce boundary condition, used only if it is selected.
	Boundary::ElasticBounce *elastic_bounce;
	/// DeLeffe boundary condition, used only if it is selected.
	Boundary::DeLeffe *de_Leffe;
	/// Ghost particles, used only if it is selected.
	Boundary::GhostParticles *ghost_particles;
	/// 0th order correction stage, used only if it is selected.
	Shepard *shepard_tool;
	/// Corrector stage.
	Corrector *corrector;
	/// Domain bounds test.
	Domain *domain;
	/// Time-step calculation for the next iteration.
	TimeStep *time_step;
	/// Density reinitialization process.
	DensityInterpolation *dens_int;
	/// Motions array.
	std::deque<Movement::Movement*> motions;
	/// Sensors measurement.
	Sensors *sensors;
	/// Energy computation.
	Energy *energy_tool;
	/// Bounds computation.
	Bounds *bounds_tool;
	/// Portals array.
	std::deque<Portal::Portal*> portals;
private:
	/** Setup the OpenCL stuff.
	 * @return false if the OpenCL environment has been succesfully built,
	 * true otherwise
	 */
	bool setupOpenCL();
	/** Prints all the available platforms and devices returned by OpenCL.
	 * @return false if the OpenCL environment can be succesfully built,
	 * true otherwise
	 */
	bool queryOpenCL();
	/** Get a platform from the available ones.
	 * @return false if a platform could be obtained, true otherwise
	 */
	bool getPlatform();
	/** Get the available devices in the selected platform.
	 * @return false if the devices have been succesfully obtained, true
	 * otherwise
	 */
	bool getDevices();

};

}}  // namespace

#endif // CALCSERVER_H_INCLUDED
