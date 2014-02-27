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

#ifndef PROBLEMSETUP_H_INCLUDED
#define PROBLEMSETUP_H_INCLUDED

/** \def __NO_OUTPUT_MODE__
 * No print output is selected.
 */
#define __NO_OUTPUT_MODE__ 0
/** \def __FPS_MODE__
 * FPS (frames per second) mode is selected
 * @remarks could not be combinable with IPF mode
 */
#define __FPS_MODE__ 1 << 0
/** \def __IPF_MODE__
 * IPF (iterations per frame) mode is selected
 * @remarks could not be combinable with FPS mode
 */
#define __IPF_MODE__ 1 << 1
/** \def __TIME_MODE__
 * Time control mode is selected
 * @remarks Is equal to __FPS_MODE__
 */
#define __TIME_MODE__ 1 << 0
/** \def __ITER_MODE__
 * Iterations control mode is selected
 * @remarks Is equal to __IPF_MODE__
 */
#define __ITER_MODE__ 1 << 1
/** \def __FRAME_MODE__
 * Frames control mode is selected
 * @remarks Is equal to __IPF_MODE__
 */
#define __FRAME_MODE__ 1 << 2

/** \def __H5Part__
 * H5Part output is selected.
 */
#define __H5Part__ 1 << 0
/** \def __VTK__
 * VTK output is selected.
 */
#define __VTK__ 1 << 1
/** \def __TECPLOT__
 * TECPLOT output is selected.
 */
#define __TECPLOT__ 1 << 2

/** \def __DT_VARIABLE__
 * Time step will be calculated every
 * iteration. This method guarantee
 * Courant number conservation.
 */
#define __DT_VARIABLE__ 0
/** \def __DT_FIXCALCULATED__
 * Time step will calculated at the start
 * of simulation, and conserved along it.
 */
#define __DT_FIXCALCULATED__ 1
/** \def __DT_FIX__
 * Time step will provided by user.
 */
#define __DT_FIX__ 2

// ----------------------------------------------------------------------------
// Include Prerequisites
// ----------------------------------------------------------------------------
#include <sphPrerequisites.h>

// ----------------------------------------------------------------------------
// Include standar libraries
// ----------------------------------------------------------------------------
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <vector>
#include <deque>

// ----------------------------------------------------------------------------
// Include OpenCL libraries
// ----------------------------------------------------------------------------
#include <CL/cl.h>

// ----------------------------------------------------------------------------
// Include Singleton abstract class
// ----------------------------------------------------------------------------
#include <Singleton.h>

namespace Aqua{ namespace InputOutput{

/** \class ProblemSetup ProblemSetup.h ProblemSetup.h
 *  Storing class filled with problem parameters (input, output, SPH, ...)
 */
class ProblemSetup : public Aqua::Singleton<Aqua::InputOutput::ProblemSetup>
{
public:
	/** Constructor.
	 */
	ProblemSetup();

	/** Destructor.
	 */
	~ProblemSetup();

	/** Method to call when all input has been parsed
	 * return false if all gone right.\n true otherwise.
	 */
	bool perform();

	/** \struct sphSettings
	 * Struct that stores general
	 * settings of the program.
	 */
	struct sphSettings
	{
		/** Se5t the default values
		 */
		void init();
		/** Clear the allocated memory
		 */
		void destroy();

	    /** Verbose level:
         *   - 2 = Slowest alternative, prints relevant data every timestep.
         *   - 1 = Prints relevant data every frame (Default value).
         *   - 0 = Fastest alternative, prints minimum data.
	     */
	    int verbose_level;

	    /// Identifier of the OpenCL platform to use.
	    unsigned int platform_id;

	    /// Identifier of the OpenCL device inside the platform.
	    unsigned int device_id;

	    /// Type of device allowed
	    cl_device_type device_type;
	}settings;

	/** \struct sphOpenCLKernels
	 * Data structure used to store the kernel file paths.
	 */
	struct sphOpenCLKernels
	{
		/** Init the kernel paths.
		 */
		void init();
		/** Delete the kernel paths.
		 */
		void destroy();

	    /// Predictor
		char *predictor;
	    /// LinkList
		char *link_list;
	    /// Rates
		char *rates;
	    /// Corrector
		char *corrector;
	    /// TimeStep
		char *time_step;

	    // Non specific tools
	    /// Reductions (It can be used for prefix sums, to compute the max, ...)
		char *reduction;
	    /// Radix sort
		char *radix_sort;

		// Optional tools
	    /// Density interpolation (formerly density reinitialization)
		char *dens_int;
	    /// 0th order correction
		char *shepard;
	    /// Elastic Bounce boundary condition
		char *elastic_bounce;
	    /// DeLeffe boundary condition (formerly Boundary Integrals method)
		char *de_Leffe;
	    /// Ghost particles boundary condition
		char *ghost;
	    /// Torque calculation
		char *torque;
	    /// Energies calculation
		char *energy;
	    /// Bounds calculation
		char *bounds;
	    /// Domain test (to exclude particles out of it)
		char *domain;
	    /// Particles teleporting models
		char *portal;
	}OpenCL_kernels;

	/** \struct sphTimingParameters
	 * Data structure used to store the time control options.
	 */
	struct sphTimingParameters
	{
	    /** Starting time.
	     * @note If a negative value is provided a stabilization period will
	     * be performed.
	     * @note During the stabilization stage motions will not be computed,
	     * and particles velocity will be vanished (infinite viscosity model).
	     */
	    float t0;

	    /** Starting time step.
	     * @remarks The starting time step is designed just for the previous
	     * simulations continuation.
	     */
	    float dt0;

	    /** Starting step.
	     * @remarks The starting step is designed just for the previous
	     * simulations continuation.
	     */
	    unsigned int step0;

	    /** Starting frame.
	     * @remarks The starting frame is designed just for the previous
	     * simulations continuation.
	     */
	    unsigned int frame0;

	    /// Provided by user time step \f$ \Delta t \f$.
	    float dt;

	    /// Minimum time step allowed.
	    float dt_min;

		/// Courant factor.
		float courant;

	    /** If a minimum time step is set, it can be additionally forced to
	     * clamp the particles velocity to preserve the Courant factor.
	     */
	    bool velocity_clamp;

		/** Type of criteria to know when the simulation has finished.
		 * @remarks Must a combination of the following options:
		 * - __TIME_MODE__ : Maximum simulated time.
		 * - __ITER_MODE__ : Maximum number os time steps.
		 * - __FRAME_MODE__ : Maximum number of output frames printed.
		 * @warning If any simulation stop cirteria is selected, the simulation
		 * will not stop
		 */
		unsigned int sim_end_mode;

		/// Maximum simulation time to be reached (in seconds).
		float sim_end_time;

		/// Maximum number of steps to be performed.
		int sim_end_step;

		/// Maximum number of frames to be performed.
		int sim_end_frame;

		/** Type of criteria to know when the log file must be printed/updated.
		 * @remarks Must a combination of the following options:
		 * - __NO_OUTPUT_MODE__ : The log file will not be never printed.
		 * - __FPS_MODE__ : Updates per second.
		 * - __IPF_MODE__ : Time steps before to update.
		 */
		unsigned int log_mode;

		/// Update rate (per second of simulation).
		float log_fps;

		/// Time steps to be done before to update the data.
		int log_ipf;

		/** Type of criteria to know when the energy report file must be
		 * printed/updated.
		 * @remarks Must a combination of the following options:
		 * - __NO_OUTPUT_MODE__ : The log file will not be never printed.
		 * - __FPS_MODE__ : Updates per second.
		 * - __IPF_MODE__ : Time steps before to update.
		 */
		unsigned int energy_mode;

		/// Update rate (per second of simulation).
		float energy_fps;

		/// Time steps to be done before to update the data.
		int energy_ipf;

		/** Type of criteria to know when the bounds of fluid file must be
		 * printed/updated.
		 * @remarks Must a combination of the following options:
		 * - __NO_OUTPUT_MODE__ : The log file will not be never printed.
		 * - __FPS_MODE__ : Updates per second.
		 * - __IPF_MODE__ : Time steps before to update.
		 */
		unsigned int bounds_mode;

		/// Update rate (per second of simulation).
		float bounds_fps;

		/// Time steps to be done before to update the data.
		int bounds_ipf;

		/** Type of criteria to know when a full fluid output must be
		 * printed/updated.
		 * @remarks Must a combination of the following options:
		 * - __NO_OUTPUT_MODE__ : The log file will not be never printed.
		 * - __FPS_MODE__ : Updates per second.
		 * - __IPF_MODE__ : Time steps before to update.
		 */
		unsigned int output_mode;

		/// Printing rate (per second of simulation).
		float output_fps;

		/// Time steps to be done before to print the data again.
		int output_ipf;

		/** Time step mode.
		 * @remarks Must be one of the following options:
		 * - __DT_VARIABLE__ If time step must be calculated each step.
		 * - __DT_FIXCALCULATED__ If time step must be calculated at the start of simulation.
		 * - __DT_FIX__ If the time step is provided by user.
		 */
	    unsigned int dt_mode;
	}time_opts;

	/** \struct sphSPHParameters
	 * Data structure used to store the SPH options.
	 */
	struct sphSPHParameters
	{
		/// Garivity aceleration \f$ \mathbf{g} \f$.
		vec g;
		/// Kernel height factor \f$ \frac{h}{\Delta r} \f$.
		float hfac;
		/// Tipical particle distance \f$ \Delta r \f$.
		vec deltar;
		/// kernel height \f$ h \f$.
		float h;
		/// Characteristic sound speed \f$ c_s \f$.
		float cs;
	    /// LinkList steps (steps before LinkList must be performed again).
	    unsigned int link_list_steps;

	    /** Density interpolation steps (steps before density interpolation
         * must be performed again).
	     * 0 if density interpolation will not be applied.
	     */
	    unsigned int dens_int_steps;
	    /** Minimum tolerated value for the density, 0  by default. If the
	     * density of a particle goes below this value it will be clamped.
	     */
	    float rho_min;
	    /** Maximum tolerated value for the density, -1 by default.
	     * If this value is lower or equal to the minimum value,
	     * maximum density value clamping will not be applied.
	     */
        float rho_max;

	    /** Boundary type:
	     *   - 0 = Elastic bounce.
	     *   - 1 = Fixed particles.
	     *   - 2 = DeLeffe.
	     * @note Ghost particles boundary condition are managed in a different
         * way.
	     */
	    unsigned int boundary_type;
	    /** Slip condition:
	     * 0 = No Slip condition.
	     * 1 = Free Slip condition.
	     */
	    unsigned int slip_condition;
	    /** Elastic factor for elastic bounce boundary condition. When a
	     * particle velocity is mirrored to avoid the wall trespass, part of
	     * the energy can be removed. Elastic factor indicates the amount of
	     * energy preserved.
	     */
	    float elastic_factor;
	    /** Minimum distance to the wall \f$ r_{min} \f$. If a particle is
	     * nearer than this distance a elastic bounce will be performed.
	     * @note Distance is provided as kernel height rate
	     * \f$ \vert \mathbf{r}_{min} \vert = r_{min} h \f$.
	     */
	    float elastic_dist;

	    /** 0th order correction (formerly Shepard correction)
	     * \f$ \gamma(\mathbf{x}) = \int_{Omega}
	     *     W(\mathbf{y} - \mathbf{x})
	     *     \mathrm{d}\mathbf{x} \f$.
	     * Must be one of the following values:
	     *   - 0 = No 0th order correction.
	     *   - 1 = It will be applied to the computed forces.
	     *   - 2 = It will be applied to computed density rate.
	     *   - 3 = It will be applied to both of them.
	     * @remarks 0th order correction (formerly Shepard correction) is
         * mandatory when Boundary Integrals based boundary condition is
         * imposed, but since Shepard term can carry some instabilities, it
         * is strongly recommended to do not apply it if it is not required by
         * the formulation.
	     */
	    unsigned int has_shepard;

	    /** Domain bounds flag. If true, domain bounds will be established such
	     * that the particles out of the domian will be replaced by zero mass
	     * fixed particles. The particle will not be deleted.
	     */
	    bool has_domain;
	    /// Minimum domain bounds coordinates.
	    vec domain_min;
	    /// Maximum domain bounds coordinates.
	    vec domain_max;
	    /// Must domain be translated following the motions.
	    bool domain_motion;
	}SPH_opts;

	/** \struct sphFluidParameters
	 * Data structure used to store the Fluids options.
	 */
	struct sphFluidParameters
	{
		/** Init the fluid options structure
		 */
		void init();

		/** Delete the allocated memory
		 */
		void destroy();

	    /// Number of particles
	    unsigned int n;

		/** Gamma exponent for the EOS \f$ \gamma \f$.
		 * \f$ p = \frac{c_s^2 \rho_0}{\gamma}
		 \left(
            \left( \frac{\rho}{\rho_0} \right)^\gamma - 1
         \right) \f$.
		 * @note The fluids which an undefined gamma value will use this one.
		 */
		float gamma;

		/// Density of reference
		float refd;

		/// Dynamic viscosity \f$ \mu \f$
		float visc_dyn;

		/// Kinematic viscosity \f$ \nu \f$
		float visc_kin;

	    /** Minimum artificial viscosity \f$\alpha\f$ value.
	     * \f$\alpha\f$ define the artificial viscosity of the fluid, such that
         * the dynamic viscosity will be: \n
	     * \f$ \mu \geq \frac{\alpha}{8} \rho c_s h \f$. \n
	     * Hence if the fluid viscosity is not larger enough to achieve the
	     * requested \f$\alpha\f$ it will be conveniently modified.
	     * @note \f$\alpha \geq 0.01\f$ is suggested if \f$\delta\f$ is null.
	     */
		float alpha;

		/** Continuity equation diffusive term factor (formerly
         * \f$\delta\f$-SPH).
		 * This term have some benefits in front of the artificial viscosity,
		 * and therefore you may consider setting \f$\alpha = 0\f$ and
		 * \f$\delta = 1\f$.
		 */
        float delta;

		/** Dynamic viscosity (corrected by the articial viscosity parameter
         * \f$\alpha\f$)
         */
		float visc_dyn_corrected;

		/// Fluid particles data file to load
		char *in_path;

		/// Fluid particles data file to load format
		char *in_format;

		/// Fluid particles data file to write
		char *out_path;

		/// Fluid particles data file to load format
		char *out_format;
	}*fluids;

	/// Number of fluids
	unsigned int n_fluids;
	/** Add a fluid to the list
	 */
	void addFluid();

	/** \class sphMoveParameters ProblemSetup.h ProblemSetup.h
	 * Data structure used to store the motions data.
	 */
	class sphMoveParameters
	{
	public:
	    /** Constructor
	     */
	    sphMoveParameters();
	    /** Destructor
	     */
	    ~sphMoveParameters();

	    /** Type of movement: \n
	     *   - -1 = No movement.
	     *   - 0 = Quaternion (Fully controlled by the OpenCL script).
	     *   - 1 = LIQuaternion (Lineary interpolation from a data file).
	     *   - 2 = C1Quaternion (Continuous C1 interpolation from a data file).
	     *   - 3 = ScriptQuaternion (Python controlled motion).
	     */
	    int type;

	    /// Motions definition XML file
	    char* path;
	};

	/// Array of motions
	std::deque<sphMoveParameters*> motions;

	/** \class sphSensorsParameters ProblemSetup.h ProblemSetup.h
	 * Data structure used to store the sensors.
	 */
	class sphSensorsParameters
	{
	public:
		/** Constructor
		 */
		sphSensorsParameters();
		/** Destructor
		 */
		~sphSensorsParameters();

	    /// Frecuency of output
	    float fps;

		/// Script
		char *script;

	    /// Array of positions
	    std::deque<vec> pos;

	    /** Method to add a sensor.
	     * @param position Position of sensor.
	     * @return false if all gone right, true otherwise.
	     */
	    bool add(vec position);
	}SensorsParameters;

	/** \struct sphPortal
	 * Data structure used to store the portals.
	 * Portals are rectangular planes that transfer particles that pass troguth them
	 * from one side to the other.
	 */
	struct sphPortal
	{
	    /** \struct Portal
	     * Specific portal storage.
	     * @note Interior normals are considered.
	     */
	    struct Portal
	    {
	        /// Start corner
	        vec corner;
	        /// Up vector
	        vec up;
	        /// Side vector
	        vec side;
	        /// Normal
	        vec normal;
	    }in,out;
	};
	/// Array of portal pairs
	std::deque<sphPortal*> portals;

	/** \class sphGhostParticles ProblemSetup.h ProblemSetup.h
	 * Data structure used to store the Ghost particles data.
	 * and walls.
	 */
	class sphGhostParticles
	{
	public:
	    #ifdef HAVE_3D
	        /** \struct Wall
	         * Wall defined by 4 points.
	         */
	        struct Wall
	        {
	            /// 1st Corner
	            vec p1;
	            /// 2nd Corner
	            vec p2;
	            /// 3rd Corner
	            vec p3;
	            /// 4th Corner
	            vec p4;
	            /// Normal to face
	            vec n;
	            /// 1st Corner velocity
	            vec v1;
	            /// 2nd Corner velocity
	            vec v2;
	            /// 3rd Corner velocity
	            vec v3;
	            /// 4th Corner velocity
	            vec v4;
	        };
	        /// Stored walls
	        std::deque<Wall*> walls;

	        /** Add a quadangular wall.
	         * @param p1 1st corner of wall.
	         * @param p2 2nd corner of wall.
	         * @param p3 3rd corner of wall.
	         * @param p4 4th corner of wall.
	         * @param v1 1st corner of wall velocity.
	         * @param v2 2nd corner of wall velocity.
	         * @param v3 3rd corner of wall velocity.
	         * @param v4 4th corner of wall velocity.
	         * @return false if all gone right, true otherwise.
	         * @remarks Normal will be computed using corners data,
	         * in order to get best results is strongly recommended
	         * try to use planar faces.
	         * @warning Program will assume corners connected as
	         * \f$p_1 \rightarrow p_2 \rightarrow p_3 \rightarrow p_4 \rightarrow p_1\f$.
	         */
	        bool add(vec p1,
                     vec p2,
                     vec p3,
                     vec p4,
                     vec v1,
                     vec v2,
                     vec v3,
                     vec v4);

	        /** Add a triangular wall.
	         * @param p1 1st corner of wall.
	         * @param p2 2nd corner of wall.
	         * @param p3 3rd corner of wall.
	         * @param v1 1st corner of wall velocity.
	         * @param v2 2nd corner of wall velocity.
	         * @param v3 3rd corner of wall velocity.
	         * @return false if all gone right, true otherwise.
	         */
	        bool add(vec p1,
                     vec p2,
                     vec p3,
                     vec v1,
                     vec v2,
                     vec v3);
	    #else
	        /** \struct Wall
	         * Wall defined by 2 points.
	         */
	        struct Wall
	        {
	            /// 1st Corner
	            vec p1;
	            /// 2nd Corner
	            vec p2;
	            /// Normal to face
	            vec n;
	            /// 1st Corner velocity
	            vec v1;
	            /// 2nd Corner velocity
	            vec v2;
	        };
	        /// Stored walls
	        std::deque<Wall*> walls;

	        /** Add wall.
	         * @param p1 1st corner of wall.
	         * @param p2 2nd corner of wall.
	         * @param v1 1st corner of wall velocity.
	         * @param v2 2nd corner of wall velocity.
	         * @return false if all gone right, true otherwise.
	         */
	        bool add(vec p1,
                     vec p2,
                     vec v1,
                     vec v2);
	    #endif

	    /** Pressure extension model. \n
	     * 0 = ASM (Antisymmetric) \n
	     * 1 = SSM (Symmetric) \n
	     * 2 = Takeda (Variation of ASM).
	     * @note Hydrostatic pressure will be added.
	     */
	    uint p_extension;

	    /** Normal velocity extension model. \n
	     * 0 = ASM (Antisymmetric) \n
	     * 1 = SSM (Symmetric) \n
	     * 2 = Takeda (Variation of ASM). \n
	     * 3 = U0M (No velocity will considered).
	     * @note Velocity extension is composed to wall velocity.
	     */
	    uint vn_extension;

	    /** Tangent velocity extension model. \n
	     * 0 = ASM (Antisymmetric) \n
	     * 1 = SSM (Symmetric) \n
	     * 2 = Takeda (Variation of ASM). \n
	     * 3 = U0M (No velocity will considered).
	     * @note Velocity extension is composed to wall velocity.
	     */
	    uint vt_extension;

	}ghost_particles;

};

}}  // namespace

#endif // PROBLEMSETUP_H_INCLUDED
