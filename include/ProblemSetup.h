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
		/** Init the kernels path
		 */
		void init();
		/** Destroy the paths
		 */
		void destroy();

	    /** Verbose level: \n
	     * 2 Slowest alternative, prints relevant data every timestep.
	     * 1 (Default value) Prints relevant data every frame.
	     * 0 Fastest alternative, prints minimum data.
	     */
	    int verbose_level;

	    /** Start mode: \n
	     * 0 Starts at t=0.
	     * 1 Starts at latest time (using H5Part file).
	     */
	    int start_mode;

	    /// Path of the file to read fluid.
	    char* fluid_file;

	    /// Identifier of the OpenCL platform to use.
	    unsigned int platform_id;

	    /// Identifier of the OpenCL device inside the platform.
	    unsigned int device_id;

	    /// Type of device allowed
	    cl_device_type device_type;
	}settings;

	/** \struct sphOpenCLKernels
	 * Struct that stores the kernels
	 * that will be used.
	 */
	struct sphOpenCLKernels
	{
		/** Init the kernels path
		 */
		void init();
		/** Destroy the paths
		 */
		void destroy();

	    /// Predictor kernel
		char *predictor;
	    /// LinkList kernel
		char *link_list;
	    /// Rates kernel
		char *rates;
	    /// Corrector kernel
		char *corrector;
	    /// TimeStep kernel
		char *time_step;

	    // Non specific tools
	    /// Reduction tool
		char *reduction;
	    /// Radix sort tool kernels
		char *radix_sort;

		// Optional tools
	    /// Density interpolation
		char *dens_int;
	    /// 0th order correction
		char *shepard;
	    /// Elastic Bounce boundary condition
		char *elastic_bounce;
	    /// DeLeffe boundary condition
		char *de_Leffe;
	    /// Torque calculation
		char *torque;
	    /// Energies calculation
		char *energy;
	    /// Bounds calculation
		char *bounds;
	    /// Domain bounds test
		char *domain;
	    /// Teleporting particles
		char *portal;
	    /// Ghost particles boundary condition
		char *ghost;
	}OpenCL_kernels;

	/** \struct sphTimingParameters
	 *  Struct that stores timing parameters
	 */
	struct sphTimingParameters
	{
		/** Type of criteria to know when the simulation reachs the final.
		 * @remarks Must a combination of this:\n
		 * __TIME_MODE__ (default option): The simulation stops when this time is reached (Time/T in XML files).\n
		 * __ITER_MODE__ : The simulation stops when this number of frames is performed (Frames/F in XML files).
		 * @warning If any simulation stop cirteria is selected, the simulation will not stop never.\n
		 * Incompatible with H5Part output.
		 */
		unsigned int SimTimingMode;
		/// Maximum simulation time to be reached (in seconds).
		float SimMaxTime;
		/// Maximum number of frames to be performed.
		int SimMaxSteps;
		/// Maximum number of frames to be performed.
		int SimMaxFrames;

		/** Type of criteria to know when a log file must be printed.
		 * @remarks Must a combination of this:\n
		 * __NO_OUTPUT_MODE__ (default option) Discard all other options (No in XML files).\n
		 * __FPS_MODE__ Log files printed per second mode (FPS in XML files).\n
		 * __IPF_MODE__ Iteration to print a Log file (IPF in XML files).
		 */
		unsigned int LogTimingMode;
		/// Log printed per second (in seconds).
		float LogFPS;
		/// Frames to print a log.
		int LogIPF;

		/** Type of criteria to know when a energy report file must be printed.
		 * @remarks Must a combination of this:\n
		 * __NO_OUTPUT_MODE__ (default option) Discard all other options (No in XML files).\n
		 * __FPS_MODE__ Energy reports printed per second mode (FPS in XML files).\n
		 * __IPF_MODE__ Iteration to print a Log file (IPF in XML files).
		 */
		unsigned int ReTimingMode;
		/// Energy report printed per second (in seconds).
		float ReFPS;
		/// Frames to print a energy report.
		int ReIPF;

		/** Type of criteria to know when a Bounds report file must be printed.
		 * @remarks Must a combination of this:\n
		 * __NO_OUTPUT_MODE__ (default option) Discard all other options (No in XML files).\n
		 * __FPS_MODE__ Bounds reports printed per second mode (FPS in XML files).\n
		 * __IPF_MODE__ Iteration to print a Log file (IPF in XML files).
		 */
		unsigned int BoundsTimingMode;
		/// Bounds report printed per second (in seconds).
		float BoundsFPS;
		/// Frames to print a Bounds report.
		int BoundsIPF;

		/** Type of criteria to know when an output file must be printed.
		 * @remarks Must a combination of this:\n
		 * __NO_OUTPUT_MODE__ (default option) Discard all other options (No in XML files).\n
		 * __FPS_MODE__ Output files printed per second mode (FPS in XML files).\n
		 * __IPF_MODE__ Iteration to print a Log file (IPF in XML files).
		 */
		unsigned int OutputMode;
		/** Format of the output.
		 * @remarks Must be a combination of this:\n
		 * __H5Part__ H5Parts format files to view into ParaView (H5Part reader plugin must be activated).\n
		 * __VTK__ VTK format files to view into ParaView.\n
		 * __TECPLOT__ Tecplot format (in ascii) files to view into Tecplot.
		 * @warning If any format is specified, any file will be printed.\n
		 * To could use this formats, the program must be compiled with support for this.
		 */
		unsigned int OutputFormat;
		/// Output files printed per second (in seconds).
		float OutputFPS;
		/// Frames to print an output file.
		int OutputIPF;

		/** Time step mode.
		 * @remarks Must be one of this: \n
		 * __DT_VARIABLE__ If time step must be calculated every iteration in order to preserve Courant number. \n
		 * __DT_FIXCALCULATED__ If time step must be calculated at the start of simulation, and conserved along it. \n
		 * __DT_FIX__ If time step is provided by user.
		 */
	    unsigned int dtMode;
	    /** Provided by user time step.
	     */
	    float dt;
	    /** Minimum time step allowed.
	     */
	    float mindt;
	    /** If a minimum time step is set, additionally can be forced to clamp the
	     * velocity of the particles in order to ensure that no more than 0.1h distance can be moved in a time step.
	     */
	    bool clampV;
	    /** Stabilization time. During this time movements will not computed, and accelerations will used as velocities,
	     * in order to reduce inertial particles movements. Outputs will not printed too.
	     */
	    float stabTime;
	}TimeParameters;

	/** \struct sphSPHParameters
	 *  Struct that stores SPH parameters
	 */
	struct sphSPHParameters
	{
		/** Default gamma
		 * @note All the fluids uses this gamma as default.
		 */
		float gamma;
		/// Garivity aceleration
		vec g;
		/// Factor between particles distance and kernel height.
		float hfac;
		/// Tipical particle distance
		vec deltar;
		/// kernel height (internal data)
		float h;
		/// Maximum density found
		float rhomax;
		/// Minimum density found
		float rhomin;
		/// Density rate
		float rhorat;
		/// Sound speed
		float cs;
		/// Time step divisor
		float DivDt;
	    /// LinkList steps (steps before LinkList must be performed).
	    unsigned int LLSteps;

	    /** Density interpolation steps (steps before density interpolation must be performed).
	     * 0 if density interpolation don't needed.
	     */
	    unsigned int DensSteps;
	    /// Minimum tolerated value for the density, 0  by default
	    float minDens;
	    /** Maximum tolerated value for the density, -1 by default.
	     * If this value is lower or equal to the minimum value,
	     * maximum density value clamping will not be applied.
	     */
        float maxDens;

	    /** Boundary type: \n
	     * 0 = Elastic bounce. \n
	     * 1 = Fixed particles. \n
	     * 2 = DeLeffe.
	     */
	    unsigned int Boundary;
	    /** Slip condition: \n
	     * 0 = No Slip condition. \n
	     * 1 = Free Slip condition.
	     */
	    unsigned int SlipCondition;
	    /** Elastic factor for Elastic Bounce boundary condition. When a particle is reflected
	     * because will pass trought boundary, part of the energy can be removed. Elastic
	     * factor indicates the amount of energy preserved.
	     */
	    float BoundElasticFactor;
	    /** Minimum distance to wall rate, if a particle is nearer than this distance
	     * velocity will forced to don't advance against the wall.
	     * @note Distance is provided as kernel height rate.
	     */
	    float BoundDist;

	    /** 0th order correction (formerly Shepard correction) application: \n
	     * 0 = No 0th order correction. \n
	     * 1 = 0th order correction will be applied over computed forces. \n
	     * 2 = 0th order correction will be applied over computed density rate. \n
	     * 3 = 0th order correction will be applied over both of them.
	     * @remarks 0th order correction (formerly Shepard correction) is mandatory when
	     * free surface is being computed, or DeLeffe boundary condition imposed,
	     * but since Shepard term can be a little bit unstable, will not executed if is
	     * disabled from options (is active by default).
	     */
	    unsigned int isShepard;

	    /** Domain bounds flag. If true, domain bounds will used.
	     * If a particle go out from domain bounds, will be corrected
	     * as zero mass fixed particle. The particle will not deleted
	     * or removed.
	     */
	    bool hasDomain;
	    /// Minimum domain bounds coordinates
	    vec minDomain;
	    /// Maximum domain bounds coordinates
	    vec maxDomain;
	    /// Must domain be translated
	    bool moveDomain;
	}SPHParameters;

	/** \struct sphFluidParameters
	 *  Struct that stores Fluid parameters
	 */
	struct sphFluidParameters
	{
	    /// Number of particles
	    unsigned int n;
		/// Gamma for this fluid
		float gamma;
		/// Referenced density
		float refd;
		/// Dynamic viscosity
		float Viscdyn;
		/// Kinematic viscosity
		float Visckin;
	    /** Minimum artificial viscosity alpha value.
	     * alpha define the artificial viscosity term
	     * that achieve the method stability. Dynamic
	     * viscosity will be: \n
	     * \f$ \mu \geq \frac{\alpha}{8} \rho v h \f$
	     */
		float alpha;
		/** Continuity equation diffusive term factor.
		 * You may use 1.0 and disable the artificial viscosity
         * (\f$\alpha=0\f$), or disable it (\f$\delta=0\f$) and use the
         * artificial viscosity to grant the stable integration
         * (\f$\alpha \geq 0.01\f$ is strongly recommended).
		 */
        float delta;
		/// Alpha corrected dynamic viscosity (internal data)
		float ViscdynCorr;
		/// Script path
		char *Path;
		/// Script
		char *Script;
		/// Load distribution file path
		char *LoadPath;

		/** Init the fluid parameters (also allocating distribution)
		 * @param P Problem setup, to take default values.
		 */
		void init(ProblemSetup *P);
		/** Destroy the fluid parameters
		 */
		void destroy();
	}*FluidParameters;
	/// Number of fluids
	unsigned int nFluids;
	/** Number of sphFluidParameters allocated (internal data).
	 * @remarks This variable is designed to avoid DYNAMIC allocation problems.
	 */
	unsigned int dimFluids;
	/** Method to add a fluid
	 * @warning Don't increment nFluids manually.
	 */
	void AddFluid();

	/** \class sphMoveParameters ProblemSetup.h ProblemSetup.h
	 * Struct with movement data.
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
	     * <ul><li>-1 = No movement.</li>
	     * <li>0 = Quaternion.</li>
	     * Manually specified quaternion. In this mode a
	     * quaternion is provided. The quaternion is
	     * specified by center and moving axis.
	     * Unless you provide custom kernel, this movement
	     * only moves fix particles (of all solids)
	     * <li>1 = LIQuaternion.</li>
	     * Linear interpolated quaternion. In this mode a
	     * quaternion file is provided. The quaternion is
	     * specified by center, x axis vector and y axis
	     * vector (x,y normalized). Unless you provide
	     * custom kernel, this movement only moves fix
	     * particles (of all solids)
	     * <li>2 = C1Quaternion.</li>
	     * C1 interpolated quaternion. In this mode a
	     * quaternion file is provided. The quaternion is
	     * specified by center, x axis vector and y axis
	     * vector (x,y normalized). Unless you provide
	     * custom kernel, this movement only moves fix
	     * particles (of all solids)
	     * <li>3 = ScriptQuaternion.</li>
	     * External scripted quaternion. In this mode,
	     * an external script provided by user controls
	     * the quaternion movement. Unless you provide
	     * custom kernel, this movement only moves fix
	     * particles (of all solids)
	     * </ul>
	     */
	    int MoveType;

	    /// Movement definition file
	    char* defFile;
	};

	/// Array of movements
	std::deque<sphMoveParameters*> MoveParameters;

	/** \class sphSensorsParameters ProblemSetup.h ProblemSetup.h
	 * Struct that stores Sensors parameters
	 */
	class sphSensorsParameters
	{
	public:
	    /// Frecuency of output
	    float fps;

		/// Script
		char *script;

	    /// Array of positions
	    std::deque<vec> pos;

	    /** Array of sensors types: \n
	     * 0 = Mean value sensor. \n
	     * 1 = Maximum value sensor.
	     */
	    std::deque<cl_ushort> mod;

		/** Constructor
		 */
		sphSensorsParameters();
		/** Destructor
		 */
		~sphSensorsParameters();
	    /** Method to add a sensor.
	     * @param position Position of sensor.
	     * @param mode Sensor type: \n
	     * 0 = Mean value sensor. \n
	     * 1 = Maximum value sensor.
	     * @return false if all gone right. \n true otherwise.
	     */
	    bool add(vec position, cl_ushort mode);
	}SensorsParameters;

	/** \struct sphPortal
	 * Struct that stores portals pair parameters.
	 * Portals are rectangular planes that transfer particles that pass troguth them
	 * from one side to the other one. Portals are clasified as in and out.
	 */
	struct sphPortal
	{
	    /** \struct Portal
	     * Specific portal storage, that get inlet or outlet.
	     * @note Interior normals.
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
	std::deque<sphPortal*> Portals;

	/** \class sphGhostParticles ProblemSetup.h ProblemSetup.h
	 * Struct that stores ghost particles boundary condition parameters
	 * and walls.
	 */
	class sphGhostParticles
	{
	public:
	    #ifdef HAVE_3D
	        /** \struct Wall
	         * Wall defined by 4 walls.
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
	         * @return false if all gone right. \n true otherwise.
	         * @remarks Normal will be computed using corners data,
	         * in order to get best results is strongly recommended
	         * try to use planar faces.
	         * @warning Program will assume corners connected as
	         * \f$p_1 \rightarrow p_2 \rightarrow p_3 \rightarrow p_4 \rightarrow p_1\f$.
	         */
	        bool add(vec p1, vec p2, vec p3, vec p4);

	        /** Add a triangular wall.
	         * @param p1 1st corner of wall.
	         * @param p2 2nd corner of wall.
	         * @param p3 3rd corner of wall.
	         * @return false if all gone right. \n true otherwise.
	         */
	        bool add(vec p1, vec p2, vec p3);
	    #else
	        /** \struct Wall
	         * Wall defined by 4 walls.
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
	         * @return false if all gone right. \n true otherwise.
	         */
	        bool add(vec p1, vec p2);
	    #endif

	    /** Pressure extension model. \n
	     * 0 = ASM (Antisymmetric) \n
	     * 1 = SSM (Symmetric) \n
	     * 2 = Takeda (Variation of ASM).
	     * @note Hydrostatic pressure will be added.
	     */
	    uint pressModel;

	    /** Normal velocity extension model. \n
	     * 0 = ASM (Antisymmetric) \n
	     * 1 = SSM (Symmetric) \n
	     * 2 = Takeda (Variation of ASM). \n
	     * 3 = U0M (No velocity will considered).
	     * @note Velocity extension is composed to wall velocity.
	     */
	    uint nVelModel;

	    /** Tangent velocity extension model. \n
	     * 0 = ASM (Antisymmetric) \n
	     * 1 = SSM (Symmetric) \n
	     * 2 = Takeda (Variation of ASM). \n
	     * 3 = U0M (No velocity will considered).
	     * @note Velocity extension is composed to wall velocity.
	     */
	    uint tVelModel;

	}GhostParticles;

};

}}  // namespace

#endif // PROBLEMSETUP_H_INCLUDED
