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

#include <limits>

// ----------------------------------------------------------------------------
// Include the main header
// ----------------------------------------------------------------------------
#include <ProblemSetup.h>

// ----------------------------------------------------------------------------
// Include the screen manager
// ----------------------------------------------------------------------------
#include <ScreenManager.h>

// ----------------------------------------------------------------------------
// Include auxiliar methods
// ----------------------------------------------------------------------------
#include <AuxiliarMethods.h>

namespace Aqua{ namespace InputOutput{

ProblemSetup::ProblemSetup()
{
	//! 1st.- General settings
	settings.init();
	//! 2nd.- OpenCL kernels
	OpenCLKernels.init();
	//! 3rd.- Init timing values
	TimeParameters.SimTimingMode = __NO_OUTPUT_MODE__;
	TimeParameters.SimMaxTime = 0.f;
	TimeParameters.SimMaxSteps  = 0;
	TimeParameters.SimMaxFrames = 0;
	TimeParameters.LogTimingMode = __NO_OUTPUT_MODE__;
	TimeParameters.LogFPS = 0.f;
	TimeParameters.LogIPF = 0;
	TimeParameters.ReTimingMode = __NO_OUTPUT_MODE__;
	TimeParameters.ReFPS = 0.f;
	TimeParameters.ReIPF = 0;
	TimeParameters.BoundsTimingMode = __NO_OUTPUT_MODE__;
	TimeParameters.BoundsFPS = 0.f;
	TimeParameters.BoundsIPF = 0;
	TimeParameters.OutputMode = __NO_OUTPUT_MODE__;
	TimeParameters.OutputFormat = __NO_OUTPUT_MODE__;
	TimeParameters.OutputFPS = 0.f;
	TimeParameters.OutputIPF = 0;
	TimeParameters.dtMode = __DT_VARIABLE__;
	TimeParameters.dt = 0.f;
	TimeParameters.mindt = 0.f;
	TimeParameters.clampV = false;
	TimeParameters.stabTime = 0.f;
	//! 4th.- SPH parameters (as Vortex problem)
	SPHParameters.gamma = 7.f;
	SPHParameters.g.x = 0.f;
	SPHParameters.g.y = 0.f;
	SPHParameters.hfac = 4.f/3.f;
	SPHParameters.deltar.x = 0.05f;
	SPHParameters.deltar.y = 0.05f;
	SPHParameters.h = SPHParameters.hfac * (SPHParameters.deltar.x + SPHParameters.deltar.y)/2.f;
	SPHParameters.rhomax = 1.f;
	SPHParameters.rhomin = 1.f;
	SPHParameters.rhorat = 1.f;
	SPHParameters.cs    = 5.f;
	SPHParameters.DivDt = 8.f;
	SPHParameters.LLSteps = 1;
	SPHParameters.DensSteps = 0;
	SPHParameters.minDens = 0.f;
	SPHParameters.maxDens = -1.f;
	SPHParameters.Boundary = 1;
	SPHParameters.SlipCondition = 0;
	SPHParameters.BoundElasticFactor = 1.f;
	SPHParameters.BoundDist = 0.1f;
	SPHParameters.isShepard = 0;
	SPHParameters.hasDomain = false;
	SPHParameters.minDomain.x = 0.f;
	SPHParameters.minDomain.y = 0.f;
	SPHParameters.maxDomain.x = 0.f;
	SPHParameters.maxDomain.y = 0.f;
	SPHParameters.moveDomain  = false;
	#ifdef HAVE_3D
	    SPHParameters.g.z = 0.f;
	    SPHParameters.g.w = 0.f;
	    SPHParameters.deltar.z = 0.05f;
	    SPHParameters.deltar.w = 0.f;
	    SPHParameters.h = SPHParameters.hfac * (SPHParameters.deltar.x + SPHParameters.deltar.y + SPHParameters.deltar.z)/3.f;
	    SPHParameters.minDomain.z = 0.f;
	    SPHParameters.minDomain.w = 0.f;
	    SPHParameters.maxDomain.z = 0.f;
	    SPHParameters.maxDomain.w = 0.f;
	#endif
	//! 5th.- Fluid parameters
	nFluids = 0;
	dimFluids = 10;
	FluidParameters = new sphFluidParameters[dimFluids];
	//! 6th.- Ghost particles parameters.
	GhostParticles.pressModel = 1;
	GhostParticles.nVelModel  = 0;
	GhostParticles.tVelModel  = 1;
}

ProblemSetup::~ProblemSetup()
{
	unsigned int i;
	settings.destroy();
	OpenCLKernels.destroy();
	for(i=0;i<nFluids;i++)
	{
		FluidParameters[i].destroy();
	}
	delete[] FluidParameters; FluidParameters=0;
	for(i=0;i<MoveParameters.size();i++){
	    delete MoveParameters.at(i);
	}
	MoveParameters.clear();
	for(i=0;i<Portals.size();i++){
	    delete Portals.at(i);
	}
	Portals.clear();
	for(i=0;i<GhostParticles.walls.size();i++){
	    delete GhostParticles.walls.at(i);
	}
	GhostParticles.walls.clear();
}

bool ProblemSetup::perform()
{
	unsigned int i;
	InputOutput::ScreenManager *S = InputOutput::ScreenManager::singleton();
	char msg[512];
	strcpy(msg, "");
	//! 1st.- Check for errors
	if(nFluids <= 0)
	{
        sprintf(msg, "(ProblemSetup::perform): There are not any fluid.\n");
        S->addMessage(3, msg);
        sprintf(msg, "\tDo you include any Fluid section at input?.\n");
        S->addMessage(0, msg);
		return true;
	}
	//! 2nd.- Perform calculations
    float K = 8.f;
	#ifndef HAVE_3D
        K = 15.f;
	    SPHParameters.h = SPHParameters.hfac * (SPHParameters.deltar.x + SPHParameters.deltar.y) / 2.f;
	#else
	    SPHParameters.h = SPHParameters.hfac * (SPHParameters.deltar.x + SPHParameters.deltar.y + SPHParameters.deltar.z) / 3.f;
	#endif
	SPHParameters.rhomax = FluidParameters[0].refd;
	SPHParameters.rhomin = FluidParameters[0].refd;
	for(i=0;i<nFluids;i++)
	{
	    FluidParameters[i].ViscdynCorr = FluidParameters[i].Viscdyn;
		if(FluidParameters[i].alpha > K * FluidParameters[i].ViscdynCorr / FluidParameters[i].refd / SPHParameters.h / SPHParameters.cs){
		    FluidParameters[i].ViscdynCorr = FluidParameters[i].alpha / K * FluidParameters[i].refd * SPHParameters.h * SPHParameters.cs;
            sprintf(msg, "(ProblemSetup::perform): fluid %u dynamic viscosity corrected\n", i);
            S->addMessage(2, msg);
            sprintf(msg, "\tValue changed from %g [Pa/s] to %g [Pa/s] (alpha = %g)\n",
	                FluidParameters[i].Viscdyn,
	                FluidParameters[i].ViscdynCorr,
	                FluidParameters[i].alpha);
            S->addMessage(0, msg);
		}
		FluidParameters[i].Visckin = FluidParameters[i].ViscdynCorr / FluidParameters[i].refd;
		FluidParameters[i].alpha = K * FluidParameters[i].ViscdynCorr / FluidParameters[i].refd / SPHParameters.h / SPHParameters.cs;
		if(FluidParameters[i].refd > SPHParameters.rhomax)
			SPHParameters.rhomax = FluidParameters[i].refd;
		if(FluidParameters[i].refd < SPHParameters.rhomin)
			SPHParameters.rhomin = FluidParameters[i].refd;
	}
	SPHParameters.rhorat = SPHParameters.rhomax / SPHParameters.rhomin;

	if(SPHParameters.maxDens <= SPHParameters.minDens){
        SPHParameters.maxDens = std::numeric_limits<float>::max();
	}
	return false;
}

void ProblemSetup::sphSettings::init()
{
	verbose_level = 1;
	start_mode = 0;
	fluid_file = new char[256];
	strcpy(fluid_file, "zParticles.h5part");
	platform_id = 0;
	device_id = 0;
	device_type = CL_DEVICE_TYPE_ALL;
}

void ProblemSetup::sphSettings::destroy()
{
	delete[] fluid_file; fluid_file=0;
}

void ProblemSetup::sphOpenCLKernels::init()
{
	//! 1st.- Alloc memory for paths
	Predictor     = new char[256];
	LinkList      = new char[256];
	Rates         = new char[256];
	Corrector     = new char[256];
	TimeStep      = new char[256];
	Reduction     = new char[256];
	RadixSort     = new char[256];
	DensInt       = new char[256];
	Shepard       = new char[256];
	ElasticBounce = new char[256];
	DeLeffe       = new char[256];
	Torque        = new char[256];
	Energy        = new char[256];
	Bounds        = new char[256];
	Domain        = new char[256];
	Portal        = new char[256];
	Ghost         = new char[256];
	//! 3rd.- Set default paths
	strcpy(Predictor,    "Input/Common/Kernels/Predictor");
	strcpy(LinkList,     "Input/Common/Kernels/LinkList");
	strcpy(Rates,        "Input/Common/Kernels/Rates");
	strcpy(Corrector,    "Input/Common/Kernels/Corrector");
	strcpy(TimeStep,     "Input/Common/Kernels/TimeStep");
	strcpy(Reduction,    "Input/Common/Kernels/Reduction");
	strcpy(RadixSort,    "Input/Common/Kernels/RadixSort");
	strcpy(DensInt,      "Input/Common/Kernels/DensInt");
	strcpy(Shepard,      "Input/Common/Kernels/Shepard");
	strcpy(ElasticBounce,"Input/Common/Kernels/Boundary/ElasticBounce");
	strcpy(DeLeffe,      "Input/Common/Kernels/Boundary/DeLeffe");
	strcpy(Ghost,        "Input/Common/Kernels/Boundary/GhostParticles");
	strcpy(Torque,       "Input/Common/Kernels/Torque");
	strcpy(Energy,       "Input/Common/Kernels/Energy");
	strcpy(Bounds,       "Input/Common/Kernels/Bounds");
	strcpy(Domain,       "Input/Common/Kernels/Domain");
	strcpy(Portal,       "Input/Common/Kernels/Portal/Portal");
}

void ProblemSetup::sphOpenCLKernels::destroy()
{
	delete[] Predictor; Predictor=0;
	delete[] LinkList; LinkList=0;
	delete[] Rates; Rates=0;
	delete[] Corrector; Corrector=0;
	delete[] TimeStep; TimeStep=0;
	delete[] Reduction; Reduction=0;
	delete[] RadixSort; RadixSort=0;
	delete[] DensInt; DensInt=0;
	delete[] Shepard; Shepard=0;
	delete[] ElasticBounce; ElasticBounce=0;
	delete[] DeLeffe; DeLeffe=0;
	delete[] Ghost; Ghost=0;
	delete[] Torque; Torque=0;
	delete[] Energy; Energy=0;
	delete[] Bounds; Bounds=0;
	delete[] Domain; Domain=0;
	delete[] Portal; Portal=0;
}

void ProblemSetup::sphFluidParameters::init(ProblemSetup *P)
{
	//! 1st.- Default values.
	n       = 0;
	gamma   = P->SPHParameters.gamma;
	refd    = 1.f;
	Viscdyn = 0.744f;
	Visckin = Viscdyn/refd;
	alpha   = 0.f;
	delta   = 0.f;
	//! 2nd.- Alloc memory for scripts
	Script = new char[256];
	Path = new char[256];
	LoadPath = new char[256];
	//! 3rd.- Script default path
	strcpy(Path, "");
	strcpy(Script, "");
	strcpy(LoadPath, "");
}

void ProblemSetup::sphFluidParameters::destroy()
{
	delete[] Script;
	delete[] Path;
	delete[] LoadPath;
}

ProblemSetup::sphMoveParameters::sphMoveParameters()
	: MoveType(0)
	, defFile(0)
{
	defFile = new char[256];
	strcpy(defFile,"");
}

ProblemSetup::sphMoveParameters::~sphMoveParameters()
{
	if(defFile) delete[] defFile; defFile=0;
}

ProblemSetup::sphSensorsParameters::sphSensorsParameters()
{
	//! 1st.- Default values.
	fps = 60.f;
	pos.clear();
	mod.clear();
	//! 2nd.- Alloc memory
	script = new char[256];
	strcpy(script, "");
}

ProblemSetup::sphSensorsParameters::~sphSensorsParameters()
{
	if(script) delete[] script; script=0;
}

bool ProblemSetup::sphSensorsParameters::add(vec position, cl_ushort mode)
{
	pos.push_back(position);
	mod.push_back(mode);
	return false;
}

void ProblemSetup::AddFluid()
{
	unsigned int i;
	//! 1st.- Add the fluid to the index.
	nFluids++;
	//! 2nd.- Analize if we need alloc more fluids
	if(nFluids >= dimFluids)
	{
		sphFluidParameters *Backup = new sphFluidParameters[dimFluids];
		for(i=0;i<dimFluids;i++)
		{
			Backup[i] = FluidParameters[i];
		}
		delete[] FluidParameters;
		FluidParameters = new sphFluidParameters[dimFluids+10];
		for(i=0;i<dimFluids;i++)
		{
			FluidParameters[i] = Backup[i];
		}
		delete[] Backup; Backup=0;
		dimFluids+=10;
	}
	//! 3rd.- Init the new fluid
	FluidParameters[nFluids-1].init(this);
}

#ifdef HAVE_3D
	bool ProblemSetup::sphGhostParticles::add(vec p1, vec p2, vec p3, vec p4)
	{
	    vec d1 = sub(p2, p1);
	    vec d2 = sub(p4, p1);
	    vec n1 = cross(d1, d2);
	    d1 = sub(p1, p2);
	    d2 = sub(p3, p2);
	    vec n2 = cross(d1, d2);
	    d1 = sub(p2, p3);
	    d2 = sub(p4, p2);
	    vec n3 = cross(d1, d2);
	    d1 = sub(p3, p4);
	    d2 = sub(p1, p4);
	    vec n4 = cross(d1, d2);
	    vec n;
	    n.x = n1.x + n2.x + n3.x + n4.x;
	    n.y = n1.y + n2.y + n3.y + n4.y;
	    n.z = n1.z + n2.z + n3.z + n4.z;
	    n.w = 0.f;
	    Wall *wall = new Wall();
	    wall->p1   = p1;
	    wall->p2   = p2;
	    wall->p3   = p3;
	    wall->p4   = p4;
	    wall->n    = normalize(n);
	    wall->v1   = Vzero();
	    wall->v2   = Vzero();
	    wall->v3   = Vzero();
	    wall->v4   = Vzero();
	    walls.push_back(wall);
	    return false;
	}
	bool ProblemSetup::sphGhostParticles::add(vec p1, vec p2, vec p3)
	{
	    vec d1 = sub(p2, p1);
	    vec d2 = sub(p3, p1);
	    vec n  = cross(d1, d2);
	    Wall *wall = new Wall();
	    wall->p1   = p1;
	    wall->p2   = p2;
	    wall->p3   = p3;
	    wall->p4   = p3;
	    wall->n    = normalize(n);
	    wall->v1   = Vzero();
	    wall->v2   = Vzero();
	    wall->v3   = Vzero();
	    wall->v4   = Vzero();
	    walls.push_back(wall);
	    return false;
	}
#else
	bool ProblemSetup::sphGhostParticles::add(vec p1, vec p2)
	{
	    vec n;
	    n.x = p1.y - p2.y;
	    n.y = p2.x - p1.x;
	    Wall *wall = new Wall();
	    wall->p1   = p1;
	    wall->p2   = p2;
	    wall->n    = normalize(n);
	    wall->v1   = Vzero();
	    wall->v2   = Vzero();
	    walls.push_back(wall);
	    return false;
	}
#endif

}}  // namespace
