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
	OpenCL_kernels.init();
	//! 3rd.- Init timing values
	time_opts.sim_end_mode = __NO_OUTPUT_MODE__;
	time_opts.sim_end_time = 0.f;
	time_opts.sim_end_step  = 0;
	time_opts.sim_end_frame = 0;
	time_opts.log_mode = __NO_OUTPUT_MODE__;
	time_opts.log_fps = 0.f;
	time_opts.log_ipf = 0;
	time_opts.energy_mode = __NO_OUTPUT_MODE__;
	time_opts.energy_fps = 0.f;
	time_opts.energy_ipf = 0;
	time_opts.bounds_mode = __NO_OUTPUT_MODE__;
	time_opts.bounds_fps = 0.f;
	time_opts.bounds_ipf = 0;
	time_opts.output_mode = __NO_OUTPUT_MODE__;
	time_opts.output_format = __NO_OUTPUT_MODE__;
	time_opts.output_fps = 0.f;
	time_opts.output_ipf = 0;
	time_opts.dt_mode = __DT_VARIABLE__;
	time_opts.dt = 0.f;
	time_opts.dt_min = 0.f;
	time_opts.velocity_clamp = false;
	time_opts.stabilization_time = 0.f;
	//! 4th.- SPH parameters (as Vortex problem)
	SPH_opts.gamma = 7.f;
	SPH_opts.g.x = 0.f;
	SPH_opts.g.y = 0.f;
	SPH_opts.hfac = 4.f/3.f;
	SPH_opts.deltar.x = 0.05f;
	SPH_opts.deltar.y = 0.05f;
	SPH_opts.h = SPH_opts.hfac * (SPH_opts.deltar.x + SPH_opts.deltar.y)/2.f;
	SPH_opts.cs    = 5.f;
	SPH_opts.dt_divisor = 8.f;
	SPH_opts.link_list_steps = 1;
	SPH_opts.dens_int_steps = 0;
	SPH_opts.rho_min = 0.f;
	SPH_opts.rho_max = -1.f;
	SPH_opts.boundary_type = 2;
	SPH_opts.slip_condition = 0;
	SPH_opts.elastic_factor = 1.f;
	SPH_opts.elastic_dist = 0.1f;
	SPH_opts.has_shepard = 0;
	SPH_opts.has_domain = false;
	SPH_opts.domain_min.x = 0.f;
	SPH_opts.domain_min.y = 0.f;
	SPH_opts.domain_max.x = 0.f;
	SPH_opts.domain_max.y = 0.f;
	SPH_opts.domain_motion  = false;
	#ifdef HAVE_3D
	    SPH_opts.g.z = 0.f;
	    SPH_opts.g.w = 0.f;
	    SPH_opts.deltar.z = 0.05f;
	    SPH_opts.deltar.w = 0.f;
	    SPH_opts.h = SPH_opts.hfac * (SPH_opts.deltar.x + SPH_opts.deltar.y + SPH_opts.deltar.z)/3.f;
	    SPH_opts.domain_min.z = 0.f;
	    SPH_opts.domain_min.w = 0.f;
	    SPH_opts.domain_max.z = 0.f;
	    SPH_opts.domain_max.w = 0.f;
	#endif
	//! 5th.- Fluid parameters
	nFluids = 0;
	dimFluids = 10;
	fluids = new sphFluidParameters[dimFluids];
	//! 6th.- Ghost particles parameters.
	GhostParticles.pressModel = 1;
	GhostParticles.nVelModel  = 0;
	GhostParticles.tVelModel  = 1;
}

ProblemSetup::~ProblemSetup()
{
	unsigned int i;
	settings.destroy();
	OpenCL_kernels.destroy();
	for(i=0;i<nFluids;i++)
	{
		fluids[i].destroy();
	}
	delete[] fluids; fluids=0;
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
	    SPH_opts.h = SPH_opts.hfac * (SPH_opts.deltar.x + SPH_opts.deltar.y) / 2.f;
	#else
	    SPH_opts.h = SPH_opts.hfac * (SPH_opts.deltar.x + SPH_opts.deltar.y + SPH_opts.deltar.z) / 3.f;
	#endif
	for(i=0;i<nFluids;i++)
	{
	    fluids[i].visc_dyn_corrected = fluids[i].visc_dyn;
		if(fluids[i].alpha > K * fluids[i].visc_dyn_corrected / fluids[i].refd / SPH_opts.h / SPH_opts.cs){
		    fluids[i].visc_dyn_corrected = fluids[i].alpha / K * fluids[i].refd * SPH_opts.h * SPH_opts.cs;
            sprintf(msg, "(ProblemSetup::perform): fluid %u dynamic viscosity corrected\n", i);
            S->addMessage(2, msg);
            sprintf(msg, "\tValue changed from %g [Pa/s] to %g [Pa/s] (alpha = %g)\n",
	                fluids[i].visc_dyn,
	                fluids[i].visc_dyn_corrected,
	                fluids[i].alpha);
            S->addMessage(0, msg);
		}
		fluids[i].visc_kin = fluids[i].visc_dyn_corrected / fluids[i].refd;
		fluids[i].alpha = K * fluids[i].visc_dyn_corrected / fluids[i].refd / SPH_opts.h / SPH_opts.cs;
	}

	if(SPH_opts.rho_max <= SPH_opts.rho_min){
        SPH_opts.rho_max = std::numeric_limits<float>::max();
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
	predictor      = new char[256];
	link_list       = new char[256];
	rates          = new char[256];
	corrector      = new char[256];
	time_step      = new char[256];
	reduction      = new char[256];
	radix_sort     = new char[256];
	dens_int       = new char[256];
	shepard        = new char[256];
	elastic_bounce = new char[256];
	de_Leffe       = new char[256];
	ghost         = new char[256];
	torque        = new char[256];
	energy        = new char[256];
	bounds        = new char[256];
	domain        = new char[256];
	portal        = new char[256];
	//! 3rd.- Set default paths
	strcpy(predictor,    "Input/Common/Kernels/Predictor");
	strcpy(link_list,     "Input/Common/Kernels/LinkList");
	strcpy(rates,        "Input/Common/Kernels/Rates");
	strcpy(corrector,    "Input/Common/Kernels/Corrector");
	strcpy(time_step,     "Input/Common/Kernels/TimeStep");
	strcpy(reduction,    "Input/Common/Kernels/Reduction");
	strcpy(radix_sort,    "Input/Common/Kernels/RadixSort");
	strcpy(dens_int,      "Input/Common/Kernels/DensInt");
	strcpy(shepard,      "Input/Common/Kernels/Shepard");
	strcpy(elastic_bounce,"Input/Common/Kernels/Boundary/ElasticBounce");
	strcpy(de_Leffe,      "Input/Common/Kernels/Boundary/DeLeffe");
	strcpy(ghost,        "Input/Common/Kernels/Boundary/GhostParticles");
	strcpy(torque,       "Input/Common/Kernels/Torque");
	strcpy(energy,       "Input/Common/Kernels/Energy");
	strcpy(bounds,       "Input/Common/Kernels/Bounds");
	strcpy(domain,       "Input/Common/Kernels/Domain");
	strcpy(portal,       "Input/Common/Kernels/Portal/Portal");
}

void ProblemSetup::sphOpenCLKernels::destroy()
{
	delete[] predictor; predictor=0;
	delete[] link_list; link_list=0;
	delete[] rates; rates=0;
	delete[] corrector; corrector=0;
	delete[] time_step; time_step=0;
	delete[] reduction; reduction=0;
	delete[] radix_sort; radix_sort=0;
	delete[] dens_int; dens_int=0;
	delete[] shepard; shepard=0;
	delete[] elastic_bounce; elastic_bounce=0;
	delete[] de_Leffe; de_Leffe=0;
	delete[] ghost; ghost=0;
	delete[] torque; torque=0;
	delete[] energy; energy=0;
	delete[] bounds; bounds=0;
	delete[] domain; domain=0;
	delete[] portal; portal=0;
}

void ProblemSetup::sphFluidParameters::init()
{
    ProblemSetup *P = ProblemSetup::singleton();
	//! 1st.- Default values.
	n       = 0;
	gamma   = P->SPH_opts.gamma;
	refd    = 1.f;
	visc_dyn = 0.744f;
	visc_kin = visc_dyn/refd;
	alpha   = 0.f;
	delta   = 0.f;
	//! 2nd.- Alloc memory for scripts
	Script = new char[256];
	Path = new char[256];
	path = new char[256];
	//! 3rd.- Script default path
	strcpy(Path, "");
	strcpy(Script, "");
	strcpy(path, "");
}

void ProblemSetup::sphFluidParameters::destroy()
{
	delete[] Script;
	delete[] Path;
	delete[] path;
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
			Backup[i] = fluids[i];
		}
		delete[] fluids;
		fluids = new sphFluidParameters[dimFluids+10];
		for(i=0;i<dimFluids;i++)
		{
			fluids[i] = Backup[i];
		}
		delete[] Backup; Backup=0;
		dimFluids+=10;
	}
	//! 3rd.- Init the new fluid
	fluids[nFluids-1].init();
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
