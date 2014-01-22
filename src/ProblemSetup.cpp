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
	time_opts.sim_end_step = 0;
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
	time_opts.t0 = 0.f;
	//! 4th.- SPH parameters (as Vortex problem)
	SPH_opts.g.x = 0.f;
	SPH_opts.g.y = 0.f;
	SPH_opts.hfac = 4.f/3.f;
	SPH_opts.deltar.x = 0.05f;
	SPH_opts.deltar.y = 0.05f;
	SPH_opts.h = SPH_opts.hfac * (SPH_opts.deltar.x + SPH_opts.deltar.y)/2.f;
	SPH_opts.cs = 15.f;
	SPH_opts.courant = 0.25f;
	SPH_opts.link_list_steps = 1;
	SPH_opts.dens_int_steps = 0;
	SPH_opts.rho_min = 0.f;
	SPH_opts.rho_max = -1.f;
	SPH_opts.boundary_type = 0;
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
	n_fluids = 0;
	fluids = NULL;
	//! 6th.- Ghost particles parameters.
	ghost_particles.p_extension = 1;
	ghost_particles.vn_extension = 0;
	ghost_particles.vt_extension = 1;
}

ProblemSetup::~ProblemSetup()
{
	unsigned int i;
	settings.destroy();
	OpenCL_kernels.destroy();
	for(i=0;i<n_fluids;i++)
	{
		fluids[i].destroy();
	}
	delete[] fluids; fluids=0;
	for(i=0;i<motions.size();i++){
	    delete motions.at(i);
	}
	motions.clear();
	for(i=0;i<portals.size();i++){
	    delete portals.at(i);
	}
	portals.clear();
	for(i=0;i<ghost_particles.walls.size();i++){
	    delete ghost_particles.walls.at(i);
	}
	ghost_particles.walls.clear();
}

bool ProblemSetup::perform()
{
	unsigned int i;
	InputOutput::ScreenManager *S = InputOutput::ScreenManager::singleton();
	char msg[512];
	strcpy(msg, "");
	//! 1st.- Check for errors
	if(n_fluids <= 0)
	{
        sprintf(msg, "There are not any fluid.\n");
        S->addMessageF(3, msg);
        sprintf(msg, "\tDo you include any Fluid section at input?.\n");
        S->addMessage(0, msg);
		return true;
	}
	//! 2nd.- Perform calculations
    float K = 8.f;
	#ifndef HAVE_3D
        K = 15.f;
	    SPH_opts.h = SPH_opts.hfac * (SPH_opts.deltar.x
                                   + SPH_opts.deltar.y) / 2.f;
	#else
	    SPH_opts.h = SPH_opts.hfac * (SPH_opts.deltar.x
                                   + SPH_opts.deltar.y
                                   + SPH_opts.deltar.z) / 3.f;
	#endif
	for(i=0;i<n_fluids;i++) {
	    fluids[i].visc_dyn_corrected = fluids[i].visc_dyn;
		if(fluids[i].alpha > K * fluids[i].visc_dyn_corrected / fluids[i].refd / SPH_opts.h / SPH_opts.cs){
		    fluids[i].visc_dyn_corrected = fluids[i].alpha / K * fluids[i].refd * SPH_opts.h * SPH_opts.cs;
            sprintf(msg, "fluid %u dynamic viscosity corrected\n", i);
            S->addMessageF(2, msg);
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
	platform_id = 0;
	device_id = 0;
	device_type = CL_DEVICE_TYPE_ALL;
}

void ProblemSetup::sphSettings::destroy()
{
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
	strcpy(predictor,     "Input/Common/Kernels/Predictor");
	strcpy(link_list,     "Input/Common/Kernels/LinkList");
	strcpy(rates,         "Input/Common/Kernels/Rates");
	strcpy(corrector,     "Input/Common/Kernels/Corrector");
	strcpy(time_step,     "Input/Common/Kernels/TimeStep");
	strcpy(reduction,     "Input/Common/Kernels/Reduction");
	strcpy(radix_sort,    "Input/Common/Kernels/RadixSort");
	strcpy(dens_int,      "Input/Common/Kernels/DensInt");
	strcpy(shepard,       "Input/Common/Kernels/Shepard");
	strcpy(elastic_bounce,"Input/Common/Kernels/Boundary/ElasticBounce");
	strcpy(de_Leffe,      "Input/Common/Kernels/Boundary/DeLeffe");
	strcpy(ghost,         "Input/Common/Kernels/Boundary/GhostParticles");
	strcpy(torque,        "Input/Common/Kernels/Torque");
	strcpy(energy,        "Input/Common/Kernels/Energy");
	strcpy(bounds,        "Input/Common/Kernels/Bounds");
	strcpy(domain,        "Input/Common/Kernels/Domain");
	strcpy(portal,        "Input/Common/Kernels/Portal/Portal");
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

	n       = 0;
	gamma   = 7.f;
	refd    = 1000.f;
	visc_dyn = 0.798e-3f;
	visc_kin = visc_dyn/refd;
	alpha   = 0.f;
	delta   = 0.f;

	in_path = new char[1024];
	in_format = new char[1024];
	strcpy(path, "");
	strcpy(in_format, "ASCII"):

	out_path = new char[1024];
	out_format = new char[1024];
	strcpy(path, "");
	strcpy(out_format, "VTK"):
}

void ProblemSetup::sphFluidParameters::destroy()
{
	delete[] path; path=NULL;

	delete[] Script; Script=NULL;
	delete[] Path; Path=NULL;
}

ProblemSetup::sphMoveParameters::sphMoveParameters()
	: type(0)
	, path(0)
{
	path = new char[256];
	strcpy(path,"");
}

ProblemSetup::sphMoveParameters::~sphMoveParameters()
{
	if(path) delete[] path; path=0;
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

void ProblemSetup::addFluid()
{
    ScreenManager *S = ScreenManager::singleton();
	n_fluids++;
    sphFluidParameters *backup = fluids;
    fluids = new sphFluidParameters[n_fluids];
    if(!fluids){
        S->addMessageF(3, "Memory cannot be allocated.\n");
        exit(255);
    }
    if(backup){
        memcpy(fluids, backup, (n_fluids-1)*sizeof(sphFluidParameters));
        free(backup);
        backup = NULL;
    }
    fluids[n_fluids-1].init();
}

#ifdef HAVE_3D
	bool ProblemSetup::sphGhostParticles::add(vec p1,
                                              vec p2,
                                              vec p3,
                                              vec p4,
                                              vec v1,
                                              vec v2,
                                              vec v3,
                                              vec v4)
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
	    wall->v1   = v1;
	    wall->v2   = v2;
	    wall->v3   = v3;
	    wall->v4   = v4;
	    walls.push_back(wall);
	    return false;
	}
	bool ProblemSetup::sphGhostParticles::add(vec p1,
                                              vec p2,
                                              vec p3,
                                              vec v1,
                                              vec v2,
                                              vec v3)
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
	    wall->v1   = v1;
	    wall->v2   = v2;
	    wall->v3   = v3;
	    wall->v4   = v3;
	    walls.push_back(wall);
	    return false;
	}
#else
	bool ProblemSetup::sphGhostParticles::add(vec p1,
                                              vec p2,
                                              vec v1,
                                              vec v2)
	{
	    vec n;
	    n.x = p1.y - p2.y;
	    n.y = p2.x - p1.x;
	    Wall *wall = new Wall();
	    wall->p1   = p1;
	    wall->p2   = p2;
	    wall->n    = normalize(n);
	    wall->v1   = v1;
	    wall->v2   = v2;
	    walls.push_back(wall);
	    return false;
	}
#endif

}}  // namespace
