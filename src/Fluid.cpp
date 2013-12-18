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

// ----------------------------------------------------------------------------
// Include the main header
// ----------------------------------------------------------------------------
#include <Fluid.h>

// ----------------------------------------------------------------------------
// Include calculation server
// ----------------------------------------------------------------------------
#include <ProblemSetup.h>

// ----------------------------------------------------------------------------
// Include calculation server
// ----------------------------------------------------------------------------
#include <CalcServer.h>

// ----------------------------------------------------------------------------
// Include the screen manager
// ----------------------------------------------------------------------------
#include <ScreenManager.h>

namespace Aqua{ namespace InputOutput{

Fluid::Fluid()
{
	unsigned int i;
	ProblemSetup *P = ProblemSetup::singleton();
	ScreenManager *S = ScreenManager::singleton();
	char msg[512];
	//! 1st.- Take particles amount
	nParticle = 0;
	for(i=0;i<P->n_fluids;i++) {
	    nParticle += P->fluids[i].n;
	}
	unsigned int num_sensors = P->SensorsParameters.pos.size();
	nParticle += num_sensors;
	if(nParticle<=0){
	    S->addMessage(3, "(Fluid::Fluid): Any particle found.\n");
		exit(1);
	}
	sprintf(msg, "(Fluid::Fluid): Number of particles = %u\n", nParticle);
	S->addMessage(1, msg);
	//! 2nd.- Set the number of fluids
	nFluid = P->n_fluids;
	//! 3rd.- Init the fluid (allocating memory)
	imove = new int[nParticle];
	ifluid = new int[nParticle];
	pos = new vec[nParticle];
	normal = new vec[nParticle];
	v = new vec[nParticle];
	dens = new float[nParticle];
	drdt = new float[nParticle];
	hp = new float[nParticle];
	f = new vec[nParticle];
	mass = new float[nParticle];
	press = new float[nParticle];
	shepard = new float[nParticle];
	shepard_gradient = new vec[nParticle];
	//! 4th.- Set sensors (That have been already specified)
	for(i=nParticle-num_sensors;i<nParticle;i++) {
	    imove[i]    = 0;
	    ifluid[i]   = 0;
	    pos[i]      = P->SensorsParameters.pos.at(i-nParticle+num_sensors);
	    normal[i].x = 0.f; normal[i].y = 0.f;
	    v[i].x      = 0.f; v[i].y = 0.f;
	    dens[i]     = 0.f;
	    drdt[i]     = 0.f;
	    hp[i]       = P->SPH_opts.h;
	    f[i].x      = 0.f; f[i].y = 0.f;
	    mass[i]     = 0.f;
	    press[i]    = 0.f;
	    shepard[i]  = 0.f;
	    shepard_gradient[i].x = 0.f;
	    shepard_gradient[i].y = 0.f;
	    #ifdef HAVE_3D
	        normal[i].z = 0.f;
	        normal[i].w = 0.f;
	        v[i].z      = 0.f;
	        v[i].w      = 0.f;
	        f[i].z      = 0.f;
	        f[i].w      = 0.f;
	        shepard_gradient[i].z = 0.f;
	        shepard_gradient[i].w = 0.f;
	    #endif
	}
    sprintf(msg, "(Fluid::Fluid): %u particles allocated OK, we can continue happily! ;-)\n", nParticle);
	S->addMessage(1, msg);
}


Fluid::~Fluid()
{
	ProblemSetup *P = ProblemSetup::singleton();
	delete[] imove; imove=0;
	delete[] ifluid; ifluid=0;
	delete[] pos; pos=0;
	delete[] normal; normal=0;
	delete[] v; v=0;
	delete[] dens; dens=0;
	delete[] drdt; drdt=0;
	delete[] hp; hp=0;
	delete[] f; f=0;
	delete[] mass; mass=0;
	delete[] press; press=0;
	delete[] shepard; shepard=0;
	delete[] shepard_gradient; shepard_gradient=0;
}

bool Fluid::retrieveData()
{
	ProblemSetup *P = ProblemSetup::singleton();
	CalcServer::CalcServer *C = CalcServer::CalcServer::singleton();
	int flag;
	flag  = C->getData((void*)imove,   C->imove,   sizeof(cl_int)*nParticle);
	flag |= C->getData((void*)ifluid,  C->ifluid,  sizeof(cl_int)*nParticle);
	flag |= C->getData((void*)pos,     C->posin,   sizeof(vec)*nParticle);
	flag |= C->getData((void*)normal,  C->normal,  sizeof(vec)*nParticle);
	flag |= C->getData((void*)v,       C->vin,     sizeof(vec)*nParticle);
	flag |= C->getData((void*)f,       C->fin,     sizeof(vec)*nParticle);
	flag |= C->getData((void*)dens,    C->densin,  sizeof(cl_float)*nParticle);
	flag |= C->getData((void*)drdt,    C->drdt,    sizeof(cl_float)*nParticle);
	flag |= C->getData((void*)press,   C->press,   sizeof(cl_float)*nParticle);
	flag |= C->getData((void*)mass,    C->mass,    sizeof(cl_float)*nParticle);
	flag |= C->getData((void*)shepard, C->shepard, sizeof(cl_float)*nParticle);
	flag |= C->getData((void*)shepard_gradient, C->shepard_gradient, sizeof(vec)*nParticle);
	flag |= C->getData((void*)hp,      C->hpin,    sizeof(cl_float)*nParticle);
	if(flag)
	    return true;
	return false;
}

}}  // namespace
