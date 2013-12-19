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

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

#include <Fluid.h>
#include <ProblemSetup.h>
#include <CalcServer.h>
#include <ScreenManager.h>

namespace Aqua{ namespace InputOutput{

Fluid::Fluid()
{
	unsigned int i, num_sensors;
	char msg[512];
	ProblemSetup *P = ProblemSetup::singleton();
	ScreenManager *S = ScreenManager::singleton();

	num_particles = 0;
	for(i=0;i<P->n_fluids;i++) {
	    num_particles += P->fluids[i].n;
	}
	num_sensors = P->SensorsParameters.pos.size();
	num_particles += num_sensors;
	if(num_particles<=0){
	    S->addMessageF(3, "Any particle found.\n");
		exit(1);
	}
	sprintf(msg, "Number of particles = %u\n", num_particles);
	S->addMessageF(1, msg);

	num_fluids = P->n_fluids;

	imove = new int[num_particles];
	ifluid = new int[num_particles];
	pos = new vec[num_particles];
	normal = new vec[num_particles];
	v = new vec[num_particles];
	dens = new float[num_particles];
	drdt = new float[num_particles];
	hp = new float[num_particles];
	f = new vec[num_particles];
	mass = new float[num_particles];
	press = new float[num_particles];
	shepard = new float[num_particles];
	shepard_gradient = new vec[num_particles];

	for(i=num_particles-num_sensors;i<num_particles;i++) {
	    imove[i]    = 0;
	    ifluid[i]   = 0;
	    pos[i]      = P->SensorsParameters.pos.at(i-num_particles+num_sensors);
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
    sprintf(msg, "%u particles allocated OK, we can continue happily! ;-)\n", num_particles);
	S->addMessageF(1, msg);
}


Fluid::~Fluid()
{
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
	CalcServer::CalcServer *C = CalcServer::CalcServer::singleton();
	int flag;
	flag  = C->getData((void*)imove,   C->imove,   sizeof(cl_int)*num_particles);
	flag |= C->getData((void*)ifluid,  C->ifluid,  sizeof(cl_int)*num_particles);
	flag |= C->getData((void*)pos,     C->posin,   sizeof(vec)*num_particles);
	flag |= C->getData((void*)normal,  C->normal,  sizeof(vec)*num_particles);
	flag |= C->getData((void*)v,       C->vin,     sizeof(vec)*num_particles);
	flag |= C->getData((void*)f,       C->fin,     sizeof(vec)*num_particles);
	flag |= C->getData((void*)dens,    C->densin,  sizeof(cl_float)*num_particles);
	flag |= C->getData((void*)drdt,    C->drdt,    sizeof(cl_float)*num_particles);
	flag |= C->getData((void*)press,   C->press,   sizeof(cl_float)*num_particles);
	flag |= C->getData((void*)mass,    C->mass,    sizeof(cl_float)*num_particles);
	flag |= C->getData((void*)hp,      C->hpin,    sizeof(cl_float)*num_particles);
	flag |= C->getData((void*)shepard, C->shepard, sizeof(cl_float)*num_particles);
	flag |= C->getData((void*)shepard_gradient, C->shepard_gradient, sizeof(vec)*num_particles);
	if(flag)
	    return true;
	return false;
}

}}  // namespace
