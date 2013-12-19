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

#include <Input.h>

#ifdef HAVE_H5PART
	#include <H5Part.h>
	#include <hdf5.h>
#endif // HAVE_H5PART

#include <Fluid.h>
#include <TimeManager.h>
#include <ProblemSetup.h>
#include <CalcServer.h>
#include <AuxiliarMethods.h>
#include <ScreenManager.h>

#include <Input/ASCII.h>
#include <Input/GiD.h>
#include <Input/XML.h>

namespace Aqua{ namespace InputOutput{

bool input()
{
	Fluid *F = Fluid::singleton();
	ProblemSetup *P = ProblemSetup::singleton();
    ScreenManager *S = ScreenManager::singleton();

	unsigned int i;
	unsigned int start=0, n=0, Start=0, N=0;
	for(i=0;i<P->n_fluids;i++){
	    // File path
	    const char* path = P->fluids[i].path;
	    if(!strlen(path)){
	        continue;
	    }
	    // Get extension
	    const char* file_type = getExtensionFromFilePath(path);
	    // Get data
	    n = P->fluids[i].n;
	    // Load file with selected reader
	    if(!strcmp(file_type,"xml")){    // XML file
	        if(Input::loadXML(path, i, start, n, P->fluids[i].refd, P->SPH_opts.h, F))
	            return true;
	    }
	    if(!strcmp(file_type,"gid")){    // GiD file
	        if(Input::loadGiD(path, i, start, n, P->fluids[i].refd, P->SPH_opts.h, F))
	            return true;
	    }
	    else{                            // Plain text tabulated file
	        if(Input::loadASCII(path, i, start, n, P->fluids[i].refd, P->SPH_opts.h, F))
	            return true;
	    }
	    start += n;
	}

	if(P->settings.start_mode == 1) {
	    #ifdef HAVE_H5PART
	        if(Input::loadH5Part())
	            return true;
	    #else
            S->addMessageF(3, "Can't read H5Part input file.\n");
            S->addMessage(0, "\tFormat not supported (Consider recompile using H5PART support).\n");
	        return true;
	    #endif // HAVE_H5PART
	}
	return false;
}

#ifdef HAVE_H5PART
bool Input::loadH5Part()
{
	unsigned int i;
	CalcServer::CalcServer *C = CalcServer::CalcServer::singleton();
	Fluid *F = Fluid::singleton();
	TimeManager *T = TimeManager::singleton();
	ProblemSetup *P = ProblemSetup::singleton();
    ScreenManager *S = ScreenManager::singleton();
    char msg[256];
	H5PartFile *H5PartFileID;
	h5part_int64_t n, num_fluids, step;
	double dTime, ddt;
	h5part_int64_t *imove, *ifluid;
	double *x, *y, *nx, *ny, *vx, *vy, *dvx, *dvy, *press, *dens, *ddens, *hp, *mass;
	#ifdef HAVE_3D
	    double *z, *nz, *vz, *dvz;
	#endif

	sprintf("Loading fluid from \"%s\"\n", P->settings.fluid_file);
	S->addMessageF(1, msg);
	#ifdef HAVE_MPI
	    H5PartFileID = H5PartOpenFileParallel(P->settings.fluid_file,H5PART_READ,MPI_COMM_WORLD);
	#else
	    H5PartFileID = H5PartOpenFile(P->settings.fluid_file,H5PART_WRITE);
	#endif

	h5part_int64_t zeroFrame;
	H5PartReadFileAttrib(H5PartFileID, "zeroFrame", &zeroFrame);
	T->startFrame(zeroFrame);
	int iFrame = H5PartGetNumSteps(H5PartFileID);
	H5PartSetStep(H5PartFileID,iFrame-1);
	n = H5PartGetNumParticles(H5PartFileID);
	if(F->n() != (unsigned int)n){
	    S->addMessageF(3, "Number of particles of the file don't match than fluid particles.\n");
        sprintf(msg, "\tFluid particles = %u\n", F->n());
        S->addMessage(0, msg);
        sprintf(msg, "\tFile particles  = %u\n", n);
        S->addMessage(0, msg);
	    return true;
	}
	T->frame(iFrame+T->startFrame());
	T->startFrame(T->frame());

	H5PartReadStepAttrib(H5PartFileID,"nfluid",&num_fluids);
	if(C->num_fluids != (unsigned int)num_fluids){
	    S->addMessageF(3, "Number of fluids of the file don't match than calcserver.\n");
        sprintf(msg, "\tFluid on server = %u\n", C->num_fluids);
        S->addMessage(0, msg);
        sprintf(msg, "\tFile fluids      = %u\n", num_fluids);
        S->addMessage(0, msg);
	    return true;
	}
	H5PartReadStepAttrib(H5PartFileID,"step",&step);
	H5PartReadStepAttrib(H5PartFileID,"time",&dTime);
	H5PartReadStepAttrib(H5PartFileID,"dt",&ddt);
	T->step(step);
	T->time(dTime);
	T->startTime(dTime);
	T->dt(ddt);
	T->outputStep(step);
	T->outputTime(dTime);

	imove = new h5part_int64_t[F->n()];
	ifluid = new h5part_int64_t[F->n()];
	x     = new double[F->n()];
	y     = new double[F->n()];
	nx    = new double[F->n()];
	ny    = new double[F->n()];
	vx    = new double[F->n()];
	vy    = new double[F->n()];
	dvx   = new double[F->n()];
	dvy   = new double[F->n()];
	press = new double[F->n()];
	dens  = new double[F->n()];
	ddens = new double[F->n()];
	hp    = new double[F->n()];
	mass  = new double[F->n()];
	#ifdef HAVE_3D
	    z     = new double[F->n()];
	    nz    = new double[F->n()];
	    vz    = new double[F->n()];
	    dvz   = new double[F->n()];
	#endif
	H5PartReadDataFloat64(H5PartFileID,"x",x);
	H5PartReadDataFloat64(H5PartFileID,"y",y);
	H5PartReadDataFloat64(H5PartFileID,"nx",nx);
	H5PartReadDataFloat64(H5PartFileID,"ny",ny);
	H5PartReadDataFloat64(H5PartFileID,"vx",vx);
	H5PartReadDataFloat64(H5PartFileID,"vy",vy);
	H5PartReadDataFloat64(H5PartFileID,"dvx",dvx);
	H5PartReadDataFloat64(H5PartFileID,"dvy",dvy);
	H5PartReadDataFloat64(H5PartFileID,"press",press);
	H5PartReadDataFloat64(H5PartFileID,"dens",dens);
	H5PartReadDataFloat64(H5PartFileID,"ddens",ddens);
	H5PartReadDataFloat64(H5PartFileID,"hp",hp);
	H5PartReadDataFloat64(H5PartFileID,"mass",mass);
	H5PartReadDataInt64(H5PartFileID,"imove",imove);
	H5PartReadDataInt64(H5PartFileID,"ifluid",ifluid);
	#ifdef HAVE_3D
	    z     = new double[F->n()];
	    vz    = new double[F->n()];
	    nz    = new double[F->n()];
	    dvz   = new double[F->n()];
	    H5PartReadDataFloat64(H5PartFileID,"z",z);
	    H5PartReadDataFloat64(H5PartFileID,"vz",vz);
	    H5PartReadDataFloat64(H5PartFileID,"nz",nz);
	    H5PartReadDataFloat64(H5PartFileID,"dvz",dvz);
	#endif

	for(i=0;i<F->n();i++)
	{
	    F->pos[i].x    = x[i];
	    F->pos[i].y    = y[i];
	    F->normal[i].x = nx[i];
	    F->normal[i].y = ny[i];
	    F->v[i].x      = vx[i];
	    F->v[i].y      = vy[i];
	    F->f[i].x      = dvx[i];
	    F->f[i].y      = dvy[i];
	    F->press[i]    = press[i];
	    F->dens[i]     = dens[i];
	    F->drdt[i]     = ddens[i];
	    F->hp[i]       = hp[i];
	    F->mass[i]     = mass[i];
	    F->imove[i]    = imove[i];
	    F->ifluid[i]   = ifluid[i];
	    #ifdef HAVE_3D
	        F->pos[i].z    = z[i];
	        F->normal[i].z = nz[i];
	        F->v[i].z      = vz[i];
	        F->f[i].z      = dvz[i];
	        F->pos[i].w    = 0.f;
	        F->normal[i].w = 0.f;
	        F->v[i].w      = 0.f;
	        F->f[i].w      = 0.f;
	    #endif
	}
	delete[] x; x=0;
	delete[] y; y=0;
	delete[] nx; nx=0;
	delete[] ny; ny=0;
	delete[] vx; vx=0;
	delete[] vy; vy=0;
	delete[] dvx; dvx=0;
	delete[] dvy; dvy=0;
	delete[] press; press=0;
	delete[] dens; dens=0;
	delete[] ddens; ddens=0;
	delete[] hp; hp=0;
	delete[] mass; mass=0;
	delete[] imove; imove=0;
	delete[] ifluid; ifluid=0;
	#ifdef HAVE_3D
	    delete[] z; z=0;
	    delete[] nz; nz=0;
	    delete[] vz; vz=0;
	    delete[] dvz; dvz=0;
	#endif

	H5PartCloseFile(H5PartFileID);
	return false;
}

#endif // HAVE_H5PART

}}  // namespaces
