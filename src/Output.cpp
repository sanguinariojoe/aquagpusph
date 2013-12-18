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
#include <limits.h>

// ----------------------------------------------------------------------------
// Include the h5part library
// ----------------------------------------------------------------------------
#ifdef HAVE_H5PART
	#include <H5Part.h>
	#include <hdf5.h>
#endif // HAVE_H5PART

// ----------------------------------------------------------------------------
// Include the vtk library
// ----------------------------------------------------------------------------
#ifdef HAVE_VTK
	#include <vtkVersion.h>
	// VTK smart pointer
	#include <vtkSmartPointer.h>
	// Data
	#include <vtkFloatArray.h>
	#include <vtkIntArray.h>
	#include <vtkUnsignedIntArray.h>
	#include <vtkPointData.h>
	#include <vtkPoints.h>
	#include <vtkPixel.h>
	#include <vtkLine.h>
	#include <vtkCellArray.h>
	#include <vtkCellData.h>
	// UnstructuredGrid
	#include <vtkXMLUnstructuredGridWriter.h>
	#include <vtkUnstructuredGrid.h>
	#include <vtkDataSetMapper.h>
#endif

// ----------------------------------------------------------------------------
// Include the main header
// ----------------------------------------------------------------------------
#include <Output.h>

// ----------------------------------------------------------------------------
// Include the Host-Server fluid transfer layer
// ----------------------------------------------------------------------------
#include <Fluid.h>

// ----------------------------------------------------------------------------
// Include the simulation time manager
// ----------------------------------------------------------------------------
#include <TimeManager.h>

// ----------------------------------------------------------------------------
// Include the calculation server
// ----------------------------------------------------------------------------
#include <CalcServer.h>

// ----------------------------------------------------------------------------
// Include the File manager header
// ----------------------------------------------------------------------------
#include <FileManager.h>

// ----------------------------------------------------------------------------
// Include the File manager header
// ----------------------------------------------------------------------------
#include <ScreenManager.h>

namespace Aqua{ namespace InputOutput{

bool output()
{
	ScreenManager *S = ScreenManager::singleton();
	FileManager *Files = FileManager::singleton();
	ProblemSetup *P = ProblemSetup::singleton();
	if(Files->isH5Part()){
	    #ifdef HAVE_H5PART
	        Output::printH5Part();
	    #else // HAVE_H5PART
	        S->addMessage(3, "(output): H5Part output called, but not supported.\n");
	        return true;
	    #endif // HAVE_H5PART
	}
	if(Files->isVTK()){
	    #ifdef HAVE_VTK
	        Output::printVTK();
	    #else // HAVE_VTK
	        S->addMessage(3, "(output): VTK output called, but not supported.\n");
	        return true;
	    #endif // HAVE_VTK
	}
	if(Files->isTecplot()){
	    #ifdef HAVE_TECPLOT
	        Output::printTecplot();
	    #else // HAVE_TECPLOT
	        S->addMessage(3, "(output): Tecplot output called, but not supported.\n");
	        return true;
		#endif // HAVE_TECPLOT
	}
	return false;
}

/** Auxiliar function to config a name for output files. The name
 * will has the same number of digits if n <= N.
 * @param fileName.
 * @param n Index of the output file.
 * @param N Maximum index of the output file considered.
 * @param pre Name prefix (tec for Tecplot and vtk for VTK).
 * @param ext Extension of file (dat for Tecplot and vtu for VTK).
 */
void _createOutputName(char* fileName, unsigned int n, unsigned int N, const char *pre, const char *ext)
{
	strcpy(fileName, pre);
	strcat(fileName, ".");
	while(n < N){
	    strcat(fileName, "0");
	    N = N/10;
	}
	sprintf(fileName, "%s%u", fileName, n);
	strcat(fileName, ".");
	strcat(fileName, ext);
}

bool Output::printTecplot()
{
	// Variables
	unsigned int i;
	ScreenManager *S = ScreenManager::singleton();
	char file_name[256];
	FILE *FileID;
	float vv;
	Fluid *F = Fluid::singleton();
	TimeManager *T = TimeManager::singleton();
	//! 1st: Take a name for file
	_createOutputName((char*)file_name, T->frame(), UINT_MAX/10, "tec", "dat");
	//! 2nd: Open file
	FileID = fopen( file_name, "w" );
	if(!FileID){
        S->addMessage(3, "(Output::printTecplot): Can't open output file.\n");
        return true;
	}
	//! 3rd: Writing file head
	fprintf(FileID,"TITLE = \"2D SPH\"\n");
	#ifndef HAVE_3D
	    fprintf(FileID,"VARIABLES = \"X\" \"Y\" \"Vx\" \"Vy\" \"VV\" \"dvx\" \"dvy\" \"press\" \"dens\" \"ddens\"\n");
	#else // HAVE_3D
	    fprintf(FileID,"VARIABLES = \"X\" \"Y\" \"Z\" \"Vx\" \"Vy\" \"Vz\" \"VV\" \"dvx\" \"dvy\" \"dvz\" \"press\" \"dens\" \"ddens\"\n");
	#endif // HAVE_3D
	//! 4th: Writing fluid particles
	if(F->n()>0) {
	    fprintf(FileID,"ZONE T = \"FP\"\n");
	    fprintf(FileID,"I=%d J=1 K=1 ZONETYPE=Ordered\n",F->n());
	    fprintf(FileID,"DATAPACKING = POINT \n");
	    for(i=0;i<F->n();i++) {
	        #ifndef HAVE_3D
	            vv  = sqrt(pow(F->v[i].x,2.f) + pow(F->v[i].y,2.f));
	            fprintf(FileID, "%f %f %f %f %f %f %f %f %f %f\n",
	                    F->pos[i].x,F->pos[i].y,F->v[i].x,F->v[i].y,vv,F->f[i].x,F->f[i].y,F->press[i],F->dens[i],F->drdt[i]);
	        #else // HAVE_3D
	            vv  = sqrt(pow(F->v[i].x,2.f) + pow(F->v[i].y,2.f) + pow(F->v[i].z,2.f));
	            fprintf(FileID, "%f %f %f %f %f %f %f %f %f %f %f %f %f\n",
	                    F->pos[i].x,F->pos[i].y,F->pos[i].z,F->v[i].x,F->v[i].y,F->v[i].z,vv,F->f[i].x,F->f[i].y,F->f[i].z,F->press[i],F->dens[i],F->drdt[i]);
	        #endif // HAVE_3D
	    }
	}
	//! 5th: Close file
	fclose(FileID);
	return false;
}

#ifdef HAVE_H5PART
static int firstPrint = 1;

bool Output::printH5Part()
{
	CalcServer::CalcServer *C = CalcServer::CalcServer::singleton();
	Fluid *F = Fluid::singleton();
	TimeManager *T = TimeManager::singleton();
	FileManager *Files = FileManager::singleton();
	unsigned int i;
	h5part_int64_t n, num_fluids, iStep;
	h5part_int64_t *imove, *ifluid, *id;
	double dTime, ddt;
	double *x, *y, *z, *nx, *ny, *vx, *vy, *vv, *dvx, *dvy, *press, *dens, *ddens, *hp, *mass, *sumW, *gradWx, *gradWy;
	#ifdef HAVE_3D
	    double *nz, *vz, *dvz, *gradWz;
	#endif // HAVE_3D
	//! 1st: Casting operation
	x      = new double[F->n()];
	y      = new double[F->n()];
	z      = new double[F->n()];
	nx     = new double[F->n()];
	ny     = new double[F->n()];
	vx     = new double[F->n()];
	vy     = new double[F->n()];
	vv     = new double[F->n()];
	dvx    = new double[F->n()];
	dvy    = new double[F->n()];
	press  = new double[F->n()];
	dens   = new double[F->n()];
	ddens  = new double[F->n()];
	hp     = new double[F->n()];
	mass   = new double[F->n()];
	sumW   = new double[F->n()];
	gradWx = new double[F->n()];
	gradWy = new double[F->n()];
	imove  = new h5part_int64_t[F->n()];
	ifluid = new h5part_int64_t[F->n()];
	id     = new h5part_int64_t[F->n()];
	#ifdef HAVE_3D
	    nz  = new double[F->n()];
	    vz  = new double[F->n()];
	    dvz = new double[F->n()];
	    gradWz = new double[F->n()];
	#endif // HAVE_3D
	for(i=0;i<F->n();i++)
	{
	    x[i]     = F->pos[i].x;
	    y[i]     = F->pos[i].y;
	    z[i]     = 0.f;
	    nx[i]    = F->normal[i].x;
	    ny[i]    = F->normal[i].y;
	    vx[i]    = F->v[i].x;
	    vy[i]    = F->v[i].y;
	    vv[i]    = sqrt(pow(F->v[i].x,2.0f) + pow(F->v[i].y,2.0f));
	    dvx[i]   = F->f[i].x;
	    dvy[i]   = F->f[i].y;
	    press[i] = F->press[i];
	    dens[i]  = F->dens[i];
	    ddens[i] = F->drdt[i];
	    hp[i]    = F->hp[i];
	    mass[i]  = F->mass[i];
	    sumW[i]  = F->shepard[i];
	    gradWx[i]= F->shepard_gradient[i].x;
	    gradWy[i]= F->shepard_gradient[i].y;
	    imove[i] = F->imove[i];
	    ifluid[i] = F->ifluid[i];
	    id[i]    = i;
	    #ifdef HAVE_3D
	        z[i]     = F->pos[i].z;
	        nz[i]    = F->normal[i].z;
	        vz[i]    = F->v[i].z;
	        vv[i]    = sqrt(pow(F->v[i].x,2.f) + pow(F->v[i].y,2.f) + pow(F->v[i].z,2.f));
	        dvz[i]   = F->f[i].z;
	        gradWz[i]= F->shepard_gradient[i].z;
	    #endif // HAVE_3D
	}
	n = F->n();
	num_fluids = C->num_fluids;
	iStep = T->step();
	dTime = T->time();
	ddt = T->dt();
	if(firstPrint){
	    firstPrint = 0;
	    h5part_int64_t zeroFrame = T->startFrame();
	    H5PartWriteFileAttrib(Files->h5File(), "zeroFrame", H5T_NATIVE_INT64, &zeroFrame,1);
	}
	//! 2st: Writing step parameters
	H5PartSetStep(Files->h5File(),T->frame()-1-T->startFrame()); // Set of step
	H5PartSetNumParticles(Files->h5File(),n); // Set number of particles
	H5PartWriteStepAttrib(Files->h5File(),"nfluid",H5T_NATIVE_INT64,&num_fluids,1);
	H5PartWriteStepAttrib(Files->h5File(),"step",H5T_NATIVE_INT64,&iStep,1);
	H5PartWriteStepAttrib(Files->h5File(),"time",H5T_NATIVE_DOUBLE,&dTime,1);
	H5PartWriteStepAttrib(Files->h5File(),"dt",H5T_NATIVE_DOUBLE,&ddt,1);
	//! 3rd: Writing datasets
	H5PartWriteDataFloat64(Files->h5File(),"x",x);
	H5PartWriteDataFloat64(Files->h5File(),"y",y);
	H5PartWriteDataFloat64(Files->h5File(),"z",z);
	if(Files->isField("nx"))
	    H5PartWriteDataFloat64(Files->h5File(),"nx",nx);
	if(Files->isField("ny"))
	    H5PartWriteDataFloat64(Files->h5File(),"ny",ny);
	#ifdef HAVE_3D
	    if(Files->isField("nz"))
	        H5PartWriteDataFloat64(Files->h5File(),"nz",nz);
	#endif // HAVE_3D
	if(Files->isField("vx"))
	    H5PartWriteDataFloat64(Files->h5File(),"vx",vx);
	if(Files->isField("vy"))
	    H5PartWriteDataFloat64(Files->h5File(),"vy",vy);
	#ifdef HAVE_3D
	    if(Files->isField("vz"))
	        H5PartWriteDataFloat64(Files->h5File(),"vz",vz);
	#endif // HAVE_3D
	if(Files->isField("vv"))
	    H5PartWriteDataFloat64(Files->h5File(),"vv",vv);
	if(Files->isField("dvx"))
	    H5PartWriteDataFloat64(Files->h5File(),"dvx",dvx);
	if(Files->isField("dvy"))
	    H5PartWriteDataFloat64(Files->h5File(),"dvy",dvy);
	#ifdef HAVE_3D
	    if(Files->isField("dvz"))
	        H5PartWriteDataFloat64(Files->h5File(),"dvz",dvz);
	#endif // HAVE_3D
	if(Files->isField("press"))
	    H5PartWriteDataFloat64(Files->h5File(),"press",press);
	if(Files->isField("dens"))
	    H5PartWriteDataFloat64(Files->h5File(),"dens",dens);
	if(Files->isField("ddens"))
	    H5PartWriteDataFloat64(Files->h5File(),"ddens",ddens);
	if(Files->isField("hp"))
	    H5PartWriteDataFloat64(Files->h5File(),"hp",hp);
	if(Files->isField("mass"))
	    H5PartWriteDataFloat64(Files->h5File(),"mass",mass);
	if(Files->isField("sumW"))
	    H5PartWriteDataFloat64(Files->h5File(),"sumW",sumW);
	if(Files->isField("gradWx"))
	    H5PartWriteDataFloat64(Files->h5File(),"gradWx",gradWx);
	if(Files->isField("gradWy"))
	    H5PartWriteDataFloat64(Files->h5File(),"gradWy",gradWy);
	#ifdef HAVE_3D
	    if(Files->isField("gradWz"))
	        H5PartWriteDataFloat64(Files->h5File(),"gradWz",gradWz);
	#endif // HAVE_3D
	if(Files->isField("imove"))
	    H5PartWriteDataInt64(Files->h5File(),"imove",imove);
	if(Files->isField("ifluid"))
	    H5PartWriteDataInt64(Files->h5File(),"ifluid",ifluid);
	if(Files->isField("id"))
	    H5PartWriteDataInt64(Files->h5File(),"id",id);
	//! 4th.- Free memory
	delete[] x;
	delete[] y;
	delete[] z;
	delete[] nx;
	delete[] ny;
	delete[] vx;
	delete[] vy;
	delete[] vv;
	delete[] dvx;
	delete[] dvy;
	delete[] press;
	delete[] dens;
	delete[] ddens;
	delete[] hp;
	delete[] mass;
	delete[] sumW;
	delete[] gradWx;
	delete[] gradWy;
	delete[] imove;
	delete[] ifluid;
	delete[] id;
	#ifdef HAVE_3D
	    delete[] nz;
	    delete[] vz;
	    delete[] dvz;
	    delete[] gradWz;
	#endif // HAVE_3D
	return false;
}

bool assemblyH5Part()
{
	InputOutput::FileManager *F = InputOutput::FileManager::singleton();
	ScreenManager *S = ScreenManager::singleton();
	int i,j, nFiles=0, iFrame=0;
	int dim=0;
	char file_name[256], msg[256];
	H5PartFile *Input, *Output;
	h5part_int64_t n, num_fluids, iStep;
	h5part_int64_t *imove=0, *ifluid=0;
	double dTime, ddt;
	double *x=0, *y=0, *z=0, *nx=0, *ny=0, *vx=0, *vy=0, *vv=0, *dvx=0, *dvy=0,
	       *press=0, *dens=0, *ddens=0, *hp=0, *mass=0;
	#ifdef HAVE_3D
	double *nz=0, *vz=0, *dvz=0;
	#endif // HAVE_3D

    S->addMessage(1, "(assemblyH5Part): Reassembling H5Part data files.\n");
	//! 1st.- Find valid files
	sprintf(file_name, "%sParticles%d.h5part", F->outputPrefix(), nFiles);
	while(isFile(file_name)){
	    nFiles++;
	    sprintf(file_name, "%sParticles%d.h5part", F->outputPrefix(), nFiles);
	}
	if(nFiles == 0){
        sprintf(msg, "(assemblyH5Part): Any valid file found!\n");
        S->addMessage(3, msg);
	    return true;
	}
	if(nFiles == 1){
        sprintf(msg, "(assemblyH5Part): Only 1 file found, reassembly is don't needed therefore\n");
        S->addMessage(1, msg);
	    return false;
	}
	//! 2nd.- Open output
	sprintf(msg, "(assemblyH5Part): Saving into %sParticles.h5part.\n", F->outputPrefix());
    S->addMessage(1, msg);
	sprintf(file_name, "%sParticles.h5part", F->outputPrefix());
	Output = H5PartOpenFileParallel(file_name,H5PART_WRITE,MPI_COMM_WORLD);
	H5PartWriteFileAttribString(Output, "rmsUnit", "m");
	h5part_int64_t zeroFrame = 0;
	H5PartWriteFileAttrib(Output, "zeroFrame", H5T_NATIVE_INT64, &zeroFrame,1);
	//! 3rd.- Reading files
	S->addMessage(1, "(assemblyH5Part): Transfering data.\n");
	for(i=0;i<nFiles;i++){
	    //! \t 3.a.- Open input file
        sprintf(msg, "\tFile %d/%d...\n", i+1, nFiles);
        S->addMessage(0, msg);
	    sprintf(file_name, "%sParticles%d.h5part", F->outputPrefix(), i);
	    #ifdef HAVE_MPI
	        Input = H5PartOpenFileParallel(file_name,H5PART_READ,MPI_COMM_WORLD);
	    #else // HAVE_MPI
	        Input = H5PartOpenFile(file_name,H5PART_READ);
	    #endif // HAVE_MPI
	    //! \t 3.b.- Reads basic data
	    int nFrames = H5PartGetNumSteps(Input);
	    //! \t 3.c.- Reads all the time steps
	    for(j=0;j<nFrames;j++){
	        // Set output step
	        H5PartSetStep(Output,iFrame);
	        iFrame++;
	        // Set input step
	        H5PartSetStep(Input,j);
	        // Reads the number of particles
	        n = H5PartGetNumParticles(Input);
	        H5PartSetNumParticles(Output,n);
	        // Reads step, time, and dt
	        H5PartReadStepAttrib(Input,"step",&iStep);
	        H5PartReadStepAttrib(Input,"time",&dTime);
	        H5PartReadStepAttrib(Input,"dt",&ddt);
	        H5PartWriteStepAttrib(Output,"step",H5T_NATIVE_INT64,&iStep,1);
	        H5PartWriteStepAttrib(Output,"time",H5T_NATIVE_DOUBLE,&dTime,1);
	        H5PartWriteStepAttrib(Output,"dt",H5T_NATIVE_DOUBLE,&ddt,1);
	        // Read fluid info
	        H5PartReadStepAttrib(Input,"nfluid",&num_fluids);
	        H5PartWriteStepAttrib(Output,"nfluid",H5T_NATIVE_INT64,&num_fluids,1);
	        // Allocate memory
	        if(dim != n){
	            delete[] x; x=0;
	            delete[] y; y=0;
	            delete[] z; z=0;
	            delete[] nx; nx=0;
	            delete[] ny; ny=0;
	            delete[] vx; vx=0;
	            delete[] vy; vy=0;
	            delete[] vv; vv=0;
	            delete[] dvx; dvx=0;
	            delete[] dvy; dvy=0;
	            delete[] press; press=0;
	            delete[] dens; dens=0;
	            delete[] ddens; ddens=0;
	            delete[] hp; hp=0;
	            delete[] mass; mass=0;
	            delete[] imove; imove=0;
	            delete[] ifluid; ifluid=0;
	            x     = new double[n];
	            y     = new double[n];
	            z     = new double[n];
	            nx    = new double[n];
	            ny    = new double[n];
	            vx    = new double[n];
	            vy    = new double[n];
	            vv    = new double[n];
	            dvx   = new double[n];
	            dvy   = new double[n];
	            press = new double[n];
	            dens  = new double[n];
	            ddens = new double[n];
	            hp    = new double[n];
	            mass  = new double[n];
	            imove = new h5part_int64_t[n];
	            ifluid = new h5part_int64_t[n];
	            dim = n;
	            #ifdef HAVE_3D
	                delete[] vz; vz=0;
	                delete[] nz; nz=0;
	                delete[] dvz; dvz=0;
	                nz  = new double[n];
	                vz  = new double[n];
	                dvz = new double[n];
	            #endif // HAVE_3D
	        }
	        // Reads Datasheets
	        H5PartReadDataFloat64(Input,"x",x);
	        H5PartReadDataFloat64(Input,"y",y);
	        H5PartReadDataFloat64(Input,"z",z);
	        H5PartReadDataFloat64(Input,"nx",nx);
	        H5PartReadDataFloat64(Input,"ny",ny);
	        H5PartReadDataFloat64(Input,"vx",vx);
	        H5PartReadDataFloat64(Input,"vy",vy);
	        H5PartReadDataFloat64(Input,"vv",vv);
	        H5PartReadDataFloat64(Input,"dvx",dvx);
	        H5PartReadDataFloat64(Input,"dvy",dvy);
	        H5PartReadDataFloat64(Input,"press",press);
	        H5PartReadDataFloat64(Input,"dens",dens);
	        H5PartReadDataFloat64(Input,"ddens",ddens);
	        H5PartReadDataFloat64(Input,"hp",hp);
	        H5PartReadDataFloat64(Input,"mass",mass);
	        H5PartReadDataInt64(Input,"imove",imove);
	        H5PartReadDataInt64(Input,"ifluid",ifluid);
	        H5PartWriteDataFloat64(Output,"x",x);
	        H5PartWriteDataFloat64(Output,"y",y);
	        H5PartWriteDataFloat64(Output,"z",z);
	        H5PartWriteDataFloat64(Output,"nx",nx);
	        H5PartWriteDataFloat64(Output,"ny",ny);
	        H5PartWriteDataFloat64(Output,"vx",vx);
	        H5PartWriteDataFloat64(Output,"vy",vy);
	        H5PartWriteDataFloat64(Output,"vv",vv);
	        H5PartWriteDataFloat64(Output,"dvx",dvx);
	        H5PartWriteDataFloat64(Output,"dvy",dvy);
	        H5PartWriteDataFloat64(Output,"press",press);
	        H5PartWriteDataFloat64(Output,"dens",dens);
	        H5PartWriteDataFloat64(Output,"ddens",ddens);
	        H5PartWriteDataFloat64(Output,"hp",hp);
	        H5PartWriteDataFloat64(Output,"mass",mass);
	        H5PartWriteDataInt64(Output,"imove",imove);
	        H5PartWriteDataInt64(Output,"ifluid",ifluid);
	        #ifdef HAVE_3D
	            H5PartReadDataFloat64(Input,"nz",nz);
	            H5PartReadDataFloat64(Input,"vz",vz);
	            H5PartReadDataFloat64(Input,"dvz",dvz);
	            H5PartWriteDataFloat64(Output,"nz",nz);
	            H5PartWriteDataFloat64(Output,"vz",vz);
	            H5PartWriteDataFloat64(Output,"dvz",dvz);
	        #endif // HAVE_3D
            sprintf(msg, "\t\t%d%% (%d%%)\n", (int)(100.f*((float)(i) + (j+1.f)/nFrames)/nFiles), (int)(100.f*(j+1.f)/nFrames));
            S->addMessage(0, msg);
	    }
	    //! \t 3.d.- Close input file
	    H5PartCloseFile(Input);
	}
	//! 4th.- Close output file and exit
	H5PartCloseFile(Output);
    sprintf(msg, "(assemblyH5Part): %sParticles.h5part written!\n", F->outputPrefix());
    S->addMessage(1, msg);
	delete[] x; x=0;
	delete[] y; y=0;
	delete[] z; z=0;
	delete[] nx; nx=0;
	delete[] ny; ny=0;
	delete[] vx; vx=0;
	delete[] vy; vy=0;
	delete[] vv; vv=0;
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
	    delete[] vz; vz=0;
	    delete[] nz; nz=0;
	    delete[] dvz; dvz=0;
	#endif // HAVE_3D
	return false;
}
#endif // HAVE_H5PART

#ifdef HAVE_VTK
	bool _printVTK(vtkSmartPointer<vtkXMLUnstructuredGridWriter> fd)
	{
	    unsigned int i;
	    Fluid *F = Fluid::singleton();
	    FileManager *Files = FileManager::singleton();
	    vtkSmartPointer<vtkPixel> pixel;
	    float vect[3] = {0.f, 0.f, 0.f};
	    // Prepare points arrays
	    vtkSmartPointer<vtkPoints> points = vtkSmartPointer<vtkPoints>::New();
	    vtkSmartPointer<vtkFloatArray> normal = vtkSmartPointer<vtkFloatArray>::New();
	    normal->SetNumberOfComponents(3);
	    normal->SetName("n");
	    vtkSmartPointer<vtkFloatArray> v = vtkSmartPointer<vtkFloatArray>::New();
	    v->SetNumberOfComponents(3);
	    v->SetName("v");
	    vtkSmartPointer<vtkFloatArray> dv = vtkSmartPointer<vtkFloatArray>::New();
	    dv->SetNumberOfComponents(3);
	    dv->SetName("dv/dt");
	    vtkSmartPointer<vtkFloatArray> p = vtkSmartPointer<vtkFloatArray>::New();
	    p->SetNumberOfComponents(1);
	    p->SetName("press");
	    vtkSmartPointer<vtkFloatArray> dens = vtkSmartPointer<vtkFloatArray>::New();
	    dens->SetNumberOfComponents(1);
	    dens->SetName("dens");
	    vtkSmartPointer<vtkFloatArray> drdt = vtkSmartPointer<vtkFloatArray>::New();
	    drdt->SetNumberOfComponents(1);
	    drdt->SetName("ddens/dt");
	    vtkSmartPointer<vtkFloatArray> h = vtkSmartPointer<vtkFloatArray>::New();
	    h->SetNumberOfComponents(1);
	    h->SetName("h");
	    vtkSmartPointer<vtkFloatArray> m = vtkSmartPointer<vtkFloatArray>::New();
	    m->SetNumberOfComponents(1);
	    m->SetName("mass");
	    vtkSmartPointer<vtkFloatArray> W = vtkSmartPointer<vtkFloatArray>::New();
	    W->SetNumberOfComponents(1);
	    W->SetName("sumW");
	    vtkSmartPointer<vtkFloatArray> dW = vtkSmartPointer<vtkFloatArray>::New();
	    dW->SetNumberOfComponents(3);
	    dW->SetName("gradW");
	    vtkSmartPointer<vtkIntArray> imove = vtkSmartPointer<vtkIntArray>::New();
	    imove->SetNumberOfComponents(1);
	    imove->SetName("imove");
	    vtkSmartPointer<vtkIntArray> ifluid = vtkSmartPointer<vtkIntArray>::New();
	    ifluid->SetNumberOfComponents(1);
	    ifluid->SetName("ifluid");
	    vtkSmartPointer<vtkUnsignedIntArray> id = vtkSmartPointer<vtkUnsignedIntArray>::New();
	    id->SetNumberOfComponents(1);
	    id->SetName("id");
	    vtkSmartPointer<vtkCellArray> cells = vtkSmartPointer<vtkCellArray>::New();
	    // Set data
	    for(i=0;i<F->n();i++){
	        #ifdef HAVE_3D
	            points->InsertNextPoint(F->pos[i].x, F->pos[i].y, F->pos[i].z);
	            vect[0] = F->normal[i].x; vect[1] = F->normal[i].y; vect[2] = F->normal[i].z;
	            normal->InsertNextTupleValue(vect);
	            vect[0] = F->v[i].x; vect[1] = F->v[i].y; vect[2] = F->v[i].z;
	            v->InsertNextTupleValue(vect);
	            vect[0] = F->f[i].x; vect[1] = F->f[i].y; vect[2] = F->f[i].z;
	            dv->InsertNextTupleValue(vect);
	            vect[0] = F->shepard_gradient[i].x; vect[1] = F->shepard_gradient[i].y; vect[2] = F->shepard_gradient[i].z;
	            dW->InsertNextTupleValue(vect);
	        #else // HAVE_3D
	            points->InsertNextPoint(F->pos[i].x, F->pos[i].y, 0.f);
	            vect[0] = F->normal[i].x; vect[1] = F->normal[i].y;
	            normal->InsertNextTupleValue(vect);
	            vect[0] = F->v[i].x; vect[1] = F->v[i].y;
	            v->InsertNextTupleValue(vect);
	            vect[0] = F->f[i].x; vect[1] = F->f[i].y;
	            dv->InsertNextTupleValue(vect);
	            vect[0] = F->shepard_gradient[i].x; vect[1] = F->shepard_gradient[i].y;
	            dW->InsertNextTupleValue(vect);
	        #endif // HAVE_3D
	        p->InsertNextValue(F->press[i]);
	        dens->InsertNextValue(F->dens[i]);
	        drdt->InsertNextValue(F->drdt[i]);
	        h->InsertNextValue(F->hp[i]);
	        m->InsertNextValue(F->mass[i]);
	        W->InsertNextValue(F->shepard[i]);
	        imove->InsertNextValue(F->imove[i]);
	        ifluid->InsertNextValue(F->ifluid[i]);
	        id->InsertNextValue(i);
	        // Store the data into the cells
	        pixel = vtkSmartPointer<vtkPixel>::New();
	        pixel->GetPointIds()->SetId(0, i);
	        cells->InsertNextCell(pixel);
	    }
	    // Setup the unstructured grid
	    vtkSmartPointer<vtkUnstructuredGrid> unstructuredGrid = vtkSmartPointer<vtkUnstructuredGrid>::New();
	    unstructuredGrid->SetPoints(points);
	    unstructuredGrid->SetCells(pixel->GetCellType(), cells);
	    if(Files->isField("n"))
	        unstructuredGrid->GetCellData()->AddArray(normal);
	    if(Files->isField("v"))
	        unstructuredGrid->GetCellData()->AddArray(v);
	    if(Files->isField("dv/dt"))
	        unstructuredGrid->GetCellData()->AddArray(dv);
	    if(Files->isField("press"))
	        unstructuredGrid->GetCellData()->AddArray(p);
	    if(Files->isField("dens"))
	        unstructuredGrid->GetCellData()->AddArray(dens);
	    if(Files->isField("ddens/dt"))
	        unstructuredGrid->GetCellData()->AddArray(drdt);
	    if(Files->isField("h"))
	        unstructuredGrid->GetCellData()->AddArray(h);
	    if(Files->isField("mass"))
	        unstructuredGrid->GetCellData()->AddArray(m);
	    if(Files->isField("sumW"))
	        unstructuredGrid->GetCellData()->AddArray(W);
	    if(Files->isField("gradW"))
	        unstructuredGrid->GetCellData()->AddArray(dW);
	    if(Files->isField("imove"))
	        unstructuredGrid->GetCellData()->AddArray(imove);
	    if(Files->isField("ifluid"))
	        unstructuredGrid->GetCellData()->AddArray(ifluid);
	    if(Files->isField("id"))
	        unstructuredGrid->GetCellData()->AddArray(id);
	    // Write file
	    #if VTK_MAJOR_VERSION <= 5
	        fd->SetInput(unstructuredGrid);
	    #else // VTK_MAJOR_VERSION
	        fd->SetInputData(unstructuredGrid);
	    #endif // VTK_MAJOR_VERSION
	    fd->Write();
	}

	bool Output::printVTK()
	{
	    char fileName[256];
	    TimeManager *T = TimeManager::singleton();
	    _createOutputName((char*)fileName, T->frame(), UINT_MAX/10, "vtk", "vtu");
	    vtkSmartPointer<vtkXMLUnstructuredGridWriter> fd = vtkSmartPointer<vtkXMLUnstructuredGridWriter>::New();
	    fd->SetFileName(fileName);
	    return _printVTK(fd);
	}
#endif // HAVE_VTK

}}  //namespace
