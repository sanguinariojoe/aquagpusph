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
#include <CalcServer/Movements/Quaternion.h>

// ----------------------------------------------------------------------------
// Include the calculation server
// ----------------------------------------------------------------------------
#include <CalcServer.h>

// ----------------------------------------------------------------------------
// Include the Time manager
// ----------------------------------------------------------------------------
#include <TimeManager.h>

// ----------------------------------------------------------------------------
// Include the Screen manager
// ----------------------------------------------------------------------------
#include <ScreenManager.h>

#ifdef xmlAttribute
	#undef xmlAttribute
#endif
#define xmlAttribute(elem, att) XMLString::transcode( elem->getAttribute(XMLString::transcode(att)) )
using namespace xercesc;

namespace Aqua{ namespace CalcServer{ namespace Movement{

Quaternion::Quaternion()
	: Movement()
	, mRelPos(0)
	, mRelNormal(0)
{
	CalcServer *C = CalcServer::singleton();
	//! Set inital values
	mCOR.x = 0.f;
	mCOR.y = 0.f;
	#ifdef HAVE_3D
	    mCOR.z = 0.f;
	    mCOR.w = 0.f;
	#endif
	mAxis[0].x = 1.f;
	mAxis[1].y = 1.f;
	#ifdef HAVE_3D
	    mAxis[2].z = 1.f;
	#endif
	//! Allocate memory at server
	unsigned int N = C->N;
	mRelPos = C->allocMemory(N*sizeof(vec));
	if(!mRelPos)
	    exit(255);
    mRelNormal = C->allocMemory(N*sizeof(vec));
	if(!mRelNormal)
	    exit(255);
	//! Compute initial relative positions
	if(set(mCOR,mAxis,true))
	    exit(255);
}

Quaternion::~Quaternion()
{
	unsigned int i;
	if(mRelPos)clReleaseMemObject(mRelPos); mRelPos=0;
	if(mRelNormal)clReleaseMemObject(mRelNormal); mRelNormal=0;
	for(i=0;i<walls.size();i++){
	    delete walls.at(i);
	}
	walls.clear();
}

bool Quaternion::set(vec cor, mat axis, bool initial){
	mCOR  = cor;
	mAxis = axis;
	if(initial){
	    mOldCOR  = cor;
	    mOldAxis = axis;
	    if(computePos())
	        return true;
	    if(computeWalls())
	        return true;
	    if(computeDomain())
	        return true;
	}
	return false;
}

bool Quaternion::execute()
{
	InputOutput::ScreenManager *S = InputOutput::ScreenManager::singleton();
	CalcServer *C = CalcServer::singleton();
	vec X, Y, Z, oldX, oldY, oldZ;
	X    = mAxis[0];
	Y    = mAxis[1];
	oldX = mOldAxis[0];
	oldY = mOldAxis[1];
	#ifdef HAVE_3D
	    Z    = mAxis[2];
	    oldZ = mOldAxis[2];
	#else
	    Z.x = 0.f; Z.y = 0.f;
	    oldZ.x = 0.f; oldZ.y = 0.f;
	#endif
	//! Send variables to kernel
	cl_int err_code=0;
	err_code |= sendArgument(_kernel,  0, sizeof(cl_mem  ), (void*)&(C->imove));
	err_code |= sendArgument(_kernel,  1, sizeof(cl_mem  ), (void*)&(C->ifluid));
	err_code |= sendArgument(_kernel,  2, sizeof(cl_mem  ), (void*)&(C->posin));
	err_code |= sendArgument(_kernel,  3, sizeof(cl_mem  ), (void*)&(C->normal));
	err_code |= sendArgument(_kernel,  4, sizeof(cl_mem  ), (void*)&(C->vin));
	err_code |= sendArgument(_kernel,  5, sizeof(cl_mem  ), (void*)&(C->densin));
	err_code |= sendArgument(_kernel,  6, sizeof(cl_mem  ), (void*)&(C->mass));
	err_code |= sendArgument(_kernel,  7, sizeof(cl_mem  ), (void*)&(C->hpin));
	err_code |= sendArgument(_kernel,  8, sizeof(cl_mem  ), (void*)&(mRelPos));
	err_code |= sendArgument(_kernel,  9, sizeof(cl_mem  ), (void*)&(mRelNormal));
	err_code |= sendArgument(_kernel, 10, sizeof(cl_uint ), (void*)&(C->N));
	err_code |= sendArgument(_kernel, 11, sizeof(cl_float), (void*)&(C->dt));
	err_code |= sendArgument(_kernel, 12, sizeof(vec     ), (void*)&(mCOR));
	err_code |= sendArgument(_kernel, 13, sizeof(vec     ), (void*)&(X));
	err_code |= sendArgument(_kernel, 14, sizeof(vec     ), (void*)&(Y));
	err_code |= sendArgument(_kernel, 15, sizeof(vec     ), (void*)&(Z));
	err_code |= sendArgument(_kernel, 16, sizeof(vec     ), (void*)&(mOldCOR));
	err_code |= sendArgument(_kernel, 17, sizeof(vec     ), (void*)&(oldX));
	err_code |= sendArgument(_kernel, 18, sizeof(vec     ), (void*)&(oldY));
	err_code |= sendArgument(_kernel, 19, sizeof(vec     ), (void*)&(oldZ));
	if(err_code != CL_SUCCESS) {
		S->addMessage(3, "(Quaternion::execute): Failure sending variables to _kernel.\n");
	    return true;
	}
	//! Execute the kernel
	#ifdef HAVE_GPUPROFILE
	    cl_event event;
	    cl_ulong end, start;
	    profileTime(0.f);
	    err_code = clEnqueueNDRangeKernel(C->command_queue, _kernel, 1, NULL, &_global_work_size, NULL, 0, NULL, &event);
	#else
	    err_code = clEnqueueNDRangeKernel(C->command_queue, _kernel, 1, NULL, &_global_work_size, NULL, 0, NULL, NULL);
	#endif
	if(err_code != CL_SUCCESS) {
		S->addMessage(3, "(Rates::Execute): I cannot execute the kernel.\n");
	    if(err_code == CL_INVALID_WORK_GROUP_SIZE)
	        S->addMessage(0, "\tInvalid local work group size.\n");
	    else if(err_code == CL_OUT_OF_RESOURCES)
	        S->addMessage(0, "\tDevice out of resources.\n");
	    else if(err_code == CL_MEM_OBJECT_ALLOCATION_FAILURE)
	        S->addMessage(0, "\tAllocation error at device.\n");
	    else if(err_code == CL_OUT_OF_HOST_MEMORY)
	        S->addMessage(0, "\tfailure to allocate resources required by the OpenCL implementation on the host.\n");
	    return true;
	}
	// Profile the kernel execution
	#ifdef HAVE_GPUPROFILE
	    err_code = clWaitForEvents(1, &event);
	    if(err_code != CL_SUCCESS) {
	        S->addMessage(3, "(Rates::Execute): Impossible to wait for the sorting kernel ends.\n");
	        return true;
	    }
	    err_code |= clGetEventProfilingInfo(event, CL_PROFILING_COMMAND_END, sizeof(cl_ulong), &end, 0);
	    err_code |= clGetEventProfilingInfo(event, CL_PROFILING_COMMAND_START, sizeof(cl_ulong), &start, 0);
	    if(err_code != CL_SUCCESS) {
	        S->addMessage(3, "(Rates::Execute): I cannot profile the sorting kernel execution.\n");
	        return true;
	    }
	    profileTime(profileTime() + (end - start)/1000.f);  // 10^-3 ms
	#endif
	//! Call to move walls (ghost particles)
	if(executeWalls())
	    return true;
	if(executeDomain())
	    return true;
	//! Backup the quaternion.
	mOldCOR  = mCOR;
	mOldAxis = mAxis;
	return false;
}

bool Quaternion::executeWalls()
{
	unsigned int i;
	unsigned int N = walls.size();
	if(!N)
	    return false;
	InputOutput::ProblemSetup *P = InputOutput::ProblemSetup::singleton();
	CalcServer *C = CalcServer::singleton();
	vec X, Y, Z, oldX, oldY, oldZ;
	X    = mAxis[0];
	Y    = mAxis[1];
	oldX = mOldAxis[0];
	oldY = mOldAxis[1];
	#ifdef HAVE_3D
	    Z    = mAxis[2];
	    oldZ = mOldAxis[2];
	#else
	    Z.x = 0.f; Z.y = 0.f;
	    oldZ.x = 0.f; oldZ.y = 0.f;
	#endif
	for(i=0;i<N;i++){
	    vec newPos, oldPos;
	    newPos = mCOR;
	    newPos = add(newPos, mult(walls.at(i)->p1.x,X));
	    newPos = add(newPos, mult(walls.at(i)->p1.y,Y));
	    #ifdef HAVE_3D
	        newPos = add(newPos, mult(walls.at(i)->p1.z,Z));
	    #endif
	    oldPos = mOldCOR;
	    oldPos = add(oldPos, mult(walls.at(i)->p1.x,oldX));
	    oldPos = add(oldPos, mult(walls.at(i)->p1.y,oldY));
	    #ifdef HAVE_3D
	        oldPos = add(oldPos, mult(walls.at(i)->p1.z,oldZ));
	    #endif
	    P->ghost_particles.walls.at(i)->p1 = newPos;
	    if(C->dt <= 0.f){
	        P->ghost_particles.walls.at(i)->v1 = Vzero();
	    }
	    else{
	        P->ghost_particles.walls.at(i)->v1 = mult(1.f/C->dt, sub(newPos, oldPos));
	    }
	    newPos = mCOR;
	    newPos = add(newPos, mult(walls.at(i)->p2.x,X));
	    newPos = add(newPos, mult(walls.at(i)->p2.y,Y));
	    #ifdef HAVE_3D
	        newPos = add(newPos, mult(walls.at(i)->p2.z,Z));
	    #endif
	    oldPos = mOldCOR;
	    oldPos = add(oldPos, mult(walls.at(i)->p2.x,oldX));
	    oldPos = add(oldPos, mult(walls.at(i)->p2.y,oldY));
	    #ifdef HAVE_3D
	        oldPos = add(oldPos, mult(walls.at(i)->p2.z,oldZ));
	    #endif
	    P->ghost_particles.walls.at(i)->p2 = newPos;
	    if(C->dt <= 0.f){
	        P->ghost_particles.walls.at(i)->v2 = Vzero();
	    }
	    else{
	        P->ghost_particles.walls.at(i)->v2 = mult(1.f/C->dt, sub(newPos, oldPos));
	    }
	    #ifdef HAVE_3D
	        newPos = mCOR;
	        newPos = add(newPos, mult(walls.at(i)->p3.x,X));
	        newPos = add(newPos, mult(walls.at(i)->p3.y,Y));
	        newPos = add(newPos, mult(walls.at(i)->p3.z,Z));
	        oldPos = mOldCOR;
	        oldPos = add(oldPos, mult(walls.at(i)->p3.x,oldX));
	        oldPos = add(oldPos, mult(walls.at(i)->p3.y,oldY));
	        oldPos = add(oldPos, mult(walls.at(i)->p3.z,oldZ));
	        P->ghost_particles.walls.at(i)->p3 = newPos;
	        if(C->dt <= 0.f){
	            P->ghost_particles.walls.at(i)->v3 = Vzero();
	        }
	        else{
	            P->ghost_particles.walls.at(i)->v3 = mult(1.f/C->dt, sub(newPos, oldPos));
	        }
	    #endif
	    #ifdef HAVE_3D
	        newPos = mCOR;
	        newPos = add(newPos, mult(walls.at(i)->p4.x,X));
	        newPos = add(newPos, mult(walls.at(i)->p4.y,Y));
	        newPos = add(newPos, mult(walls.at(i)->p4.z,Z));
	        oldPos = mOldCOR;
	        oldPos = add(oldPos, mult(walls.at(i)->p4.x,oldX));
	        oldPos = add(oldPos, mult(walls.at(i)->p4.y,oldY));
	        oldPos = add(oldPos, mult(walls.at(i)->p4.z,oldZ));
	        P->ghost_particles.walls.at(i)->p4 = newPos;
	        if(C->dt <= 0.f){
	            P->ghost_particles.walls.at(i)->v4 = Vzero();
	        }
	        else{
	            P->ghost_particles.walls.at(i)->v4 = mult(1.f/C->dt, sub(newPos, oldPos));
	        }
	    #endif
	    newPos = mult(walls.at(i)->n.x,X);
	    newPos = add(newPos, mult(walls.at(i)->n.y,Y));
	    #ifdef HAVE_3D
	        newPos = add(newPos, mult(walls.at(i)->n.z,Z));
	    #endif
	    P->ghost_particles.walls.at(i)->n = newPos;
	}
	return false;
}

bool Quaternion::executeDomain()
{
	InputOutput::ProblemSetup *P = InputOutput::ProblemSetup::singleton();
	if(!P->SPH_opts.domain_motion)
        return false;
	P->SPH_opts.domain_min = add(domain_min, mCOR);
	P->SPH_opts.domain_max = add(domain_max, mCOR);
	return false;
}

bool Quaternion::_parse(xercesc::DOMElement *root)
{
	return false;
}

bool Quaternion::computePos()
{
	InputOutput::ScreenManager *S = InputOutput::ScreenManager::singleton();
	InputOutput::ProblemSetup *P = InputOutput::ProblemSetup::singleton();
	CalcServer *C = CalcServer::singleton();
	char msg[1024];
	unsigned int i;
	unsigned int N = C->N;
	//! Get original positions and normals
	vec *hRelPos = new vec[N];
	if(!hRelPos){
	    S->addMessage(3, "(Quaternion::computePos): Can't allocate temporal memory at host.\n");
	    sprintf(msg, "\t%lu bytes required.\n", N*sizeof(vec));
	    S->addMessage(0, msg);
	    return true;
	}
	if(C->getData(hRelPos, C->pos, N*sizeof(vec)))
	    return true;
	vec *hRelNormal = new vec[N];
	if(!hRelNormal){
	    S->addMessage(3, "(Quaternion::computePos): Can't allocate temporal memory at host.\n");
	    sprintf(msg, "\t%lu bytes required.\n", N*sizeof(vec));
	    S->addMessage(0, msg);
	    return true;
	}
	if(C->getData(hRelNormal, C->normal, N*sizeof(vec)))
	    return true;
	//! Project over original quaterrnion
	for(i=0;i<N;i++){
	    vec point  = hRelPos[i];
	    vec normal = hRelNormal[i];
	    vec x      = mAxis[0];
	    vec y      = mAxis[1];
	    point.x   -= mCOR.x;
	    point.y   -= mCOR.y;
	    #ifdef HAVE_3D
	        vec z    = mAxis[2];
	        point.z -= mCOR.z;
	        point.w  = 0.f;
	    #endif
	    hRelPos[i].x = dot(point,x);
	    hRelPos[i].y = dot(point,y);
	    hRelNormal[i].x = dot(normal,x);
	    hRelNormal[i].y = dot(normal,y);
	    #ifdef HAVE_3D
	        hRelPos[i].z = dot(point,z);
	        hRelPos[i].w = 0.f;
	        hRelNormal[i].z = dot(normal,z);
	        hRelNormal[i].w = 0.f;
	    #endif
	}
	//! Send data to server.
	if(C->sendData(mRelPos, hRelPos, N*sizeof(vec)))
	    return true;
	if(C->sendData(mRelNormal, hRelNormal, N*sizeof(vec)))
	    return true;
	delete[] hRelPos; hRelPos=0;
	return false;
}

bool Quaternion::computeWalls()
{
	CalcServer *C = CalcServer::singleton();
	InputOutput::ProblemSetup *P = InputOutput::ProblemSetup::singleton();
	unsigned int i;
	unsigned int N = P->ghost_particles.walls.size();
	for(i=0;i<walls.size();i++){
	    delete walls.at(i);
	}
	walls.clear();
	for(i=0;i<N;i++){
	    //! Copy wall
	    InputOutput::ProblemSetup::sphGhostParticles::Wall *wall =
	                new InputOutput::ProblemSetup::sphGhostParticles::Wall();
	    wall->p1 = P->ghost_particles.walls.at(i)->p1;
	    wall->p2 = P->ghost_particles.walls.at(i)->p2;
	    wall->n  = P->ghost_particles.walls.at(i)->n;
	    wall->v1 = P->ghost_particles.walls.at(i)->v1;
	    wall->v2 = P->ghost_particles.walls.at(i)->v2;
	    #ifdef HAVE_3D
	        wall->p3 = P->ghost_particles.walls.at(i)->p3;
	        wall->p4 = P->ghost_particles.walls.at(i)->p4;
	        wall->v3 = P->ghost_particles.walls.at(i)->v3;
	        wall->v4 = P->ghost_particles.walls.at(i)->v4;
	    #endif
	    walls.push_back(wall);
	    //! Transform vertexes and normal to COR coordinates
	    vec x      = mAxis[0];
	    vec y      = mAxis[1];
	    #ifdef HAVE_3D
	        vec z    = mAxis[2];
	    #endif
	    vec aux;
	    aux        = wall->p1;
	    wall->p1.x = dot(sub(aux, mCOR),x);
	    wall->p1.y = dot(sub(aux, mCOR),y);
	    #ifdef HAVE_3D
	        wall->p1.z = dot(sub(aux, mCOR),z);
	        wall->p1.w = 0.f;
	    #endif
	    aux        = wall->p2;
	    wall->p2.x = dot(sub(aux, mCOR),x);
	    wall->p2.y = dot(sub(aux, mCOR),y);
	    #ifdef HAVE_3D
	        wall->p2.z = dot(sub(aux, mCOR),z);
	        wall->p2.w = 0.f;
	    #endif
	    #ifdef HAVE_3D
	        aux        = wall->p3;
	        wall->p3.x = dot(sub(aux, mCOR),x);
	        wall->p3.y = dot(sub(aux, mCOR),y);
	        wall->p3.z = dot(sub(aux, mCOR),z);
	        wall->p3.w = 0.f;
	    #endif
	    #ifdef HAVE_3D
	        aux        = wall->p4;
	        wall->p4.x = dot(sub(aux, mCOR),x);
	        wall->p4.y = dot(sub(aux, mCOR),y);
	        wall->p4.z = dot(sub(aux, mCOR),z);
	        wall->p4.w = 0.f;
	    #endif
	    aux        = wall->n;
	    wall->n.x  = dot(aux,x);
	    wall->n.y  = dot(aux,y);
	    #ifdef HAVE_3D
	        wall->n.z = dot(aux,z);
	        wall->n.w = 0.f;
	    #endif
	}
	return false;
}

bool Quaternion::computeDomain()
{
	InputOutput::ProblemSetup *P = InputOutput::ProblemSetup::singleton();
	if(!P->SPH_opts.domain_motion)
        return false;
	domain_min = sub(P->SPH_opts.domain_min, mCOR);
	domain_max = sub(P->SPH_opts.domain_max, mCOR);
	return false;
}

}}} // namespaces
