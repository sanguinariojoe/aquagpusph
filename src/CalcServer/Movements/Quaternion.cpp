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

#include <CalcServer/Movements/Quaternion.h>
#include <CalcServer.h>
#include <TimeManager.h>
#include <ScreenManager.h>

#ifdef xmlAttribute
	#undef xmlAttribute
#endif
#define xmlAttribute(elem, att) XMLString::transcode( elem->getAttribute(XMLString::transcode(att)) )
using namespace xercesc;

namespace Aqua{ namespace CalcServer{ namespace Movement{

Quaternion::Quaternion()
	: Movement()
	, _pos(0)
	, _normal(0)
{
	CalcServer *C = CalcServer::singleton();
	//! Set inital values
	_cor.x = 0.f;
	_cor.y = 0.f;
	#ifdef HAVE_3D
	    _cor.z = 0.f;
	    _cor.w = 0.f;
	#endif
	_axis[0].x = 1.f;
	_axis[1].y = 1.f;
	#ifdef HAVE_3D
	    _axis[2].z = 1.f;
	#endif

	unsigned int N = C->N;
	_pos = C->allocMemory(N*sizeof(vec));
	if(!_pos)
	    exit(255);
    _normal = C->allocMemory(N*sizeof(vec));
	if(!_normal)
	    exit(255);

	if(set(_cor,_axis,true))
	    exit(255);
}

Quaternion::~Quaternion()
{
	unsigned int i;
	if(_pos)clReleaseMemObject(_pos); _pos=0;
	if(_normal)clReleaseMemObject(_normal); _normal=0;
	for(i=0;i<_walls.size();i++){
	    delete _walls.at(i);
	}
	_walls.clear();
}

bool Quaternion::set(vec cor, mat axis, bool initial){
	_cor  = cor;
	_axis = axis;
	if(initial){
	    _old_cor  = cor;
	    _old_axis = axis;
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
	X    = _axis[0];
	Y    = _axis[1];
	oldX = _old_axis[0];
	oldY = _old_axis[1];
	#ifdef HAVE_3D
	    Z    = _axis[2];
	    oldZ = _old_axis[2];
	#else
	    Z.x = 0.f; Z.y = 0.f;
	    oldZ.x = 0.f; oldZ.y = 0.f;
	#endif

	cl_int err_code=0;
	err_code |= sendArgument(_kernel,
                             0,
                             sizeof(cl_mem),
                             (void*)&(C->imove));
	err_code |= sendArgument(_kernel,
                             1,
                             sizeof(cl_mem),
                             (void*)&(C->ifluid));
	err_code |= sendArgument(_kernel,
                             2,
                             sizeof(cl_mem),
                             (void*)&(C->posin));
	err_code |= sendArgument(_kernel,
                             3,
                             sizeof(cl_mem),
                             (void*)&(C->normal));
	err_code |= sendArgument(_kernel,
                             4,
                             sizeof(cl_mem),
                             (void*)&(C->vin));
	err_code |= sendArgument(_kernel,
                             5,
                             sizeof(cl_mem),
                             (void*)&(C->densin));
	err_code |= sendArgument(_kernel,
                             6,
                             sizeof(cl_mem),
                             (void*)&(C->mass));
	err_code |= sendArgument(_kernel,
                             7,
                             sizeof(cl_mem),
                             (void*)&(C->hpin));
	err_code |= sendArgument(_kernel,
                             8,
                             sizeof(cl_mem),
                             (void*)&(_pos));
	err_code |= sendArgument(_kernel,
                             9,
                             sizeof(cl_mem),
                             (void*)&(_normal));
	err_code |= sendArgument(_kernel,
                             10,
                             sizeof(cl_uint),
                             (void*)&(C->N));
	err_code |= sendArgument(_kernel,
                             11,
                             sizeof(cl_float),
                             (void*)&(C->dt));
	err_code |= sendArgument(_kernel,
                             12,
                             sizeof(vec),
                             (void*)&(_cor));
	err_code |= sendArgument(_kernel,
                             13,
                             sizeof(vec),
                             (void*)&(X));
	err_code |= sendArgument(_kernel,
                             14,
                             sizeof(vec),
                             (void*)&(Y));
	err_code |= sendArgument(_kernel,
                             15,
                             sizeof(vec),
                             (void*)&(Z));
	err_code |= sendArgument(_kernel,
                             16,
                             sizeof(vec),
                             (void*)&(_old_cor));
	err_code |= sendArgument(_kernel,
                             17,
                             sizeof(vec),
                             (void*)&(oldX));
	err_code |= sendArgument(_kernel,
                             18,
                             sizeof(vec),
                             (void*)&(oldY));
	err_code |= sendArgument(_kernel,
                             19,
                             sizeof(vec),
                             (void*)&(oldZ));
	if(err_code != CL_SUCCESS) {
		S->addMessageF(3, "Failure sending variables to the kernel.\n");
	    return true;
	}
	//! Execute the kernel
	#ifdef HAVE_GPUPROFILE
	    cl_event event;
	    cl_ulong end, start;
	    profileTime(0.f);
	    err_code = clEnqueueNDRangeKernel(C->command_queue,
                                          _kernel,
                                          1,
                                          NULL,
                                          &_global_work_size,
                                          NULL,
                                          0,
                                          NULL,
                                          &event);
	#else
	    err_code = clEnqueueNDRangeKernel(C->command_queue,
                                          _kernel,
                                          1,
                                          NULL,
                                          &_global_work_size,
                                          NULL,
                                          0,
                                          NULL,
                                          NULL);
	#endif
	if(err_code != CL_SUCCESS) {
		S->addMessage(3, "I cannot execute the kernel.\n");
        S->printOpenCLError(err_code);
        return true;
    }
    #ifdef HAVE_GPUPROFILE
        err_code = clWaitForEvents(1, &event);
        if(err_code != CL_SUCCESS) {
            S->addMessage(3, "Impossible to wait for the velocity clamping kernel end.\n");
            S->printOpenCLError(err_code);
            return true;
        }
        err_code = clGetEventProfilingInfo(event,
                                           CL_PROFILING_COMMAND_END,
                                           sizeof(cl_ulong),
                                           &end,
                                           0);
        if(err_code != CL_SUCCESS) {
            S->addMessage(3, "I cannot profile the velocity clamping kernel execution.\n");
            S->printOpenCLError(err_code);
            return true;
        }
        err_code = clGetEventProfilingInfo(event,
                                           CL_PROFILING_COMMAND_START,
                                           sizeof(cl_ulong),
                                           &start,
                                           0);
        if(err_code != CL_SUCCESS) {
            S->addMessage(3, "I cannot profile the velocity clamping kernel execution.\n");
            S->printOpenCLError(err_code);
            return true;
        }
	    profileTime(profileTime() + (end - start)/1000.f);  // 10^-3 ms
	#endif

	if(executeWalls())
	    return true;
	if(executeDomain())
	    return true;

	_old_cor  = _cor;
	_old_axis = _axis;
	return false;
}

bool Quaternion::executeWalls()
{
	unsigned int i;
	unsigned int N = _walls.size();
	if(!N)
	    return false;
	InputOutput::ProblemSetup *P = InputOutput::ProblemSetup::singleton();
	CalcServer *C = CalcServer::singleton();
	vec X, Y, Z, oldX, oldY, oldZ;
	X    = _axis[0];
	Y    = _axis[1];
	oldX = _old_axis[0];
	oldY = _old_axis[1];
	#ifdef HAVE_3D
	    Z    = _axis[2];
	    oldZ = _old_axis[2];
	#else
	    Z.x = 0.f; Z.y = 0.f;
	    oldZ.x = 0.f; oldZ.y = 0.f;
	#endif
	for(i=0;i<N;i++){
	    vec newPos, oldPos;
	    newPos = _cor;
	    newPos = add(newPos, mult(_walls.at(i)->p1.x,X));
	    newPos = add(newPos, mult(_walls.at(i)->p1.y,Y));
	    #ifdef HAVE_3D
	        newPos = add(newPos, mult(_walls.at(i)->p1.z,Z));
	    #endif
	    oldPos = _old_cor;
	    oldPos = add(oldPos, mult(_walls.at(i)->p1.x,oldX));
	    oldPos = add(oldPos, mult(_walls.at(i)->p1.y,oldY));
	    #ifdef HAVE_3D
	        oldPos = add(oldPos, mult(_walls.at(i)->p1.z,oldZ));
	    #endif
	    P->ghost_particles.walls.at(i)->p1 = newPos;
	    if(C->dt <= 0.f){
	        P->ghost_particles.walls.at(i)->v1 = Vzero();
	    }
	    else{
	        P->ghost_particles.walls.at(i)->v1 = mult(1.f/C->dt, sub(newPos, oldPos));
	    }
	    newPos = _cor;
	    newPos = add(newPos, mult(_walls.at(i)->p2.x,X));
	    newPos = add(newPos, mult(_walls.at(i)->p2.y,Y));
	    #ifdef HAVE_3D
	        newPos = add(newPos, mult(_walls.at(i)->p2.z,Z));
	    #endif
	    oldPos = _old_cor;
	    oldPos = add(oldPos, mult(_walls.at(i)->p2.x,oldX));
	    oldPos = add(oldPos, mult(_walls.at(i)->p2.y,oldY));
	    #ifdef HAVE_3D
	        oldPos = add(oldPos, mult(_walls.at(i)->p2.z,oldZ));
	    #endif
	    P->ghost_particles.walls.at(i)->p2 = newPos;
	    if(C->dt <= 0.f){
	        P->ghost_particles.walls.at(i)->v2 = Vzero();
	    }
	    else{
	        P->ghost_particles.walls.at(i)->v2 = mult(1.f/C->dt, sub(newPos, oldPos));
	    }
	    #ifdef HAVE_3D
	        newPos = _cor;
	        newPos = add(newPos, mult(_walls.at(i)->p3.x,X));
	        newPos = add(newPos, mult(_walls.at(i)->p3.y,Y));
	        newPos = add(newPos, mult(_walls.at(i)->p3.z,Z));
	        oldPos = _old_cor;
	        oldPos = add(oldPos, mult(_walls.at(i)->p3.x,oldX));
	        oldPos = add(oldPos, mult(_walls.at(i)->p3.y,oldY));
	        oldPos = add(oldPos, mult(_walls.at(i)->p3.z,oldZ));
	        P->ghost_particles.walls.at(i)->p3 = newPos;
	        if(C->dt <= 0.f){
	            P->ghost_particles.walls.at(i)->v3 = Vzero();
	        }
	        else{
	            P->ghost_particles.walls.at(i)->v3 = mult(1.f/C->dt, sub(newPos, oldPos));
	        }
	    #endif
	    #ifdef HAVE_3D
	        newPos = _cor;
	        newPos = add(newPos, mult(_walls.at(i)->p4.x,X));
	        newPos = add(newPos, mult(_walls.at(i)->p4.y,Y));
	        newPos = add(newPos, mult(_walls.at(i)->p4.z,Z));
	        oldPos = _old_cor;
	        oldPos = add(oldPos, mult(_walls.at(i)->p4.x,oldX));
	        oldPos = add(oldPos, mult(_walls.at(i)->p4.y,oldY));
	        oldPos = add(oldPos, mult(_walls.at(i)->p4.z,oldZ));
	        P->ghost_particles.walls.at(i)->p4 = newPos;
	        if(C->dt <= 0.f){
	            P->ghost_particles.walls.at(i)->v4 = Vzero();
	        }
	        else{
	            P->ghost_particles.walls.at(i)->v4 = mult(1.f/C->dt, sub(newPos, oldPos));
	        }
	    #endif
	    newPos = mult(_walls.at(i)->n.x,X);
	    newPos = add(newPos, mult(_walls.at(i)->n.y,Y));
	    #ifdef HAVE_3D
	        newPos = add(newPos, mult(_walls.at(i)->n.z,Z));
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
	P->SPH_opts.domain_min = add(_domain_min, _cor);
	P->SPH_opts.domain_max = add(_domain_max, _cor);
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

	vec *hRelPos = new vec[N];
	if(!hRelPos){
	    S->addMessageF(3, "I Cannot allocate temporal memory at host.\n");
	    sprintf(msg, "\t%lu bytes required.\n", N*sizeof(vec));
	    S->addMessage(0, msg);
	    return true;
	}
	if(C->getData(hRelPos, C->pos, N*sizeof(vec)))
	    return true;
	vec *hRelNormal = new vec[N];
	if(!hRelNormal){
	    S->addMessageF(3, "I Cannot allocate temporal memory at host.\n");
	    sprintf(msg, "\t%lu bytes required.\n", N*sizeof(vec));
	    S->addMessage(0, msg);
	    return true;
	}
	if(C->getData(hRelNormal, C->normal, N*sizeof(vec)))
	    return true;

	for(i=0;i<N;i++){
	    vec point  = hRelPos[i];
	    vec normal = hRelNormal[i];
	    vec x      = _axis[0];
	    vec y      = _axis[1];
	    point.x   -= _cor.x;
	    point.y   -= _cor.y;
	    #ifdef HAVE_3D
	        vec z    = _axis[2];
	        point.z -= _cor.z;
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

	if(C->sendData(_pos, hRelPos, N*sizeof(vec)))
	    return true;
	if(C->sendData(_normal, hRelNormal, N*sizeof(vec)))
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
	for(i=0;i<_walls.size();i++){
	    delete _walls.at(i);
	}
	_walls.clear();
	for(i=0;i<N;i++){

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
	    _walls.push_back(wall);

	    vec x      = _axis[0];
	    vec y      = _axis[1];
	    #ifdef HAVE_3D
	        vec z    = _axis[2];
	    #endif
	    vec aux;
	    aux        = wall->p1;
	    wall->p1.x = dot(sub(aux, _cor),x);
	    wall->p1.y = dot(sub(aux, _cor),y);
	    #ifdef HAVE_3D
	        wall->p1.z = dot(sub(aux, _cor),z);
	        wall->p1.w = 0.f;
	    #endif
	    aux        = wall->p2;
	    wall->p2.x = dot(sub(aux, _cor),x);
	    wall->p2.y = dot(sub(aux, _cor),y);
	    #ifdef HAVE_3D
	        wall->p2.z = dot(sub(aux, _cor),z);
	        wall->p2.w = 0.f;
	    #endif
	    #ifdef HAVE_3D
	        aux        = wall->p3;
	        wall->p3.x = dot(sub(aux, _cor),x);
	        wall->p3.y = dot(sub(aux, _cor),y);
	        wall->p3.z = dot(sub(aux, _cor),z);
	        wall->p3.w = 0.f;
	    #endif
	    #ifdef HAVE_3D
	        aux        = wall->p4;
	        wall->p4.x = dot(sub(aux, _cor),x);
	        wall->p4.y = dot(sub(aux, _cor),y);
	        wall->p4.z = dot(sub(aux, _cor),z);
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
	_domain_min = sub(P->SPH_opts.domain_min, _cor);
	_domain_max = sub(P->SPH_opts.domain_max, _cor);
	return false;
}

}}} // namespaces
