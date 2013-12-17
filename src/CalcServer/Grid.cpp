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
// Include the Problem setup manager header
// ----------------------------------------------------------------------------
#include <ProblemSetup.h>

// ----------------------------------------------------------------------------
// Include the Problem setup manager header
// ----------------------------------------------------------------------------
#include <ScreenManager.h>

// ----------------------------------------------------------------------------
// Include the main header
// ----------------------------------------------------------------------------
#include <CalcServer/Grid.h>

// ----------------------------------------------------------------------------
// Include the calculation server
// ----------------------------------------------------------------------------
#include <CalcServer.h>

namespace Aqua{ namespace CalcServer{

Grid::Grid()
    : Kernel("Grid")
	, maximum(NULL)
	, minimum(NULL)
{
	InputOutput::ScreenManager *S = InputOutput::ScreenManager::singleton();
	CalcServer *C = CalcServer::singleton();
    char operation[512];
    strcpy(operation, "");
    strcat(operation, "c.x = (a.x > b.x) ? a.x : b.x;\n");
    strcat(operation, "\tc.y = (a.y > b.y) ? a.y : b.y;\n");
    #ifdef HAVE_3D
        strcat(operation, "\tc.z = (a.z > b.z) ? a.z : b.z;\n");
        strcat(operation, "\tc.w = 0.f;\n");
        maximum = new Reduction(C->pos, C->N, "vec", "(vec)(-INFINITY,-INFINITY,-INFINITY,0.f)", operation);
    #else
        maximum = new Reduction(C->pos, C->N, "vec", "(vec)(-INFINITY,-INFINITY)", operation);
    #endif
    strcpy(operation, "");
    strcat(operation, "c.x = (a.x < b.x) ? a.x : b.x;\n");
    strcat(operation, "\tc.y = (a.y < b.y) ? a.y : b.y;\n");
    #ifdef HAVE_3D
        strcat(operation, "\tc.z = (a.z < b.z) ? a.z : b.z;\n");
        strcat(operation, "\tc.w = 0.f;\n");
        minimum = new Reduction(C->pos, C->N, "vec", "(vec)(INFINITY,INFINITY,INFINITY,0.f)", operation);
    #else
        minimum = new Reduction(C->pos, C->N, "vec", "(vec)(INFINITY,INFINITY)", operation);
    #endif
	S->addMessage(1, "(Grid::Grid): Grid ready to work!\n");
}

Grid::~Grid()
{
	InputOutput::ScreenManager *S = InputOutput::ScreenManager::singleton();
	S->addMessage(1, "(Grid::~Grid): Destroying grid maximum position reduction processor...\n");
	if(maximum)delete maximum; maximum=NULL;
	S->addMessage(1, "(Grid::~Grid): Destroying grid minimum position reduction processor...\n");
	if(minimum)delete minimum; minimum=NULL;
}

bool Grid::execute()
{
	InputOutput::ScreenManager *S = InputOutput::ScreenManager::singleton();
	CalcServer *C = CalcServer::singleton();
	unsigned int i;
	int clFlag=0;
	float sep=2.f;
	#if defined(__GAUSS_KERNEL_TYPE__)
		sep = 3.f;
	#endif
    cl_mem reduced = maximum->execute();
    if(!reduced)
        return true;
	if(C->getData((void *)&C->posmax, reduced, sizeof(vec)))
		return true;
    reduced = minimum->execute();
    if(!reduced)
        return true;
	if(C->getData((void *)&C->posmin, reduced, sizeof(vec)))
		return true;
	C->rdist = 1.0f / (C->CellFac * sep * C->h);
	C->lvec.x   = (unsigned int)((C->posmax.x - C->posmin.x) * C->rdist)+(unsigned int)6;
	C->lvec.y   = (unsigned int)((C->posmax.y - C->posmin.y) * C->rdist)+(unsigned int)6;
	C->lxy	  = C->lvec.x * C->lvec.y;
	#ifdef HAVE_3D
		C->lvec.z  = (unsigned int)((C->posmax.z - C->posmin.z) * C->rdist)+(unsigned int)6;
		C->lvec.w  = 0;
		C->lxy	*= C->lvec.z;
	#endif
	return false;
}

}}  // namespace
