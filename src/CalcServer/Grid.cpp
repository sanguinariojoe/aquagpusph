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

#include <ProblemSetup.h>
#include <ScreenManager.h>
#include <CalcServer/Grid.h>
#include <CalcServer.h>

namespace Aqua{ namespace CalcServer{

Grid::Grid()
    : Kernel("Grid")
	, _maximum(NULL)
	, _minimum(NULL)
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
        _maximum = new Reduction(C->pos, C->N, "vec",
                                 "(vec)(-INFINITY,-INFINITY,-INFINITY,0.f)",
                                 operation);
    #else
        _maximum = new Reduction(C->pos, C->N, "vec",
                                 "(vec)(-INFINITY,-INFINITY)",
                                 operation);
    #endif
    strcpy(operation, "");
    strcat(operation, "c.x = (a.x < b.x) ? a.x : b.x;\n");
    strcat(operation, "\tc.y = (a.y < b.y) ? a.y : b.y;\n");
    #ifdef HAVE_3D
        strcat(operation, "\tc.z = (a.z < b.z) ? a.z : b.z;\n");
        strcat(operation, "\tc.w = 0.f;\n");
        _minimum = new Reduction(C->pos, C->N, "vec",
                                 "(vec)(INFINITY,INFINITY,INFINITY,0.f)",
                                 operation);
    #else
        _minimum = new Reduction(C->pos, C->N, "vec",
                                 "(vec)(INFINITY,INFINITY)",
                                 operation);
    #endif
	S->addMessageF(1, "Grid ready to work!\n");
}

Grid::~Grid()
{
	InputOutput::ScreenManager *S = InputOutput::ScreenManager::singleton();
	S->addMessageF(1, "Destroying grid maximum position reduction processor...\n");
	if(_maximum)delete _maximum; _maximum=NULL;
	S->addMessageF(1, "Destroying grid minimum position reduction processor...\n");
	if(_minimum)delete _minimum; _minimum=NULL;
}

bool Grid::execute()
{
	InputOutput::ScreenManager *S = InputOutput::ScreenManager::singleton();
	CalcServer *C = CalcServer::singleton();
	unsigned int i;
	int err_code=0;
	float sep=2.f;
	#if defined(__GAUSS_KERNEL_TYPE__)
		sep = 3.f;
	#endif
    cl_mem reduced = _maximum->execute();
    if(!reduced)
        return true;
	if(C->getData((void *)&C->pos_max, reduced, sizeof(vec)))
		return true;
    reduced = _minimum->execute();
    if(!reduced)
        return true;
	if(C->getData((void *)&C->pos_min, reduced, sizeof(vec)))
		return true;
	C->cell_length = 1.f / (C->cell_length_factor * sep * C->h);
	C->num_cells_vec.x = (unsigned int)((C->pos_max.x - C->pos_min.x) * C->cell_length)+(unsigned int)6;
	C->num_cells_vec.y = (unsigned int)((C->pos_max.y - C->pos_min.y) * C->cell_length)+(unsigned int)6;
	C->num_cells = C->num_cells_vec.x * C->num_cells_vec.y;
	#ifdef HAVE_3D
		C->num_cells_vec.z = (unsigned int)((C->pos_max.z - C->pos_min.z) * C->cell_length)+(unsigned int)6;
		C->num_cells_vec.w = 0;
		C->num_cells *= C->num_cells_vec.z;
	#endif
	return false;
}

}}  // namespace
