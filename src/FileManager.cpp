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
#include <string.h>

#include <FileManager.h>
#include <ScreenManager.h>
#include <ProblemSetup.h>
#include <InputOutput/ASCII.h>

namespace Aqua{ namespace InputOutput{

FileManager::FileManager()
	: _in_file(NULL)
	, _log(NULL)
	, _energy(NULL)
	, _bounds(NULL)
{
    inputFile("Input.xml");
    _state = new State();
    _log = new Log();
}

FileManager::~FileManager()
{
    unsigned int i;

    if(_in_file)
        delete[] _in_file;
    _in_file = NULL;
    if(_state)
        delete _state;
    _state = NULL;
    if(_log)
        delete _log;
    _log = NULL;
    if(_energy)
        delete _energy;
    _energy = NULL;
    if(_bounds)
        delete _bounds;
    _bounds = NULL;

    for(i=0; i<_loaders.size(); i++){
        delete _loaders.at(i);
        _loaders.at(i) = NULL;
    }
    _loaders.clear();

    for(i=0; i<_savers.size(); i++){
        delete _savers.at(i);
        _savers.at(i) = NULL;
    }
    _savers.clear();
}

void FileManager::inputFile(const char* path)
{
    size_t len;

    if(_in_file)
        delete[] _in_file;
    _in_file = NULL;

    if(!path)
        return;

    len = strlen(path) + 1;
    _in_file = new char[len];
    strcpy(_in_file, path);
}

FILE* FileManager::energyFile()
{
    if(_energy)
        return _energy->fileHandler();
}

FILE* FileManager::boundsFile()
{
    if(_bounds)
        return _bounds->fileHandler();
}

Fluid* FileManager::load()
{
    unsigned int i, n=0;
    char msg[1024];
    ProblemSetup *P = ProblemSetup::singleton();
    ScreenManager *S = ScreenManager::singleton();

    // Load the XML definition file
    if(_state->load()){
        return NULL;
    }

    // Setup the problem setup
	if(P->perform()) {
		return NULL;
	}

    // Build the additional reporters requested
    if(P->time_opts.energy_mode != __NO_OUTPUT_MODE__){
        _energy = new Energy();
        if(!_energy)
            return NULL;
    }
    if(P->time_opts.bounds_mode != __NO_OUTPUT_MODE__){
        _bounds = new Bounds();
        if(!_bounds)
            return NULL;
    }

    // Now we can build the fluid manager and the particles loaders/savers
    Fluid *F = new Fluid();
    if(!F)
        return NULL;
    for(i=0; i<P->n_fluids; i++){
        if(!strcmp(P->fluids[i].in_format, "ASCII")){
            ASCII *loader = new ASCII(n, P->fluids[i].n, i);
            if(!loader)
                return NULL;
            n += P->fluids[i].n;
            _loaders.push_back((Particles*)loader);
        }
        else{
            sprintf(msg,
                    "Unknow \"%s\" file format.\n",
                    P->fluids[i].in_format);
            S->addMessageF(3, msg);
            return NULL;
        }
        if(!strcmp(P->fluids[i].out_format, "ASCII")){
            ASCII *saver = new ASCII(n, P->fluids[i].n, i);
            if(!saver)
                return NULL;
            n += P->fluids[i].n;
            _savers.push_back((Particles*)saver);
        }
        else{
            sprintf(msg,
                    "Unknow \"%s\" file format.\n",
                    P->fluids[i].out_format);
            S->addMessageF(3, msg);
            return NULL;
        }
    }

    // Execute the loaders
    for(i=0; i<P->n_fluids; i++){
        if(_loaders.at(i)->load())
            return NULL;
    }

    return F;
}

}}  // namespace
