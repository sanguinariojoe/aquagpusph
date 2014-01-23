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
#include <iostream>
#include <sstream>
#include <stdexcept>
#include <math.h>
#include <string>
#include <string.h>
#include <sys/time.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <list>
#include <unistd.h>
#include <errno.h>

#include <FileManager.h>
#include <ScreenManager.h>
#include <ProblemSetup.h>

namespace Aqua{ namespace InputOutput{

FileManager::FileManager()
	: _in_file(NULL)
	, _log(NULL)
{
    inputFile("Input.xml");
    _state = new State();
    _log = new Log();
}

FileManager::~FileManager()
{
    if(_in_file)
        delete[] _in_file;
    _in_file = NULL;
    if(_state)
        delete _state;
    _state = NULL;
    if(_log)
        delete _log;
    _log = NULL;
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

}}  // namespace
