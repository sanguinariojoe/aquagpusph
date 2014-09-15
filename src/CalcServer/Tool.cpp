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

#include <CalcServer/Tool.h>
#include <CalcServer.h>

namespace Aqua{ namespace CalcServer{

Tool::Tool(const char* tool_name)
	: _name(0)
{
	name(tool_name);
}

Tool::~Tool()
{
	if(_name) delete[] _name; _name=0;
}

void Tool::name(const char* tool_name)
{
	if(_name) delete[] _name; _name=0;
	_name = new char[strlen(tool_name)+1];
	strcpy(_name, tool_name);
}

}}  // namespace
