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

#ifndef TOOL_H_INCLUDED
#define TOOL_H_INCLUDED

#include <sphPrerequisites.h>

namespace Aqua{ namespace CalcServer{

/** @class Tool Tool.h CalcServer/Tool.h
 * @brief Tools base class. The way that AQUAgpusph compute each problem is set
 * through a set of tools that are computed sequentially. Several tools can be
 * considered, for instance:
 *   -# A single OpenCL kernel
 *   -# A more complex OpenCL tool like Reductions or LinkList
 *   -# Python scripts
 *   -# Variables set
 */
class Tool
{
public:
	/** Constructor.
	 * @param tool_name Name of the tool. Useful to identify errors.
	 */
	Tool(const char* tool_name);

	/** Destructor
	 */
	~Tool();

	/** Set the tool name.
	 * @param tool_name Tool name.
	 */
	void name(const char* tool_name);

	/** Get the tool name.
	 * @return Tool name.
	 */
	const char* name(){return (const char*)_name;}

    /** Execute the tool.
     * @return false if all gone right, true otherwise.
     */
    virtual bool execute()=0;

private:
	/// Kernel name
	char* _name;
};

}}  // namespace

#endif // TOOL_H_INCLUDED
