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

#ifndef INPUT_H_INCLUDED
#define INPUT_H_INCLUDED

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

// ----------------------------------------------------------------------------
// Include the h5part library
// ----------------------------------------------------------------------------
#ifdef HAVE_H5PART
	#include <H5Part.h>
	#include <hdf5.h>
#endif // HAVE_H5PART

namespace Aqua{ namespace InputOutput{

/** Read files to load the fluid state.
 * @return true if all gone right. \n False otherwise.
 */
bool input();

/// @namespace Input namespace to hide the methods included.
namespace Input{

	#ifdef HAVE_H5PART
	    /** Reads H5Part to set the fluid.
	     * @return true if all gone right. \n False otherwise.
	     */
	    bool loadH5Part();
	#endif // HAVE_H5PART

}   // namespace Input

}}  // namespaces

#endif // INPUT_H_INCLUDED
