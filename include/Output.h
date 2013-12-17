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

#ifndef OUTPUT_H_INCLUDED
#define OUTPUT_H_INCLUDED

// ----------------------------------------------------------------------------
// Include Prerequisites
// ----------------------------------------------------------------------------
#include <sphPrerequisites.h>

// ----------------------------------------------------------------------------
// Include some auxiliar methods
// ----------------------------------------------------------------------------
#include <AuxiliarMethods.h>

namespace Aqua{ namespace InputOutput{

/** Prints the output, with the formats selected at input.
 * @return false if all gone right. \n true otherwise.
 */
bool output();

#ifdef HAVE_H5PART
	/** Reassembly output H5Part data files.
	 * @return false if all gone right. \n false otherwise
	 */
	bool assemblyH5Part();
#endif // HAVE_H5PART

/// @namespace Output namespace to hide the methods included.
namespace Output{
	/** Prints a tecplot file.
	 * @return false if all gone right. \n false otherwise
	 */
	bool printTecplot();

	#ifdef HAVE_H5PART
	    /** Fills the H5Part file.
	     * @return false if all gone right. \n false otherwise
	     */
	    bool printH5Part();
	#endif // HAVE_H5PART

	#ifdef HAVE_VTK
	    /** Writes VTK output file.
	     * @return false if all gone right. \n false otherwise
	     */
	    bool printVTK();
	#endif
}   // namespace Output

}}  // namespace

#endif // OUTPUT_H_INCLUDED
