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

#ifndef INPUTOUTPUT_H_INCLUDED
#define INPUTOUTPUT_H_INCLUDED

namespace Aqua{
namespace InputOutput{

/** \class InputOutput InputOutput.h InputOutput/InputOutput.h
 *  Input/Output base class.
 */
class InputOutput
{
public:
	/** Constructor
	 */
	InputOutput(){};

	/** Destructor
	 */
	virtual ~InputOutput(){};

    /** Save the data
     * @return false if all gone right, true otherwise.
     */
    virtual bool save() = 0;

    /** Load the data
     * @return false if all gone right, true otherwise.
     */
    virtual bool load() = 0;

};  // class InputOutput

}}  // namespaces

#endif // INPUTOUTPUT_H_INCLUDED
